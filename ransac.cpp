//#include "stdafx.h"
#include <Windows.h>
#include "ransac.h"
#include "publicfunction.h"

#include "math.h"


#define LM_CIRCUL 200
#define RANSAC_CIRCUL 500
#define MAX_DATA 100000000
#define RAND_BLOCK 8
#define LEASTTHREAD 0.01


/************************************************************************/
/* 函数声明                                                             */
/************************************************************************/
double** cal_trans(double left_point[4][2],double right_point[4][2]);
double cal_mid_error(MatchPoint* pot,double** trans_h,double* dv,int num);
MatchPoints LMedS(MatchPoints pots,int imgwidth,int imgheight,int matnum,int& inline_num);
MatchPoints sift_ransac(MatchPoints pots,int imgwidth,int imgheight,int matnum,int& inline_num);
int randpoint(MatchPoint* pot,double left_point[4][2],double right_point[4][2],int block_width,int block_height,int matnum);
void randpoint1(MatchPoints** pots,double left_point[4][2],double right_point[4][2]);
void randpoint2(MatchPoints pots,double left_point[4][2],double right_point[4][2]);
MatchPoints** reognizepoint(MatchPoints pots, int imgwidth, int imgheight);
MatchPoint* savpoint(MatchPoint* oldpoint,int* weight,int matnum,int inline_num);


//流程参考尚明珠论文《一种基于特征的全自动图像拼接算法》
/****************************************************** 
* 函数名称： 
* de_gross_error()
*  
* 函数参数：
* MatchPoints matchpoints       - 剔除粗差前得点
* MatchPoints &new_matchpoints  - 剔除粗差后的点
* int imgwidth                  - 匹配两张影像中左影像的宽(因为是将左影像分为RAND_BLOCK*RAND_BLOCK个块进行随机抽样的)
* int imgheight                 - 匹配两张影像中左影像的高
* int &inline_num               - 剔除粗差后点的个数
*
* 返回值：
* 图像间的投影变换系数
*  
*说明：利用RANSAC和LMedS算法剔除粗差，模型为图像的投影变换系数
********************************************************/
MatchPoints del_gross_error(MatchPoints matchpoints,int imgwidth,int imgheight) 
{
	MatchPoints new_matchpoints;
	int inline_num=0;
	
	/************************************************************************/
	/* RANSAC算法剔除粗差                                                   */
	/************************************************************************/
	new_matchpoints=sift_ransac(matchpoints,imgwidth,imgheight,matchpoints.mat_key_num,inline_num);
	new_matchpoints.mat_key_num=inline_num;
	
	/************************************************************************/
	/* LMedS最小中值平方法剔除粗差                                          */
	/************************************************************************/
	//inline_num=0;
	//new_matchpoints=LMedS(new_matchpoints,imgwidth,imgheight,new_matchpoints.mat_key_num,inline_num);
	//new_matchpoints.mat_key_num=inline_num;
	
	return new_matchpoints;
}
/****************************************************** 
* 函数名称： 
* savpoint()
*  
* 函数参数：
* MatchPoint* oldpoint         - 两张图像的匹配点
* int* weight                  - 点的权重，1为内点，否则为外点
* int matnum                   - 剔除粗差前匹配点的个数
* int inline_num               - 内点的个数
* 
* 返回值：
* 剔除粗差后匹配点
*  
*说明：计算剔除粗差后的匹配点
********************************************************/
MatchPoint* savpoint(MatchPoint* oldpoint,int* weight,int matnum,int inline_num)
{
	MatchPoint* newpoint=new MatchPoint[matnum];
	int num=0;
	for (int j=0;j<matnum;j++)
	{
		if (weight[j]==1)
		{	
			newpoint[num].imageX1=oldpoint[j].imageX1;
			newpoint[num].imageY1=oldpoint[j].imageY1;
			newpoint[num].imageX2=oldpoint[j].imageX2;
			newpoint[num].imageY2=oldpoint[j].imageY2;
			newpoint[num].id1=oldpoint[j].id1;
			newpoint[num].id2=oldpoint[j].id2;

			num++;
		}
	}
	return newpoint;
}
/****************************************************** 
* 函数名称： 
* LMedS()
*  
* 函数参数：
* MatchPoints* pot             - 两张图像的匹配点
* ImageFile imagefile          - 图像信息
* int imgwidth                 - 匹配两张影像中左影像的宽
* int imgheight                - 匹配两张影像中左影像的高
* int& inline_num              - 剔除粗差后的内点数
* 
* 返回值：
* LMedS剔除粗差后的匹配点
*  
*说明：LMedS(最小中值平方法)剔除粗差
********************************************************/
MatchPoints LMedS(MatchPoints pots,int imgwidth,int imgheight,int matnum,int& inline_num)
{	
	double** dv=new double*[LM_CIRCUL];   //投影变换误差
	for (int k=0;k<LM_CIRCUL;k++)
	{
		dv[k]=new double[matnum];
		memset(dv[k],0,sizeof(double)*matnum);
	}
	int block_width=imgwidth/RAND_BLOCK;     //图像分为LM_CIRCUL个区域后，每个区域列的宽
	int block_height=imgheight/RAND_BLOCK;   //图像分为LM_CIRCUL个区域后，每个区域行的高

	/************************************************************************/
	/* 随机抽样，计算投影系数，根据投影系数计算误差和中值                   */
	/************************************************************************/
	double* med_dv=new double[LM_CIRCUL];//投影变换误差的中值
	double left_point[4][2];//用于计算单映矩阵的4对匹配点左影像上的点
	double right_point[4][2];//用于计算单映矩阵的4对匹配点右影像上的点

	for (int k=0;k<LM_CIRCUL;k++)//循环LM_CIRCUL次，抽取样本
	{
		/******************************************************************************/
		/*将图像均匀分为BLOCK行BLOCK列个区域,随机抽取4个区，并从每个区中随即抽取一个角点*/
		/******************************************************************************/
		if (!randpoint(pots.matpoint,left_point,right_point,block_width,block_height,matnum))
		{	
			med_dv[k]=MAX_DATA;
			continue;//找四个匹配点失败,则退出此次循环
		}
		/**********************************************************************************************/
		/* 由4对匹配点计算投影变换矩阵（即单映性矩阵）,利用该矩阵进行投影变换，计算每个点的投影误差   */
		/**********************************************************************************************/
		double** trans_h=cal_trans(left_point,right_point);
		if (trans_h)
		{
			med_dv[k]=cal_mid_error(pots.matpoint,trans_h,dv[k],matnum);//利用旋转矩阵计算所有点旋转后的坐标残差，并计算误差距离
			
			for (int i=0;i<3;i++)//释放内存
			{
				delete []trans_h[i];
			}
			delete []trans_h;
		} 
		else
		{
			med_dv[k]=MAX_DATA;//计算旋转矩阵H失败
		}
	}//LM循环

	/************************************************************************/
	/* 寻找所有中值的最小值，以及其对应的旋转矩阵                           */
	/************************************************************************/
	int minnum=0;
	double* temp_dev=new double[LM_CIRCUL];
	memcpy(temp_dev,med_dv,sizeof(double)*LM_CIRCUL);
	BubbleSort(temp_dev,LM_CIRCUL,minnum);
	//释放内存
	delete []temp_dev;
	//cout<<"最小中值："<<med_dv[minnum]<<endl;
	if (med_dv[minnum]==MAX_DATA)
	{//所有点的中值计算均失败,则RANSAC剔除粗差失败，返回原始的匹配点
		delete []med_dv;//释放内存
		return pots; 
	}
	/************************************************************************/
	/* 计算中值的阈值(公式来自王琳论文) 并判断内点                          */
	/************************************************************************/
	double sigma=1.4826*(1+5.0/(matnum-4))*sqrt(/*sqrt(*/med_dv[minnum]/*)*/);//当matnum>=2*4时可以采用此公式！！！！！！！！
	delete []med_dv;

	int* weight=new int[matnum];
	memset(weight,0,sizeof(int)*matnum);
	inline_num=0;//内点的个数

	//根据阈值判断内点
	for (int j=0;j<matnum;j++)
	{
		if (dv[minnum][j]<2.5*sigma*2.5*sigma)
		{
			weight[j]=1;//是内点
			inline_num++;
		}
	}
	//释放内存
	for (int k=0;k<LM_CIRCUL;k++)
	{
		delete []dv[k];
	}
	delete []dv; 
	/************************************************************************/
	/* 保存剔除粗差后的点                                                   */
	/************************************************************************/
	MatchPoints new_pot;
	new_pot.matpoint=savpoint(pots.matpoint,weight,matnum,inline_num);
	//释放内存
	delete []weight;
	return new_pot;

}
/****************************************************** 
* 函数名称： 
* sift_ransac()
*  
* 函数参数：
* MatchPoints* pot             - 两张图像的匹配点
* int imgwidth                 - 匹配两张影像中左影像的宽
* int imgheight                - 匹配两张影像中左影像的高
* int matnum                   - 匹配点的个数
* int& inline_num              - 剔除粗差后的内点数
* 
* 返回值：
* RANSAC剔除粗差后的匹配点
*  
*说明：RANSAC剔除粗差
********************************************************/
MatchPoints sift_ransac(MatchPoints pots,int imgwidth,int imgheight,int matnum,int& inline_num)
{
	double left_point[4][2];       //用于计算单映矩阵的4对匹配点左影像上的点
	double right_point[4][2];      //用于计算单映矩阵的4对匹配点右影像上的点
	double* dv=new double[matnum]; //投影变换误差
	memset(dv,0,sizeof(double)*matnum);
	
	/************************************************************************/
	/* 以格网的形式重新组织匹配点                                           */
	/************************************************************************/
	/*
	MatchPoints** re_pots=reognizepoint(pots,imgwidth,imgheight);
	//块中的点数不为0的块数
	int zero_num=0;
	for (int i=0;i<RAND_BLOCK;i++)
	{
		for (int j=0;j<RAND_BLOCK;j++)
		{
			if (re_pots[i][j].mat_key_num>0)
			{
				zero_num++;
			}
		}
	}
	*/

	/************************************************************************/
	/* 随机抽样，计算投影系数，根据投影系数计算误差和中值                   */
	/************************************************************************/
	int** weight=new int*[RANSAC_CIRCUL];//内点的权重
	int ransac_cir=0;                    //总的循环次数
	double sample_num=(double)MAX_DATA;  //样本数
	int max_inline_num=0;                //最大内点集合中包含的点数，也是最重确定的内点个数
	int max_inline_index=0;              //最大内点集合的索引
	//最大循环抽样RANSAC_CIRCUL次
	do 
	{
		weight[ransac_cir]=new int[matnum];//权重
		memset(weight[ransac_cir],0,sizeof(int)*matnum);

		/******************************************************************************/
		/*随机选点                                                                    */
		/******************************************************************************/
		/*if (zero_num>4)
		{		
			randpoint1(re_pots,left_point,right_point);
		} 
		else*/
		{	
			randpoint2(pots,left_point,right_point);
		}
		/**********************************************************************************************/
		/* 由4对匹配点计算投影变换矩阵（即单映性矩阵）,利用该矩阵进行投影变换，计算每个点的投影误差   */
		/**********************************************************************************************/
		double** trans_h=cal_trans(left_point,right_point);	
		double med_dv;//投影变换误差的中值
		if (trans_h)
		{
			med_dv=cal_mid_error(pots.matpoint,trans_h,dv,matnum);
			//释放内存
			for (int i=0;i<3;i++)
			{
				delete []trans_h[i];
			}
			delete []trans_h;
		} 
		else
		{
			//计算旋转矩阵H失败
			cout<<"计算旋转矩阵H失败！"<<endl;
			ransac_cir++;
			continue;
		}
		/************************************************************************/
		/* 计算中值的阈值(公式来自王琳论文) 并判断内点                          */
		/************************************************************************/
		double sigma=1.4826*(1+5.0/(matnum-4))*sqrt(sqrt(med_dv));
		//根据阈值判断内点
		int inline_nums=0;//每次循环的内点个数
		for (int j=0;j<matnum;j++)
		{
			if (dv[j]<2.5*sigma*2.5*sigma)
			{
				weight[ransac_cir][j]=1;
				inline_nums++;
			}
		}
		/************************************************************************/
		/* 计算误匹配率，和样本数值                                             */
		/************************************************************************/	
		if (inline_nums>max_inline_num)
		{
			max_inline_index=ransac_cir;
			max_inline_num=inline_nums;
			double error_ratio=(double)inline_nums/matnum;
			sample_num=log(LEASTTHREAD)/(log(1.0-pow(error_ratio,4)));	
			//double error_ratio=1-(double)inline_nums/matnum;
			//sample_num=log(LEASTTHREAD)/(log(1.0-pow(1.0-pow(error_ratio,4),ransac_cir)));
		}
		ransac_cir++;
	} while (sample_num>=ransac_cir&&ransac_cir<500);

	//cout<<"循环"<<ransac_cir<<"次"<<endl;
	//cout<<"最大点数"<<max_inline_num<<endl;
	//释放内存
	delete []dv; 
	
	/*
	for (int i=0;i<RAND_BLOCK;i++)
	{
		for (int j=0;j<RAND_BLOCK;j++)
		{	
			delete []re_pots[i][j].matpoint;
		}
		delete []re_pots[i];
	}
	delete []re_pots;
	*/

	/************************************************************************/
	/* 保存剔除粗差后的点                                                   */
	/************************************************************************/
	if (max_inline_num!=0)
	{
		MatchPoints new_pot;
		new_pot.matpoint=savpoint(pots.matpoint,weight[max_inline_index],matnum,max_inline_num);
		inline_num=max_inline_num;
		//释放内存
		for (int i=0;i<ransac_cir;i++)
		{
			delete []weight[i];
		}
		delete []weight;
		return new_pot;
	} 
	else
	{
		inline_num=0;
		//释放内存
		for (int i=0;i<ransac_cir;i++)
		{
			delete []weight[i];
		}
		delete []weight;
		return pots;
	}
}

//重新组织匹配点
MatchPoints** reognizepoint(MatchPoints pots, int imgwidth, int imgheight)
{
	int** p_num=new int*[RAND_BLOCK];
	for (int i=0;i<RAND_BLOCK;i++)
	{
		p_num[i]=new int[RAND_BLOCK];
		memset(p_num[i],0,sizeof(int)*RAND_BLOCK);
	}
	int block_width=imgwidth/RAND_BLOCK;//图像分为LM_CIRCUL个区域后，每个区域列的宽
	int block_height=imgheight/RAND_BLOCK;//图像分为LM_CIRCUL个区域后，每个区域行的高
	//遍历所有点,记录每块内的点数和每个点所在的块的行列号	

	int* r=new int[pots.mat_key_num];
	int* c=new int[pots.mat_key_num];
	for (int i=0;i<pots.mat_key_num;i++)
	{
		r[i]=(pots.matpoint[i].imageX1/block_height)<RAND_BLOCK?(pots.matpoint[i].imageX1/block_height):(RAND_BLOCK-1);
		c[i]=(pots.matpoint[i].imageY1/block_width)<RAND_BLOCK?(pots.matpoint[i].imageY1/block_width):(RAND_BLOCK-1);
		p_num[r[i]][c[i]]++;
	}
	//为每个块开辟内存	
	MatchPoints** re_pots=new MatchPoints*[RAND_BLOCK];
	for (int i=0;i<RAND_BLOCK;i++)
	{		
		re_pots[i]=new MatchPoints[RAND_BLOCK];
		for (int j=0;j<RAND_BLOCK;j++)
		{
			re_pots[i][j].matpoint=new MatchPoint[p_num[i][j]]; 
			re_pots[i][j].mat_key_num=p_num[i][j];
		}
	}
	//为每一个块赋值
	for (int i=0;i<pots.mat_key_num;i++)
	{
		re_pots[r[i]][c[i]].matpoint[--p_num[r[i]][c[i]]]=pots.matpoint[i];
	}

	//释放内存
	delete []r;
	delete []c;
	for (int i=0;i<RAND_BLOCK;i++)
	{
		delete []p_num[i];
	}
	delete []p_num;
	return re_pots;
}


/****************************************************** 
* 函数名称： 
* randpoint()
*  
* 函数参数：
* MatchPoint* pot              - 两张图像的匹配点
* double left_point[4][2]      - 左图像上的点
* double right_point[4][2]     - 右图像上的点
* int block_width              - 图像分为RANSAC_BLOCK个区域后，每个区域列的宽
* int block_height             - 图像分为RANSAC_BLOCK个区域后，每个区域行的高
* int matnum                   - 匹配点的个数
* 
* 返回值：
* 1成功，0失败
*  
*说明：将图像分为RANSAC_BLOCK*RANSAC_BLOCK区域，随机选择4个区域，每个区域中随机选择一个点
********************************************************/
int randpoint(MatchPoint* pot,double left_point[4][2],double right_point[4][2],int block_width,int block_height,int matnum)
{
	int indnum=0;//记录随机选择点的循环次数
	int cirthred=0;//记录随机选择块的循环次数
	int ransac_cir=RAND_BLOCK*RAND_BLOCK*2;
	int inde[4];//随机选取的4个点在匹配点中的索引
	do 
	{
		int rect[4];
		rect[0]=RAND_BLOCK*RAND_BLOCK*(rand()/(RAND_MAX+1.0));
		do 
		{
			rect[1]=RAND_BLOCK*RAND_BLOCK*(rand()/(RAND_MAX+1.0));
		} while (rect[1]==rect[0]);	
		do 
		{
			rect[2]=RAND_BLOCK*RAND_BLOCK*(rand()/(RAND_MAX+1.0));
		} while ((rect[2]==rect[0])||(rect[2]==rect[1]));	
		do 
		{
			rect[3]=RAND_BLOCK*RAND_BLOCK*(rand()/(RAND_MAX+1.0));
		} while ((rect[3]==rect[0])||(rect[3]==rect[1])||(rect[3]==rect[2]));	

		indnum=0;
		int cirnum=0;//记录循环的次数
		do 
		{
			cirnum++;
			int row=rect[indnum]/RAND_BLOCK;
			int col=rect[indnum]%RAND_BLOCK;
			inde[indnum]=matnum*(rand()/(RAND_MAX+1.0));
			if (pot[inde[indnum]].imageY1>=col*block_width &&pot[inde[indnum]].imageY1<(col+1)*block_width &&   //imageY表示列，imageX表示行
				pot[inde[indnum]].imageX1>=(row)*block_height &&pot[inde[indnum]].imageX1<(row+1)*block_height)
			{
				left_point[indnum][0]=pot[inde[indnum]].imageY1;left_point[indnum][1]=pot[inde[indnum]].imageX1;
				right_point[indnum][0]=pot[inde[indnum]].imageY2;right_point[indnum][1]=pot[inde[indnum]].imageX2;
				indnum++;
			}
		} while ((indnum<4)&&(cirnum<matnum));//没找到四个匹配点，或者总次数已经大于总点数则重新选择块
		cirthred++;

	} while ( cirthred<ransac_cir && indnum<4 );

	if (cirthred>=ransac_cir)
	{	
		return 0;//失败
	}
	return 1;
}
//图像分块选点
void randpoint1(MatchPoints** pots,double left_point[4][2],double right_point[4][2])
{
	int indnum=0;
	int rect[4];
	do
	{	
		int eq=0;
		rect[indnum]=RAND_BLOCK*RAND_BLOCK*(rand()/(RAND_MAX+1.0));
		int row=rect[indnum]/RAND_BLOCK;
		int col=rect[indnum]%RAND_BLOCK;
		if (pots[row][col].mat_key_num>0)
		{
			int inde_num=pots[row][col].mat_key_num*(rand()/(RAND_MAX+1.0));
			left_point[indnum][0]=pots[row][col].matpoint[inde_num].imageY1;
			left_point[indnum][1]=pots[row][col].matpoint[inde_num].imageX1;
			right_point[indnum][0]=pots[row][col].matpoint[inde_num].imageY2;
			right_point[indnum][1]=pots[row][col].matpoint[inde_num].imageX2;
			//四个不同的块
			for (int i=0;i<indnum;i++)
			{
				eq+=(rect[i]==rect[indnum]);
			}
			if (eq==0)
			{	
				indnum++;
			}
		}
	}while ( indnum<4 );
}
//不分块选点
void randpoint2(MatchPoints pots,double left_point[4][2],double right_point[4][2])
{
	int indnum=0;
	//int rect[4];	
	while (indnum<4)
	{
		int inde_num=pots.mat_key_num*(rand()/(RAND_MAX+1.0));	
		left_point[indnum][0]=pots.matpoint[inde_num].imageY1;
		left_point[indnum][1]=pots.matpoint[inde_num].imageX1;
		right_point[indnum][0]=pots.matpoint[inde_num].imageY2;
		right_point[indnum][1]=pots.matpoint[inde_num].imageX2;
		indnum++;
	}
}
/****************************************************** 
* 函数名称： 
* cal_mid_error()
*  
* 函数参数：
* MatchPoints* pot      - 两张图像的匹配点
* double** trans_h      - 投影变换矩阵
* int num               - 匹配点的个数
*
* 返回值：
* 投影误差的中值
*  
*说明：利用投影变换矩阵，计算投影误差的中值
********************************************************/
double cal_mid_error(MatchPoint* pot,double** trans_h,double* dv,int num)
{
	double** temp_h=new double*[3];
	temp_h[0]=new double[3];
	memcpy(temp_h[0],trans_h[0],sizeof(double)*3);
	temp_h[1]=new double[3];
	memcpy(temp_h[1],trans_h[1],sizeof(double)*3);
	temp_h[2]=new double[3];
	memcpy(temp_h[2],trans_h[2],sizeof(double)*3);

	//利用旋转矩阵计算所有点投影后的坐标，并计算误差距离
	double** inv_trans_h=new double*[3];
	for (int j=0;j<3;j++)
	{
		inv_trans_h[j]=new double[3];
		memset(inv_trans_h[j],0,sizeof(double)*3);
	}

	//计算旋转矩阵H的逆矩阵
	int inv_result=Inv(temp_h,inv_trans_h,3);
	
	
	if (inv_result)
	{
		for (int j=0;j<num;j++)
		{
			//公式x1=(h0*x+h1*y+h2)/(h6*x+h7*y+1)
			//    y1=(h3*x+h4*y+h5)/(h6*x+h7*y+1)          
			double dx=pot[j].imageY1-(trans_h[0][0]*pot[j].imageY2+trans_h[0][1]*pot[j].imageX2+trans_h[0][2])/(trans_h[2][0]*pot[j].imageY2+trans_h[2][1]*pot[j].imageX2+trans_h[2][2]);
			double dy=pot[j].imageX1-(trans_h[1][0]*pot[j].imageY2+trans_h[1][1]*pot[j].imageX2+trans_h[1][2])/(trans_h[2][0]*pot[j].imageY2+trans_h[2][1]*pot[j].imageX2+trans_h[2][2]);
			double dxt=pot[j].imageY2-(inv_trans_h[0][0]*pot[j].imageY1+inv_trans_h[0][1]*pot[j].imageX1+inv_trans_h[0][2])/(inv_trans_h[2][0]*pot[j].imageY1+inv_trans_h[2][1]*pot[j].imageX1+inv_trans_h[2][2]);
			double dyt=pot[j].imageX2-(inv_trans_h[1][0]*pot[j].imageY1+inv_trans_h[1][1]*pot[j].imageX1+inv_trans_h[1][2])/(inv_trans_h[2][0]*pot[j].imageY1+inv_trans_h[2][1]*pot[j].imageX1+inv_trans_h[2][2]);
			dv[j]=dx*dx+dy*dy+dxt*dxt+dyt*dyt;//残差平方
		}
		//释放内存
		for (int i=0;i<3;i++)
		{
			delete []temp_h[i];
			delete []inv_trans_h[i];
		}
		delete []temp_h;
		delete []inv_trans_h;
		//查找中值
		return Find_Mid(dv,num);
	} 
	else
	{//旋转矩阵H求逆失败
		cout<<"旋转矩阵H求逆失败！"<<endl;
		//释放内存
		for (int i=0;i<3;i++)
		{
			delete []temp_h[i];
			delete []inv_trans_h[i];
		}
		delete []temp_h;
		delete []inv_trans_h;
		return (double)MAX_DATA;
	}
}
/****************************************************** 
* 函数名称： 
* inline_cal_trans()
*  
* 函数参数：
* MatchPoint* pot       - 两张图像的匹配点
* int num               - 匹配点的个数
*
* 返回值：
* 投影变换矩阵
*  
*说明：利用所有内点计算投影矩阵
********************************************************/
double** inline_cal_trans(MatchPoint* pot,int &num)
{
	////测试后100个数据的精度
	//num-=100;
	//系数矩阵XY
	double** coe_xy=new double*[2*num];
	for (int j=0;j<2*num;j++)
	{
		coe_xy[j]=new double[8];
	}
	for (int j=0;j<num;j++)
	{
		coe_xy[2*j][0]=pot[j].imageY2;
		coe_xy[2*j][1]=pot[j].imageX2;
		coe_xy[2*j][2]=1;
		coe_xy[2*j][3]=0;
		coe_xy[2*j][4]=0;
		coe_xy[2*j][5]=0;
		coe_xy[2*j][6]=-pot[j].imageY2*pot[j].imageY1;
		coe_xy[2*j][7]=-pot[j].imageX2*pot[j].imageY1;

		coe_xy[2*j+1][0]=0;
		coe_xy[2*j+1][1]=0;
		coe_xy[2*j+1][2]=0;
		coe_xy[2*j+1][3]=pot[j].imageY2;
		coe_xy[2*j+1][4]=pot[j].imageX2;
		coe_xy[2*j+1][5]=1;
		coe_xy[2*j+1][6]=-pot[j].imageY2*pot[j].imageX1;
		coe_xy[2*j+1][7]=-pot[j].imageX2*pot[j].imageX1;
	}

	//常数项xy
	double** xy=new double*[2*num];
	for (int j=0;j<num;j++)
	{	
		xy[2*j]=new double[1];
		xy[2*j+1]=new double[1];
		xy[2*j][0]=pot[j].imageY1;
		xy[2*j+1][0]=pot[j].imageX1;
	}

	//最小二乘
	//long st=GetTickCount();
	int num2=2*num;
	double** deta_h=Least_Square(coe_xy,xy,num2,8);

	num=num2/2;
	for (int j=0;j<num;j++)
	{
		pot[j].imageY2=coe_xy[2*j][0];
		pot[j].imageX2=coe_xy[2*j][1];
		pot[j].imageY1=coe_xy[2*j][6]/(-pot[j].imageY2);
		pot[j].imageX1=coe_xy[2*j+1][6]/(-pot[j].imageY2);
	}
	//cout<<deta_h[0][0]<<"  "<<deta_h[1][0]<<"  "<<deta_h[2][0]<<endl;
	//cout<<deta_h[3][0]<<"  "<<deta_h[4][0]<<"  "<<deta_h[5][0]<<endl;
	//cout<<deta_h[6][0]<<"  "<<deta_h[7][0]<<endl;
	//long et=GetTickCount();
	//cout<<"最小二乘时间"<<et-st<<endl;
	//释放内存
	for (int j=0;j<2*num;j++)
	{
		delete []xy[j];
		delete []coe_xy[j];
	}
	delete []xy;
	delete []coe_xy;

	return deta_h;
}


/****************************************************** 
* 函数名称： 
* cal_trans()
*  
* 函数参数：
* double left_point[4][2]    - 左影像4个点
* double right_point[4][2]   - 右影像4个点（左影像4个点的匹配点）
*
* 返回值：
*  投影变换矩阵
*  
*说明：计算投影变换矩阵
********************************************************/
double** cal_trans(double left_point[4][2],double right_point[4][2])
{
	double** trans_coe=new double*[3];
	for (int i=0;i<3;i++)
	{
		trans_coe[i]=new double[3];
	}
	double** trans_xy=new double*[8];
	double** inv_trans_xy=new double*[8];

	for (int i=0;i<8;i++)
	{
		trans_xy[i]=new double[8];
		inv_trans_xy[i]=new double[8];
		memset(inv_trans_xy[i],0,sizeof(double)*8);
	}
	//系数
	for (int i=0;i<4;i++)
	{
		trans_xy[2*i][0]=right_point[i][0];
		trans_xy[2*i][1]=right_point[i][1];
		trans_xy[2*i][2]=1;
		trans_xy[2*i][3]=0;
		trans_xy[2*i][4]=0;
		trans_xy[2*i][5]=0;
		trans_xy[2*i][6]=-right_point[i][0]*left_point[i][0];
		trans_xy[2*i][7]=-right_point[i][1]*left_point[i][0];

		trans_xy[2*i+1][0]=0;
		trans_xy[2*i+1][1]=0;
		trans_xy[2*i+1][2]=0;
		trans_xy[2*i+1][3]=right_point[i][0];
		trans_xy[2*i+1][4]=right_point[i][1];
		trans_xy[2*i+1][5]=1;
		trans_xy[2*i+1][6]=-right_point[i][0]*left_point[i][1];
		trans_xy[2*i+1][7]=-right_point[i][1]*left_point[i][1];

	}

	//常数项
	double* trans_x=new double[8];
	trans_x[0] =  left_point[0][0];
	trans_x[1] =  left_point[0][1];
	trans_x[2] =  left_point[1][0];
	trans_x[3] =  left_point[1][1];
	trans_x[4] =  left_point[2][0];
	trans_x[5] =  left_point[2][1];
	trans_x[6] =  left_point[3][0];
	trans_x[7] =  left_point[3][1];

	//系数求逆
	if (Inv(trans_xy,inv_trans_xy,8))
	{
		for (int i=0;i<8;i++)
		{
			double tempsum=0.0;
			for (int j=0;j<8;j++)
			{
				tempsum+=inv_trans_xy[i][j]*trans_x[j];
			}
			trans_coe[i/3][i%3]=tempsum;
		}
		trans_coe[2][2]=1;
		//释放内存
		for (int i=0;i<8;i++)
		{
			delete []trans_xy[i];
			delete []inv_trans_xy[i];
		}
		delete []inv_trans_xy;
		delete []trans_xy;
		delete []trans_x;

		return trans_coe;
	} 
	else
	{
		return NULL;
	}
}

