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
/* ��������                                                             */
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


//���̲ο����������ġ�һ�ֻ���������ȫ�Զ�ͼ��ƴ���㷨��
/****************************************************** 
* �������ƣ� 
* de_gross_error()
*  
* ����������
* MatchPoints matchpoints       - �޳��ֲ�ǰ�õ�
* MatchPoints &new_matchpoints  - �޳��ֲ��ĵ�
* int imgwidth                  - ƥ������Ӱ������Ӱ��Ŀ�(��Ϊ�ǽ���Ӱ���ΪRAND_BLOCK*RAND_BLOCK����������������)
* int imgheight                 - ƥ������Ӱ������Ӱ��ĸ�
* int &inline_num               - �޳��ֲ���ĸ���
*
* ����ֵ��
* ͼ����ͶӰ�任ϵ��
*  
*˵��������RANSAC��LMedS�㷨�޳��ֲģ��Ϊͼ���ͶӰ�任ϵ��
********************************************************/
MatchPoints del_gross_error(MatchPoints matchpoints,int imgwidth,int imgheight) 
{
	MatchPoints new_matchpoints;
	int inline_num=0;
	
	/************************************************************************/
	/* RANSAC�㷨�޳��ֲ�                                                   */
	/************************************************************************/
	new_matchpoints=sift_ransac(matchpoints,imgwidth,imgheight,matchpoints.mat_key_num,inline_num);
	new_matchpoints.mat_key_num=inline_num;
	
	/************************************************************************/
	/* LMedS��С��ֵƽ�����޳��ֲ�                                          */
	/************************************************************************/
	//inline_num=0;
	//new_matchpoints=LMedS(new_matchpoints,imgwidth,imgheight,new_matchpoints.mat_key_num,inline_num);
	//new_matchpoints.mat_key_num=inline_num;
	
	return new_matchpoints;
}
/****************************************************** 
* �������ƣ� 
* savpoint()
*  
* ����������
* MatchPoint* oldpoint         - ����ͼ���ƥ���
* int* weight                  - ���Ȩ�أ�1Ϊ�ڵ㣬����Ϊ���
* int matnum                   - �޳��ֲ�ǰƥ���ĸ���
* int inline_num               - �ڵ�ĸ���
* 
* ����ֵ��
* �޳��ֲ��ƥ���
*  
*˵���������޳��ֲ���ƥ���
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
* �������ƣ� 
* LMedS()
*  
* ����������
* MatchPoints* pot             - ����ͼ���ƥ���
* ImageFile imagefile          - ͼ����Ϣ
* int imgwidth                 - ƥ������Ӱ������Ӱ��Ŀ�
* int imgheight                - ƥ������Ӱ������Ӱ��ĸ�
* int& inline_num              - �޳��ֲ����ڵ���
* 
* ����ֵ��
* LMedS�޳��ֲ���ƥ���
*  
*˵����LMedS(��С��ֵƽ����)�޳��ֲ�
********************************************************/
MatchPoints LMedS(MatchPoints pots,int imgwidth,int imgheight,int matnum,int& inline_num)
{	
	double** dv=new double*[LM_CIRCUL];   //ͶӰ�任���
	for (int k=0;k<LM_CIRCUL;k++)
	{
		dv[k]=new double[matnum];
		memset(dv[k],0,sizeof(double)*matnum);
	}
	int block_width=imgwidth/RAND_BLOCK;     //ͼ���ΪLM_CIRCUL�������ÿ�������еĿ�
	int block_height=imgheight/RAND_BLOCK;   //ͼ���ΪLM_CIRCUL�������ÿ�������еĸ�

	/************************************************************************/
	/* �������������ͶӰϵ��������ͶӰϵ������������ֵ                   */
	/************************************************************************/
	double* med_dv=new double[LM_CIRCUL];//ͶӰ�任������ֵ
	double left_point[4][2];//���ڼ��㵥ӳ�����4��ƥ�����Ӱ���ϵĵ�
	double right_point[4][2];//���ڼ��㵥ӳ�����4��ƥ�����Ӱ���ϵĵ�

	for (int k=0;k<LM_CIRCUL;k++)//ѭ��LM_CIRCUL�Σ���ȡ����
	{
		/******************************************************************************/
		/*��ͼ����ȷ�ΪBLOCK��BLOCK�и�����,�����ȡ4����������ÿ�������漴��ȡһ���ǵ�*/
		/******************************************************************************/
		if (!randpoint(pots.matpoint,left_point,right_point,block_width,block_height,matnum))
		{	
			med_dv[k]=MAX_DATA;
			continue;//���ĸ�ƥ���ʧ��,���˳��˴�ѭ��
		}
		/**********************************************************************************************/
		/* ��4��ƥ������ͶӰ�任���󣨼���ӳ�Ծ���,���øþ������ͶӰ�任������ÿ�����ͶӰ���   */
		/**********************************************************************************************/
		double** trans_h=cal_trans(left_point,right_point);
		if (trans_h)
		{
			med_dv[k]=cal_mid_error(pots.matpoint,trans_h,dv[k],matnum);//������ת����������е���ת�������в������������
			
			for (int i=0;i<3;i++)//�ͷ��ڴ�
			{
				delete []trans_h[i];
			}
			delete []trans_h;
		} 
		else
		{
			med_dv[k]=MAX_DATA;//������ת����Hʧ��
		}
	}//LMѭ��

	/************************************************************************/
	/* Ѱ��������ֵ����Сֵ���Լ����Ӧ����ת����                           */
	/************************************************************************/
	int minnum=0;
	double* temp_dev=new double[LM_CIRCUL];
	memcpy(temp_dev,med_dv,sizeof(double)*LM_CIRCUL);
	BubbleSort(temp_dev,LM_CIRCUL,minnum);
	//�ͷ��ڴ�
	delete []temp_dev;
	//cout<<"��С��ֵ��"<<med_dv[minnum]<<endl;
	if (med_dv[minnum]==MAX_DATA)
	{//���е����ֵ�����ʧ��,��RANSAC�޳��ֲ�ʧ�ܣ�����ԭʼ��ƥ���
		delete []med_dv;//�ͷ��ڴ�
		return pots; 
	}
	/************************************************************************/
	/* ������ֵ����ֵ(��ʽ������������) ���ж��ڵ�                          */
	/************************************************************************/
	double sigma=1.4826*(1+5.0/(matnum-4))*sqrt(/*sqrt(*/med_dv[minnum]/*)*/);//��matnum>=2*4ʱ���Բ��ô˹�ʽ����������������
	delete []med_dv;

	int* weight=new int[matnum];
	memset(weight,0,sizeof(int)*matnum);
	inline_num=0;//�ڵ�ĸ���

	//������ֵ�ж��ڵ�
	for (int j=0;j<matnum;j++)
	{
		if (dv[minnum][j]<2.5*sigma*2.5*sigma)
		{
			weight[j]=1;//���ڵ�
			inline_num++;
		}
	}
	//�ͷ��ڴ�
	for (int k=0;k<LM_CIRCUL;k++)
	{
		delete []dv[k];
	}
	delete []dv; 
	/************************************************************************/
	/* �����޳��ֲ��ĵ�                                                   */
	/************************************************************************/
	MatchPoints new_pot;
	new_pot.matpoint=savpoint(pots.matpoint,weight,matnum,inline_num);
	//�ͷ��ڴ�
	delete []weight;
	return new_pot;

}
/****************************************************** 
* �������ƣ� 
* sift_ransac()
*  
* ����������
* MatchPoints* pot             - ����ͼ���ƥ���
* int imgwidth                 - ƥ������Ӱ������Ӱ��Ŀ�
* int imgheight                - ƥ������Ӱ������Ӱ��ĸ�
* int matnum                   - ƥ���ĸ���
* int& inline_num              - �޳��ֲ����ڵ���
* 
* ����ֵ��
* RANSAC�޳��ֲ���ƥ���
*  
*˵����RANSAC�޳��ֲ�
********************************************************/
MatchPoints sift_ransac(MatchPoints pots,int imgwidth,int imgheight,int matnum,int& inline_num)
{
	double left_point[4][2];       //���ڼ��㵥ӳ�����4��ƥ�����Ӱ���ϵĵ�
	double right_point[4][2];      //���ڼ��㵥ӳ�����4��ƥ�����Ӱ���ϵĵ�
	double* dv=new double[matnum]; //ͶӰ�任���
	memset(dv,0,sizeof(double)*matnum);
	
	/************************************************************************/
	/* �Ը�������ʽ������֯ƥ���                                           */
	/************************************************************************/
	/*
	MatchPoints** re_pots=reognizepoint(pots,imgwidth,imgheight);
	//���еĵ�����Ϊ0�Ŀ���
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
	/* �������������ͶӰϵ��������ͶӰϵ������������ֵ                   */
	/************************************************************************/
	int** weight=new int*[RANSAC_CIRCUL];//�ڵ��Ȩ��
	int ransac_cir=0;                    //�ܵ�ѭ������
	double sample_num=(double)MAX_DATA;  //������
	int max_inline_num=0;                //����ڵ㼯���а����ĵ�����Ҳ������ȷ�����ڵ����
	int max_inline_index=0;              //����ڵ㼯�ϵ�����
	//���ѭ������RANSAC_CIRCUL��
	do 
	{
		weight[ransac_cir]=new int[matnum];//Ȩ��
		memset(weight[ransac_cir],0,sizeof(int)*matnum);

		/******************************************************************************/
		/*���ѡ��                                                                    */
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
		/* ��4��ƥ������ͶӰ�任���󣨼���ӳ�Ծ���,���øþ������ͶӰ�任������ÿ�����ͶӰ���   */
		/**********************************************************************************************/
		double** trans_h=cal_trans(left_point,right_point);	
		double med_dv;//ͶӰ�任������ֵ
		if (trans_h)
		{
			med_dv=cal_mid_error(pots.matpoint,trans_h,dv,matnum);
			//�ͷ��ڴ�
			for (int i=0;i<3;i++)
			{
				delete []trans_h[i];
			}
			delete []trans_h;
		} 
		else
		{
			//������ת����Hʧ��
			cout<<"������ת����Hʧ�ܣ�"<<endl;
			ransac_cir++;
			continue;
		}
		/************************************************************************/
		/* ������ֵ����ֵ(��ʽ������������) ���ж��ڵ�                          */
		/************************************************************************/
		double sigma=1.4826*(1+5.0/(matnum-4))*sqrt(sqrt(med_dv));
		//������ֵ�ж��ڵ�
		int inline_nums=0;//ÿ��ѭ�����ڵ����
		for (int j=0;j<matnum;j++)
		{
			if (dv[j]<2.5*sigma*2.5*sigma)
			{
				weight[ransac_cir][j]=1;
				inline_nums++;
			}
		}
		/************************************************************************/
		/* ������ƥ���ʣ���������ֵ                                             */
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

	//cout<<"ѭ��"<<ransac_cir<<"��"<<endl;
	//cout<<"������"<<max_inline_num<<endl;
	//�ͷ��ڴ�
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
	/* �����޳��ֲ��ĵ�                                                   */
	/************************************************************************/
	if (max_inline_num!=0)
	{
		MatchPoints new_pot;
		new_pot.matpoint=savpoint(pots.matpoint,weight[max_inline_index],matnum,max_inline_num);
		inline_num=max_inline_num;
		//�ͷ��ڴ�
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
		//�ͷ��ڴ�
		for (int i=0;i<ransac_cir;i++)
		{
			delete []weight[i];
		}
		delete []weight;
		return pots;
	}
}

//������֯ƥ���
MatchPoints** reognizepoint(MatchPoints pots, int imgwidth, int imgheight)
{
	int** p_num=new int*[RAND_BLOCK];
	for (int i=0;i<RAND_BLOCK;i++)
	{
		p_num[i]=new int[RAND_BLOCK];
		memset(p_num[i],0,sizeof(int)*RAND_BLOCK);
	}
	int block_width=imgwidth/RAND_BLOCK;//ͼ���ΪLM_CIRCUL�������ÿ�������еĿ�
	int block_height=imgheight/RAND_BLOCK;//ͼ���ΪLM_CIRCUL�������ÿ�������еĸ�
	//�������е�,��¼ÿ���ڵĵ�����ÿ�������ڵĿ�����к�	

	int* r=new int[pots.mat_key_num];
	int* c=new int[pots.mat_key_num];
	for (int i=0;i<pots.mat_key_num;i++)
	{
		r[i]=(pots.matpoint[i].imageX1/block_height)<RAND_BLOCK?(pots.matpoint[i].imageX1/block_height):(RAND_BLOCK-1);
		c[i]=(pots.matpoint[i].imageY1/block_width)<RAND_BLOCK?(pots.matpoint[i].imageY1/block_width):(RAND_BLOCK-1);
		p_num[r[i]][c[i]]++;
	}
	//Ϊÿ���鿪���ڴ�	
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
	//Ϊÿһ���鸳ֵ
	for (int i=0;i<pots.mat_key_num;i++)
	{
		re_pots[r[i]][c[i]].matpoint[--p_num[r[i]][c[i]]]=pots.matpoint[i];
	}

	//�ͷ��ڴ�
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
* �������ƣ� 
* randpoint()
*  
* ����������
* MatchPoint* pot              - ����ͼ���ƥ���
* double left_point[4][2]      - ��ͼ���ϵĵ�
* double right_point[4][2]     - ��ͼ���ϵĵ�
* int block_width              - ͼ���ΪRANSAC_BLOCK�������ÿ�������еĿ�
* int block_height             - ͼ���ΪRANSAC_BLOCK�������ÿ�������еĸ�
* int matnum                   - ƥ���ĸ���
* 
* ����ֵ��
* 1�ɹ���0ʧ��
*  
*˵������ͼ���ΪRANSAC_BLOCK*RANSAC_BLOCK�������ѡ��4������ÿ�����������ѡ��һ����
********************************************************/
int randpoint(MatchPoint* pot,double left_point[4][2],double right_point[4][2],int block_width,int block_height,int matnum)
{
	int indnum=0;//��¼���ѡ����ѭ������
	int cirthred=0;//��¼���ѡ����ѭ������
	int ransac_cir=RAND_BLOCK*RAND_BLOCK*2;
	int inde[4];//���ѡȡ��4������ƥ����е�����
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
		int cirnum=0;//��¼ѭ���Ĵ���
		do 
		{
			cirnum++;
			int row=rect[indnum]/RAND_BLOCK;
			int col=rect[indnum]%RAND_BLOCK;
			inde[indnum]=matnum*(rand()/(RAND_MAX+1.0));
			if (pot[inde[indnum]].imageY1>=col*block_width &&pot[inde[indnum]].imageY1<(col+1)*block_width &&   //imageY��ʾ�У�imageX��ʾ��
				pot[inde[indnum]].imageX1>=(row)*block_height &&pot[inde[indnum]].imageX1<(row+1)*block_height)
			{
				left_point[indnum][0]=pot[inde[indnum]].imageY1;left_point[indnum][1]=pot[inde[indnum]].imageX1;
				right_point[indnum][0]=pot[inde[indnum]].imageY2;right_point[indnum][1]=pot[inde[indnum]].imageX2;
				indnum++;
			}
		} while ((indnum<4)&&(cirnum<matnum));//û�ҵ��ĸ�ƥ��㣬�����ܴ����Ѿ������ܵ���������ѡ���
		cirthred++;

	} while ( cirthred<ransac_cir && indnum<4 );

	if (cirthred>=ransac_cir)
	{	
		return 0;//ʧ��
	}
	return 1;
}
//ͼ��ֿ�ѡ��
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
			//�ĸ���ͬ�Ŀ�
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
//���ֿ�ѡ��
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
* �������ƣ� 
* cal_mid_error()
*  
* ����������
* MatchPoints* pot      - ����ͼ���ƥ���
* double** trans_h      - ͶӰ�任����
* int num               - ƥ���ĸ���
*
* ����ֵ��
* ͶӰ������ֵ
*  
*˵��������ͶӰ�任���󣬼���ͶӰ������ֵ
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

	//������ת����������е�ͶӰ������꣬������������
	double** inv_trans_h=new double*[3];
	for (int j=0;j<3;j++)
	{
		inv_trans_h[j]=new double[3];
		memset(inv_trans_h[j],0,sizeof(double)*3);
	}

	//������ת����H�������
	int inv_result=Inv(temp_h,inv_trans_h,3);
	
	
	if (inv_result)
	{
		for (int j=0;j<num;j++)
		{
			//��ʽx1=(h0*x+h1*y+h2)/(h6*x+h7*y+1)
			//    y1=(h3*x+h4*y+h5)/(h6*x+h7*y+1)          
			double dx=pot[j].imageY1-(trans_h[0][0]*pot[j].imageY2+trans_h[0][1]*pot[j].imageX2+trans_h[0][2])/(trans_h[2][0]*pot[j].imageY2+trans_h[2][1]*pot[j].imageX2+trans_h[2][2]);
			double dy=pot[j].imageX1-(trans_h[1][0]*pot[j].imageY2+trans_h[1][1]*pot[j].imageX2+trans_h[1][2])/(trans_h[2][0]*pot[j].imageY2+trans_h[2][1]*pot[j].imageX2+trans_h[2][2]);
			double dxt=pot[j].imageY2-(inv_trans_h[0][0]*pot[j].imageY1+inv_trans_h[0][1]*pot[j].imageX1+inv_trans_h[0][2])/(inv_trans_h[2][0]*pot[j].imageY1+inv_trans_h[2][1]*pot[j].imageX1+inv_trans_h[2][2]);
			double dyt=pot[j].imageX2-(inv_trans_h[1][0]*pot[j].imageY1+inv_trans_h[1][1]*pot[j].imageX1+inv_trans_h[1][2])/(inv_trans_h[2][0]*pot[j].imageY1+inv_trans_h[2][1]*pot[j].imageX1+inv_trans_h[2][2]);
			dv[j]=dx*dx+dy*dy+dxt*dxt+dyt*dyt;//�в�ƽ��
		}
		//�ͷ��ڴ�
		for (int i=0;i<3;i++)
		{
			delete []temp_h[i];
			delete []inv_trans_h[i];
		}
		delete []temp_h;
		delete []inv_trans_h;
		//������ֵ
		return Find_Mid(dv,num);
	} 
	else
	{//��ת����H����ʧ��
		cout<<"��ת����H����ʧ�ܣ�"<<endl;
		//�ͷ��ڴ�
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
* �������ƣ� 
* inline_cal_trans()
*  
* ����������
* MatchPoint* pot       - ����ͼ���ƥ���
* int num               - ƥ���ĸ���
*
* ����ֵ��
* ͶӰ�任����
*  
*˵�������������ڵ����ͶӰ����
********************************************************/
double** inline_cal_trans(MatchPoint* pot,int &num)
{
	////���Ժ�100�����ݵľ���
	//num-=100;
	//ϵ������XY
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

	//������xy
	double** xy=new double*[2*num];
	for (int j=0;j<num;j++)
	{	
		xy[2*j]=new double[1];
		xy[2*j+1]=new double[1];
		xy[2*j][0]=pot[j].imageY1;
		xy[2*j+1][0]=pot[j].imageX1;
	}

	//��С����
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
	//cout<<"��С����ʱ��"<<et-st<<endl;
	//�ͷ��ڴ�
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
* �������ƣ� 
* cal_trans()
*  
* ����������
* double left_point[4][2]    - ��Ӱ��4����
* double right_point[4][2]   - ��Ӱ��4���㣨��Ӱ��4�����ƥ��㣩
*
* ����ֵ��
*  ͶӰ�任����
*  
*˵��������ͶӰ�任����
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
	//ϵ��
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

	//������
	double* trans_x=new double[8];
	trans_x[0] =  left_point[0][0];
	trans_x[1] =  left_point[0][1];
	trans_x[2] =  left_point[1][0];
	trans_x[3] =  left_point[1][1];
	trans_x[4] =  left_point[2][0];
	trans_x[5] =  left_point[2][1];
	trans_x[6] =  left_point[3][0];
	trans_x[7] =  left_point[3][1];

	//ϵ������
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
		//�ͷ��ڴ�
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

