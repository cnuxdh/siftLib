/************************************************************************/
/*提取SIFT特征                                  */
/************************************************************************/

/**************************************头文件*****************************************/
//#include "stdafx.h"
#include <iostream>
#include <math.h>

#ifdef WIN32
#include "omp.h"
#endif


#include "siftfeature.h"

using namespace std;




/**************************************声明图像处理函数*****************************************/
Key_Point* _sift_features( float *, int , int , int & );//sift特征提取
float* GetGridData(float* ,int ,int ,int ,int ,int ,int );//从数组中读取子块的数据
float* Down_Sample(float *,int &,int &);//降采样
SCALESPACE**  Build_Guss_Space(float *,int,int,int);//建立高斯空间
float* Guss_Smooth(float *,int,int,float);//高斯平滑
void Build_Dog_Space(SCALESPACE**,int,int);//建立差分高斯空间
Key_Point* DetectKeypoint(SCALESPACE**,int,int,float,float,int &);//检测极值点
bool is_extremum(SCALESPACE**,int,int,int,int,int);//判断是否为极值
bool interp_step(SCALESPACE **,int ,int ,int ,int,int ,float &,float &,float &);//三维二次拟合，求取极值点的改正值
bool interp_extremum(SCALESPACE **,int,int&,int&,int&,int,int,float,float &,float &,float &);//插入极值点，精确定位关键点
float interp_contr(SCALESPACE **,int,int,int,int,int,float,float,float);//计算精确定位的极值点对应的灰度值
bool is_too_edge_like(SCALESPACE **,int,int,int,int,int,float);
float inter_contr(int , int , int , int ,float ,float ,float);//计算精确定位后关键点的对比度
void calc_scales(Key_Point *,int);//计算关键点所在的尺度
Key_Point* calc_oritation(Key_Point* ,SCALESPACE** ,int &);//计算关键点的主方向
bool calc_grad_mag_ori(float*, int, int, int, int, float &g, float &);//计算点的方向和梯度
float dominant_ori( float*, int );//计算关键点的主方向
float* ori_hist(float* , int , int , int ,int , int , int , float);//计算方向直方图
void smooth_ori_hist( float* , int );//对直方图进行平滑
float* descr_hist(float* ,int ,int ,int ,int ,float ,float ,int ,int ,int);//计算计算特征点的描述直方图
void normalize_descr(float* ,int );//对128维描述子进行归一化
void interp_hist_entry(float* ,float ,float ,float ,float ,int ,int );//计算每个特征点邻域的点对描述子直方图的贡献
float interp_hist_peak(float ,float ,float/* ,float &*/);//对直方图进行抛物线差值
void compute_descriptors( Key_Point* , SCALESPACE **, int );//计算特征点描述子
bool adjust_for_img_down(Key_Point* ,int );//对原始影像进行降采样时，需要对关键点的坐标进行调整
int sift_Inv(double** a,double** inva,int dim );//对矩阵求逆
void adjust_for_img_grid(Key_Point* feat,int keynumber,int col,int row,int gridwidth,int imgheight);//对原始影像进行分块时，需要对关键点的坐标进行调整
//void free_keypoint(Key_Point* keyhead);
/**************************************图像处理函数实体*****************************************/

/****************************************************** 
* 函数名称： 
* sift_features()
*  
* 函数参数：
* unsigned char* imagedata- GDAL读取的图像数据
* int colum      - 图像列数
* int row        - 图像行数
* int bands      - 图像波段数
* int &keynumber - 关键点的个数
* int &gridcol   - 列方向块的个数
* int &gridrow   - 行方向块的个数
*
* 返回值：
* 提取的特征点
*  
*说明：提取图像的sift特征点，判断是否对图像进行分块
********************************************************/
Feature** sift_features( float* imgdata,int colum,int row,int &keynumber,int &gridcol,int &gridrow )
{
	/*********************************图像数据分块***************************************/
	int gridwidth = SIFT_IMG_GRID_THRESHOLD; //块的宽
	
	//确定块的高和宽方向的块个数
	gridcol=(colum+gridwidth-1)/gridwidth;       //列向块的个数
	gridrow=(row+gridwidth-1)/gridwidth;         //行向块的个数
	int *gridc=new int [gridcol]; 
	int *gridr=new int [gridrow];       
	for (int i=0;i<gridcol-1;i++)
	{
		gridc[i]=gridwidth;
	}
	for (int i=0;i<gridrow-1;i++)
	{
		gridr[i]=gridwidth;
	}
	gridr[gridrow-1]=((row-1) % gridwidth)+1; //行列末尾小块 
	gridc[gridcol-1]=((colum-1) % gridwidth)+1;

	float*** griddata=new float**[gridrow];
	
#ifdef WIN32
	int pros_num=omp_get_num_procs();
	omp_set_num_threads(pros_num*3/4);
#pragma omp parallel for
#endif
	for (int i=0;i<gridrow;i++)
	{	
		griddata[i]=new float* [gridcol];
		for (int j=0;j<gridcol;j++)
		{
			griddata[i][j]=GetGridData(imgdata,i,j,colum,gridwidth,gridr[i],gridc[j]);
		}
	}
	delete []imgdata;

	//分块提取特征点
	Feature **feat=new Feature *[gridrow];
	#ifdef WIN32
	omp_set_num_threads(pros_num*3/4);
#pragma omp parallel for
#endif

	for (int i=0;i<gridrow;i++)
	{
		feat[i]=new Feature [gridcol];
		for (int j=0;j<gridcol;j++)
		{
			int tmpnum=0;
			feat[i][j].feature=_sift_features( griddata[i][j], gridc[j], gridr[i], tmpnum );
			keynumber+=tmpnum;
			feat[i][j].row = i;
			feat[i][j].col = j;
			feat[i][j].num = tmpnum;
		}
	}

	//将分块中的坐标，转换为原始影像上的坐标
	for (int k=0;k<gridrow;k++)
	{
		for (int m=0;m<gridcol;m++)
		{
			adjust_for_img_grid(feat[k][m].feature,feat[k][m].num,
				feat[k][m].col,feat[k][m].row,SIFT_IMG_GRID_THRESHOLD,row);
		}
	}

	delete []gridr;
	delete []gridc;
	return feat;
}

/****************************************************** 
* 函数名称： 
* GetGridData()
*  
* 函数参数：
* float *imgdata  - 图像数据
* int rownum      - 块的行号
* int colnum      - 块的列号
* int imgwidth    - 原始图像的宽
* int gridwidth   - 标准块的宽
* int row         - 块的行数
* int col         - 块的列数
* 
* 返回值：
* 块的数据
*  
*说明：提取分块的数据
********************************************************/
float* GetGridData(float* imgdata,int rownum,int colnum,int imgwidth,int gridwidth,int row,int col)
{
	float* griddata=new float[row*col];
	int tmprow=rownum*gridwidth;
	int tmpcolum=colnum*gridwidth;
	for (int i=0;i<row;i++)
	{
		for (int j=0;j<col;j++)
		{
			griddata[i*col+j]=imgdata[(tmprow+i)*imgwidth+tmpcolum+j];
		}
	}
	return griddata;
}

/****************************************************** 
* 函数名称： 
* _sift_features()
*  
* 函数参数：
* float *imgdata  - 图像数据
* int colum       - 数据的列数
* int row         - 数据的行数
* int &keynumber  - 关键点的个数
* 
* 返回值：
* 提取的特征点
*  
*说明：提取图像的sift特征点
********************************************************/
Key_Point* _sift_features( float *imgdata, int colum, int row, int &keynumber )
{	
	//如果行和列中最小值小于阈值则此块不计算特征点
	if (min(colum,row)<GRID_REM_PIXEL_THRESHOLD)
	{
		keynumber=0;
		return NULL;
	}
	//首先进行降采样，减少数据量
	float *idata;
	if (SIFT_IMG_DOWN)
	{
		idata=Down_Sample(imgdata,colum,row);
		//delete []imgdata;
	}
	else
	{
		idata=imgdata;
	}

	/************************************************************************/
	/*1 建立高斯空间                                                        */
	/************************************************************************/
	//计算尺度空间最大组数，防止出现采样后图像为0的情况

	//尺度空间的组数
	int octnum;
	if (OCTAVESTYE)
	{
		octnum=1;
	} 
	else
	{
		octnum=int(log(double(min(colum, row)))/log(2.0)-6);
		//如果行和列的最小值为64则将组数设为1
		if (min(colum, row)<64)
		{
			octnum=1;
		}
		octnum=max(1,octnum);
	}
	
	SCALESPACE** scalespace;
	scalespace=Build_Guss_Space(idata,octnum,colum,row);
	//delete []idata;

	/************************************************************************/
	/*2 建立差分高斯空间                                                    */
	/************************************************************************/

	Build_Dog_Space(scalespace,octnum,SCALESPEROCTAVE);

	/************************************************************************/
	/*3 检测极值点                                                          */
	/************************************************************************/

	Key_Point* keyhead=new Key_Point[1];
	keyhead->next=DetectKeypoint(scalespace,octnum,SCALESPEROCTAVE,
		SIFT_CONTR_THRESHOLD,SIFT_CURV_THRESHOLD,keynumber);

	/************************************************************************/
	/*4 计算特征点的尺度                                                    */
	/************************************************************************/
	calc_scales(keyhead,keynumber);

	/************************************************************************/
	/*5 确定主方向                                                          */
	/************************************************************************/

	Key_Point* feathead=new Key_Point[1];
	feathead->next=calc_oritation(keyhead,scalespace,keynumber);

	//将以链表形式存储的特征点改为以指针形式存储（占用时间较多）
	Key_Point *feat=new Key_Point[keynumber];
	Key_Point *tempfeat;
	tempfeat=feathead->next;
	for (int i=0;i<keynumber;i++)
	{
		feat[i]=*tempfeat;
		feat[i].index = i;
		tempfeat=tempfeat->next;
	}
	delete []tempfeat;

	Key_Point* pnode=feathead;
	Key_Point* tempkey;
	while(pnode!=NULL)
	{		
		tempkey=pnode;
		pnode=pnode->next;
		delete []tempkey;
	}
	feathead=NULL;

	/************************************************************************/
	/*6 计算特征点的描述子                                                    */
	/************************************************************************/

	compute_descriptors(feat,scalespace,keynumber);

	//如果对原始图像进行降采样，则需要对特征点的位置进行调整
	if (SIFT_IMG_DOWN)
	{
		adjust_for_img_down(feat,keynumber);
	}

	//删除尺度空间数据
	for (int i=0;i<octnum;i++)
	{
		for (int j=0;j<SCALESPEROCTAVE+2;j++)
		{
			delete []scalespace[i][j].dog_space;
			delete []scalespace[i][j].gauss_space;
		}
		delete []scalespace[i][SCALESPEROCTAVE+2].gauss_space;
		delete []scalespace[i];
	}
	delete []scalespace;
	return feat; 
}


/****************************************************** 
* 函数名称： 
* Down_sample()
*  
* 函数参数：
* float *img-指针指向进行降采样的数据
* int &col-整型数据的引用，为图像的列数即宽度
* int &row-整型数据的引用，为图像的行数即高度
*
* 返回值：
* 缩小后的图像
*  
*说明：对图像进行隔点采样
********************************************************/
float* Down_Sample(float *img,int &col,int &row)
{
	int row1=row/2;//行缩小为原来的二分之一
	int col1=col/2;//列缩小为原来的二分之一
	
	float *temp=new float[row1*col1];//用于临时存储降采样图像数据。
	for (int i=0;i<row1;i++)
	{
		for (int j=0;j<col1;j++)
		{
			temp[i*col1+j]=img[i*2*col+2*j];//注意第二个是col不是col1,表示降采样前每行的个数，它不一定等于2*col1，可能是2*col+1
			//temp[i*col1+j]=img[i*2*col+2*j]+img[i*2*col+2*j+1]+img[(i*2+1)*col+2*j]+img[(i*2+1)*col+2*j+1];
		}
	}
	row=row1;col=col1;
	return temp;
}
/****************************************************** 
* 函数名称： 
* Build_Scale_Space()
*  
* 函数参数：
* float *pimage-容器引用指向进行降采样的图像数据
* int oct_num-高斯总组数
* int col-为图像的列数即宽度
* int row-为图像的行数即高度
*
*返回值：
* SCALESPACE** 高斯空间
*  
*说明：建立高斯空间
********************************************************/
SCALESPACE** Build_Guss_Space(float* pimage,int oct_num,int col,int row)
{
	/***********************************图像和高斯函数做卷积生成高斯空间 ***************************************/
	SCALESPACE** gs=new SCALESPACE* [oct_num];//存储高斯空间
	for (int i=0;i<oct_num;i++)
	{
		gs[i]=new SCALESPACE [SCALESPEROCTAVE+3];
	}
	
	float k_sig= float(sqrt(SIFT_SIGMA * SIFT_SIGMA - SIFT_INIT_SIGMA *SIFT_INIT_SIGMA ));//计算最底层平滑的sigma值
	float intvl_k=pow(2.0f,1.0f/SCALESPEROCTAVE);//每层之间的尺度比例系数
	//当SCALESPEROCTAVE等于1时，为避免k值太大从而造成平滑模板太大的情况，如下处理
	if (SCALESPEROCTAVE==1)
	{
		intvl_k=sqrt(3.0f);
	}
	float *intvl_sigma=new float[SCALESPEROCTAVE+3];//存储每层高斯核尺度
	float pre_sigma,total_sigma;

	/**********************************开始计算采用高斯函数的核***********************************/
	intvl_sigma[0]=SIFT_SIGMA;
	for(int j=1;j<SCALESPEROCTAVE+3;j++)
	{
		pre_sigma=pow(intvl_k,j-1)*SIFT_SIGMA;
		total_sigma=pre_sigma*intvl_k;//按照下一层的影像由上一层影像计算得出
		intvl_sigma[j]=sqrt(total_sigma*total_sigma-pre_sigma*pre_sigma);//每一层采用的高斯平滑函数的核
	}

	/**********************************高斯卷积生成高斯空间***********************************/
#ifdef WIN32
	int pros_num=omp_get_num_procs();
	omp_set_num_threads(pros_num*3/4);
#pragma omp parallel for
#endif
	for (int i=0;i<oct_num;i++)
	{
		for (int j=0;j<SCALESPEROCTAVE+3;j++)//每组有SCALESPEROCTAVE+3层
		{	
			if (i==0 && j==0)
			{	//对原始影像进行平滑，产生第一组第一层的高斯影像
				gs[0][0].gauss_space=Guss_Smooth(pimage,col,row,k_sig);
				gs[0][0].gus_column=col;//高斯图像列数
				gs[0][0].gus_row=row;//高斯图像行数
			}
			else if (j==0)
			{//从第二组开始,对于上一组中二倍于底层尺度的高斯层进行降采样
				gs[i][0].gauss_space=Down_Sample(gs[i-1][SCALESPEROCTAVE].gauss_space,col,row);	
				gs[i][0].gus_column=col;//高斯图像列数
				gs[i][0].gus_row=row;//高斯图像行数
			}
			else
			{//从第二层开始平滑
				gs[i][j].gauss_space=Guss_Smooth(gs[i][j-1].gauss_space,col,row,intvl_sigma[j]);//下一层影像是在上一层影像平滑的基础上进行平滑的
				gs[i][j].gus_column=col;//高斯图像列数
				gs[i][j].gus_row=row;//高斯图像行数
			}
		}
	}
	delete []intvl_sigma;

	return gs;

}
/****************************************************** 
* 函数名称： 
*  Gauss_Smooth()
*  
* 函数参数：
* float *pre_data-存储待平滑数据的容器
* int g_col-为待平滑数据列的数即宽度
* int g_row-为待平滑数据的行数即高度
* int g_sigma-高斯平滑的核
*
* 返回值：
* float *，存储平滑后的高斯图像
*  
* 说明：高斯平滑
********************************************************/
float* Guss_Smooth(float *pre_data,int g_col,int g_row,float g_sigma)
{
	/*********************************** 计算高斯模板的大小 ***************************************/
	int dim =1+2*(int)(2.0f*g_sigma);
	int c=(dim+1)/2;  //模板的半径加1
	int r=dim/2;       //模板的半径
	float s2=g_sigma*g_sigma;
	float *g_mat=new float[c];//用于存储卷积模板
	float temp_sum=0;
	//float g_coe=1.0f/sqrt(2.0f*PI*g_sigma);
	g_mat[0]=1.0f/*g_coe*/;//只采用高斯核，没有高斯函数前面的系数，因此在中间处的值为exp(-(1.0*0*0)/(2.0*s2))=1
	temp_sum+=g_mat[0];//首先加上中间数值

	for(int i=1;i<c;i++) 
	{
		g_mat[i]=/*g_coe**/exp(-(1.0f*i*i)/(2.0f*s2));//高斯核;
		temp_sum+=2*g_mat[i];
	}
	//模板归一化
	for (int i=0;i<c;i++)
	    g_mat[i]=g_mat[i]/temp_sum;
	float *temp_row=new float[g_row*g_col];//存行平滑后的数据
	float *nex_data=new float[g_row*g_col];//存储高斯平滑后的数据

	//将二维高斯卷积分解为水平和竖直方向的两个一维卷积，先进行一维行卷积，利用其结果再进行列卷积
    /*********************************** 行方向高斯平滑***************************************/
	#ifdef WIN32
	int pros_num=omp_get_num_procs();
	omp_set_num_threads(pros_num*3/4);
#pragma omp parallel for
#endif

	for (int i=0;i<g_row;i++)
	{
		//每一行的前dim/2个像素没有平滑
		for (int m=0;m<r;m++)	
			temp_row[i*g_col+m]=pre_data[i*g_col+m];

		//对于中间数据，直接处理，不用考虑越界问题
		for (int j=r;j<g_col-r;j++)
		{
			float sum=0;//存储模板平滑数据
			sum+=pre_data[i*g_col+j]*g_mat[0];//模板中间数据单独计算
			for (int k=1;k<c;k++)           //宽度为7的模板（-3,-2,-1,0,1,2,3）
			{
				sum+=pre_data[i*g_col+j-k]*g_mat[k];
				sum+=pre_data[i*g_col+j+k]*g_mat[k];
			}
			temp_row[i*g_col+j]=sum;//存储行平滑后的数据
		}
		//第1种方法：每一行的后dim/2个像素没有平滑
		for (int m=g_col-r;m<g_col;m++)
			temp_row[i*g_col+m]=pre_data[i*g_col+m];
	}
	/***********************************列方向高斯平滑 ***************************************/

	//前dim/2行没有平滑（可以改进只用一个循环！！！！）
	for (int m=0;m<r;m++)
		for (int j=0;j<g_col;j++)
			nex_data[m*g_col+j]=temp_row[m*g_col+j];//直接赋值
	
	//对于中间数据，直接处理，不用考虑越界问题
	#ifdef WIN32
	omp_set_num_threads(pros_num*3/4);
#pragma omp parallel for
#endif

	for (int i=r;i<g_row-r;i++)
	{
		for (int j=0;j<g_col;j++)
		{
			float sum=0;//存储模板平滑数据
			sum+=temp_row[i*g_col+j]*g_mat[0];//模板中间数据单独计算
			for (int k=1;k<c;k++)
			{
				sum+=temp_row[(i-k)*g_col+j]*g_mat[k];
				sum+=temp_row[(i+k)*g_col+j]*g_mat[k];
			}
			nex_data[i*g_col+j]=sum;//存储列平滑后行的数据
		}
	}
	//第1种方法：后dim/2行没有平滑
	for (int m=g_row-r;m<g_row;m++)
		for (int j=0;j<g_col;j++)
			nex_data[m*g_col+j]=temp_row[m*g_col+j];//直接赋值
	
	delete []g_mat;
	delete []temp_row;

	return nex_data;
}
/****************************************************** 
* 函数名称： 
*  Build_Dog_Space()
*  
* 函数参数：
* SCALESPACE **ss-高斯空间
* int octaves-高斯空间的组数
* int intvls-高斯空间的层数
*
*
* 返回值：
* 0为失败，1为成功
*  
* 说明：生成差分高斯空间
********************************************************/
void Build_Dog_Space(SCALESPACE **ss,int octaves,int intvls)
{
	#ifdef WIN32
	int pros_num=omp_get_num_procs();
	omp_set_num_threads(pros_num*3/4);
#pragma omp parallel for
#endif

	for (int o=0;o<octaves;o++)
	{
		for (int s=0;s<intvls+2;s++)
		{
			int temp=ss[o][s].gus_row*ss[o][s].gus_column;
			float *tempdog=new float[temp];
			for (int i=0;i<temp;i++)
			{
				tempdog[i]=ss[o][s+1].gauss_space[i]-ss[o][s].gauss_space[i];
			}
			ss[o][s].dog_space=tempdog;
		}
	}
}

/****************************************************** 
* 函数名称： 
*  DetectKeypoint()
*  
* 函数参数：
* SCALESPACE **ss - 尺度空间
* int octaves     - 尺度空间组数
* int intvls      - 每组尺度空间的高斯层数（总数-3）
* float contr_thr - 对比度阈值
* float curv_thr  - 边缘检测的阈值
* int &keynum     - 关键点的个数
*
* 返回值：
* 检测到的关键点，剔除了边缘点和低对比度的点
*  
* 说明：极值点检测
********************************************************/
Key_Point* DetectKeypoint(SCALESPACE **ss,int octaves,int intvls,float contr_thr,float curv_thr,int &keynum)
{
	keynum=0;
	Key_Point *keypoints= NULL;
	//根据opencv中的公式将这阈值除以每组空间的层数-3
	float prelim_contr_thr=0.5f*contr_thr/intvls;
	#ifdef WIN32
	int pros_num=omp_get_num_procs();
	omp_set_num_threads(pros_num*3/4);
#pragma omp parallel for
#endif

	for (int o=0;o<octaves;o++)
	{
		for (int s=1;s<=intvls;s++)
		{
			int row=ss[o][s].gus_row;
			int col=ss[o][s].gus_column;
			for (int r=SIFT_IMG_BORDER;r<row-SIFT_IMG_BORDER;r++)
			{
				for (int c=SIFT_IMG_BORDER;c<col-SIFT_IMG_BORDER;c++)
				{
					if(fabs(ss[o][s].dog_space[r*col+c])>prelim_contr_thr)
					{
						if (is_extremum(ss,o,s,r,c,col))
						{
							float xs=0, xr=0, xc=0;
							int ms=s, mr=r, mc=c;
							if (interp_extremum(ss,o,ms,mr,mc,intvls,col,contr_thr,xs,xr,xc))
							{	
								if(!is_too_edge_like(ss, o, ms, mr, mc,col,curv_thr))
								{	
									keynum++;
									Key_Point* k=new Key_Point[1];//采用链表存储关键点
									k->next=keypoints;
									keypoints=k;
									k->key_row=mr; //行
									k->key_column=mc;  //列
									k->initl_row=(mr+xr)*pow(2.0f,o);//在原始影像上的行位置
									k->initl_column=(mc+xc)*pow(2.0f,o);//在原始影像上的列位置
									k->key_octave=o;//组
									k->key_intvl=ms;//层
									k->sub_intvl=xs;//精确定位层的改正数
								}
							}
						}
					}
				}
			}
		}
	}
	return keypoints;
}
/****************************************************** 
* 函数名称： 
* is_extremum()
*  
* 函数参数：
* SCALESPACE ss
* int o 尺度空间的组数
* int s 尺度空间的层数
* int r 差分高斯图像的行
* int c 差分高斯图像的列
* int col 差分高斯每行的像素个数
* 
* 返回值：
* 0为失败，1为成功
*  
* 说明：判断是否为极值点
********************************************************/
bool is_extremum(SCALESPACE **ss,int o,int s,int r,int c,int col)
{
	float val=ss[o][s].dog_space[r*col+c];

	//大于零，检测极大值
	if( val > 0 )
	{
		for(int i = -1; i <= 1; i++ )
			for( int j = -1; j <= 1; j++ )
				for(int k = -1; k <= 1; k++ )
					if( val<ss[o][s+i].dog_space[(r+j)*col+c+k])//不是极大值
						if (!(i==0&&j==0&&k==0)) //排除正在检测的点，因为（o,s,r,c）肯定=val
						    return 0;
	}

//小于零，检测极小值
	else
	{
		for(int i = -1; i <= 1; i++ )
			for(int j = -1; j <= 1; j++ )
				for(int k = -1; k <= 1; k++ )
					if( val>ss[o][s+i].dog_space[(r+j)*col+c+k] )//不是极小值
						if (!(i==0&&j==0&&k==0))
							return 0;
	}
	return 1;
}
/****************************************************** 
* 函数名称： 
* interp_extremum()
*  
* 函数参数：
* SCALESPACE ss
* int o 尺度空间的组数
* int &ms 尺度空间的层数,并返回精确定位后尺度空间的层数
* int &mr 差分高斯图像的行，并返回精确定位后差分高斯图像的行
* int &mc 差分高斯图像的列，并返回精确定位后差分高斯图像的列
* int intvls 每组高斯层数-3
* int col 差分高斯图像的总列数
* int contr_thr 对比度阈值
* float &xs 尺度层的改正数
* float &xr 行方向改正数
* float &xc 列方向改正数
*
* 返回值：
* 0为失败，1为成功
*  
* 说明：精确定位极值点
********************************************************/
bool interp_extremum(SCALESPACE **ss,int o,int &ms,int &mr,int &mc,int intvls,int col,float contr_thr,float &xs,float &xr,float &xc)
{
	int i=0;
	while(i<SIFT_MAX_INTERP_STEPS)
	{
		if (!interp_step( ss, o, ms, mr, mc,col, xs, xr, xc))//精确差值是否成功
		{
			return 0;
		}
		
		if( fabs(xs)<0.5 && fabs(xr)<0.5 && fabs(xc)<0.5)
			break;

		mc += floor( xc+0.5f );
		mr += floor( xr+0.5f );
		ms += floor( xs+0.5f );
		if( ms<1||ms>intvls||//尺度超出范围   
			mc < SIFT_IMG_BORDER||mr < SIFT_IMG_BORDER ||//行或列超出边界
			mr >= ss[o][ms].gus_row - SIFT_IMG_BORDER ||mc >= ss[o][ms].gus_column - SIFT_IMG_BORDER )
		{
			return 0;
		}
		i++;
	}
	if( i >= SIFT_MAX_INTERP_STEPS )
		return 0;//循环迭代的次数超过了阈值

	//计算差值后的dog空间中的对比度的值，即在Lowe论文中公式（3）带入公式（2）所求得的值
	float contr;
	contr=interp_contr(ss,o,ms,mr,mc,col,xs,xr,xc);
	if( fabs(contr)<contr_thr/intvls)//关键点的灰度是否大于阈值
		return 0;
	return 1;
}

/****************************************************** 
* 函数名称： 
* interp_contr()
*  
* 函数参数：
* SCALESPACE ss 尺度空间
* int o 尺度空间的组数
* int s 尺度空间的层数
* int r 差分高斯图像的行
* int c 差分高斯图像的列
* int col 尺度空间总列数
* float xs 层的改正数
* float xr 行的改正数
* float xc 列的改正数
*
* 返回值：
* float contr 精确定位极值点后关键点对应的灰度值
*  
* 说明：精确定位极值点后关键点对应的灰度值
********************************************************/
float interp_contr(SCALESPACE **ss,int o,int s,int r,int c,int col,float xs,float xr,float xc)
{
	float Dx,Dy,Ds,contr;
	/*一阶偏导矩阵
	| Dx  Dy  Ds |
	*/
	Dx= (ss[o][s].dog_space[r*col+c+1]-ss[o][s].dog_space[r*col+c-1])*0.5f;
	Dy= (ss[o][s].dog_space[(r-1)*col+c]-ss[o][s].dog_space[(r+1)*col+c])*0.5f;
	Ds= (ss[o][s+1].dog_space[r*col+c]-ss[o][s-1].dog_space[r*col+c])*0.5f;
	contr=ss[o][s].dog_space[r*col+c]+(Dx*xc+Dy*xr+Ds*xs)*0.5f;//精确定位关键点后，关键点所对应的灰度值
	return contr;
}
/****************************************************** 
* 函数名称： 
* is_too_edge_like()
*  
* 函数参数：
* SCALESPACE** ss 尺度空间
* int o 尺度空间的组数
* int s 尺度空间的层数
* int r 差分高斯图像的行
* int c 差分高斯图像的列
* int col 差分高斯图像的总列数
* float curv_thr 剔除边缘点的阈值
*
* 返回值：
* 0为失败，1为成功
*  
* 说明：剔除边缘点
********************************************************/
bool is_too_edge_like(SCALESPACE** ss, int o, int s,int r,int c,int col,float curv_thr)
{
	//计算关键点的曲率，从而去除边缘点
	float val, dxx, dyy, dxy, tr, det;
	val=ss[o][s].dog_space[r*col+c];
	dxx=ss[o][s].dog_space[r*col+c-1] + ss[o][s].dog_space[r*col+c+1]-2*val;
	dyy = ss[o][s].dog_space[(r-1)*col+c] + ss[o][s].dog_space[(r+1)*col+c]-2* val;
	dxy = ( ss[o][s].dog_space[(r+1)*col+c-1]+ss[o][s].dog_space[(r-1)*col+c+1]-ss[o][s].dog_space[(r+1)*col+c+1]-ss[o][s].dog_space[(r-1)*col+c-1])*0.25f;
	tr = dxx + dyy;//矩阵的迹
	det = dxx*dyy - dxy*dxy;//矩阵的行列式
	//二阶导数行列式的值小于零不存在极值
	if( det <= 0 )
		return 1;
	// 计算比率，根据Lowe文章中公式（4）计算curv_thr = (tr*tr)/det;
	if(tr*tr/det>=(curv_thr+1.0f)*(curv_thr+1.0f)/curv_thr)
	{
		return 1;
	}
	return 0;
}

/****************************************************** 
* 函数名称： 
*  inter_extremum()
*  
* 函数参数：
* SCALESPACE **ss 尺度空间
* int o     尺度空间的组数
* float s   尺度空间的层数
* float r   差分高斯图像的行数
* float c   差分高斯图像的列数
* int col   差分高斯图像列的总数
* float &xs 层的改正数
* float &xr 行的改正数
* float &xc 列的改正数
*
* 返回值：
* 0为失败，1为成功
*  
* 说明：精确定位极值点时，计算x，y，s的改正值
********************************************************/
bool interp_step(SCALESPACE **ss,int o, int s, int r, int c,int col,float &xs,float &xr,float &xc)
{
	float Dx,Dy,Ds,Dxx,Dyy,Dss,Dxy,Dxs,Dys;
	
	double** a=new double*[3];//用于存储二阶偏导
	double** inv_a=new double*[3];//用于存储二阶偏导的逆矩阵
	for (int i=0;i<3;i++)
	{
		a[i]=new double[3];
		inv_a[i]=new double[3];
	}
	float val=ss[o][s].dog_space[r*col+c];
	/*一阶偏导矩阵
	| Dx  Dy  Ds |
	*/
	Dx= (ss[o][s].dog_space[r*col+c+1]-ss[o][s].dog_space[r*col+c-1])*0.5f;
	Dy= (ss[o][s].dog_space[(r-1)*col+c]-ss[o][s].dog_space[(r+1)*col+c])*0.5f;
	Ds= (ss[o][s+1].dog_space[r*col+c]-ss[o][s-1].dog_space[r*col+c])*0.5f;

	/*二阶偏导矩阵
	/ Dxx  Dxy  Dxs \
	| Dxy  Dyy  Dys |
	\ Dxs  Dys  Dss /
	*/
	Dxx = ss[o][s].dog_space[r*col+c-1] + ss[o][s].dog_space[r*col+c+1]-2*val;
	Dyy = ss[o][s].dog_space[(r-1)*col+c] + ss[o][s].dog_space[(r+1)*col+c]-2*val;
	Dss = ss[o][s-1].dog_space[r*col+c] + ss[o][s+1].dog_space[r*col+c]-2*val;
	Dxy = (ss[o][s].dog_space[(r-1)*col+c+1] + ss[o][s].dog_space[(r+1)*col+c-1] - ss[o][s].dog_space[(r+1)*col+c+1] - ss[o][s].dog_space[(r-1)*col+c-1])*0.25f;
	Dxs = (ss[o][s-1].dog_space[r*col+c-1] + ss[o][s+1].dog_space[r*col+c+1] - ss[o][s+1].dog_space[r*col+c-1] - ss[o][s-1].dog_space[r*col+c+1])*0.25f;
	Dys = (ss[o][s-1].dog_space[(r+1)*col+c] + ss[o][s+1].dog_space[(r-1)*col+c] - ss[o][s+1].dog_space[(r+1)*col+c] - ss[o][s-1].dog_space[(r-1)*col+c])*0.25f;
								
	a[0][0]=Dxx;
	a[0][1]=Dxy;
	a[0][2]=Dxs;
	a[1][0]=Dxy;
	a[1][1]=Dyy;
	a[1][2]=Dys;
	a[2][0]=Dxs;
	a[2][1]=Dys;
	a[2][2]=Dss;

	//求二阶偏导矩阵的逆
	if (sift_Inv(a,inv_a,3))
	{
		// 计算x，y，s的改正值
		xr=-(Dx*inv_a[0][0]+Dy*inv_a[0][1]+Ds*inv_a[0][2]);
		xc=-(Dx*inv_a[1][0]+Dy*inv_a[1][1]+Ds*inv_a[1][2]);
		xs=-(Dx*inv_a[2][0]+Dy*inv_a[2][1]+Ds*inv_a[2][2]);
		for (int i=0;i<3;i++)
		{
			delete []a[i];
			delete inv_a[i];
		}
		delete []a;
		delete []inv_a;

		return 1;
	} 
	else
	{
		for (int i=0;i<3;i++)
		{
			delete []a[i];
			delete inv_a[i];
		}
		delete []a;
		delete []inv_a;

		return 0;
	}	
}

/****************************************************** 
* 函数名称： 
* calc_scales()
*  
* 函数参数：
* Key_Point *keypoint 关键点
* int keynumber 关键点的个数
*
* 返回值：
* 0为失败，1为成功
*  
* 说明：求关键点的尺度
********************************************************/
void calc_scales(Key_Point *keyhead,int keynumber)
{
	int i;
	double intvl;
	Key_Point *keypoint=keyhead->next;
	for( i = 0; i < keynumber; i++ )
	{
		intvl = keypoint->key_intvl + keypoint->sub_intvl;
		keypoint->scl =SIFT_SIGMA * pow( 2.0, keypoint->key_octave + intvl / SCALESPEROCTAVE );//keypoint->scl只为显示使用
		keypoint->scl_octave = SIFT_SIGMA * pow( 2.0, intvl / SCALESPEROCTAVE );
		keypoint=keypoint->next;
	}
	delete []keypoint;
}
/****************************************************** 
* 函数名称： 
* calc_oritation()
*  
* 函数参数：
* Key_Point* keypoint 关键点
* SCALESPACE** ss 尺度空间
* int &keynumber 关键点的个数
* 
* 返回值：
* keypoints 关键点
*  
* 说明：计算主方向
********************************************************/
Key_Point* calc_oritation(Key_Point* keyhead,SCALESPACE** ss,int &keynumber)
{
	float *hist;
	float omax,bin,/*mag_peak,*/PI2=PI*2.0f;//mag_peak表示差值后的主方向对应的梯度值
	int addnum=0, n=SIFT_ORI_HIST_BINS;   //addnum是由于多个主方向而引起关键点个数增加后总点数
	Key_Point *keypoints=NULL;
	Key_Point *k=keyhead->next;
	for(int i = 0; i < keynumber; i++ )
	{ 
		hist = ori_hist(ss[k->key_octave][k->key_intvl].gauss_space,k->key_row,k->key_column,ss[k->key_octave][k->key_intvl].gus_row,
			ss[k->key_octave][k->key_intvl].gus_column,SIFT_ORI_HIST_BINS,int(floor(SIFT_ORI_RADIUS*k->scl_octave+0.5f)),SIFT_ORI_SIG_FCTR*k->scl_octave);
		for(int j=0;j<SIFT_ORI_SMOOTH_PASSES;j++)
			smooth_ori_hist(hist,n);
		omax = dominant_ori(hist,n);

		//根据直方图中的最大值，计算出直方图中所有的主方向即大于0.8倍的最大值的所有方向
		float mag_thr=omax*SIFT_ORI_PEAK_RATIO;
		int l, r;
		for(int j=0;j<n;j++ )
		{
			l=(j==0)?n-1:j-1;
			r=(j+1)% n;
			if( hist[j]>hist[l] && hist[j]>hist[r] && hist[j]>=mag_thr )
			{
				bin =j+interp_hist_peak(hist[l], hist[j], hist[r]/*,mag_peak*/);
				bin =(bin<0)?n+bin:((bin>=n)?bin-n:bin);
				addnum++;
				Key_Point* kk=new Key_Point[1];//存储关键点
				kk->next=keypoints;
				keypoints=kk;
				kk[0].initl_column=k->initl_column;
				kk[0].key_column=k->key_column;
				kk[0].initl_row=k->initl_row;
				kk[0].key_intvl=k->key_intvl;
				kk[0].key_octave=k->key_octave;
				kk[0].key_row=k->key_row;
				kk[0].scl=k->scl;
				kk[0].scl_octave=k->scl_octave;
				kk[0].sub_intvl=k->sub_intvl;
				kk[0].ori=((PI2 * bin )/n)/*-PI*/;//ori的取值范围（0,2*PI）
				//kk[0].mag=mag_peak;
			}
		}
		k=k->next;
		delete []hist;
	}
	keynumber=addnum;
	delete []k;

	Key_Point* pnode=keyhead;
	Key_Point* tempkey;
	while(pnode!=NULL)
	{		
		tempkey=pnode;
		pnode=pnode->next;
		delete []tempkey;
	}
	keyhead=NULL;
	//delete []keyhead;
	return keypoints;
}
/****************************************************** 
* 函数名称： 
* ori_hist()
*  
* 函数参数：
* float* guss  高斯空间数据
* int r, int c 行号和列号
* int col      总列数
* int row      总行数
* int n        方向直方图采用的总方向数（36）
* int rad      特征点的梯度值
* floa sigma   高斯权重采用的标准差（等于1.5倍的窗口宽度）      
*
* 返回值：
* float* hist  特征点的直方图
*  
* 说明：主方向的方向直方图
********************************************************/
float* ori_hist(float* guss, int r, int c, int row,int col, int n, int rad, float sigma)
{
	
	float mag=0,ori=0,w,PI2=PI*2.0f;//mag梯度 ori主方向 w高斯权重 
	int bin;//方向值
	float sigma2 = 2.0f * sigma * sigma;
	float* hist=new float[n];//36个方向的统计直方图
	memset(hist,0,n*sizeof(float));

	for(int i = -rad; i <= rad; i++ )
		for(int j = -rad; j <= rad; j++ )
			if(calc_grad_mag_ori(guss,r+i,c+j,row,col,mag,ori))
			{
				w = exp(-(i*i+j*j)/sigma2);
				bin = floor(n*(ori/*+PI*/)/PI2+0.5f);///
				bin = (bin<n)?bin:0;
				hist[bin]+=w*mag;
			}
	return hist;
}
/****************************************************** 
* 函数名称： 
* calc_grad_mag_ori()
*  
* float* guss  高斯空间数据
* int r, int c 行号和列号
* int col      总列数
* int row      总行数
* float &mag   特征点的梯度值
* float &ori   特征点的方向值 
* 
* 返回值：
* 0为失败，1为成功
*  
* 说明：计算点的梯度和方向
********************************************************/

bool calc_grad_mag_ori(float* guss, int r, int c,int row,int col, float &mag, float &ori )
{
	float dx, dy;
	if(r>0 && r<row-1 && c>0 && c<col-1)
	{
		dx=guss[r*col+c+1]-guss[r*col+c-1];     //dx
		dy=guss[(r-1)*col+c]-guss[(r+1)*col+c]; //dy
		mag = sqrt( dx*dx + dy*dy );
		ori = atan2(dy, dx);      //求得主方向取值范围为（-PI，PI）
		ori=(ori>0)?ori:(ori+PI+PI);
		return 1;
	}
	else
		return 0;
}
/****************************************************** 
* 函数名称： 
* smooth_ori_hist()
*  
* float* hist  方向直方图
* int n        方向直方图采用的总方向数（36）
* 
* 返回值：
* 
*  
* 说明：对直方图进行平滑，模板为（0.25 0.5 0.25）
********************************************************/
void smooth_ori_hist( float* hist, int n )
{
	float prev,tmp;
	prev = hist[n-1];
	for(int i = 0; i < n; i++ )
	{
		tmp = hist[i];
		hist[i]=0.25f*prev+0.5f*hist[i]+0.25f*((i+1==n)?hist[0]:hist[i+1] );
		prev = tmp;
	}
}
/****************************************************** 
* 函数名称： 
* dominant_ori()
*  
* float* hist 方向直方图
* int n       方向直方图采用的总方向数（36）
* 
* 返回值：
* 直方图中的最大值
*  
* 说明：计算直方图中的最大值
********************************************************/
float dominant_ori( float* hist, int n )
{
	float omax;
	omax = hist[0];
	for(int i = 1; i < n; i++ )
		omax=( hist[i]>omax )?hist[i]:omax;
	return omax;
}
/****************************************************** 
* 函数名称： 
* interp_hist_peak()
*  
* 函数参数：
* float lef 主方向左边的梯度值
* float cent  主方向的梯度值
* float rigt  主方向右边的梯度值
* float &mag_peak 主方向的抛物线差值后的梯度值
*
* 返回值：
* 抛物线差值后主方向的值
*  
* 说明：对主方向进行抛物线差值
********************************************************/
float interp_hist_peak(float lef,float cent,float rigt/*,float &mag_peak*/)
{
	//mag_peak=cent-0.125f*(lef-rigt)*(lef-rigt)/(lef-2.0f*cent+rigt);
	return 0.5f*(lef-rigt)/(lef-2.0f*cent+rigt);
}
/****************************************************** 
* 函数名称： 
* compute_descriptors()
*  
* 函数参数：
* Key_Point* feat 特征点
* SCALESPACE **ss 尺度空间
* int keynumber 关键点的个数
* 
* 返回值：
*  
* 说明：计算特征点的128维向量即描述子
********************************************************/
void compute_descriptors( Key_Point* feat, SCALESPACE **ss, int keynum)
{
	int d=SIFT_DESCR_WIDTH,n=SIFT_DESCR_HIST_BINS;
	for(int i=0;i<keynum;i++)
	{		
		int num=SIFT_DESCR_BINS;//描述子的维数 128
		//计算直方图
		float *hist= descr_hist(ss[feat[i].key_octave][feat[i].key_intvl].gauss_space,feat[i].key_row,
			feat[i].key_column,ss[feat[i].key_octave][feat[i].key_intvl].gus_row,
			ss[feat[i].key_octave][feat[i].key_intvl].gus_column,feat[i].ori,feat[i].scl_octave,d,n,num);
		//归一化

		normalize_descr(hist,num);
		//为了减小较大梯度值对特征向量的影响，将归一化后的梯度值限制在SIFT_DESCR_MAG_THR=0.2内。
		for(int j = 0; j < num; j++ )
		{
			if(hist[j]>SIFT_DESCR_MAG_THR)
				hist[j]=SIFT_DESCR_MAG_THR;
		}
		//再次归一化
		normalize_descr(hist,num);
		//特征向量乘以系数512，从而将特征向量化为整数
		for(int j=0; j<num; j++)
		{
			feat[i].descriptor[j] =float( min(255, int(SIFT_INT_DESCR_FCTR * hist[j])));
		}
		delete []hist;
	}
}
/****************************************************** 
* 函数名称： 
* descr_hist()
*  
* 函数参数：
*
* float* guss  高斯空间数据
* int r, int c 行号和列号
* int row      总行数
* int col      总列数
* float &ori   特征点的方向值
* float scl    只与层相关计算所得的尺度
* int d        描述子采用的区域参数  4
* int n        描述子直方图的列数  8
* int num      描述子的维数
* 返回值：
* float *hist 128维的直方图描述
*  
* 说明：计算特征点的直方图描述
********************************************************/
float* descr_hist(float* guss,int r,int c,int row,int col,float ori,float scl,int d,int n,int num)
{
	
	float cos_t,sin_t, hist_width, exp_denom, r_rot, c_rot, grad_mag,
		grad_ori, w, rbin, cbin, obin, bins_per_rad, PI2=2.0f*PI;
	int radius;
	//定义并初始化描述子直方图
	float* hist=new float[num];
	memset(hist,0,num*sizeof(float));

	cos_t=cos(ori);
	sin_t=sin(ori);
	bins_per_rad=n/PI2;//用于弧度和8个方向之间转换
	exp_denom=d*d*0.5f;//用于计算高斯权重
	hist_width=SIFT_DESCR_SCL_FCTR*scl;//4×4块中每块的边长	
	radius=int(hist_width*sqrt(2.0f)*(d+1.0f)*0.5f+0.5f);//统计描述子的区域半径
	for(int i = -radius; i <= radius; i++ )
		for(int j = -radius; j <= radius; j++ )
		{
			c_rot=(j*cos_t-i*sin_t)/hist_width;//像点旋转后在以特征点为原点的块（4X4）坐标系中的列号
			r_rot=(j*sin_t+i*cos_t)/hist_width;
			rbin=r_rot+d/2-0.5f;
			cbin=c_rot+d/2-0.5f;//像点旋转后在以左上角为原点的块（4X4）坐标系中的列号

			if(rbin>-1.0f && rbin<d && cbin>-1.0f && cbin<d)//判断旋转后的点在统计区域内
				if(calc_grad_mag_ori(guss, r + i, c + j, row, col, grad_mag, grad_ori))
				{
					grad_ori -= ori;  //方向归一化
					while( grad_ori < 0 )//因为grad-ori和ori的取值范围都是（-PI，PI），所以需要对相减的结果进行调整
						grad_ori += PI2;
					while( grad_ori >= PI2 )
						grad_ori -= PI2;

					obin = grad_ori * bins_per_rad;//方向所属于的直方图的列
					w = exp(-(c_rot*c_rot+r_rot*r_rot)/exp_denom);//以块的长度为基本单位来计算高斯权重
					interp_hist_entry(hist, rbin, cbin, obin, grad_mag*w, d, n);
				}
		}
	return hist;
}
/****************************************************** 
* 函数名称： 
* interp_hist_entry()
*  
* 函数参数：
*
* float* hist  128维的直方图
* float rbin   点在4个行块中的行号
* float cbin   点在4个列块中的列号
* float obin   点的方向在8个方向值中的序号
* float mag    特征点的方向值
* int d        描述子采用的区域参数  4
* int n        描述子直方图的列数  8
* 
* 返回值：
* 
*  
* 说明：对于特征点邻域内的每个像点，根据其梯度值和权重（1-d），计算它对与其相邻
的行（如rbin=2.3，则其对于第二行和第三行有贡献）、列和方向的梯度值的贡献
********************************************************/
void interp_hist_entry(float* hist,float rbin,float cbin,float obin,float mag,int d,int n)
{
	float d_r, d_c, d_o, v_r,v_c , v_o;
	int r0, c0, o0, rb, cb, ob;

	r0 = int(floor( rbin ));
	c0 = int(floor( cbin ));
	o0 = int(floor( obin ));
	d_r = rbin - r0;
	d_c = cbin - c0;
	d_o = obin - o0;

	for(int r=0;r<=1;r++)
	{
		rb = r0 + r;
		if(rb>=0 && rb<d)
		{
			v_r=mag*((r==0)?1.0f-d_r:d_r);
			for(int c=0;c<=1;c++)
			{
				cb = c0 + c;
				if(cb>=0 && cb<d)
				{
					v_c =v_r*((c==0)?1.0f-d_c:d_c );
					for(int o=0;o<=1;o++)
					{
						ob =(o0+o)%n;
						v_o =v_c*((o==0)?1.0f-d_o:d_o);
						hist[(rb*d+cb)*n+ob]+= v_o;//将贡献值加入到128维直方图中
					}
				}
			}
		}
	}
}
/****************************************************** 
* 函数名称： 
* normalize_descr()
*  
* 函数参数：
*
* float* hist  描述子直方图
* int d        描述子维数 128
* 
* 返回值：
* 
*  
* 说明：对于128维描述子进行归一化
********************************************************/
void normalize_descr(float* hist,int d)
{
	float cur,len_inv,len_sq=0.0f;

	for(int i=0;i<d;i++)
	{
		cur=hist[i];
		len_sq+=cur*cur;
	}
	len_inv=1.0f/sqrt(len_sq);

	for(int i=0;i<d;i++)
		hist[i]*=len_inv;
}

/****************************************************** 
* 函数名称： 
* adjust_for_img_down()
*  
* 函数参数：
* Key_Point* feat - 关键点
* int keynumber   - 关键点的个数
*
* 返回值：
* 0为失败，1为成功
*  
* 说明：对初始影像进行降采样时，需要对关键点的坐标进行调整
********************************************************/
bool adjust_for_img_down(Key_Point* feat,int keynumber)
{
	for(int i = 0; i < keynumber; i++ )
	{
		feat[i].initl_row *= 2.0f;
		feat[i].initl_column *= 2.0f;
		feat[i].scl *= 2.0f;
	}
	return 1;
}

/****************************************************** 
* 函数名称： 
* adjust_for_img_grid()
*  
* 函数参数：
* Key_Point* feat - 子块关键点
* int keynumber   - 子块关键点的个数
* int col         - 子块的列号
* int row         - 子块的行号
* int gridwidth   - 子块宽
* int imgheight   - 整张影像的高
*
* 返回值：
* 0为失败，1为成功
*  
* 说明：分块的图像，需要对关键点的坐标进行调整，恢复到原始图像上的坐标
********************************************************/
void adjust_for_img_grid(Key_Point* feat,int keynumber,int col,int row,int gridwidth,int imgheight)
{
	for(int i = 0; i < keynumber; i++ )
	{
		feat[i].initl_row =row*gridwidth+feat[i].initl_row ;//关键点坐标的行initl_row，是以左下角为坐标原点
		//feat[i].initl_row =imgheight-(row*gridwidth+feat[i].initl_row) ;//关键点坐标的行initl_row是以左上角为坐标原点
		feat[i].initl_column = col*gridwidth+feat[i].initl_column ;
	}
}
/****************************************************** 
* 函数名称： 
* Inv()
*  
* 函数参数：
* double** a         - 待求逆的矩阵
* double** inva      - 求逆后的矩阵
* int dim            - 矩阵的维数
*
* 返回值：
* 矩阵的逆矩阵
*  
*说明：矩阵求逆
********************************************************/

int sift_Inv(double** a,double** inva,int dim )
{
	int i,j,k,p;//记录列主元的列数p及数值c_max
	double temp,c_max=0;
	//单位矩阵赋值
	double **b=new double *[dim];
	for (i=0;i<dim;i++)
	{
		b[i]=new double [dim];
		memset(b[i],0,sizeof(double)*dim);
		b[i][i]=1;
	}

	//高斯-约当消元
	for(i=0;i<dim;i++)
	{  
		/*******************选列主元**********************/
		for(j=i+1,c_max=a[i][i],p=i;j<dim;j++)//c_max=a[i][i]不能放在循环体内，否则没有找最后一列的最大值
		{  
			if(fabs(a[j][i])>fabs(c_max))  
			{p=j; c_max=a[j][i]; }     
		}  
		//奇异矩阵不存在逆矩阵
		if(fabs(c_max)<1e-10) 
		{return 0;}
		/*****************换行、归一*********************/
		if (i!=p)
		{
			for(j=0;j<dim;j++)
			{   
				temp=a[i][j];
				a[i][j]=a[p][j]/c_max; 
				a[p][j]=temp;
				//单位矩阵
				temp=b[i][j];
				b[i][j]=b[p][j]/c_max; 
				b[p][j]=temp;
			}
		}
		else
		{
			for (j=0;j<dim;j++)
			{
				a[i][j]/=c_max;
				b[i][j]/=c_max;
			}
		}
		/**********************消元**********************/
		for( k=0;k<dim;k++)  
		{
			if(k==i) continue;
			for(j=0,temp=a[k][i];j<dim;j++)
			{
				a[k][j]-=temp*a[i][j];
				b[k][j]-=temp*b[i][j];
			}
		}
	}
	//将逆矩阵的值赋值给输入函数参数
	for (i=0;i<dim;i++)
	{
		for (j=0;j<dim;j++)
		{
			inva[i][j]=b[i][j];
		}
	}
	//释放内存
	for (int i=0;i<dim;i++)
	{
		delete[]b[i];
	}
	delete []b;
	return 1;
}


DLL_EXPORT Key_Point* SiftFeaturesFloat( float *imgdata, int colum, int row, int &keynumber )
{
	return _sift_features(imgdata, colum, row, keynumber);
}


//删除链表
//void free_keypoint(Key_Point* keyhead)
//{
//	Key_Point* tempkey;
//	while(tempkey!=NULL)
//	{
//		tempkey=keyhead;
//		keyhead=keyhead->next;
//		delete []tempkey;
//	}
//	keyhead=NULL;
//}