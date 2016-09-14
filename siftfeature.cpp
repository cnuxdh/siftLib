/************************************************************************/
/*��ȡSIFT����                                  */
/************************************************************************/

/**************************************ͷ�ļ�*****************************************/
//#include "stdafx.h"
#include <iostream>
#include <math.h>

#ifdef WIN32
#include "omp.h"
#endif


#include "siftfeature.h"

using namespace std;




/**************************************����ͼ������*****************************************/
Key_Point* _sift_features( float *, int , int , int & );//sift������ȡ
float* GetGridData(float* ,int ,int ,int ,int ,int ,int );//�������ж�ȡ�ӿ������
float* Down_Sample(float *,int &,int &);//������
SCALESPACE**  Build_Guss_Space(float *,int,int,int);//������˹�ռ�
float* Guss_Smooth(float *,int,int,float);//��˹ƽ��
void Build_Dog_Space(SCALESPACE**,int,int);//������ָ�˹�ռ�
Key_Point* DetectKeypoint(SCALESPACE**,int,int,float,float,int &);//��⼫ֵ��
bool is_extremum(SCALESPACE**,int,int,int,int,int);//�ж��Ƿ�Ϊ��ֵ
bool interp_step(SCALESPACE **,int ,int ,int ,int,int ,float &,float &,float &);//��ά������ϣ���ȡ��ֵ��ĸ���ֵ
bool interp_extremum(SCALESPACE **,int,int&,int&,int&,int,int,float,float &,float &,float &);//���뼫ֵ�㣬��ȷ��λ�ؼ���
float interp_contr(SCALESPACE **,int,int,int,int,int,float,float,float);//���㾫ȷ��λ�ļ�ֵ���Ӧ�ĻҶ�ֵ
bool is_too_edge_like(SCALESPACE **,int,int,int,int,int,float);
float inter_contr(int , int , int , int ,float ,float ,float);//���㾫ȷ��λ��ؼ���ĶԱȶ�
void calc_scales(Key_Point *,int);//����ؼ������ڵĳ߶�
Key_Point* calc_oritation(Key_Point* ,SCALESPACE** ,int &);//����ؼ����������
bool calc_grad_mag_ori(float*, int, int, int, int, float &g, float &);//�����ķ�����ݶ�
float dominant_ori( float*, int );//����ؼ����������
float* ori_hist(float* , int , int , int ,int , int , int , float);//���㷽��ֱ��ͼ
void smooth_ori_hist( float* , int );//��ֱ��ͼ����ƽ��
float* descr_hist(float* ,int ,int ,int ,int ,float ,float ,int ,int ,int);//������������������ֱ��ͼ
void normalize_descr(float* ,int );//��128ά�����ӽ��й�һ��
void interp_hist_entry(float* ,float ,float ,float ,float ,int ,int );//����ÿ������������ĵ��������ֱ��ͼ�Ĺ���
float interp_hist_peak(float ,float ,float/* ,float &*/);//��ֱ��ͼ���������߲�ֵ
void compute_descriptors( Key_Point* , SCALESPACE **, int );//����������������
bool adjust_for_img_down(Key_Point* ,int );//��ԭʼӰ����н�����ʱ����Ҫ�Թؼ����������е���
int sift_Inv(double** a,double** inva,int dim );//�Ծ�������
void adjust_for_img_grid(Key_Point* feat,int keynumber,int col,int row,int gridwidth,int imgheight);//��ԭʼӰ����зֿ�ʱ����Ҫ�Թؼ����������е���
//void free_keypoint(Key_Point* keyhead);
/**************************************ͼ������ʵ��*****************************************/

/****************************************************** 
* �������ƣ� 
* sift_features()
*  
* ����������
* unsigned char* imagedata- GDAL��ȡ��ͼ������
* int colum      - ͼ������
* int row        - ͼ������
* int bands      - ͼ�񲨶���
* int &keynumber - �ؼ���ĸ���
* int &gridcol   - �з����ĸ���
* int &gridrow   - �з����ĸ���
*
* ����ֵ��
* ��ȡ��������
*  
*˵������ȡͼ���sift�����㣬�ж��Ƿ��ͼ����зֿ�
********************************************************/
Feature** sift_features( float* imgdata,int colum,int row,int &keynumber,int &gridcol,int &gridrow )
{
	/*********************************ͼ�����ݷֿ�***************************************/
	int gridwidth = SIFT_IMG_GRID_THRESHOLD; //��Ŀ�
	
	//ȷ����ĸߺͿ���Ŀ����
	gridcol=(colum+gridwidth-1)/gridwidth;       //�����ĸ���
	gridrow=(row+gridwidth-1)/gridwidth;         //�����ĸ���
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
	gridr[gridrow-1]=((row-1) % gridwidth)+1; //����ĩβС�� 
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

	//�ֿ���ȡ������
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

	//���ֿ��е����꣬ת��ΪԭʼӰ���ϵ�����
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
* �������ƣ� 
* GetGridData()
*  
* ����������
* float *imgdata  - ͼ������
* int rownum      - ����к�
* int colnum      - ����к�
* int imgwidth    - ԭʼͼ��Ŀ�
* int gridwidth   - ��׼��Ŀ�
* int row         - �������
* int col         - �������
* 
* ����ֵ��
* �������
*  
*˵������ȡ�ֿ������
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
* �������ƣ� 
* _sift_features()
*  
* ����������
* float *imgdata  - ͼ������
* int colum       - ���ݵ�����
* int row         - ���ݵ�����
* int &keynumber  - �ؼ���ĸ���
* 
* ����ֵ��
* ��ȡ��������
*  
*˵������ȡͼ���sift������
********************************************************/
Key_Point* _sift_features( float *imgdata, int colum, int row, int &keynumber )
{	
	//����к�������СֵС����ֵ��˿鲻����������
	if (min(colum,row)<GRID_REM_PIXEL_THRESHOLD)
	{
		keynumber=0;
		return NULL;
	}
	//���Ƚ��н�����������������
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
	/*1 ������˹�ռ�                                                        */
	/************************************************************************/
	//����߶ȿռ������������ֹ���ֲ�����ͼ��Ϊ0�����

	//�߶ȿռ������
	int octnum;
	if (OCTAVESTYE)
	{
		octnum=1;
	} 
	else
	{
		octnum=int(log(double(min(colum, row)))/log(2.0)-6);
		//����к��е���СֵΪ64��������Ϊ1
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
	/*2 ������ָ�˹�ռ�                                                    */
	/************************************************************************/

	Build_Dog_Space(scalespace,octnum,SCALESPEROCTAVE);

	/************************************************************************/
	/*3 ��⼫ֵ��                                                          */
	/************************************************************************/

	Key_Point* keyhead=new Key_Point[1];
	keyhead->next=DetectKeypoint(scalespace,octnum,SCALESPEROCTAVE,
		SIFT_CONTR_THRESHOLD,SIFT_CURV_THRESHOLD,keynumber);

	/************************************************************************/
	/*4 ����������ĳ߶�                                                    */
	/************************************************************************/
	calc_scales(keyhead,keynumber);

	/************************************************************************/
	/*5 ȷ��������                                                          */
	/************************************************************************/

	Key_Point* feathead=new Key_Point[1];
	feathead->next=calc_oritation(keyhead,scalespace,keynumber);

	//����������ʽ�洢���������Ϊ��ָ����ʽ�洢��ռ��ʱ��϶ࣩ
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
	/*6 �����������������                                                    */
	/************************************************************************/

	compute_descriptors(feat,scalespace,keynumber);

	//�����ԭʼͼ����н�����������Ҫ���������λ�ý��е���
	if (SIFT_IMG_DOWN)
	{
		adjust_for_img_down(feat,keynumber);
	}

	//ɾ���߶ȿռ�����
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
* �������ƣ� 
* Down_sample()
*  
* ����������
* float *img-ָ��ָ����н�����������
* int &col-�������ݵ����ã�Ϊͼ������������
* int &row-�������ݵ����ã�Ϊͼ����������߶�
*
* ����ֵ��
* ��С���ͼ��
*  
*˵������ͼ����и������
********************************************************/
float* Down_Sample(float *img,int &col,int &row)
{
	int row1=row/2;//����СΪԭ���Ķ���֮һ
	int col1=col/2;//����СΪԭ���Ķ���֮һ
	
	float *temp=new float[row1*col1];//������ʱ�洢������ͼ�����ݡ�
	for (int i=0;i<row1;i++)
	{
		for (int j=0;j<col1;j++)
		{
			temp[i*col1+j]=img[i*2*col+2*j];//ע��ڶ�����col����col1,��ʾ������ǰÿ�еĸ���������һ������2*col1��������2*col+1
			//temp[i*col1+j]=img[i*2*col+2*j]+img[i*2*col+2*j+1]+img[(i*2+1)*col+2*j]+img[(i*2+1)*col+2*j+1];
		}
	}
	row=row1;col=col1;
	return temp;
}
/****************************************************** 
* �������ƣ� 
* Build_Scale_Space()
*  
* ����������
* float *pimage-��������ָ����н�������ͼ������
* int oct_num-��˹������
* int col-Ϊͼ������������
* int row-Ϊͼ����������߶�
*
*����ֵ��
* SCALESPACE** ��˹�ռ�
*  
*˵����������˹�ռ�
********************************************************/
SCALESPACE** Build_Guss_Space(float* pimage,int oct_num,int col,int row)
{
	/***********************************ͼ��͸�˹������������ɸ�˹�ռ� ***************************************/
	SCALESPACE** gs=new SCALESPACE* [oct_num];//�洢��˹�ռ�
	for (int i=0;i<oct_num;i++)
	{
		gs[i]=new SCALESPACE [SCALESPEROCTAVE+3];
	}
	
	float k_sig= float(sqrt(SIFT_SIGMA * SIFT_SIGMA - SIFT_INIT_SIGMA *SIFT_INIT_SIGMA ));//������ײ�ƽ����sigmaֵ
	float intvl_k=pow(2.0f,1.0f/SCALESPEROCTAVE);//ÿ��֮��ĳ߶ȱ���ϵ��
	//��SCALESPEROCTAVE����1ʱ��Ϊ����kֵ̫��Ӷ����ƽ��ģ��̫�����������´���
	if (SCALESPEROCTAVE==1)
	{
		intvl_k=sqrt(3.0f);
	}
	float *intvl_sigma=new float[SCALESPEROCTAVE+3];//�洢ÿ���˹�˳߶�
	float pre_sigma,total_sigma;

	/**********************************��ʼ������ø�˹�����ĺ�***********************************/
	intvl_sigma[0]=SIFT_SIGMA;
	for(int j=1;j<SCALESPEROCTAVE+3;j++)
	{
		pre_sigma=pow(intvl_k,j-1)*SIFT_SIGMA;
		total_sigma=pre_sigma*intvl_k;//������һ���Ӱ������һ��Ӱ�����ó�
		intvl_sigma[j]=sqrt(total_sigma*total_sigma-pre_sigma*pre_sigma);//ÿһ����õĸ�˹ƽ�������ĺ�
	}

	/**********************************��˹������ɸ�˹�ռ�***********************************/
#ifdef WIN32
	int pros_num=omp_get_num_procs();
	omp_set_num_threads(pros_num*3/4);
#pragma omp parallel for
#endif
	for (int i=0;i<oct_num;i++)
	{
		for (int j=0;j<SCALESPEROCTAVE+3;j++)//ÿ����SCALESPEROCTAVE+3��
		{	
			if (i==0 && j==0)
			{	//��ԭʼӰ�����ƽ����������һ���һ��ĸ�˹Ӱ��
				gs[0][0].gauss_space=Guss_Smooth(pimage,col,row,k_sig);
				gs[0][0].gus_column=col;//��˹ͼ������
				gs[0][0].gus_row=row;//��˹ͼ������
			}
			else if (j==0)
			{//�ӵڶ��鿪ʼ,������һ���ж����ڵײ�߶ȵĸ�˹����н�����
				gs[i][0].gauss_space=Down_Sample(gs[i-1][SCALESPEROCTAVE].gauss_space,col,row);	
				gs[i][0].gus_column=col;//��˹ͼ������
				gs[i][0].gus_row=row;//��˹ͼ������
			}
			else
			{//�ӵڶ��㿪ʼƽ��
				gs[i][j].gauss_space=Guss_Smooth(gs[i][j-1].gauss_space,col,row,intvl_sigma[j]);//��һ��Ӱ��������һ��Ӱ��ƽ���Ļ����Ͻ���ƽ����
				gs[i][j].gus_column=col;//��˹ͼ������
				gs[i][j].gus_row=row;//��˹ͼ������
			}
		}
	}
	delete []intvl_sigma;

	return gs;

}
/****************************************************** 
* �������ƣ� 
*  Gauss_Smooth()
*  
* ����������
* float *pre_data-�洢��ƽ�����ݵ�����
* int g_col-Ϊ��ƽ�������е��������
* int g_row-Ϊ��ƽ�����ݵ��������߶�
* int g_sigma-��˹ƽ���ĺ�
*
* ����ֵ��
* float *���洢ƽ����ĸ�˹ͼ��
*  
* ˵������˹ƽ��
********************************************************/
float* Guss_Smooth(float *pre_data,int g_col,int g_row,float g_sigma)
{
	/*********************************** �����˹ģ��Ĵ�С ***************************************/
	int dim =1+2*(int)(2.0f*g_sigma);
	int c=(dim+1)/2;  //ģ��İ뾶��1
	int r=dim/2;       //ģ��İ뾶
	float s2=g_sigma*g_sigma;
	float *g_mat=new float[c];//���ڴ洢���ģ��
	float temp_sum=0;
	//float g_coe=1.0f/sqrt(2.0f*PI*g_sigma);
	g_mat[0]=1.0f/*g_coe*/;//ֻ���ø�˹�ˣ�û�и�˹����ǰ���ϵ����������м䴦��ֵΪexp(-(1.0*0*0)/(2.0*s2))=1
	temp_sum+=g_mat[0];//���ȼ����м���ֵ

	for(int i=1;i<c;i++) 
	{
		g_mat[i]=/*g_coe**/exp(-(1.0f*i*i)/(2.0f*s2));//��˹��;
		temp_sum+=2*g_mat[i];
	}
	//ģ���һ��
	for (int i=0;i<c;i++)
	    g_mat[i]=g_mat[i]/temp_sum;
	float *temp_row=new float[g_row*g_col];//����ƽ���������
	float *nex_data=new float[g_row*g_col];//�洢��˹ƽ���������

	//����ά��˹����ֽ�Ϊˮƽ����ֱ���������һά������Ƚ���һά�о�������������ٽ����о��
    /*********************************** �з����˹ƽ��***************************************/
	#ifdef WIN32
	int pros_num=omp_get_num_procs();
	omp_set_num_threads(pros_num*3/4);
#pragma omp parallel for
#endif

	for (int i=0;i<g_row;i++)
	{
		//ÿһ�е�ǰdim/2������û��ƽ��
		for (int m=0;m<r;m++)	
			temp_row[i*g_col+m]=pre_data[i*g_col+m];

		//�����м����ݣ�ֱ�Ӵ������ÿ���Խ������
		for (int j=r;j<g_col-r;j++)
		{
			float sum=0;//�洢ģ��ƽ������
			sum+=pre_data[i*g_col+j]*g_mat[0];//ģ���м����ݵ�������
			for (int k=1;k<c;k++)           //���Ϊ7��ģ�壨-3,-2,-1,0,1,2,3��
			{
				sum+=pre_data[i*g_col+j-k]*g_mat[k];
				sum+=pre_data[i*g_col+j+k]*g_mat[k];
			}
			temp_row[i*g_col+j]=sum;//�洢��ƽ���������
		}
		//��1�ַ�����ÿһ�еĺ�dim/2������û��ƽ��
		for (int m=g_col-r;m<g_col;m++)
			temp_row[i*g_col+m]=pre_data[i*g_col+m];
	}
	/***********************************�з����˹ƽ�� ***************************************/

	//ǰdim/2��û��ƽ�������ԸĽ�ֻ��һ��ѭ������������
	for (int m=0;m<r;m++)
		for (int j=0;j<g_col;j++)
			nex_data[m*g_col+j]=temp_row[m*g_col+j];//ֱ�Ӹ�ֵ
	
	//�����м����ݣ�ֱ�Ӵ������ÿ���Խ������
	#ifdef WIN32
	omp_set_num_threads(pros_num*3/4);
#pragma omp parallel for
#endif

	for (int i=r;i<g_row-r;i++)
	{
		for (int j=0;j<g_col;j++)
		{
			float sum=0;//�洢ģ��ƽ������
			sum+=temp_row[i*g_col+j]*g_mat[0];//ģ���м����ݵ�������
			for (int k=1;k<c;k++)
			{
				sum+=temp_row[(i-k)*g_col+j]*g_mat[k];
				sum+=temp_row[(i+k)*g_col+j]*g_mat[k];
			}
			nex_data[i*g_col+j]=sum;//�洢��ƽ�����е�����
		}
	}
	//��1�ַ�������dim/2��û��ƽ��
	for (int m=g_row-r;m<g_row;m++)
		for (int j=0;j<g_col;j++)
			nex_data[m*g_col+j]=temp_row[m*g_col+j];//ֱ�Ӹ�ֵ
	
	delete []g_mat;
	delete []temp_row;

	return nex_data;
}
/****************************************************** 
* �������ƣ� 
*  Build_Dog_Space()
*  
* ����������
* SCALESPACE **ss-��˹�ռ�
* int octaves-��˹�ռ������
* int intvls-��˹�ռ�Ĳ���
*
*
* ����ֵ��
* 0Ϊʧ�ܣ�1Ϊ�ɹ�
*  
* ˵�������ɲ�ָ�˹�ռ�
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
* �������ƣ� 
*  DetectKeypoint()
*  
* ����������
* SCALESPACE **ss - �߶ȿռ�
* int octaves     - �߶ȿռ�����
* int intvls      - ÿ��߶ȿռ�ĸ�˹����������-3��
* float contr_thr - �Աȶ���ֵ
* float curv_thr  - ��Ե������ֵ
* int &keynum     - �ؼ���ĸ���
*
* ����ֵ��
* ��⵽�Ĺؼ��㣬�޳��˱�Ե��͵ͶԱȶȵĵ�
*  
* ˵������ֵ����
********************************************************/
Key_Point* DetectKeypoint(SCALESPACE **ss,int octaves,int intvls,float contr_thr,float curv_thr,int &keynum)
{
	keynum=0;
	Key_Point *keypoints= NULL;
	//����opencv�еĹ�ʽ������ֵ����ÿ��ռ�Ĳ���-3
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
									Key_Point* k=new Key_Point[1];//��������洢�ؼ���
									k->next=keypoints;
									keypoints=k;
									k->key_row=mr; //��
									k->key_column=mc;  //��
									k->initl_row=(mr+xr)*pow(2.0f,o);//��ԭʼӰ���ϵ���λ��
									k->initl_column=(mc+xc)*pow(2.0f,o);//��ԭʼӰ���ϵ���λ��
									k->key_octave=o;//��
									k->key_intvl=ms;//��
									k->sub_intvl=xs;//��ȷ��λ��ĸ�����
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
* �������ƣ� 
* is_extremum()
*  
* ����������
* SCALESPACE ss
* int o �߶ȿռ������
* int s �߶ȿռ�Ĳ���
* int r ��ָ�˹ͼ�����
* int c ��ָ�˹ͼ�����
* int col ��ָ�˹ÿ�е����ظ���
* 
* ����ֵ��
* 0Ϊʧ�ܣ�1Ϊ�ɹ�
*  
* ˵�����ж��Ƿ�Ϊ��ֵ��
********************************************************/
bool is_extremum(SCALESPACE **ss,int o,int s,int r,int c,int col)
{
	float val=ss[o][s].dog_space[r*col+c];

	//�����㣬��⼫��ֵ
	if( val > 0 )
	{
		for(int i = -1; i <= 1; i++ )
			for( int j = -1; j <= 1; j++ )
				for(int k = -1; k <= 1; k++ )
					if( val<ss[o][s+i].dog_space[(r+j)*col+c+k])//���Ǽ���ֵ
						if (!(i==0&&j==0&&k==0)) //�ų����ڼ��ĵ㣬��Ϊ��o,s,r,c���϶�=val
						    return 0;
	}

//С���㣬��⼫Сֵ
	else
	{
		for(int i = -1; i <= 1; i++ )
			for(int j = -1; j <= 1; j++ )
				for(int k = -1; k <= 1; k++ )
					if( val>ss[o][s+i].dog_space[(r+j)*col+c+k] )//���Ǽ�Сֵ
						if (!(i==0&&j==0&&k==0))
							return 0;
	}
	return 1;
}
/****************************************************** 
* �������ƣ� 
* interp_extremum()
*  
* ����������
* SCALESPACE ss
* int o �߶ȿռ������
* int &ms �߶ȿռ�Ĳ���,�����ؾ�ȷ��λ��߶ȿռ�Ĳ���
* int &mr ��ָ�˹ͼ����У������ؾ�ȷ��λ���ָ�˹ͼ�����
* int &mc ��ָ�˹ͼ����У������ؾ�ȷ��λ���ָ�˹ͼ�����
* int intvls ÿ���˹����-3
* int col ��ָ�˹ͼ���������
* int contr_thr �Աȶ���ֵ
* float &xs �߶Ȳ�ĸ�����
* float &xr �з��������
* float &xc �з��������
*
* ����ֵ��
* 0Ϊʧ�ܣ�1Ϊ�ɹ�
*  
* ˵������ȷ��λ��ֵ��
********************************************************/
bool interp_extremum(SCALESPACE **ss,int o,int &ms,int &mr,int &mc,int intvls,int col,float contr_thr,float &xs,float &xr,float &xc)
{
	int i=0;
	while(i<SIFT_MAX_INTERP_STEPS)
	{
		if (!interp_step( ss, o, ms, mr, mc,col, xs, xr, xc))//��ȷ��ֵ�Ƿ�ɹ�
		{
			return 0;
		}
		
		if( fabs(xs)<0.5 && fabs(xr)<0.5 && fabs(xc)<0.5)
			break;

		mc += floor( xc+0.5f );
		mr += floor( xr+0.5f );
		ms += floor( xs+0.5f );
		if( ms<1||ms>intvls||//�߶ȳ�����Χ   
			mc < SIFT_IMG_BORDER||mr < SIFT_IMG_BORDER ||//�л��г����߽�
			mr >= ss[o][ms].gus_row - SIFT_IMG_BORDER ||mc >= ss[o][ms].gus_column - SIFT_IMG_BORDER )
		{
			return 0;
		}
		i++;
	}
	if( i >= SIFT_MAX_INTERP_STEPS )
		return 0;//ѭ�������Ĵ�����������ֵ

	//�����ֵ���dog�ռ��еĶԱȶȵ�ֵ������Lowe�����й�ʽ��3�����빫ʽ��2������õ�ֵ
	float contr;
	contr=interp_contr(ss,o,ms,mr,mc,col,xs,xr,xc);
	if( fabs(contr)<contr_thr/intvls)//�ؼ���ĻҶ��Ƿ������ֵ
		return 0;
	return 1;
}

/****************************************************** 
* �������ƣ� 
* interp_contr()
*  
* ����������
* SCALESPACE ss �߶ȿռ�
* int o �߶ȿռ������
* int s �߶ȿռ�Ĳ���
* int r ��ָ�˹ͼ�����
* int c ��ָ�˹ͼ�����
* int col �߶ȿռ�������
* float xs ��ĸ�����
* float xr �еĸ�����
* float xc �еĸ�����
*
* ����ֵ��
* float contr ��ȷ��λ��ֵ���ؼ����Ӧ�ĻҶ�ֵ
*  
* ˵������ȷ��λ��ֵ���ؼ����Ӧ�ĻҶ�ֵ
********************************************************/
float interp_contr(SCALESPACE **ss,int o,int s,int r,int c,int col,float xs,float xr,float xc)
{
	float Dx,Dy,Ds,contr;
	/*һ��ƫ������
	| Dx  Dy  Ds |
	*/
	Dx= (ss[o][s].dog_space[r*col+c+1]-ss[o][s].dog_space[r*col+c-1])*0.5f;
	Dy= (ss[o][s].dog_space[(r-1)*col+c]-ss[o][s].dog_space[(r+1)*col+c])*0.5f;
	Ds= (ss[o][s+1].dog_space[r*col+c]-ss[o][s-1].dog_space[r*col+c])*0.5f;
	contr=ss[o][s].dog_space[r*col+c]+(Dx*xc+Dy*xr+Ds*xs)*0.5f;//��ȷ��λ�ؼ���󣬹ؼ�������Ӧ�ĻҶ�ֵ
	return contr;
}
/****************************************************** 
* �������ƣ� 
* is_too_edge_like()
*  
* ����������
* SCALESPACE** ss �߶ȿռ�
* int o �߶ȿռ������
* int s �߶ȿռ�Ĳ���
* int r ��ָ�˹ͼ�����
* int c ��ָ�˹ͼ�����
* int col ��ָ�˹ͼ���������
* float curv_thr �޳���Ե�����ֵ
*
* ����ֵ��
* 0Ϊʧ�ܣ�1Ϊ�ɹ�
*  
* ˵�����޳���Ե��
********************************************************/
bool is_too_edge_like(SCALESPACE** ss, int o, int s,int r,int c,int col,float curv_thr)
{
	//����ؼ�������ʣ��Ӷ�ȥ����Ե��
	float val, dxx, dyy, dxy, tr, det;
	val=ss[o][s].dog_space[r*col+c];
	dxx=ss[o][s].dog_space[r*col+c-1] + ss[o][s].dog_space[r*col+c+1]-2*val;
	dyy = ss[o][s].dog_space[(r-1)*col+c] + ss[o][s].dog_space[(r+1)*col+c]-2* val;
	dxy = ( ss[o][s].dog_space[(r+1)*col+c-1]+ss[o][s].dog_space[(r-1)*col+c+1]-ss[o][s].dog_space[(r+1)*col+c+1]-ss[o][s].dog_space[(r-1)*col+c-1])*0.25f;
	tr = dxx + dyy;//����ļ�
	det = dxx*dyy - dxy*dxy;//���������ʽ
	//���׵�������ʽ��ֵС���㲻���ڼ�ֵ
	if( det <= 0 )
		return 1;
	// ������ʣ�����Lowe�����й�ʽ��4������curv_thr = (tr*tr)/det;
	if(tr*tr/det>=(curv_thr+1.0f)*(curv_thr+1.0f)/curv_thr)
	{
		return 1;
	}
	return 0;
}

/****************************************************** 
* �������ƣ� 
*  inter_extremum()
*  
* ����������
* SCALESPACE **ss �߶ȿռ�
* int o     �߶ȿռ������
* float s   �߶ȿռ�Ĳ���
* float r   ��ָ�˹ͼ�������
* float c   ��ָ�˹ͼ�������
* int col   ��ָ�˹ͼ���е�����
* float &xs ��ĸ�����
* float &xr �еĸ�����
* float &xc �еĸ�����
*
* ����ֵ��
* 0Ϊʧ�ܣ�1Ϊ�ɹ�
*  
* ˵������ȷ��λ��ֵ��ʱ������x��y��s�ĸ���ֵ
********************************************************/
bool interp_step(SCALESPACE **ss,int o, int s, int r, int c,int col,float &xs,float &xr,float &xc)
{
	float Dx,Dy,Ds,Dxx,Dyy,Dss,Dxy,Dxs,Dys;
	
	double** a=new double*[3];//���ڴ洢����ƫ��
	double** inv_a=new double*[3];//���ڴ洢����ƫ���������
	for (int i=0;i<3;i++)
	{
		a[i]=new double[3];
		inv_a[i]=new double[3];
	}
	float val=ss[o][s].dog_space[r*col+c];
	/*һ��ƫ������
	| Dx  Dy  Ds |
	*/
	Dx= (ss[o][s].dog_space[r*col+c+1]-ss[o][s].dog_space[r*col+c-1])*0.5f;
	Dy= (ss[o][s].dog_space[(r-1)*col+c]-ss[o][s].dog_space[(r+1)*col+c])*0.5f;
	Ds= (ss[o][s+1].dog_space[r*col+c]-ss[o][s-1].dog_space[r*col+c])*0.5f;

	/*����ƫ������
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

	//�����ƫ���������
	if (sift_Inv(a,inv_a,3))
	{
		// ����x��y��s�ĸ���ֵ
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
* �������ƣ� 
* calc_scales()
*  
* ����������
* Key_Point *keypoint �ؼ���
* int keynumber �ؼ���ĸ���
*
* ����ֵ��
* 0Ϊʧ�ܣ�1Ϊ�ɹ�
*  
* ˵������ؼ���ĳ߶�
********************************************************/
void calc_scales(Key_Point *keyhead,int keynumber)
{
	int i;
	double intvl;
	Key_Point *keypoint=keyhead->next;
	for( i = 0; i < keynumber; i++ )
	{
		intvl = keypoint->key_intvl + keypoint->sub_intvl;
		keypoint->scl =SIFT_SIGMA * pow( 2.0, keypoint->key_octave + intvl / SCALESPEROCTAVE );//keypoint->sclֻΪ��ʾʹ��
		keypoint->scl_octave = SIFT_SIGMA * pow( 2.0, intvl / SCALESPEROCTAVE );
		keypoint=keypoint->next;
	}
	delete []keypoint;
}
/****************************************************** 
* �������ƣ� 
* calc_oritation()
*  
* ����������
* Key_Point* keypoint �ؼ���
* SCALESPACE** ss �߶ȿռ�
* int &keynumber �ؼ���ĸ���
* 
* ����ֵ��
* keypoints �ؼ���
*  
* ˵��������������
********************************************************/
Key_Point* calc_oritation(Key_Point* keyhead,SCALESPACE** ss,int &keynumber)
{
	float *hist;
	float omax,bin,/*mag_peak,*/PI2=PI*2.0f;//mag_peak��ʾ��ֵ����������Ӧ���ݶ�ֵ
	int addnum=0, n=SIFT_ORI_HIST_BINS;   //addnum�����ڶ�������������ؼ���������Ӻ��ܵ���
	Key_Point *keypoints=NULL;
	Key_Point *k=keyhead->next;
	for(int i = 0; i < keynumber; i++ )
	{ 
		hist = ori_hist(ss[k->key_octave][k->key_intvl].gauss_space,k->key_row,k->key_column,ss[k->key_octave][k->key_intvl].gus_row,
			ss[k->key_octave][k->key_intvl].gus_column,SIFT_ORI_HIST_BINS,int(floor(SIFT_ORI_RADIUS*k->scl_octave+0.5f)),SIFT_ORI_SIG_FCTR*k->scl_octave);
		for(int j=0;j<SIFT_ORI_SMOOTH_PASSES;j++)
			smooth_ori_hist(hist,n);
		omax = dominant_ori(hist,n);

		//����ֱ��ͼ�е����ֵ�������ֱ��ͼ�����е������򼴴���0.8�������ֵ�����з���
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
				Key_Point* kk=new Key_Point[1];//�洢�ؼ���
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
				kk[0].ori=((PI2 * bin )/n)/*-PI*/;//ori��ȡֵ��Χ��0,2*PI��
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
* �������ƣ� 
* ori_hist()
*  
* ����������
* float* guss  ��˹�ռ�����
* int r, int c �кź��к�
* int col      ������
* int row      ������
* int n        ����ֱ��ͼ���õ��ܷ�������36��
* int rad      ��������ݶ�ֵ
* floa sigma   ��˹Ȩ�ز��õı�׼�����1.5���Ĵ��ڿ�ȣ�      
*
* ����ֵ��
* float* hist  �������ֱ��ͼ
*  
* ˵����������ķ���ֱ��ͼ
********************************************************/
float* ori_hist(float* guss, int r, int c, int row,int col, int n, int rad, float sigma)
{
	
	float mag=0,ori=0,w,PI2=PI*2.0f;//mag�ݶ� ori������ w��˹Ȩ�� 
	int bin;//����ֵ
	float sigma2 = 2.0f * sigma * sigma;
	float* hist=new float[n];//36�������ͳ��ֱ��ͼ
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
* �������ƣ� 
* calc_grad_mag_ori()
*  
* float* guss  ��˹�ռ�����
* int r, int c �кź��к�
* int col      ������
* int row      ������
* float &mag   ��������ݶ�ֵ
* float &ori   ������ķ���ֵ 
* 
* ����ֵ��
* 0Ϊʧ�ܣ�1Ϊ�ɹ�
*  
* ˵�����������ݶȺͷ���
********************************************************/

bool calc_grad_mag_ori(float* guss, int r, int c,int row,int col, float &mag, float &ori )
{
	float dx, dy;
	if(r>0 && r<row-1 && c>0 && c<col-1)
	{
		dx=guss[r*col+c+1]-guss[r*col+c-1];     //dx
		dy=guss[(r-1)*col+c]-guss[(r+1)*col+c]; //dy
		mag = sqrt( dx*dx + dy*dy );
		ori = atan2(dy, dx);      //���������ȡֵ��ΧΪ��-PI��PI��
		ori=(ori>0)?ori:(ori+PI+PI);
		return 1;
	}
	else
		return 0;
}
/****************************************************** 
* �������ƣ� 
* smooth_ori_hist()
*  
* float* hist  ����ֱ��ͼ
* int n        ����ֱ��ͼ���õ��ܷ�������36��
* 
* ����ֵ��
* 
*  
* ˵������ֱ��ͼ����ƽ����ģ��Ϊ��0.25 0.5 0.25��
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
* �������ƣ� 
* dominant_ori()
*  
* float* hist ����ֱ��ͼ
* int n       ����ֱ��ͼ���õ��ܷ�������36��
* 
* ����ֵ��
* ֱ��ͼ�е����ֵ
*  
* ˵��������ֱ��ͼ�е����ֵ
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
* �������ƣ� 
* interp_hist_peak()
*  
* ����������
* float lef ��������ߵ��ݶ�ֵ
* float cent  ��������ݶ�ֵ
* float rigt  �������ұߵ��ݶ�ֵ
* float &mag_peak ������������߲�ֵ����ݶ�ֵ
*
* ����ֵ��
* �����߲�ֵ���������ֵ
*  
* ˵��������������������߲�ֵ
********************************************************/
float interp_hist_peak(float lef,float cent,float rigt/*,float &mag_peak*/)
{
	//mag_peak=cent-0.125f*(lef-rigt)*(lef-rigt)/(lef-2.0f*cent+rigt);
	return 0.5f*(lef-rigt)/(lef-2.0f*cent+rigt);
}
/****************************************************** 
* �������ƣ� 
* compute_descriptors()
*  
* ����������
* Key_Point* feat ������
* SCALESPACE **ss �߶ȿռ�
* int keynumber �ؼ���ĸ���
* 
* ����ֵ��
*  
* ˵���������������128ά������������
********************************************************/
void compute_descriptors( Key_Point* feat, SCALESPACE **ss, int keynum)
{
	int d=SIFT_DESCR_WIDTH,n=SIFT_DESCR_HIST_BINS;
	for(int i=0;i<keynum;i++)
	{		
		int num=SIFT_DESCR_BINS;//�����ӵ�ά�� 128
		//����ֱ��ͼ
		float *hist= descr_hist(ss[feat[i].key_octave][feat[i].key_intvl].gauss_space,feat[i].key_row,
			feat[i].key_column,ss[feat[i].key_octave][feat[i].key_intvl].gus_row,
			ss[feat[i].key_octave][feat[i].key_intvl].gus_column,feat[i].ori,feat[i].scl_octave,d,n,num);
		//��һ��

		normalize_descr(hist,num);
		//Ϊ�˼�С�ϴ��ݶ�ֵ������������Ӱ�죬����һ������ݶ�ֵ������SIFT_DESCR_MAG_THR=0.2�ڡ�
		for(int j = 0; j < num; j++ )
		{
			if(hist[j]>SIFT_DESCR_MAG_THR)
				hist[j]=SIFT_DESCR_MAG_THR;
		}
		//�ٴι�һ��
		normalize_descr(hist,num);
		//������������ϵ��512���Ӷ�������������Ϊ����
		for(int j=0; j<num; j++)
		{
			feat[i].descriptor[j] =float( min(255, int(SIFT_INT_DESCR_FCTR * hist[j])));
		}
		delete []hist;
	}
}
/****************************************************** 
* �������ƣ� 
* descr_hist()
*  
* ����������
*
* float* guss  ��˹�ռ�����
* int r, int c �кź��к�
* int row      ������
* int col      ������
* float &ori   ������ķ���ֵ
* float scl    ֻ�����ؼ������õĳ߶�
* int d        �����Ӳ��õ��������  4
* int n        ������ֱ��ͼ������  8
* int num      �����ӵ�ά��
* ����ֵ��
* float *hist 128ά��ֱ��ͼ����
*  
* ˵���������������ֱ��ͼ����
********************************************************/
float* descr_hist(float* guss,int r,int c,int row,int col,float ori,float scl,int d,int n,int num)
{
	
	float cos_t,sin_t, hist_width, exp_denom, r_rot, c_rot, grad_mag,
		grad_ori, w, rbin, cbin, obin, bins_per_rad, PI2=2.0f*PI;
	int radius;
	//���岢��ʼ��������ֱ��ͼ
	float* hist=new float[num];
	memset(hist,0,num*sizeof(float));

	cos_t=cos(ori);
	sin_t=sin(ori);
	bins_per_rad=n/PI2;//���ڻ��Ⱥ�8������֮��ת��
	exp_denom=d*d*0.5f;//���ڼ����˹Ȩ��
	hist_width=SIFT_DESCR_SCL_FCTR*scl;//4��4����ÿ��ı߳�	
	radius=int(hist_width*sqrt(2.0f)*(d+1.0f)*0.5f+0.5f);//ͳ�������ӵ�����뾶
	for(int i = -radius; i <= radius; i++ )
		for(int j = -radius; j <= radius; j++ )
		{
			c_rot=(j*cos_t-i*sin_t)/hist_width;//�����ת������������Ϊԭ��Ŀ飨4X4������ϵ�е��к�
			r_rot=(j*sin_t+i*cos_t)/hist_width;
			rbin=r_rot+d/2-0.5f;
			cbin=c_rot+d/2-0.5f;//�����ת���������Ͻ�Ϊԭ��Ŀ飨4X4������ϵ�е��к�

			if(rbin>-1.0f && rbin<d && cbin>-1.0f && cbin<d)//�ж���ת��ĵ���ͳ��������
				if(calc_grad_mag_ori(guss, r + i, c + j, row, col, grad_mag, grad_ori))
				{
					grad_ori -= ori;  //�����һ��
					while( grad_ori < 0 )//��Ϊgrad-ori��ori��ȡֵ��Χ���ǣ�-PI��PI����������Ҫ������Ľ�����е���
						grad_ori += PI2;
					while( grad_ori >= PI2 )
						grad_ori -= PI2;

					obin = grad_ori * bins_per_rad;//���������ڵ�ֱ��ͼ����
					w = exp(-(c_rot*c_rot+r_rot*r_rot)/exp_denom);//�Կ�ĳ���Ϊ������λ�������˹Ȩ��
					interp_hist_entry(hist, rbin, cbin, obin, grad_mag*w, d, n);
				}
		}
	return hist;
}
/****************************************************** 
* �������ƣ� 
* interp_hist_entry()
*  
* ����������
*
* float* hist  128ά��ֱ��ͼ
* float rbin   ����4���п��е��к�
* float cbin   ����4���п��е��к�
* float obin   ��ķ�����8������ֵ�е����
* float mag    ������ķ���ֵ
* int d        �����Ӳ��õ��������  4
* int n        ������ֱ��ͼ������  8
* 
* ����ֵ��
* 
*  
* ˵�������������������ڵ�ÿ����㣬�������ݶ�ֵ��Ȩ�أ�1-d��������������������
���У���rbin=2.3��������ڵڶ��к͵������й��ף����кͷ�����ݶ�ֵ�Ĺ���
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
						hist[(rb*d+cb)*n+ob]+= v_o;//������ֵ���뵽128άֱ��ͼ��
					}
				}
			}
		}
	}
}
/****************************************************** 
* �������ƣ� 
* normalize_descr()
*  
* ����������
*
* float* hist  ������ֱ��ͼ
* int d        ������ά�� 128
* 
* ����ֵ��
* 
*  
* ˵��������128ά�����ӽ��й�һ��
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
* �������ƣ� 
* adjust_for_img_down()
*  
* ����������
* Key_Point* feat - �ؼ���
* int keynumber   - �ؼ���ĸ���
*
* ����ֵ��
* 0Ϊʧ�ܣ�1Ϊ�ɹ�
*  
* ˵�����Գ�ʼӰ����н�����ʱ����Ҫ�Թؼ����������е���
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
* �������ƣ� 
* adjust_for_img_grid()
*  
* ����������
* Key_Point* feat - �ӿ�ؼ���
* int keynumber   - �ӿ�ؼ���ĸ���
* int col         - �ӿ���к�
* int row         - �ӿ���к�
* int gridwidth   - �ӿ��
* int imgheight   - ����Ӱ��ĸ�
*
* ����ֵ��
* 0Ϊʧ�ܣ�1Ϊ�ɹ�
*  
* ˵�����ֿ��ͼ����Ҫ�Թؼ����������е������ָ���ԭʼͼ���ϵ�����
********************************************************/
void adjust_for_img_grid(Key_Point* feat,int keynumber,int col,int row,int gridwidth,int imgheight)
{
	for(int i = 0; i < keynumber; i++ )
	{
		feat[i].initl_row =row*gridwidth+feat[i].initl_row ;//�ؼ����������initl_row���������½�Ϊ����ԭ��
		//feat[i].initl_row =imgheight-(row*gridwidth+feat[i].initl_row) ;//�ؼ����������initl_row�������Ͻ�Ϊ����ԭ��
		feat[i].initl_column = col*gridwidth+feat[i].initl_column ;
	}
}
/****************************************************** 
* �������ƣ� 
* Inv()
*  
* ����������
* double** a         - ������ľ���
* double** inva      - �����ľ���
* int dim            - �����ά��
*
* ����ֵ��
* ����������
*  
*˵������������
********************************************************/

int sift_Inv(double** a,double** inva,int dim )
{
	int i,j,k,p;//��¼����Ԫ������p����ֵc_max
	double temp,c_max=0;
	//��λ����ֵ
	double **b=new double *[dim];
	for (i=0;i<dim;i++)
	{
		b[i]=new double [dim];
		memset(b[i],0,sizeof(double)*dim);
		b[i][i]=1;
	}

	//��˹-Լ����Ԫ
	for(i=0;i<dim;i++)
	{  
		/*******************ѡ����Ԫ**********************/
		for(j=i+1,c_max=a[i][i],p=i;j<dim;j++)//c_max=a[i][i]���ܷ���ѭ�����ڣ�����û�������һ�е����ֵ
		{  
			if(fabs(a[j][i])>fabs(c_max))  
			{p=j; c_max=a[j][i]; }     
		}  
		//������󲻴��������
		if(fabs(c_max)<1e-10) 
		{return 0;}
		/*****************���С���һ*********************/
		if (i!=p)
		{
			for(j=0;j<dim;j++)
			{   
				temp=a[i][j];
				a[i][j]=a[p][j]/c_max; 
				a[p][j]=temp;
				//��λ����
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
		/**********************��Ԫ**********************/
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
	//��������ֵ��ֵ�����뺯������
	for (i=0;i<dim;i++)
	{
		for (j=0;j<dim;j++)
		{
			inva[i][j]=b[i][j];
		}
	}
	//�ͷ��ڴ�
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


//ɾ������
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