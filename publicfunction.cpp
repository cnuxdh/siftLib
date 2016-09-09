//#include "stdafx.h"
#include "string.h"
#include "publicfunction.h"
#include <math.h>
//#include "CXimageHeader\ximadef.h"

/****************************************************** 
* �������ƣ� 
* BubbleSort()
*  
* ����������
* double* Array       - �����������
* int n               - ����
* int &minnum         - ��Сֵ���ڣ�����ǰ�ģ�λ��
* 
* ����ֵ��
*
* ��
*  
*˵����ð�ݷ�����(��С����)
********************************************************/
void BubbleSort(double* Array,int n,int &minnum)
{
	double mindata;
	int* minn=new int[n-1];
	memset(minn,0,sizeof(int)*(n-1));

	for(int i=0;i<n-1;i++) 
	{ 
		mindata=Array[i];
		minn[i]=i;
		for(int j=i+1;j<n;j++) 
		{ 
			if(mindata>Array[j])//�ȽϽ�������Ԫ�� 
			{ 
				mindata=Array[j];minn[i]=j;
			} 
		} 
		double temp; 
		temp=Array[minn[i]]; Array[minn[i]]=Array[i]; Array[i]=temp; 
	} 
	minnum=minn[0];
	delete []minn;
}
/****************************************************** 
* �������ƣ� 
* Find_Mid()
*  
* ����������
* double* array    - ��������ֵ������
* int nunber       - �ܸ���
* 
* ����ֵ��
* �������ֵ
*  
*˵��������һ��������ֵ
********************************************************/
double Find_Mid(double* Array,int number)
{
	if (number<=1)
	{
		return Array[0];
	}
	//ԭ���ݵ�˳���ܸı䣬��˿����µ�һ���������
	double* sort_array=new double[number];
	memcpy(sort_array,Array,sizeof(double)*number);
	//����	
	int minnum=0;
	BubbleSort(sort_array,number,minnum);
	//������ֵ	
	double mid;
	if (number%2==0)
	{
		mid=(sort_array[number/2]+sort_array[number/2+1])/2.0;
	} 
	else
	{	
		mid=sort_array[number/2];
	}
	delete []sort_array;
	return mid;
}


/****************************************************** 
* �������ƣ� 
* Biliner()
*  
* ����������
* unsigned char data11        - ԭʼͼ�������ֵ(����)
* unsigned char data12        - ԭʼͼ�������ֵ(����)
* unsigned char data21        - ԭʼͼ�������ֵ(����)
* unsigned char data22        - ԭʼͼ�������ֵ(����)
* double oldw                 - ����õ���ԭʼͼ����������
* double oldh                 - ����õ���ԭʼͼ����������
* 
* ����ֵ��
* unsigned char ˫���Բ�ֵ�������ֵ
*  
*˵����˫���Բ�ֵ
********************************************************/
unsigned char Bilinear( unsigned char data11, unsigned char data12,unsigned char data21,unsigned char data22,double oldw,double oldh)
{
	double detaw=oldw -(int)oldw;  //��y
	double detah=oldh -(int)oldh;  //��x
	//��ֵ��ʽ
	int data=((data22-data21-data12+data11)*detah*detaw+(data12-data11)*detaw+(data21-data11)*detah+data11)+0.5;
	//�ų�С��0�ʹ���255�����
	data=(data>255)?255:data;
	data=(data<0)?0:data;
	return (unsigned char)data;
}
/****************************************************** 
* �������ƣ� 
* Matrix_trans()
*  
* ����������
* double** matrix       - ����ת�õľ���
* int row               - ԭ���������
* int col               - ԭ���������
*
* ����ֵ��
* ת�þ���
*  
*˵���������ת��
********************************************************/

double** Matrix_trans(double** matrix,int row,int col)
{
	double** matrix_t=new double*[col];
	for (int i=0;i<col;i++)
	{
		matrix_t[i]=new double[row];
	}

	for (int i=0;i<col;i++)
	{
		for (int j=0;j<row;j++)
		{
			matrix_t[i][j]=matrix[j][i];
		}
	}
	return matrix_t;
}
/****************************************************** 
* �������ƣ� 
* Matrix_mul()
*  
* ����������
* double** left_matrix       - �����
* double** right_matrix      - �Ҿ���
* int left_row               - �������
* int left_col               - ������У��Ҿ�����У�
* int right_col              - �Ҿ������
*
* ����ֵ��
* ����ĳ˻�
*  
*˵�����������
********************************************************/

double** Matrix_mul(double** left_matrix,double** right_matrix,int left_row,int left_col,int right_col)
{
	double** mul_matrix=new double*[left_row];
	for (int i=0;i<left_row;i++)
	{
		mul_matrix[i]=new double[right_col];
	}

	for (int i=0;i<left_row;i++)
	{
		for (int j=0;j<right_col;j++)
		{
			double temp_sum=0.0;
			for (int k=0;k<left_col;k++)
			{
				temp_sum+=left_matrix[i][k]*right_matrix[k][j];
			}
			mul_matrix[i][j]=temp_sum;
		}
	}
	return mul_matrix;
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

int Inv(double** a,double** inva,int dim )
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

/************************************************************************
* �������ƣ� 
* Least_Square()
*  
* ����������
* double** A               - ϵ������
* double** L               - ������
* int r                    - ����A������
* int c                    - ����A������
*
* ����ֵ��
* ���������ֵ
*  
*˵������С����    ��ʽV=AX-L
*      
************************************************************************/
double** Least_Square(double** A,double** L,int &r,int c)
{
	double** X;     //������С��������Ľ��
	double sigma=0; //����
	double pre_sigma=0;
	//���
	double** V=new double*[r];
	for (int j=0;j<r;j++)
	{
		V[j]=new double[1];
		V[j][0]=0;
	}
	//ÿ�ε����޳������A
	double** tempA=new double*[r];
	for (int j=0;j<r;j++)
	{
		tempA[j]=new double[c];
		memset(tempA[j],0,sizeof(double)*c);
	}	
	//ÿ�ε����޳������L(c*1)c��1��
	double** tempL=new double*[r];
	for (int j=0;j<r;j++)
	{
		tempL[j]=new double[1];
		tempL[j][0]=0;
	}

	////////////////////////////////////////////////////////////////////
	do
	{
		pre_sigma=sigma;
		//��ϵ�������ת��
		double** AT = Matrix_trans(A,r,c);

		//ϵ����ת�ó���ϵ��
		double** ATA=Matrix_mul(AT,A,c,r,c);

		//��˽������
		double** inv_ATA=new double*[c];
		for (int j=0;j<c;j++)
		{
			inv_ATA[j]=new double[c];
		}
		if (Inv(ATA,inv_ATA,c))
		{
			//������ٳ���ϵ����ת��
			double** mul_xym=Matrix_mul(inv_ATA,AT,c,c,r);
			
			//�����Գ�����
			X=Matrix_mul(mul_xym,L,c,r,1);

			//�����
			double** AX=Matrix_mul(A,X,r,c,1);	
			for (int j=0;j<r;j++)
			{
				V[j][0]=AX[j][0]-L[j][0];
			}

			//��������
			double** VT=Matrix_trans(V,r,1);
			double** VTV=Matrix_mul(VT,V,1,r,1);
			sigma=sqrt(VTV[0][0]/(r-c));

			int tempr=0;//�޳������ĸ���
			for (int j=0;j<r;j+=2)
			{
				if ( abs(V[j][0])<3*sigma && abs(V[j+1][0])<3*sigma )
				{			
					memcpy(tempA[tempr],A[j],sizeof(double)*c);
					tempL[tempr++]=L[j];
					memcpy(tempA[tempr],A[j+1],sizeof(double)*c);
					tempL[tempr++]=L[j+1];
				}
			}
			//���޳��ֲ��ĵ㸳ֵ��A������ѭ��
			for (int j=0;j<tempr;j++)
			{
				memcpy(A[j],tempA[j],sizeof(double)*c);
				L[j][0]=tempL[j][0];
			}
			//�ͷ��ڴ�
			for (int i=0;i<c;i++)
			{
				delete[]mul_xym[i];
				delete[]inv_ATA[i];
				delete[]ATA[i];
				delete[]AT[i];
			}
			delete[]mul_xym;
			delete[]inv_ATA;
			delete[]ATA;
			delete[]AT;

			for (int i=0;i<r;i++)
			{
				delete[]AX[i];
			}
			delete[]AX;
			delete[]VT[0];
			delete[]VT;
			delete[]VTV[0];
			delete[]VTV;
			
			//�ı��ܵ���
			r=tempr;
		}
		else
		{
			//�ͷ��ڴ�
			for (int i=0;i<c;i++)
			{
				delete[]inv_ATA[i];
				delete[]ATA[i];
				delete[]AT[i];
			}
			delete[]inv_ATA;
			delete[]ATA;
			delete []AT;

			for (int i=0;i<r;i++)
			{
				delete[]V[i];
			}
			delete[]V;

			return NULL;
		}
	} while ( abs(sigma-pre_sigma)>0.01 && r>=c);

	for (int i=0;i<r;i++)
	{
		delete[]V[i];
	}
	delete[]V;
	//�ͷ��ڴ�
	return X;
}

/****************************************************** 
* �������ƣ� 
* GetGradData()
*  
* ����������
* unsigned char* imgdata - ��ȡ��ͼ������
* int colum     - ��ȡ��ͼ����������
* int row       - ��ȡ��ͼ����������
* int BitCount  - ��ȡ��ͼ��Ĳ�����
*
* ����ֵ��
* ������ȡ�������ͼ������
*  
*˵������ȡͼ��ĻҶ�����
********************************************************/
unsigned char *GetGrayData( unsigned char* imgdata,int colum,int row,int BitCount ) 
{	
	unsigned char* imgraydata=new unsigned char[colum*row];
	
	if (BitCount==8)
	{
		int perline=(colum+3)/4*4;
		for (int i=0;i<row;i++)
			for (int j=0;j<colum;j++)
				imgraydata[i*colum+j]=imgdata[i*perline+j];
	} 
	else if(BitCount==24)
	{
		int perline=(colum*3+3)/4*4;
		for (int i=0;i<row;i++)
			for (int j=0;j<colum;j++)
				imgraydata[i*colum+j]=unsigned char(0.114*imgdata[i*perline+3*j]+
				0.587*imgdata[i*perline+3*j+1]+	0.299*imgdata[i*perline+3*j+2]);
	}
	return imgraydata;
}

/****************************************************** 
* �������ƣ� 
* ImgNomalize()
*  
* ����������
* unsigned char* imgraydata - ��ȡ��ͼ������
* int colum     - ��ȡ��ͼ����������
* int row       - ��ȡ��ͼ����������
*
* ����ֵ��
* ������ȡ�������ͼ������
*  
*˵����ͼ��ĻҶ�ֱֵ�ӳ���255��һ��
********************************************************/
float *ImgNomalize( unsigned char* imgraydata,int colum,int row ) 
{	
	//float temp=1.0f/255;
	float temp = 1.0;
	float* imagedata=new float[colum*row];
	for (int i=0;i<colum*row;i++)
		imagedata[i]=float(imgraydata[i])*temp;
	
	return imagedata;
}
/****************************************************** 
* �������ƣ� 
* descr_dist_sq()
*  
* ����������
* Key_Point* f1 - ��һ������
* Key_Point* f2 - �ڶ�������
* 
*
* ����ֵ��
* ��ֵ�������е����
*  
*˵�������������������ŷʽ����
********************************************************/
float descr_dist_sq( Key_Point* f1, Key_Point* f2, int nDim)
{
	float diff, dsq = 0,d= nDim/*SIFT_DESCR_BINS*/;
	float* descr1, * descr2;
	/*if( f2->d != d )
		return DBL_MAX;*/
	descr1 = f1->descriptor;
	descr2 = f2->descriptor;

	for(int i = 0; i < d; i++ )
	{
		diff = descr1[i] - descr2[i];
		dsq += diff*diff;
	}
	return dsq;
}