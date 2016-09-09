//#include "stdafx.h"
#include "string.h"
#include "publicfunction.h"
#include <math.h>
//#include "CXimageHeader\ximadef.h"

/****************************************************** 
* 函数名称： 
* BubbleSort()
*  
* 函数参数：
* double* Array       - 待排序的数组
* int n               - 总数
* int &minnum         - 最小值所在（排序前的）位置
* 
* 返回值：
*
* 无
*  
*说明：冒泡法排序(由小到大)
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
			if(mindata>Array[j])//比较交换相邻元素 
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
* 函数名称： 
* Find_Mid()
*  
* 函数参数：
* double* array    - 待查找中值的数组
* int nunber       - 总个数
* 
* 返回值：
* 数组的中值
*  
*说明：查找一组数的中值
********************************************************/
double Find_Mid(double* Array,int number)
{
	if (number<=1)
	{
		return Array[0];
	}
	//原数据的顺序不能改变，因此拷贝新的一组进行排序
	double* sort_array=new double[number];
	memcpy(sort_array,Array,sizeof(double)*number);
	//排序	
	int minnum=0;
	BubbleSort(sort_array,number,minnum);
	//返回中值	
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
* 函数名称： 
* Biliner()
*  
* 函数参数：
* unsigned char data11        - 原始图像的像素值(左下)
* unsigned char data12        - 原始图像的像素值(右下)
* unsigned char data21        - 原始图像的像素值(左上)
* unsigned char data22        - 原始图像的像素值(右上)
* double oldw                 - 计算得到的原始图像点的列坐标
* double oldh                 - 计算得到的原始图像点的行坐标
* 
* 返回值：
* unsigned char 双线性插值后的像素值
*  
*说明：双线性插值
********************************************************/
unsigned char Bilinear( unsigned char data11, unsigned char data12,unsigned char data21,unsigned char data22,double oldw,double oldh)
{
	double detaw=oldw -(int)oldw;  //Δy
	double detah=oldh -(int)oldh;  //Δx
	//插值公式
	int data=((data22-data21-data12+data11)*detah*detaw+(data12-data11)*detaw+(data21-data11)*detah+data11)+0.5;
	//排除小于0和大于255的情况
	data=(data>255)?255:data;
	data=(data<0)?0:data;
	return (unsigned char)data;
}
/****************************************************** 
* 函数名称： 
* Matrix_trans()
*  
* 函数参数：
* double** matrix       - 待求转置的矩阵
* int row               - 原矩阵的行数
* int col               - 原矩阵的列数
*
* 返回值：
* 转置矩阵
*  
*说明：求矩阵转置
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
* 函数名称： 
* Matrix_mul()
*  
* 函数参数：
* double** left_matrix       - 左矩阵
* double** right_matrix      - 右矩阵
* int left_row               - 左矩阵行
* int left_col               - 左矩阵列（右矩阵的行）
* int right_col              - 右矩阵的列
*
* 返回值：
* 矩阵的乘积
*  
*说明：矩阵相乘
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

int Inv(double** a,double** inva,int dim )
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

/************************************************************************
* 函数名称： 
* Least_Square()
*  
* 函数参数：
* double** A               - 系数矩阵
* double** L               - 常数项
* int r                    - 矩阵A的行数
* int c                    - 矩阵A的列数
*
* 返回值：
* 所求变量的值
*  
*说明：最小二乘    公式V=AX-L
*      
************************************************************************/
double** Least_Square(double** A,double** L,int &r,int c)
{
	double** X;     //保存最小二乘求出的结果
	double sigma=0; //精度
	double pre_sigma=0;
	//误差
	double** V=new double*[r];
	for (int j=0;j<r;j++)
	{
		V[j]=new double[1];
		V[j][0]=0;
	}
	//每次迭代剔除误差后的A
	double** tempA=new double*[r];
	for (int j=0;j<r;j++)
	{
		tempA[j]=new double[c];
		memset(tempA[j],0,sizeof(double)*c);
	}	
	//每次迭代剔除误差后的L(c*1)c行1列
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
		//求系数矩阵的转置
		double** AT = Matrix_trans(A,r,c);

		//系数的转置乘以系数
		double** ATA=Matrix_mul(AT,A,c,r,c);

		//相乘结果求逆
		double** inv_ATA=new double*[c];
		for (int j=0;j<c;j++)
		{
			inv_ATA[j]=new double[c];
		}
		if (Inv(ATA,inv_ATA,c))
		{
			//求逆后再乘以系数的转置
			double** mul_xym=Matrix_mul(inv_ATA,AT,c,c,r);
			
			//最后乘以常数项
			X=Matrix_mul(mul_xym,L,c,r,1);

			//求误差
			double** AX=Matrix_mul(A,X,r,c,1);	
			for (int j=0;j<r;j++)
			{
				V[j][0]=AX[j][0]-L[j][0];
			}

			//评定精度
			double** VT=Matrix_trans(V,r,1);
			double** VTV=Matrix_mul(VT,V,1,r,1);
			sigma=sqrt(VTV[0][0]/(r-c));

			int tempr=0;//剔除误差后点的个数
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
			//将剔除粗差后的点赋值给A，用于循环
			for (int j=0;j<tempr;j++)
			{
				memcpy(A[j],tempA[j],sizeof(double)*c);
				L[j][0]=tempL[j][0];
			}
			//释放内存
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
			
			//改变总点数
			r=tempr;
		}
		else
		{
			//释放内存
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
	//释放内存
	return X;
}

/****************************************************** 
* 函数名称： 
* GetGradData()
*  
* 函数参数：
* unsigned char* imgdata - 读取的图像数据
* int colum     - 读取的图像数据列数
* int row       - 读取的图像数据行数
* int BitCount  - 读取的图像的波段数
*
* 返回值：
* 用于提取特征点的图像数据
*  
*说明：读取图像的灰度数据
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
* 函数名称： 
* ImgNomalize()
*  
* 函数参数：
* unsigned char* imgraydata - 读取的图像数据
* int colum     - 读取的图像数据列数
* int row       - 读取的图像数据行数
*
* 返回值：
* 用于提取特征点的图像数据
*  
*说明：图像的灰度值直接除以255归一化
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
* 函数名称： 
* descr_dist_sq()
*  
* 函数参数：
* Key_Point* f1 - 第一个特征
* Key_Point* f2 - 第二个特征
* 
*
* 返回值：
* 中值在数组中的序号
*  
*说明：计算两个特征间的欧式距离
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