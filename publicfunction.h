#ifndef PUBLICFUNCTION_H_
#define PUBLICFUNCTION_H_
#include "siftfeature.h"

extern void BubbleSort(double* Array,int n,int &minnum);
extern double Find_Mid(double* dv,int number);
extern unsigned char Bilinear( unsigned char data11, unsigned char data12,unsigned char data21,unsigned char data22,double oldw,double oldh);
extern double** Matrix_trans(double** matrix,int row,int col);
extern double** Matrix_mul(double** left_matrix,double** right_matrix,int left_row,int left_col,int right_col);
extern int Inv(double** a,double** inva,int dim );
extern double** Least_Square(double** A,double** L,int &r,int c);
extern unsigned char* GetGrayData( unsigned char* ,int ,int ,int );//读取图像数据
extern float* ImgNomalize( unsigned char* imgraydata,int colum,int row ); //图像归一化 
extern float descr_dist_sq( Key_Point* f1, Key_Point* f2, int nDim=128 );
#endif 