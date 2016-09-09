/************************************************************************/
/* 提取SIFT特征                                                         */
/************************************************************************/
#ifndef SIFTFEATURE_H_FILE
#define SIFTFEATURE_H_FILE


#include "siftExports.h"



/******************************* 预定义宏 *****************************/
#define PI 3.141592653589793f
#define SIFT_IMG_DOWN 0               //是否进行降采样，0为不进行降采样
#define GRID_REM_PIXEL_THRESHOLD 10   //对于不满足整块的子块直接舍去
#define SIFT_IMG_GRID_THRESHOLD 1024  //图像分块的大小
#define OCTAVESTYE 0                  //高斯空间分组的方式，1表示高斯空间只有一组，其他为根据图像大小分组，
#define SCALESPEROCTAVE 3             //每组层数S
#define SIFT_SIGMA 1.6f               //
#define SIFT_INIT_SIGMA 0.5f        //初始sigma值
#define SIFT_CONTR_THRESHOLD 0.04f  //剔除低对比度的点的阈值,在检测极值点时采用
#define SIFT_CURV_THRESHOLD 10      //剔除边缘响应点时的阈值
#define SIFT_IMG_BORDER 5           //距离影像边缘小于SIFT_IMG_BORDER个像素的区域不计算极值
#define SIFT_MAX_INTERP_STEPS 5     //精确定位极值点时循环的最大次数
#define SIFT_ORI_SIG_FCTR 1.5f      //计算主方向时，确定统计梯度区域大小的参数
#define SIFT_ORI_RADIUS 3.0f * SIFT_ORI_SIG_FCTR //计算主方向时的区域半径
#define SIFT_ORI_HIST_BINS 36       //关键点方向描述的直方图个数
#define SIFT_ORI_SMOOTH_PASSES 2    //主方向直方图平滑的次数
#define SIFT_ORI_PEAK_RATIO 0.8f    //第二主方向为峰值的0.8
#define SIFT_DESCR_SCL_FCTR 3       //确定描述子统计邻域的参数
#define SIFT_DIST_THRESHOLD 0.6f    //最近距离与次近距离之比
#define SIFT_DESCR_MAG_THR  0.2f    //特征向量归一化的阈值
#define SIFT_INT_DESCR_FCTR 512.0f  //最后将特征向量小数值改为整数的系数
#define SIFT_DESCR_WIDTH 4          //描述子区域的宽度4×4
#define SIFT_DESCR_HIST_BINS 8      //描述子方向直方图个数
#define SIFT_DESCR_BINS   SIFT_DESCR_WIDTH*SIFT_DESCR_WIDTH*SIFT_DESCR_HIST_BINS//描述子的维数 128


/**************************************全局变量*****************************************/

//保存记录相应组层的高斯空间
struct SCALESPACE
{
	float* gauss_space;//用于存储高斯空间的图像数据
	float* dog_space;  //存储差分高斯
	float gus_sigma;   //存储高斯尺度，即影像高斯平滑后的尺度
	int gus_column;    //尺度层的列数
	int gus_row;       //尺度层的行数
};

//用链表，保存记录关键点的组数、层数、坐标和所在尺度
struct Key_Point
{
	int index;                 //特征点的点号
	int key_octave,key_intvl;  //特征点所在的组和层 
	float sub_intvl;           //精确定位后特征点所在层的改正数
	int key_row;               //特征点的位置行数
	int key_column;            //特征点的位置列数
	float ori/*,mag*/;         //主方向oritation以及梯度幅值
	float initl_row;           //在原始图像上的行号
	float initl_column;        //在原始图像上的列号
	float scl;                 //关键点的尺度
	float scl_octave;          //
	//float descriptor[SIFT_DESCR_BINS];     //特征向量
	float descriptor[128];     //特征向量
	struct Key_Point *next;    //链表
	void* feature_data;        //用户自定义的数据类型,在BBF算法里使用
};

struct Feature
{
	int num;           //每个子块的特征点个数
	int row;           //子块所在的行号
	int col;           //子块所在的列号
	Key_Point *feature;//子块包含的特征
};

//extern "C" Feature** sift_features( float* imagedata,int colum,int row,int &keynumber,int &gridcol,int &gridrow );

DLL_EXPORT Key_Point* SiftFeaturesFloat( float *imgdata, int colum, int row, int &keynumber );



#endif