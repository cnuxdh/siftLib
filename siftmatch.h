#ifndef SIFTMATCH_H
#define SIFTMATCH_H


//#include "opencv/cxtypes.h"
#include <string>
#include "iostream"
#include "siftfeature.h"




#define KDTREE_BBF_MAX_NN_CHKS 200 //最大循环次数
#define max_overlap 2        //计算5度重叠的点
//
using namespace std;


//存储穷举匹配的点
struct All_MatchPoint
{
	int leftrow;    //匹配点在左图像中的行号
	int leftcolumn; //匹配点在左图像中的列号
	int rightrow;   //匹配点在右图像中的行号
	int rightcolumn;//匹配点在右图像中的列号
};

//匹配点
struct MatchPoint
{
	int   id1,id2;
	float imageX1,imageY1;               //匹配点在左影像上（X，Y)（行,列）的值（左上角为原点）
	float imageX2,imageY2;               //匹配点在右影像上（X，Y)（行,列）的值（左上角为原点）
	struct MatchPoint* next;             //链表指针，指向下一个匹配点
};
//两张影像上的匹配点
struct MatchPoints
{
	string imageID1,imageID2;            //匹配影像的编号
	int mat_key_num;                     //匹配点的个数
	MatchPoint* matpoint;                //匹配点
};

//存储影像的信息
struct ImageFile
{
	string imageID;  //影像的ID
	string imagename;//影像的名字，包括扩展名
	double* posdata;
	int imgWidth;    //影像的宽
	int imgHeight;   //影像的高
	int bands;       //影像波段数       
	int gridcol;     //影像宽所分的块数
	int gridrow;     //影像高所分的块数
	int keynumber;   //影像提取的特征点数
	string filepath; //影像的路径
	Feature** feat;  //影像的特征点
};


//extern MatchPoints* sift_match(ImageFile* imagefile,int image_num,int &pointnum);

//extern void all_match(Key_Point* feat1,Key_Point* feat2,IplImage* stacked, int width1,int width2 ,int keynum1,int keynum2,int &matchnum);


DLL_EXPORT void SiftPairMatch(Key_Point* feat1, int nFeat1, Key_Point* feat2, int nFeat2,
							  MatchPoint** pMatchs, int* nMatch, int nDim=128);



#endif 