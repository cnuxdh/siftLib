// match.cpp : 定义控制台应用程序的入口点。
/************************************************************************/
/* 提取SIFT特征                                                         */
/************************************************************************/
#include "stdafx.h"
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include "omp.h"

//#include "opencv/cv.h"
//#include "opencv/cxcore.h"
//#include "opencv/highgui.h"

//#include "Cimage.h"
#include "Wallisfilter.h"
#include "kdtree.h"
#include "siftfeature.h"
#include "imgfeatures.h"
#include "siftmatch.h"
#include "publicfunction.h"
#include "ransac.h"
#include "iofile.h"

using namespace std;
#define INPUT_GRID_THRED 6000 //分块读入影像时，每次读入的行数
/******************************** 全局变量 ************************************/
string abspath = "data\\";

/**********************************主函数*************************************/
int _tmain(int argc, _TCHAR* argv[])
{

	//////////////////////////////////////////////////////////////////////////
	////读取pos数据
	ImageFile* imagefile;//
	int image_num;
	int pos_num;
	vector<string> pos_vect;
	if (readtxt(abspath,pos_num,pos_vect))
	{
		image_num=pos_num/8;           //根据pos数据计算影像的个数
		imagefile=new ImageFile[image_num];//
		for (int i=0;i<image_num;i++)
		{
			imagefile[i].posdata=new double [6];
			stringstream ss;
			for (int j=0;j<6;j++)
			{	
				ss.clear();
				ss<<pos_vect[i*8+j+2];
				ss>>imagefile[i].posdata[j];
			}
			ss.clear();
			ss<<pos_vect[i*8];
			ss>>imagefile[i].imageID;
			imagefile[i].imagename=pos_vect[i*8+1];
		}
	} 
	else
	{
		exit(0);
	}

	/****************************读取图像数据********************************/	
	for (int i=0;i<image_num;i++)
	{
		//读取文件名
		string strpath;
		strpath=abspath+imagefile[i].imagename;
		const char* cfilepath=strpath.c_str();
		
		//读取图像
		Cimage image;
		image.Load(cfilepath); 
		if( image.info.pImage==NULL )
		{
			cout<<"打开"<<cfilepath<<"影像失败！"<<endl;
			return 0;
		}
		cout<<"提取"<<cfilepath<<"的特征点";
		int imgWidth=image.head.biWidth; 
		int imgHeight=image.head.biHeight;
		int BitCount=image.head.biBitCount;
		printf(" ht:%d wd:%d \n", imgHeight, imgWidth);


		/**********************************读取图像数据**************************************/
		//获取灰度数据
		unsigned char *imgraydata=GetGrayData( image.info.pImage,imgWidth,imgHeight,BitCount );
		image.~Cimage();
		
		/*
		//Wallis滤波
		unsigned char *imgwaldata=wallisfilter(imgraydata,imgWidth,imgHeight);
		delete []imgraydata;imgraydata=NULL;		
		//灰度图像归一化
		float* imgnomdata=ImgNomalize(imgwaldata,imgWidth,imgHeight);
		delete []imgwaldata;
		imgwaldata=NULL;
		*/

		float* imgnomdata=ImgNomalize(imgraydata,imgWidth,imgHeight);

		//提取sift点
		int keynumber=0,gridrow=0,gridcol=0;
		Feature **feat;
		feat=sift_features(imgnomdata,imgWidth,imgHeight,keynumber,gridcol,gridrow);
		cout<<keynumber<<endl;


		//存储影像信息
		//imagefile[i].imageID=imagefile[i].imageID;   //影像ID
		imagefile[i].filepath=strpath;               //路径
		imagefile[i].feat=feat;                      //特征点，二维指针
		imagefile[i].imgWidth=imgWidth;              //影像的宽
		imagefile[i].imgHeight=imgHeight;            //影像的高
		imagefile[i].bands=BitCount/8;               //影像的波段数
		imagefile[i].gridcol=gridcol;                //影像宽方向的块数
		imagefile[i].gridrow=gridrow;                //影像高方向的快数
		imagefile[i].keynumber=keynumber;            //影像上的特征点个数

		//显示影像特征点
		/*IplImage* img;

		img = cvLoadImage(cfilepath,1);
		if( ! img )
			cout<<"unable to load image from"<<cfilepath<<endl;
		
		for (int k=0;k<gridrow;k++)
		{
			for (int j=0;j<gridcol;j++)
			{
				draw_features( img, feat[k][j].feature, feat[k][j].num );
			}
		}
		cvNamedWindow( cfilepath, 0 );
		cvShowImage( cfilepath, img );
		cvWaitKey( 0 );
		string temppath=abspath+imagefile[i].imagename+"out.tif";
		cvSaveImage(temppath.c_str(),img); 
		cvDestroyWindow(cfilepath); 
		cvReleaseImage( &img );*/
	}

	/************************************************************************/
	/*航带内相邻影像匹配（模型匹配）                                        */
	/************************************************************************/
	int pointnum=0;        //总点数
	MatchPoints* matchpoints=sift_match( imagefile, image_num, pointnum );
	//去除点位相同的同名点
	pointnum=0;
	for (int i=0;i<image_num-1;i++)
	{
		int tempnum=matchpoints[i].mat_key_num;
		int newtempnum=tempnum;
		int temp=0;
		for (int j=0;j<tempnum;j++)
		{
			temp=0;
			for (int k=j+1;k<tempnum-1;k++)
			{
				if ((matchpoints[i].matpoint[j].imageX1==matchpoints[i].matpoint[k].imageX1 &&
				     matchpoints[i].matpoint[j].imageY1==matchpoints[i].matpoint[k].imageY1)||
					 (matchpoints[i].matpoint[j].imageX2==matchpoints[i].matpoint[k].imageX2 &&
					 matchpoints[i].matpoint[j].imageY2==matchpoints[i].matpoint[k].imageY2 ))
				{
					temp++;
					break;
				}
			}
			if (temp)
			{
				newtempnum--;
			}
		}
		int inum=0;
		MatchPoint *matpoint=new MatchPoint[newtempnum];
		for (int j=0;j<tempnum;j++)
		{
			temp=0;
			for (int k=j+1;k<tempnum;k++)
			{
				if ((matchpoints[i].matpoint[j].imageX1==matchpoints[i].matpoint[k].imageX1 &&
					matchpoints[i].matpoint[j].imageY1==matchpoints[i].matpoint[k].imageY1)||
					(matchpoints[i].matpoint[j].imageX2==matchpoints[i].matpoint[k].imageX2 &&
					matchpoints[i].matpoint[j].imageY2==matchpoints[i].matpoint[k].imageY2 ))
				{
					temp++;
					break;
				}
			}
			if (!temp)
			{
				matpoint[inum].imageX1=matchpoints[i].matpoint[j].imageX1;
				matpoint[inum].imageY1=matchpoints[i].matpoint[j].imageY1;
				matpoint[inum].imageX2=matchpoints[i].matpoint[j].imageX2;
				matpoint[inum].imageY2=matchpoints[i].matpoint[j].imageY2;
				inum++;
			}
		}			
		delete []matchpoints[i].matpoint;matchpoints[i].matpoint=NULL;
		matchpoints[i].matpoint=matpoint;
		matchpoints[i].mat_key_num=inum;
		pointnum+=newtempnum;
	}
	//删除特征点
	for (int i=0;i<image_num;i++)
	{
		for (int j=0;j<imagefile[i].gridrow;j++)
		{	
			for (int k=0;k<imagefile[i].gridcol;k++)
			{
				if (imagefile[i].feat[j][k].num!=0)
				{	
					delete []imagefile[i].feat[j][k].feature;
				}
			}
		}
	}

	/************************************************************************/
	/* 剔除粗差点                                                           */
	/************************************************************************/
	MatchPoints* new_matchpoints=new MatchPoints[image_num-1];
	double** trans_h=new double*[image_num-1];
	int new_pointnum=0;
	for (int i=0;i<image_num-1;i++) 
	{
		new_matchpoints[i]=del_gross_error(matchpoints[i],imagefile[i].imgWidth,imagefile[i].imgHeight);
		
		//cout<<"匹配点数："<<matchpoints[i].mat_key_num<<endl;
		/************************************************************************/
		/* 利用最小二乘计算投影变换系数                                         */
		/************************************************************************/
		double** deta_h=inline_cal_trans(new_matchpoints[i].matpoint,new_matchpoints[i].mat_key_num);	

		//为两图像间的投影变换矩阵开辟内存
		trans_h[i] = new double[9];

		//给每张影像的投影矩阵赋值
		for (int j=0;j<8;j++)
		{
			trans_h[i][j]=deta_h[j][0];
			//cout<<trans_h[i][j]<<endl;
		}
		trans_h[i][8]=1;

		new_matchpoints[i].imageID1=matchpoints[i].imageID1;
		new_matchpoints[i].imageID2=matchpoints[i].imageID2;
		new_pointnum+=new_matchpoints[i].mat_key_num;

		//释放内存
		for (int j=0;j<8;j++)
		{
			delete []deta_h[j];
		}
		delete []deta_h;
	}

	/************************************************************************/
	/*输出并显示匹配结果                                                    */
	/************************************************************************/
	cout<<"匹配点数："<<new_pointnum<<endl;
	if (output_point( new_matchpoints,imagefile,image_num,new_pointnum,trans_h )!=1)
	{
		cout<<"输出匹配点文件失败！"<<endl;
	}
	return 0;
}
