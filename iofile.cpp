#include "stdafx.h"
#include "iofile.h"
#include "imgfeatures.h"
#include "publicfunction.h"
#include <sstream>
#include <fstream>
#include "omp.h"
//#include "Cimage.h"
//#include "opencv/highgui.h"
#include <iomanip>


void PixelToWorld(double  a[6], float P, float L, double &x, double &y)		//像素坐标向地理坐标转换函数
{
	x = a[0] + a[1] * P + a[2] * L;
	y = a[3] + a[4] * P + a[5] * L;
}
/****************************************************** 
* 函数名称： 
* output_point()
*  
* 函数参数：
* Image_MatchPoints* img_match_point          - 匹配点
* ImageFile* imagefile                        - 图像信息
* int image_num                               - 图像个数
* int pointnum                                - 匹配点数
* 
* 返回值：
* 1成功，0失败
*  
*说明：保存匹配点，并显示
********************************************************/
/*
int output_point( MatchPoints* matchpoints,ImageFile* imagefile,int image_num,int pointnum,double **trans_h )
{
	
	ofstream outfile("data\\matchpoint.txt",ios::out);

	if (!outfile)
	{
		cout<<"输出结果失败！"<<endl;
		return 0;
	}
	
	int potnum=0;//输出数据时的PointID	
	
	outfile<<"点号 左宽 左高 右宽 右高"<<endl;
	for (int i=0;i<image_num-1;i++)
	{	
		for (int j=0;j<matchpoints[i].mat_key_num;j++)
		{	
			outfile<<setiosflags(ios::fixed)<<setprecision(0);
			outfile<<setw(10)<<j;
			outfile<<setiosflags(ios::fixed)<<setprecision(7);
			outfile<<setw(20)<<matchpoints[i].matpoint[j].imageY1;
			outfile<<setw(20)<<matchpoints[i].matpoint[j].imageX1;
			outfile<<setw(20)<<matchpoints[i].matpoint[j].imageY2;
			outfile<<setw(20)<<matchpoints[i].matpoint[j].imageX2;
			outfile<<endl;

		}	
	}
	outfile.close();
	
	for (int i=0;i<image_num-1;i++)
	{	
		//显示匹配的影像1
		string winname=imagefile[i].imagename+" and "+imagefile[i+1].imagename;

		IplImage* img1, * img2, * stacked;
		const char* cfilepath1=imagefile[i].filepath.c_str();
		img1 = cvLoadImage(cfilepath1,1);
		if( ! img1 )
			cout<<"unable to load image from"<<cfilepath1<<endl;
		const char* cfilepath2=imagefile[i+1].filepath.c_str();
		img2 = cvLoadImage(cfilepath2,1);
		if( ! img2 )
			cout<<"unable to load image from"<<cfilepath2<<endl;
		stacked = stack_imgs( img1, img2 );

		CvPoint pt1, pt2;
		for (int j=0;j<matchpoints[i].mat_key_num;j++)
		{
			pt1 = cvPoint( cvRound( matchpoints[i].matpoint[j].imageY1 ),cvRound(matchpoints[i].matpoint[j].imageX1 ));//imageY是列，imageX是行
			pt2 = cvPoint( cvRound( matchpoints[i].matpoint[j].imageY2 ),cvRound(matchpoints[i].matpoint[j].imageX2 ));
			pt2.x += imagefile[i].imgWidth;
			cvLine( stacked, pt1, pt2, CV_RGB(255,0,255), 1,8, 0 );
		}

		cvNamedWindow( winname.c_str(), 0 );
		cvShowImage( winname.c_str(), stacked );
		cvWaitKey( 0 );
		string temppath="data//匹配图像.tif";
		cvSaveImage(temppath.c_str(),stacked); 

		cvReleaseImage( &stacked );
		cvReleaseImage( &img1 );
		cvReleaseImage( &img2 );

		//////////////////////////////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////
		//显示匹配影像2
		IplImage* img3, * img4, * stacked1;
		const char* cfilepath3=imagefile[i].filepath.c_str();
		img3 = cvLoadImage(cfilepath3,1);
		if( ! img3 )
			cout<<"unable to load image from"<<cfilepath3<<endl;
		const char* cfilepath4=imagefile[i+1].filepath.c_str();
		img4 = cvLoadImage(cfilepath4,1);
		if( ! img4 )
			cout<<"unable to load image from"<<cfilepath4<<endl;
		stacked1 = stack_imgs( img3, img4 );

		for (int j=0;j<matchpoints[i].mat_key_num;j++)
		{
			//先画左图像
			int center_x = cvRound( matchpoints[i].matpoint[j].imageY1 );
			int center_y = cvRound(matchpoints[i].matpoint[j].imageX1 );//imageY是列，imageX是行
			int start_x = center_x - 5;
			int start_y = center_y - 5;
			int end_x = center_x + 5;
			int end_y = center_y + 5;
			//四个点
			CvPoint start1 = cvPoint( start_x, center_y );
			CvPoint end1 = cvPoint( end_x, center_y );
			CvPoint start2 = cvPoint( center_x, start_y );
			CvPoint end2 = cvPoint( center_x, end_y );
			//画线
			cvLine(stacked1,start1,end1, CV_RGB(255,0,255), 1,8, 0 );
			cvLine(stacked1,start2,end2, CV_RGB(255,0,255), 1,8, 0 );
			//再画右图像
			center_x = cvRound( matchpoints[i].matpoint[j].imageY2 )+imagefile[i].imgWidth;
			center_y = cvRound(matchpoints[i].matpoint[j].imageX2 );//imageY是列，imageX是行
			start_x = center_x - 5;
			start_y = center_y - 5;
			end_x = center_x + 5;
			end_y = center_y + 5;
			//四个点
			start1 = cvPoint( start_x, center_y );
			end1 = cvPoint( end_x, center_y );
			start2 = cvPoint( center_x, start_y );
			end2 = cvPoint( center_x, end_y );
			//画线
			cvLine(stacked1,start1,end1, CV_RGB(255,0,255), 1,8, 0 );
			cvLine(stacked1,start2,end2, CV_RGB(255,0,255), 1,8, 0 );
		}
		cvNamedWindow( winname.c_str(), 0 );
		cvShowImage( winname.c_str(), stacked1 );
		cvWaitKey( 0 );
		cvSaveImage("data//匹配图像1.tif",stacked1); 
		cvReleaseImage( &stacked1 );
		cvReleaseImage( &img3 );
		cvReleaseImage( &img4 );

	}
	return 1;
}
*/


/****************************************************** 
* 函数名称： 
* readtxt()
*  
* 函数参数：
* string abspath           - 文件路径
* int &num                 - 文件中数据的个数
* vector<string> &vect     - 文件中的数据
* 
* 返回值：
* 1成功，0失败
*  
*说明：读取POS数据
********************************************************/
int readtxt(string abspath,int &num,vector<string> &vect)
{
	//读取pos数据
	string path=abspath+"pos.txt";
	ifstream inputfile(path.c_str(),ios::in|ios::_Nocreate);
	if (!inputfile)
	{
		cout<<"打开外方位元素文件失败"<<endl;
		return 0;
	}
        //用于存储从pos文件中读取的数据
	string temp;
	num=0;
	while (!inputfile.eof())
	{	
		inputfile>>temp;
		vect.push_back(temp);
		num++;
	}
	inputfile.close();
	return 1;
}

