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


void PixelToWorld(double  a[6], float P, float L, double &x, double &y)		//�����������������ת������
{
	x = a[0] + a[1] * P + a[2] * L;
	y = a[3] + a[4] * P + a[5] * L;
}
/****************************************************** 
* �������ƣ� 
* output_point()
*  
* ����������
* Image_MatchPoints* img_match_point          - ƥ���
* ImageFile* imagefile                        - ͼ����Ϣ
* int image_num                               - ͼ�����
* int pointnum                                - ƥ�����
* 
* ����ֵ��
* 1�ɹ���0ʧ��
*  
*˵��������ƥ��㣬����ʾ
********************************************************/
/*
int output_point( MatchPoints* matchpoints,ImageFile* imagefile,int image_num,int pointnum,double **trans_h )
{
	
	ofstream outfile("data\\matchpoint.txt",ios::out);

	if (!outfile)
	{
		cout<<"������ʧ�ܣ�"<<endl;
		return 0;
	}
	
	int potnum=0;//�������ʱ��PointID	
	
	outfile<<"��� ��� ��� �ҿ� �Ҹ�"<<endl;
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
		//��ʾƥ���Ӱ��1
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
			pt1 = cvPoint( cvRound( matchpoints[i].matpoint[j].imageY1 ),cvRound(matchpoints[i].matpoint[j].imageX1 ));//imageY���У�imageX����
			pt2 = cvPoint( cvRound( matchpoints[i].matpoint[j].imageY2 ),cvRound(matchpoints[i].matpoint[j].imageX2 ));
			pt2.x += imagefile[i].imgWidth;
			cvLine( stacked, pt1, pt2, CV_RGB(255,0,255), 1,8, 0 );
		}

		cvNamedWindow( winname.c_str(), 0 );
		cvShowImage( winname.c_str(), stacked );
		cvWaitKey( 0 );
		string temppath="data//ƥ��ͼ��.tif";
		cvSaveImage(temppath.c_str(),stacked); 

		cvReleaseImage( &stacked );
		cvReleaseImage( &img1 );
		cvReleaseImage( &img2 );

		//////////////////////////////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////
		//��ʾƥ��Ӱ��2
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
			//�Ȼ���ͼ��
			int center_x = cvRound( matchpoints[i].matpoint[j].imageY1 );
			int center_y = cvRound(matchpoints[i].matpoint[j].imageX1 );//imageY���У�imageX����
			int start_x = center_x - 5;
			int start_y = center_y - 5;
			int end_x = center_x + 5;
			int end_y = center_y + 5;
			//�ĸ���
			CvPoint start1 = cvPoint( start_x, center_y );
			CvPoint end1 = cvPoint( end_x, center_y );
			CvPoint start2 = cvPoint( center_x, start_y );
			CvPoint end2 = cvPoint( center_x, end_y );
			//����
			cvLine(stacked1,start1,end1, CV_RGB(255,0,255), 1,8, 0 );
			cvLine(stacked1,start2,end2, CV_RGB(255,0,255), 1,8, 0 );
			//�ٻ���ͼ��
			center_x = cvRound( matchpoints[i].matpoint[j].imageY2 )+imagefile[i].imgWidth;
			center_y = cvRound(matchpoints[i].matpoint[j].imageX2 );//imageY���У�imageX����
			start_x = center_x - 5;
			start_y = center_y - 5;
			end_x = center_x + 5;
			end_y = center_y + 5;
			//�ĸ���
			start1 = cvPoint( start_x, center_y );
			end1 = cvPoint( end_x, center_y );
			start2 = cvPoint( center_x, start_y );
			end2 = cvPoint( center_x, end_y );
			//����
			cvLine(stacked1,start1,end1, CV_RGB(255,0,255), 1,8, 0 );
			cvLine(stacked1,start2,end2, CV_RGB(255,0,255), 1,8, 0 );
		}
		cvNamedWindow( winname.c_str(), 0 );
		cvShowImage( winname.c_str(), stacked1 );
		cvWaitKey( 0 );
		cvSaveImage("data//ƥ��ͼ��1.tif",stacked1); 
		cvReleaseImage( &stacked1 );
		cvReleaseImage( &img3 );
		cvReleaseImage( &img4 );

	}
	return 1;
}
*/


/****************************************************** 
* �������ƣ� 
* readtxt()
*  
* ����������
* string abspath           - �ļ�·��
* int &num                 - �ļ������ݵĸ���
* vector<string> &vect     - �ļ��е�����
* 
* ����ֵ��
* 1�ɹ���0ʧ��
*  
*˵������ȡPOS����
********************************************************/
int readtxt(string abspath,int &num,vector<string> &vect)
{
	//��ȡpos����
	string path=abspath+"pos.txt";
	ifstream inputfile(path.c_str(),ios::in|ios::_Nocreate);
	if (!inputfile)
	{
		cout<<"���ⷽλԪ���ļ�ʧ��"<<endl;
		return 0;
	}
        //���ڴ洢��pos�ļ��ж�ȡ������
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

