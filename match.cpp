// match.cpp : �������̨Ӧ�ó������ڵ㡣
/************************************************************************/
/* ��ȡSIFT����                                                         */
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
#define INPUT_GRID_THRED 6000 //�ֿ����Ӱ��ʱ��ÿ�ζ��������
/******************************** ȫ�ֱ��� ************************************/
string abspath = "data\\";

/**********************************������*************************************/
int _tmain(int argc, _TCHAR* argv[])
{

	//////////////////////////////////////////////////////////////////////////
	////��ȡpos����
	ImageFile* imagefile;//
	int image_num;
	int pos_num;
	vector<string> pos_vect;
	if (readtxt(abspath,pos_num,pos_vect))
	{
		image_num=pos_num/8;           //����pos���ݼ���Ӱ��ĸ���
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

	/****************************��ȡͼ������********************************/	
	for (int i=0;i<image_num;i++)
	{
		//��ȡ�ļ���
		string strpath;
		strpath=abspath+imagefile[i].imagename;
		const char* cfilepath=strpath.c_str();
		
		//��ȡͼ��
		Cimage image;
		image.Load(cfilepath); 
		if( image.info.pImage==NULL )
		{
			cout<<"��"<<cfilepath<<"Ӱ��ʧ�ܣ�"<<endl;
			return 0;
		}
		cout<<"��ȡ"<<cfilepath<<"��������";
		int imgWidth=image.head.biWidth; 
		int imgHeight=image.head.biHeight;
		int BitCount=image.head.biBitCount;
		printf(" ht:%d wd:%d \n", imgHeight, imgWidth);


		/**********************************��ȡͼ������**************************************/
		//��ȡ�Ҷ�����
		unsigned char *imgraydata=GetGrayData( image.info.pImage,imgWidth,imgHeight,BitCount );
		image.~Cimage();
		
		/*
		//Wallis�˲�
		unsigned char *imgwaldata=wallisfilter(imgraydata,imgWidth,imgHeight);
		delete []imgraydata;imgraydata=NULL;		
		//�Ҷ�ͼ���һ��
		float* imgnomdata=ImgNomalize(imgwaldata,imgWidth,imgHeight);
		delete []imgwaldata;
		imgwaldata=NULL;
		*/

		float* imgnomdata=ImgNomalize(imgraydata,imgWidth,imgHeight);

		//��ȡsift��
		int keynumber=0,gridrow=0,gridcol=0;
		Feature **feat;
		feat=sift_features(imgnomdata,imgWidth,imgHeight,keynumber,gridcol,gridrow);
		cout<<keynumber<<endl;


		//�洢Ӱ����Ϣ
		//imagefile[i].imageID=imagefile[i].imageID;   //Ӱ��ID
		imagefile[i].filepath=strpath;               //·��
		imagefile[i].feat=feat;                      //�����㣬��άָ��
		imagefile[i].imgWidth=imgWidth;              //Ӱ��Ŀ�
		imagefile[i].imgHeight=imgHeight;            //Ӱ��ĸ�
		imagefile[i].bands=BitCount/8;               //Ӱ��Ĳ�����
		imagefile[i].gridcol=gridcol;                //Ӱ�����Ŀ���
		imagefile[i].gridrow=gridrow;                //Ӱ��߷���Ŀ���
		imagefile[i].keynumber=keynumber;            //Ӱ���ϵ����������

		//��ʾӰ��������
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
	/*����������Ӱ��ƥ�䣨ģ��ƥ�䣩                                        */
	/************************************************************************/
	int pointnum=0;        //�ܵ���
	MatchPoints* matchpoints=sift_match( imagefile, image_num, pointnum );
	//ȥ����λ��ͬ��ͬ����
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
	//ɾ��������
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
	/* �޳��ֲ��                                                           */
	/************************************************************************/
	MatchPoints* new_matchpoints=new MatchPoints[image_num-1];
	double** trans_h=new double*[image_num-1];
	int new_pointnum=0;
	for (int i=0;i<image_num-1;i++) 
	{
		new_matchpoints[i]=del_gross_error(matchpoints[i],imagefile[i].imgWidth,imagefile[i].imgHeight);
		
		//cout<<"ƥ�������"<<matchpoints[i].mat_key_num<<endl;
		/************************************************************************/
		/* ������С���˼���ͶӰ�任ϵ��                                         */
		/************************************************************************/
		double** deta_h=inline_cal_trans(new_matchpoints[i].matpoint,new_matchpoints[i].mat_key_num);	

		//Ϊ��ͼ����ͶӰ�任���󿪱��ڴ�
		trans_h[i] = new double[9];

		//��ÿ��Ӱ���ͶӰ����ֵ
		for (int j=0;j<8;j++)
		{
			trans_h[i][j]=deta_h[j][0];
			//cout<<trans_h[i][j]<<endl;
		}
		trans_h[i][8]=1;

		new_matchpoints[i].imageID1=matchpoints[i].imageID1;
		new_matchpoints[i].imageID2=matchpoints[i].imageID2;
		new_pointnum+=new_matchpoints[i].mat_key_num;

		//�ͷ��ڴ�
		for (int j=0;j<8;j++)
		{
			delete []deta_h[j];
		}
		delete []deta_h;
	}

	/************************************************************************/
	/*�������ʾƥ����                                                    */
	/************************************************************************/
	cout<<"ƥ�������"<<new_pointnum<<endl;
	if (output_point( new_matchpoints,imagefile,image_num,new_pointnum,trans_h )!=1)
	{
		cout<<"���ƥ����ļ�ʧ�ܣ�"<<endl;
	}
	return 0;
}
