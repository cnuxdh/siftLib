#include "stdafx.h"
#include "imgfeatures.h"


/*
void draw_lowe_feature( IplImage* img, Key_Point* feat,CvScalar color );
void draw_my_feature( IplImage* img, Key_Point* feat,CvScalar color );

//显示特征点
void draw_features( IplImage* img, Key_Point* feat, int n )
{
	CvScalar color = CV_RGB( 255, 0, 255 );
	

	for(int i = 0; i < n; i++ )
		draw_my_feature( img, feat + i, color );
}
//显示特征点（不显示主方向）
void draw_my_feature( IplImage* img, Key_Point* feat,CvScalar color )
{
	int center_x = cvRound( feat->initl_column );
	int center_y = cvRound( img->height-feat->initl_row );
	double scl = feat->scl;
	int len = cvRound(5.0*scl)/2;

	int start_x = center_x-len;
	int start_y = center_y-len;
	int end_x = center_x+len;
	int end_y = center_y+len;

	CvPoint start1 = cvPoint( start_x, center_y );
	CvPoint end1 = cvPoint( end_x, center_y );
	CvPoint start2 = cvPoint( center_x, start_y );
	CvPoint end2 = cvPoint( center_x, end_y );

	cvLine( img, start1, end1, color, 1, 8, 0 );
	cvLine( img, start2, end2, color, 1, 8, 0 );
}
//显示特征点（显示主方向）
void draw_lowe_feature( IplImage* img, Key_Point* feat,CvScalar color )
{
	int len, hlen, blen, start_x, start_y, end_x, end_y, h1_x, h1_y, h2_x, h2_y;
	double scl, ori;
	double scale = 5.0;
	double hscale = 0.75;
	CvPoint start, end, h1, h2;

	
	start_x = cvRound( feat->initl_column );
	start_y = cvRound( img->height-feat->initl_row );
	scl = feat->scl;
	ori = feat->ori;
	len = cvRound( scl * scale );
	hlen = cvRound( scl * hscale );
	blen = len - hlen;
	end_x = cvRound( len *  cos( ori ) ) + start_x;
	end_y = cvRound( len * sin( ori ) ) + start_y;
	h1_x = cvRound( blen *  cos( ori + CV_PI / 18.0 ) ) + start_x;
	h1_y = cvRound( blen * sin( ori + CV_PI / 18.0 ) ) + start_y;
	h2_x = cvRound( blen *  cos( ori - CV_PI / 18.0 ) ) + start_x;
	h2_y = cvRound( blen * sin( ori - CV_PI / 18.0 ) ) + start_y;
	start = cvPoint( start_x, start_y );
	end = cvPoint( end_x, end_y );
	h1 = cvPoint( h1_x, h1_y );
	h2 = cvPoint( h2_x, h2_y );

	cvLine( img, start, end, color, 1, 8, 0 );
	cvLine( img, end, h1, color, 1, 8, 0 );
	cvLine( img, end, h2, color, 1, 8, 0 );
}

//显示匹配图像

IplImage* stack_imgs( IplImage* img1, IplImage* img2 )
{
	IplImage* stacked = cvCreateImage(cvSize((img1->width+img2->width),
		 MAX(img1->height,img2->height)),
		IPL_DEPTH_8U, 3 );

	cvZero( stacked );
	cvSetImageROI( stacked, cvRect( 0, 0, img1->width, img1->height ) );
	cvAdd( img1, stacked, stacked, NULL );
	cvSetImageROI( stacked, cvRect( img1->width, 0, img2->width, img2->height) );
	cvAdd( img2, stacked, stacked, NULL );
	cvResetImageROI( stacked );

	return stacked;
}*/

