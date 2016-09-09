#include "stdafx.h"
#include "Wallisfilter.h"
#include <math.h>


unsigned char* wallisfilter(unsigned char* grayimgdata,int width,int height)
{
	unsigned char* wallisimgdata=new unsigned char[width*height];
	int filterwindowsize=FILTERWINDOWSIZE*2+1;
	int gridwindowsize=GRIDWINDOWSIZE;
	float valueB=VALUEB;
	float valueC=VALUEC;
	float valuemean=VALUEMEAN;
	float valuesigma=VALUESIGMA;
	int gridwidthnum=(width-filterwindowsize)/gridwindowsize+1;
	int gridheightnum=(height-filterwindowsize)/gridwindowsize+1;
	float *gridr0=new float[gridwidthnum*gridheightnum];memset(gridr0,0,sizeof(float)*gridwidthnum*gridheightnum);
	float *gridr1=new float[gridwidthnum*gridheightnum];memset(gridr1,0,sizeof(float)*gridwidthnum*gridheightnum);	
	int filterwindowsize2=filterwindowsize*filterwindowsize;
	unsigned char *griddata=new unsigned char[filterwindowsize2];
	unsigned char *ptemp1,*ptemp2;
	float rmean=0.0f,rsigma=0.0f;	

	for (int i=0;i<gridheightnum;i++)
	{
		int bh=i*gridwindowsize;
		int igridwidthnum=i*gridwidthnum;
		for (int j=0;j<gridwidthnum;j++)
		{
			int bw=j*gridwindowsize;
			 ptemp2=grayimgdata+bh*width+bw;
			 ptemp1=griddata;
			//从图像中读取滤波窗口数据
			for (int k=0;k<filterwindowsize;k++)
			{
				memcpy(ptemp1,ptemp2,filterwindowsize);
				ptemp1+=filterwindowsize;
				ptemp2+=width;
			}
			//计算wallis参数
			float gridmean=0.0f;
			float gridsigma=0.0f;
			for (int k=0;k<filterwindowsize2;k++)
			{
				gridmean+=griddata[k];
				gridsigma+=griddata[k]*griddata[k];
			}
			gridmean/=filterwindowsize2;
			gridsigma=gridsigma/filterwindowsize2-gridmean*gridmean;
			gridsigma=gridsigma<0?0:gridsigma;
			gridsigma=sqrt(gridsigma);
			float r0,r1;
			if (gridsigma==0.0f)
			{
				r1=1.0f;
				r0=valueB*valuemean-valueB*gridmean;
			}
			else
			{
				r1=valueC*valuesigma/(valueC*gridsigma+(1.0f-valueC)*valuesigma);
				r0=valueB*valuemean+(1.0f-valueB-r1)*gridmean;
			}
			//////////////////////////////////////////////////////////////////////////
			gridr0[igridwidthnum+j]=r0;
			gridr1[igridwidthnum+j]=r1;
			float tempvalue=griddata[filterwindowsize2/2];
			float rc=tempvalue*r1+r0;
			rmean+=rc;
			//////////////////////////////////////////////////////////////////////////
			rsigma+=rc*rc;
		}
	}
	//////////////////////////////////////////////////////////////////////////
	rmean/=(gridheightnum*gridwidthnum);
	rsigma=rsigma/(gridheightnum*gridwidthnum)-rmean*rmean;
	rsigma=rsigma<0?0:rsigma;
	rsigma=sqrt(rsigma);
	//////////////////////////////////////////////////////////////////////////
	//对wallis参数插值，计算wallis滤波后的值
	for (int i=0;i<height;i++)
	{
		int ih=i-FILTERWINDOWSIZE;
		for (int j=0;j<width;j++)
		{
			int iw=j-FILTERWINDOWSIZE;
			int grid_h=ih/gridwindowsize;
			int grid_w=iw/gridwindowsize;

			grid_w=grid_w>gridwidthnum-2?gridwidthnum-2:grid_w;
			grid_h=grid_h>gridheightnum-2?gridheightnum-2:grid_h;
			grid_w=grid_w<0?0:grid_w;
			grid_h=grid_h<0?0:grid_h;
			float detaw=(float)(iw%gridwindowsize)/(float)gridwindowsize;
			float detah=(float)(ih%gridwindowsize)/(float)gridwindowsize;

			float data11=gridr0[grid_h*gridwidthnum+grid_w];
			float data12=gridr0[grid_h*gridwidthnum+grid_w+1];
			float data21=gridr0[(grid_h+1)*gridwidthnum+grid_w];
			float data22=gridr0[(grid_h+1)*gridwidthnum+grid_w+1];
			float r0=((data22-data21-data12+data11)*detah*detaw+(data12-data11)*detaw+(data21-data11)*detah+data11);
			
			data11=gridr1[grid_h*gridwidthnum+grid_w];
			data12=gridr1[grid_h*gridwidthnum+grid_w+1];
			data21=gridr1[(grid_h+1)*gridwidthnum+grid_w];
			data22=gridr1[(grid_h+1)*gridwidthnum+grid_w+1];
			float r1=((data22-data21-data12+data11)*detah*detaw+(data12-data11)*detaw+(data21-data11)*detah+data11);
	
			int tempdata=(grayimgdata[i*width+j]*r1+r0-rmean)*valuesigma/rsigma+valuemean+0.5f;
			tempdata=tempdata>255?255:tempdata;
			tempdata=tempdata<0?0:tempdata;
			wallisimgdata[i*width+j]=tempdata;
		}
	}
	delete []gridr0;gridr0=NULL;
	delete []gridr1;gridr1=NULL;
	delete []griddata;griddata=NULL;
	return wallisimgdata;
}