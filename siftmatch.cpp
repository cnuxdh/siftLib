//#include "stdafx.h"

#include "math.h"

//#include "opencv/cv.h"
#include "kdtree.h"
#include "siftfeature.h"
#include "siftmatch.h"
#include "publicfunction.h"
#include "omp.h"





/**********************************函数声明**************************************/
All_MatchPoint* match(Key_Point *feat1,int keynum1,Key_Point *feat2,int keynum2,int &match_num,int d);
Key_Point* get_overlap_pot(Feature** feat,int gridcol,int gridrow,int keynum);



/****************************************************** 
* 函数名称： 
* get_overlap_grid()
*  
* 函数参数：
* Feature** feat - 图像数据
* int gridcol    - 列方向块的个数
* int gridrow    - 行方向块的个数
* int keynum     - gridcol和gridrow决定的重叠区域的关键点的总个数
*
* 返回值：
* gridcol和gridrow决定的重叠区域的特征点数
*  
*说明：对于左右影像重叠的块，将其特征点合并为一个数组中，这样匹配时计算量较少
********************************************************/
Key_Point* get_overlap_pot(Feature** feat,int gridcol,int gridrow,int keynum)
{
	int n=0;
	//为重叠区域特征点开辟内存
	Key_Point *keypoint=new Key_Point[keynum];
	//存储重叠区域特征点
	for (int i=0;i<gridrow;i++)
	{
		for (int j=0;j<gridcol;j++)
		{
			for (int k=0;k<feat[i][j].num;k++)
			{
				keypoint[n++]=feat[i][j].feature[k];
			}
		}
	}
	return keypoint;	
}

/****************************************************** 
* Name:
* PairMatch()
* 
* Parameters:
* 
*
*
* Explanation:
* Matching between two feature sets based on SIFT features
********************************************************/
void SiftPairMatch(Key_Point* feat1, int nFeat1, Key_Point* feat2, int nFeat2, MatchPoint** pMatchs, int* nMatch, int nDim)
{
	MatchPoint* mp1=NULL;//存储匹配上的点

	/************************************************************************/
	//建立kd树
	Kd_Node* kd_root = kdtree_build( feat2, nFeat2, /*SIFT_DESCR_BINS*/nDim );
	int matchnum=0;
	
	for (int k=0; k<nFeat1; k++)
	{
		Key_Point** nbrs=NULL;
		Key_Point* feat=&feat1[k];
		//BBF搜索
		int kd = kdtree_bbf_knn( kd_root, feat, 2, &nbrs, KDTREE_BBF_MAX_NN_CHKS, /*SIFT_DESCR_BINS*/nDim );
		if( kd == 2 )
		{
			float d0 = descr_dist_sq( feat, nbrs[0], nDim );
			float d1 = descr_dist_sq( feat, nbrs[1], nDim );
			if( d0 < d1 * SIFT_DIST_THRESHOLD )
			{		
				MatchPoint* mp2=new MatchPoint[1];   //采用链表存储匹配的点
				mp2->next=mp1;                       //链表赋值
				mp1=mp2; 

				mp2->id1 = feat->index;
				mp2->id2 = nbrs[0]->index;
				mp2->imageY1 = feat->initl_row;    //imagefile[i].imgHeight-feat->initl_row;//左影像
				mp2->imageX1 = feat->initl_column;	
				mp2->imageY2 = nbrs[0]->initl_row; //imagefile[i+1].imgHeight-nbrs[0]->initl_row;//右影像
				mp2->imageX2 = nbrs[0]->initl_column;	
				matchnum++;        //记录点的个数
			}//最近距离和次近距离之比满足阈值
		}//是否找到两个近邻
		if (nbrs)
		{
			delete []nbrs;
		}
	}

	/************************************************************************/
	/* 将链表存储改为数组存储                                               */
	/************************************************************************/
	MatchPoint* mphead=mp1;
	*pMatchs = new MatchPoint[matchnum];
	for (int j=0;j<matchnum;j++)
	{
		(*pMatchs)[j].id1 = mp1->id1;
		(*pMatchs)[j].id2 = mp1->id2;
		(*pMatchs)[j].imageX1=mp1->imageX1;  //列
		(*pMatchs)[j].imageY1=mp1->imageY1;  //行
		(*pMatchs)[j].imageX2=mp1->imageX2;  //列
		(*pMatchs)[j].imageY2=mp1->imageY2;  //行
		
		mp1=mp1->next;
	}
	*nMatch = matchnum;

	//remove the bad matches
	//去除点位相同的同名点
	MatchPoint* matpoint = *pMatchs;
	int pointnum=0;
	int tempnum=matchnum;
	int newtempnum=tempnum;
	int temp=0;
	for (int j=0;j<tempnum;j++)
	{
		temp=0;
		for (int k=j+1; k<tempnum-1; k++)
		{
			if ((matpoint[j].imageX1==matpoint[k].imageX1 &&
				 matpoint[j].imageY1==matpoint[k].imageY1) ||
				(matpoint[j].imageX2==matpoint[k].imageX2 &&
				 matpoint[j].imageY2==matpoint[k].imageY2 ))
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
	MatchPoint *newmatpoint=new MatchPoint[newtempnum];
	for (int j=0;j<tempnum;j++)
	{
		temp=0;
		for (int k=j+1;k<tempnum;k++)
		{
			if ((matpoint[j].imageX1==matpoint[k].imageX1 &&
				 matpoint[j].imageY1==matpoint[k].imageY1)||
				(matpoint[j].imageX2==matpoint[k].imageX2 &&
				 matpoint[j].imageY2==matpoint[k].imageY2 ))
			{
				temp++;
				break;
			}
		}
		if (!temp)
		{
			newmatpoint[inum].imageX1=matpoint[j].imageX1;
			newmatpoint[inum].imageY1=matpoint[j].imageY1;
			newmatpoint[inum].imageX2=matpoint[j].imageX2;
			newmatpoint[inum].imageY2=matpoint[j].imageY2;
			newmatpoint[inum].id1=matpoint[j].id1;
			newmatpoint[inum].id2=matpoint[j].id2;
			inum++;
		}
	}			


	delete[] *pMatchs;
	*pMatchs = newmatpoint;
	*nMatch = inum;

	if (mphead)
	{
		MatchPoint* pnode=mphead;
		MatchPoint* temp;
		while(pnode!=NULL)
		{		
			temp=pnode;
			pnode=pnode->next;
			delete []temp;
		}
		mphead=NULL;
	}
}


/****************************************************** 
* 函数名称： 
* sift_match()
*  
* 函数参数：
* ImageFile* imagefile - 所有的影像
* int image_num        - 图像的个数
* int &pointnum        - 整个区域匹配点的总个数
*
* 返回值：
* 相邻影像的匹配点
*  
*说明：用BBF算法进行特征点匹配。首先对相邻影像进行匹配，
再利用匹配上的点对两度以上重叠度的点进行匹配
********************************************************/
MatchPoints* sift_match(ImageFile* imagefile,int image_num,int &pointnum)
{
	MatchPoints* mps=new MatchPoints[image_num-1];
	pointnum=0; //所有影像匹配的总点数
	for (int i=0;i<image_num-1;i++)
	{	
		MatchPoint* mp1=NULL;//存储匹配上的点
		/************************************************************************/
		/*获取相邻的两张影像的匹配点                                            */
		/************************************************************************/
		Key_Point* feat1=get_overlap_pot(imagefile[i].feat,imagefile[i].gridcol,imagefile[i].gridrow,imagefile[i].keynumber);
		Key_Point* feat2=get_overlap_pot(imagefile[i+1].feat,imagefile[i+1].gridcol,imagefile[i+1].gridrow,imagefile[i+1].keynumber);
		//建立kd树
		Kd_Node* kd_root = kdtree_build( feat2, imagefile[i+1].keynumber, SIFT_DESCR_BINS );
		int matchnum=0;

		for (int k=0;k<imagefile[i].keynumber;k++)
		{
			Key_Point** nbrs=NULL;
			Key_Point* feat=&feat1[k];
			//BBF搜索
			int kd = kdtree_bbf_knn( kd_root, feat, 2, &nbrs, KDTREE_BBF_MAX_NN_CHKS,SIFT_DESCR_BINS );
			if( kd == 2 )
			{
				float d0 = descr_dist_sq( feat, nbrs[0] );
				float d1 = descr_dist_sq( feat, nbrs[1] );
				if( d0 < d1 * SIFT_DIST_THRESHOLD )
				{		
					MatchPoint* mp2=new MatchPoint[1];   //采用链表存储匹配的点
					mp2->next=mp1;                        //链表赋值
					mp1=mp2; 
					mp2->imageX1=imagefile[i].imgHeight-feat->initl_row;//左影像
					mp2->imageY1=feat->initl_column;	
					mp2->imageX2=imagefile[i+1].imgHeight-nbrs[0]->initl_row;//右影像
					mp2->imageY2=nbrs[0]->initl_column;	
					matchnum++;        //记录点的个数
				}//最近距离和次近距离之比满足阈值
			}//是否找到两个近邻
			if (nbrs)
			{
				delete []nbrs;
			}
		}

		/************************************************************************/
		/* 将链表存储改为数组存储                                               */
		/************************************************************************/
		MatchPoint* mphead=mp1;
		mps[i].imageID1=imagefile[i].imageID;        //左影像
		mps[i].imageID2=imagefile[i+1].imageID;      //右影像
		mps[i].mat_key_num=matchnum;                 //匹配点数
		MatchPoint* mp=new MatchPoint[matchnum];
		for (int j=0;j<matchnum;j++)
		{
			mp[j].imageX1=mp1->imageX1;  //列
			mp[j].imageY1=mp1->imageY1;  //行
			mp[j].imageX2=mp1->imageX2;  //列
			mp[j].imageY2=mp1->imageY2;  //行
			mp1=mp1->next;
		}
		mps[i].matpoint=mp;
		pointnum+=matchnum;//整个区域总匹配点数	

		//释放内存
		delete []feat1;
		delete []feat2;
		if (mphead)
		{
			MatchPoint* pnode=mphead;
			MatchPoint* temp;
			while(pnode!=NULL)
			{		
				temp=pnode;
				pnode=pnode->next;
				delete []temp;
			}
			mphead=NULL;
		}
		//delete []kd_root;
		/*if (kd_root)
		{
			Kd_Node* pnode=kd_root;
			Kd_Node* temp;
			while(pnode!=NULL)
			{		
				temp=pnode;
				pnode=pnode->next;
				delete []temp;
			}
			kd_root=NULL;
		}*/

	}
	
	return mps;
}

/****************************************************** 
* 函数名称： 
* match()
*  
* 函数参数：
* Key_Point *feat1 - 第一幅影像的特征点
* int keynum1      - 第一幅影像的特征点个数
* Key_Point *feat2 - 第二幅影像的特征点
* int keynum2      - 第二幅影像的特征点个数
* int &match_num   - 匹配上的点的个数
* int d            - 特征向量的维数
* 
* 返回值：
* 匹配的点
*  
* 说明：穷举搜索方法匹配点
********************************************************/
All_MatchPoint* match(Key_Point *feat1,int keynum1,Key_Point *feat2,int keynum2,int &match_num,int d)
{
	float dist;
	int key;
	float mindist,sub_mindist;
	All_MatchPoint* matchpoint=new All_MatchPoint[keynum1];//因为keynum1<keynum2
	float* distan=new float[keynum2];//因为keynum1<keynum2

	for (int i=0;i<keynum1;i++)
	{
		//mindist=0;
		for (int j=0;j<keynum2;j++)
		{
			dist=0;
			for (int k=0;k<d;k++)
			{
				dist=dist+(feat1[i].descriptor[k]-feat2[j].descriptor[k])*
					(feat1[i].descriptor[k]-feat2[j].descriptor[k]);
			}
			distan[j]=sqrt(float(dist));
			//mindist=min(distan[j],mindist);
		}
		//计算最近距离
		mindist=distan[0];key=0;
		for (int j=1;j<keynum2;j++)
		{
			if (distan[j]<mindist)
			{
				mindist=distan[j];
				key=j;
			}
		}

		//if (mindist<300)//阈值法首先排除最近距离很大的点
		//{
		//	//计算次近距离
		if (key==0)
		{
			sub_mindist=distan[1];
		} 
		else
		{
			sub_mindist=distan[0];
		}
		for (int j=0;j<keynum2;j++)
		{
			if (j!=key)
			{
				if (distan[j]<sub_mindist)
				{
					sub_mindist=distan[j];
				}
			}
		}
		if (mindist<sub_mindist*SIFT_DIST_THRESHOLD)//记录的左影像表示特征点少的影像！！！！！！！！！！！！！！
		{
			matchpoint[match_num].leftrow = int (floor(feat1[i].initl_row+0.5f));
			matchpoint[match_num].leftcolumn = int (floor(feat1[i].initl_column+0.5f));
			matchpoint[match_num].rightrow = int (floor(feat2[key].initl_row+0.5f));
			matchpoint[match_num].rightcolumn = int (floor(feat2[key].initl_column+0.5f));
			match_num++;
		}
		/*}*/	
	}
	//清空缓存
	delete []distan;
	return matchpoint;
}

/****************************************************** 
* 函数名称： 
* all_match()
*  
* 函数参数：
* Feature** feat1     - 左影像重叠区的特征点
* Feature** feat2     - 右影像重叠区的特征点
* IplImage* stacked   - 用于显示图像
* int width1          - 左影像的宽，用于显示图像
* int height2         - 右影像的高，用于显示图像
* int keynum1         - 左影像重叠区特征点个数
* int keynum2         - 右影像重叠区特征点个数
* int &matchnum       - 左右影像重叠区匹配点的个数
*
* 返回值：
* 
*  
*说明：穷举法进行特征点匹配
********************************************************/
/*	
void all_match(Key_Point* feat1,Key_Point* feat2,IplImage* stacked,
	int width1,int width2 ,int keynum1,int keynum2,int &matchnum)
{
	CvPoint pt1, pt2;
	All_MatchPoint *matchpoint;
	int d=SIFT_DESCR_BINS;

	if (keynum1<=keynum2)//保证传入的数据中前边参数的特征点比后边的少，从而避免出现匹配数比两者最小的数大
	{
		matchpoint=match(feat1,keynum1,feat2,keynum2,matchnum,d);
	} 
	else
	{
		matchpoint=match(feat2,keynum2,feat1,keynum1,matchnum,d);
	}

	if (keynum1<=keynum2)
	{
		for (int i=0;i<matchnum;i++)
		{
			pt1 = cvPoint(cvRound( matchpoint[i].leftcolumn ) ,cvRound(matchpoint[i].leftrow ) );
			pt2 = cvPoint(cvRound( matchpoint[i].rightcolumn ),cvRound(matchpoint[i].rightrow ) );
			pt2.x += width1;
			cvLine( stacked, pt1, pt2, CV_RGB(255,0,255), 1, 8, 0 );
		}
	} 
	else
	{
		for (int i=0;i<matchnum;i++)
		{
			pt1 = cvPoint( cvRound( matchpoint[i].rightcolumn ) ,cvRound(matchpoint[i].rightrow ) );
			pt2 = cvPoint( cvRound( matchpoint[i].leftcolumn ) ,cvRound(matchpoint[i].leftrow ) );
			pt2.x += width2;
			cvLine( stacked, pt1, pt2, CV_RGB(255,0,255), 1, 8, 0 );
		}
	}

}
*/