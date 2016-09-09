//#include "stdafx.h"

#include "math.h"

//#include "opencv/cv.h"
#include "kdtree.h"
#include "siftfeature.h"
#include "siftmatch.h"
#include "publicfunction.h"
#include "omp.h"





/**********************************��������**************************************/
All_MatchPoint* match(Key_Point *feat1,int keynum1,Key_Point *feat2,int keynum2,int &match_num,int d);
Key_Point* get_overlap_pot(Feature** feat,int gridcol,int gridrow,int keynum);



/****************************************************** 
* �������ƣ� 
* get_overlap_grid()
*  
* ����������
* Feature** feat - ͼ������
* int gridcol    - �з����ĸ���
* int gridrow    - �з����ĸ���
* int keynum     - gridcol��gridrow�������ص�����Ĺؼ�����ܸ���
*
* ����ֵ��
* gridcol��gridrow�������ص��������������
*  
*˵������������Ӱ���ص��Ŀ飬����������ϲ�Ϊһ�������У�����ƥ��ʱ����������
********************************************************/
Key_Point* get_overlap_pot(Feature** feat,int gridcol,int gridrow,int keynum)
{
	int n=0;
	//Ϊ�ص����������㿪���ڴ�
	Key_Point *keypoint=new Key_Point[keynum];
	//�洢�ص�����������
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
	MatchPoint* mp1=NULL;//�洢ƥ���ϵĵ�

	/************************************************************************/
	//����kd��
	Kd_Node* kd_root = kdtree_build( feat2, nFeat2, /*SIFT_DESCR_BINS*/nDim );
	int matchnum=0;
	
	for (int k=0; k<nFeat1; k++)
	{
		Key_Point** nbrs=NULL;
		Key_Point* feat=&feat1[k];
		//BBF����
		int kd = kdtree_bbf_knn( kd_root, feat, 2, &nbrs, KDTREE_BBF_MAX_NN_CHKS, /*SIFT_DESCR_BINS*/nDim );
		if( kd == 2 )
		{
			float d0 = descr_dist_sq( feat, nbrs[0], nDim );
			float d1 = descr_dist_sq( feat, nbrs[1], nDim );
			if( d0 < d1 * SIFT_DIST_THRESHOLD )
			{		
				MatchPoint* mp2=new MatchPoint[1];   //��������洢ƥ��ĵ�
				mp2->next=mp1;                       //����ֵ
				mp1=mp2; 

				mp2->id1 = feat->index;
				mp2->id2 = nbrs[0]->index;
				mp2->imageY1 = feat->initl_row;    //imagefile[i].imgHeight-feat->initl_row;//��Ӱ��
				mp2->imageX1 = feat->initl_column;	
				mp2->imageY2 = nbrs[0]->initl_row; //imagefile[i+1].imgHeight-nbrs[0]->initl_row;//��Ӱ��
				mp2->imageX2 = nbrs[0]->initl_column;	
				matchnum++;        //��¼��ĸ���
			}//�������ʹν�����֮��������ֵ
		}//�Ƿ��ҵ���������
		if (nbrs)
		{
			delete []nbrs;
		}
	}

	/************************************************************************/
	/* ������洢��Ϊ����洢                                               */
	/************************************************************************/
	MatchPoint* mphead=mp1;
	*pMatchs = new MatchPoint[matchnum];
	for (int j=0;j<matchnum;j++)
	{
		(*pMatchs)[j].id1 = mp1->id1;
		(*pMatchs)[j].id2 = mp1->id2;
		(*pMatchs)[j].imageX1=mp1->imageX1;  //��
		(*pMatchs)[j].imageY1=mp1->imageY1;  //��
		(*pMatchs)[j].imageX2=mp1->imageX2;  //��
		(*pMatchs)[j].imageY2=mp1->imageY2;  //��
		
		mp1=mp1->next;
	}
	*nMatch = matchnum;

	//remove the bad matches
	//ȥ����λ��ͬ��ͬ����
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
* �������ƣ� 
* sift_match()
*  
* ����������
* ImageFile* imagefile - ���е�Ӱ��
* int image_num        - ͼ��ĸ���
* int &pointnum        - ��������ƥ�����ܸ���
*
* ����ֵ��
* ����Ӱ���ƥ���
*  
*˵������BBF�㷨����������ƥ�䡣���ȶ�����Ӱ�����ƥ�䣬
������ƥ���ϵĵ�����������ص��ȵĵ����ƥ��
********************************************************/
MatchPoints* sift_match(ImageFile* imagefile,int image_num,int &pointnum)
{
	MatchPoints* mps=new MatchPoints[image_num-1];
	pointnum=0; //����Ӱ��ƥ����ܵ���
	for (int i=0;i<image_num-1;i++)
	{	
		MatchPoint* mp1=NULL;//�洢ƥ���ϵĵ�
		/************************************************************************/
		/*��ȡ���ڵ�����Ӱ���ƥ���                                            */
		/************************************************************************/
		Key_Point* feat1=get_overlap_pot(imagefile[i].feat,imagefile[i].gridcol,imagefile[i].gridrow,imagefile[i].keynumber);
		Key_Point* feat2=get_overlap_pot(imagefile[i+1].feat,imagefile[i+1].gridcol,imagefile[i+1].gridrow,imagefile[i+1].keynumber);
		//����kd��
		Kd_Node* kd_root = kdtree_build( feat2, imagefile[i+1].keynumber, SIFT_DESCR_BINS );
		int matchnum=0;

		for (int k=0;k<imagefile[i].keynumber;k++)
		{
			Key_Point** nbrs=NULL;
			Key_Point* feat=&feat1[k];
			//BBF����
			int kd = kdtree_bbf_knn( kd_root, feat, 2, &nbrs, KDTREE_BBF_MAX_NN_CHKS,SIFT_DESCR_BINS );
			if( kd == 2 )
			{
				float d0 = descr_dist_sq( feat, nbrs[0] );
				float d1 = descr_dist_sq( feat, nbrs[1] );
				if( d0 < d1 * SIFT_DIST_THRESHOLD )
				{		
					MatchPoint* mp2=new MatchPoint[1];   //��������洢ƥ��ĵ�
					mp2->next=mp1;                        //����ֵ
					mp1=mp2; 
					mp2->imageX1=imagefile[i].imgHeight-feat->initl_row;//��Ӱ��
					mp2->imageY1=feat->initl_column;	
					mp2->imageX2=imagefile[i+1].imgHeight-nbrs[0]->initl_row;//��Ӱ��
					mp2->imageY2=nbrs[0]->initl_column;	
					matchnum++;        //��¼��ĸ���
				}//�������ʹν�����֮��������ֵ
			}//�Ƿ��ҵ���������
			if (nbrs)
			{
				delete []nbrs;
			}
		}

		/************************************************************************/
		/* ������洢��Ϊ����洢                                               */
		/************************************************************************/
		MatchPoint* mphead=mp1;
		mps[i].imageID1=imagefile[i].imageID;        //��Ӱ��
		mps[i].imageID2=imagefile[i+1].imageID;      //��Ӱ��
		mps[i].mat_key_num=matchnum;                 //ƥ�����
		MatchPoint* mp=new MatchPoint[matchnum];
		for (int j=0;j<matchnum;j++)
		{
			mp[j].imageX1=mp1->imageX1;  //��
			mp[j].imageY1=mp1->imageY1;  //��
			mp[j].imageX2=mp1->imageX2;  //��
			mp[j].imageY2=mp1->imageY2;  //��
			mp1=mp1->next;
		}
		mps[i].matpoint=mp;
		pointnum+=matchnum;//����������ƥ�����	

		//�ͷ��ڴ�
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
* �������ƣ� 
* match()
*  
* ����������
* Key_Point *feat1 - ��һ��Ӱ���������
* int keynum1      - ��һ��Ӱ������������
* Key_Point *feat2 - �ڶ���Ӱ���������
* int keynum2      - �ڶ���Ӱ������������
* int &match_num   - ƥ���ϵĵ�ĸ���
* int d            - ����������ά��
* 
* ����ֵ��
* ƥ��ĵ�
*  
* ˵���������������ƥ���
********************************************************/
All_MatchPoint* match(Key_Point *feat1,int keynum1,Key_Point *feat2,int keynum2,int &match_num,int d)
{
	float dist;
	int key;
	float mindist,sub_mindist;
	All_MatchPoint* matchpoint=new All_MatchPoint[keynum1];//��Ϊkeynum1<keynum2
	float* distan=new float[keynum2];//��Ϊkeynum1<keynum2

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
		//�����������
		mindist=distan[0];key=0;
		for (int j=1;j<keynum2;j++)
		{
			if (distan[j]<mindist)
			{
				mindist=distan[j];
				key=j;
			}
		}

		//if (mindist<300)//��ֵ�������ų��������ܴ�ĵ�
		//{
		//	//����ν�����
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
		if (mindist<sub_mindist*SIFT_DIST_THRESHOLD)//��¼����Ӱ���ʾ�������ٵ�Ӱ�񣡣�������������������������
		{
			matchpoint[match_num].leftrow = int (floor(feat1[i].initl_row+0.5f));
			matchpoint[match_num].leftcolumn = int (floor(feat1[i].initl_column+0.5f));
			matchpoint[match_num].rightrow = int (floor(feat2[key].initl_row+0.5f));
			matchpoint[match_num].rightcolumn = int (floor(feat2[key].initl_column+0.5f));
			match_num++;
		}
		/*}*/	
	}
	//��ջ���
	delete []distan;
	return matchpoint;
}

/****************************************************** 
* �������ƣ� 
* all_match()
*  
* ����������
* Feature** feat1     - ��Ӱ���ص�����������
* Feature** feat2     - ��Ӱ���ص�����������
* IplImage* stacked   - ������ʾͼ��
* int width1          - ��Ӱ��Ŀ�������ʾͼ��
* int height2         - ��Ӱ��ĸߣ�������ʾͼ��
* int keynum1         - ��Ӱ���ص������������
* int keynum2         - ��Ӱ���ص������������
* int &matchnum       - ����Ӱ���ص���ƥ���ĸ���
*
* ����ֵ��
* 
*  
*˵������ٷ�����������ƥ��
********************************************************/
/*	
void all_match(Key_Point* feat1,Key_Point* feat2,IplImage* stacked,
	int width1,int width2 ,int keynum1,int keynum2,int &matchnum)
{
	CvPoint pt1, pt2;
	All_MatchPoint *matchpoint;
	int d=SIFT_DESCR_BINS;

	if (keynum1<=keynum2)//��֤�����������ǰ�߲�����������Ⱥ�ߵ��٣��Ӷ��������ƥ������������С������
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