#ifndef SIFTMATCH_H
#define SIFTMATCH_H


//#include "opencv/cxtypes.h"
#include <string>
#include "iostream"
#include "siftfeature.h"




#define KDTREE_BBF_MAX_NN_CHKS 200 //���ѭ������
#define max_overlap 2        //����5���ص��ĵ�
//
using namespace std;


//�洢���ƥ��ĵ�
struct All_MatchPoint
{
	int leftrow;    //ƥ�������ͼ���е��к�
	int leftcolumn; //ƥ�������ͼ���е��к�
	int rightrow;   //ƥ�������ͼ���е��к�
	int rightcolumn;//ƥ�������ͼ���е��к�
};

//ƥ���
struct MatchPoint
{
	int   id1,id2;
	float imageX1,imageY1;               //ƥ�������Ӱ���ϣ�X��Y)����,�У���ֵ�����Ͻ�Ϊԭ�㣩
	float imageX2,imageY2;               //ƥ�������Ӱ���ϣ�X��Y)����,�У���ֵ�����Ͻ�Ϊԭ�㣩
	struct MatchPoint* next;             //����ָ�룬ָ����һ��ƥ���
};
//����Ӱ���ϵ�ƥ���
struct MatchPoints
{
	string imageID1,imageID2;            //ƥ��Ӱ��ı��
	int mat_key_num;                     //ƥ���ĸ���
	MatchPoint* matpoint;                //ƥ���
};

//�洢Ӱ�����Ϣ
struct ImageFile
{
	string imageID;  //Ӱ���ID
	string imagename;//Ӱ������֣�������չ��
	double* posdata;
	int imgWidth;    //Ӱ��Ŀ�
	int imgHeight;   //Ӱ��ĸ�
	int bands;       //Ӱ�񲨶���       
	int gridcol;     //Ӱ������ֵĿ���
	int gridrow;     //Ӱ������ֵĿ���
	int keynumber;   //Ӱ����ȡ����������
	string filepath; //Ӱ���·��
	Feature** feat;  //Ӱ���������
};


//extern MatchPoints* sift_match(ImageFile* imagefile,int image_num,int &pointnum);

//extern void all_match(Key_Point* feat1,Key_Point* feat2,IplImage* stacked, int width1,int width2 ,int keynum1,int keynum2,int &matchnum);


DLL_EXPORT void SiftPairMatch(Key_Point* feat1, int nFeat1, Key_Point* feat2, int nFeat2,
							  MatchPoint** pMatchs, int* nMatch, int nDim=128);



#endif 