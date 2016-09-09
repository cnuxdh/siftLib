/************************************************************************/
/* ��ȡSIFT����                                                         */
/************************************************************************/
#ifndef SIFTFEATURE_H_FILE
#define SIFTFEATURE_H_FILE


#include "siftExports.h"



/******************************* Ԥ����� *****************************/
#define PI 3.141592653589793f
#define SIFT_IMG_DOWN 0               //�Ƿ���н�������0Ϊ�����н�����
#define GRID_REM_PIXEL_THRESHOLD 10   //���ڲ�����������ӿ�ֱ����ȥ
#define SIFT_IMG_GRID_THRESHOLD 1024  //ͼ��ֿ�Ĵ�С
#define OCTAVESTYE 0                  //��˹�ռ����ķ�ʽ��1��ʾ��˹�ռ�ֻ��һ�飬����Ϊ����ͼ���С���飬
#define SCALESPEROCTAVE 3             //ÿ�����S
#define SIFT_SIGMA 1.6f               //
#define SIFT_INIT_SIGMA 0.5f        //��ʼsigmaֵ
#define SIFT_CONTR_THRESHOLD 0.04f  //�޳��ͶԱȶȵĵ����ֵ,�ڼ�⼫ֵ��ʱ����
#define SIFT_CURV_THRESHOLD 10      //�޳���Ե��Ӧ��ʱ����ֵ
#define SIFT_IMG_BORDER 5           //����Ӱ���ԵС��SIFT_IMG_BORDER�����ص����򲻼��㼫ֵ
#define SIFT_MAX_INTERP_STEPS 5     //��ȷ��λ��ֵ��ʱѭ����������
#define SIFT_ORI_SIG_FCTR 1.5f      //����������ʱ��ȷ��ͳ���ݶ������С�Ĳ���
#define SIFT_ORI_RADIUS 3.0f * SIFT_ORI_SIG_FCTR //����������ʱ������뾶
#define SIFT_ORI_HIST_BINS 36       //�ؼ��㷽��������ֱ��ͼ����
#define SIFT_ORI_SMOOTH_PASSES 2    //������ֱ��ͼƽ���Ĵ���
#define SIFT_ORI_PEAK_RATIO 0.8f    //�ڶ�������Ϊ��ֵ��0.8
#define SIFT_DESCR_SCL_FCTR 3       //ȷ��������ͳ������Ĳ���
#define SIFT_DIST_THRESHOLD 0.6f    //���������ν�����֮��
#define SIFT_DESCR_MAG_THR  0.2f    //����������һ������ֵ
#define SIFT_INT_DESCR_FCTR 512.0f  //�����������С��ֵ��Ϊ������ϵ��
#define SIFT_DESCR_WIDTH 4          //����������Ŀ��4��4
#define SIFT_DESCR_HIST_BINS 8      //�����ӷ���ֱ��ͼ����
#define SIFT_DESCR_BINS   SIFT_DESCR_WIDTH*SIFT_DESCR_WIDTH*SIFT_DESCR_HIST_BINS//�����ӵ�ά�� 128


/**************************************ȫ�ֱ���*****************************************/

//�����¼��Ӧ���ĸ�˹�ռ�
struct SCALESPACE
{
	float* gauss_space;//���ڴ洢��˹�ռ��ͼ������
	float* dog_space;  //�洢��ָ�˹
	float gus_sigma;   //�洢��˹�߶ȣ���Ӱ���˹ƽ����ĳ߶�
	int gus_column;    //�߶Ȳ������
	int gus_row;       //�߶Ȳ������
};

//�����������¼�ؼ������������������������ڳ߶�
struct Key_Point
{
	int index;                 //������ĵ��
	int key_octave,key_intvl;  //���������ڵ���Ͳ� 
	float sub_intvl;           //��ȷ��λ�����������ڲ�ĸ�����
	int key_row;               //�������λ������
	int key_column;            //�������λ������
	float ori/*,mag*/;         //������oritation�Լ��ݶȷ�ֵ
	float initl_row;           //��ԭʼͼ���ϵ��к�
	float initl_column;        //��ԭʼͼ���ϵ��к�
	float scl;                 //�ؼ���ĳ߶�
	float scl_octave;          //
	//float descriptor[SIFT_DESCR_BINS];     //��������
	float descriptor[128];     //��������
	struct Key_Point *next;    //����
	void* feature_data;        //�û��Զ������������,��BBF�㷨��ʹ��
};

struct Feature
{
	int num;           //ÿ���ӿ�����������
	int row;           //�ӿ����ڵ��к�
	int col;           //�ӿ����ڵ��к�
	Key_Point *feature;//�ӿ����������
};

//extern "C" Feature** sift_features( float* imagedata,int colum,int row,int &keynumber,int &gridcol,int &gridrow );

DLL_EXPORT Key_Point* SiftFeaturesFloat( float *imgdata, int colum, int row, int &keynumber );



#endif