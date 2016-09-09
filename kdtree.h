#ifndef KDTREE_H
#define KDTREE_H

#include "siftfeature.h"
/**k-d����һ���ڵ� */
struct Kd_Node
{
	int ki;                /**< partition key index �ָ�������*/
	float kv;              /**< partition key value �ָ���ֵ������õ���ֵ*/
	int leaf;              /**< 1 if node is a leaf, 0 otherwise �Ƿ�ΪҶ�ӽڵ�*/
	Key_Point* features;   /**< features at this node �ڵ��������*/
	int n;                 /**< number of features �Ըýڵ�Ϊ���׽ڵ������������*/
	Kd_Node* kd_left;      /**< left child ��ڵ�*/
	Kd_Node* kd_right;     /**< right child �ҽڵ�*/
};

/*************************** �����ӿ� *****************************/

/****************************************************** 
* �������ƣ� 
* kdtree_build()
*  
* ����������
* struct Key_Point* features - �����㣻
* int n                      - ������ĸ���
* int d                      - ����������ά��
*
* ����ֵ��
* �������kd��
*  
*˵���������������kd��
********************************************************/
extern Kd_Node* kdtree_build( Key_Point* features, int n,int d);

/****************************************************** 
* �������ƣ� 
* kdtree_bbf_knn()
*  
* ����������
* 
* struct Key_Point* feat - �����㣻
* 
*
* ����ֵ��
* �������Ľ��ڸ���
*  
*˵������BBF��best bin first���㷨��������kd������������ںʹν���
********************************************************/
extern int kdtree_bbf_knn( Kd_Node* kd_root, Key_Point* feat,
						  int k, Key_Point*** nbrs, int max_nn_chks, int d);

/****************************************************** 
* �������ƣ� 
* kdtree_release()
*  
* ����������
* 
* struct kd_node* kd_root - ��
* 
*
* ����ֵ��
* 
*  
*˵����
********************************************************/
extern void kdtree_release( Kd_Node* kd_root );


#endif