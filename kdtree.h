#ifndef KDTREE_H
#define KDTREE_H

#include "siftfeature.h"
/**k-d树的一个节点 */
struct Kd_Node
{
	int ki;                /**< partition key index 分割点的索引*/
	float kv;              /**< partition key value 分割点的值，即求得的中值*/
	int leaf;              /**< 1 if node is a leaf, 0 otherwise 是否为叶子节点*/
	Key_Point* features;   /**< features at this node 节点的特征点*/
	int n;                 /**< number of features 以该节点为父亲节点的特征的总数*/
	Kd_Node* kd_left;      /**< left child 左节点*/
	Kd_Node* kd_right;     /**< right child 右节点*/
};

/*************************** 函数接口 *****************************/

/****************************************************** 
* 函数名称： 
* kdtree_build()
*  
* 函数参数：
* struct Key_Point* features - 特征点；
* int n                      - 特征点的个数
* int d                      - 特征向量的维数
*
* 返回值：
* 特征点的kd树
*  
*说明：建立特征点的kd树
********************************************************/
extern Kd_Node* kdtree_build( Key_Point* features, int n,int d);

/****************************************************** 
* 函数名称： 
* kdtree_bbf_knn()
*  
* 函数参数：
* 
* struct Key_Point* feat - 特征点；
* 
*
* 返回值：
* 搜索到的近邻个数
*  
*说明：用BBF（best bin first）算法在特征点kd树中搜索最近邻和次近邻
********************************************************/
extern int kdtree_bbf_knn( Kd_Node* kd_root, Key_Point* feat,
						  int k, Key_Point*** nbrs, int max_nn_chks, int d);

/****************************************************** 
* 函数名称： 
* kdtree_release()
*  
* 函数参数：
* 
* struct kd_node* kd_root - ；
* 
*
* 返回值：
* 
*  
*说明：
********************************************************/
extern void kdtree_release( Kd_Node* kd_root );


#endif