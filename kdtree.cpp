//#include "stdafx.h"

#include "stdlib.h"
#include "stdio.h"
#include "string.h"

#include "kdtree.h"
#include "minpq.h"
#include "publicfunction.h"
#include <math.h>



struct BBF_Data
{
	float d;
	void* old_data;
};

/******************************* 函数声明 *******************************/

Kd_Node* kd_node_init( Key_Point*, int );
void expand_kd_node_subtree( Kd_Node*, int );
void assign_part_key( Kd_Node* ,int);
float median_select( float*, int);
float rank_select( float*, int, int );
void insertion_sort( float*, int );
int partition_array( float*, int, float );
void partition_features( Kd_Node* );
Kd_Node* explore_to_leaf( Kd_Node*, Key_Point*, Min_Pq* ,int);
int insert_into_nbr_array( Key_Point*, Key_Point**, int, int );

/******************************** 函数体 **********************************/

/****************************************************** 
* 函数名称： 
* kdtree_build()
*  
* 函数参数：
* Key_Point* features - 特征点
* int n                      - 特征点的个数
* int d                      - 特征向量的维数
*
* 返回值：
* 特征点的kd树
*  
*说明：建立特征点的kd树
********************************************************/
Kd_Node* kdtree_build( Key_Point* features, int n, int d)
{
	Kd_Node* kd_root;//存储根节点

	if( !features || n<= 0 )
	{
		/*cout<<"没有匹配点！"<<endl;*/
		return NULL;
	}

	//初始化根节点
	kd_root = kd_node_init( features, n);
	//以kd_root为根节点建立kd树
	expand_kd_node_subtree( kd_root,d );

	return kd_root;
}
/****************************************************** 
* 函数名称： 
* kdtree_bbf_knn()
*  
* 函数参数：
* kd_node* kd_root - 根节点
* Key_Point* feat  - 特征点
* int k            - 所求近邻数
* Key_Point*** nbrs- 存储特征点的近邻，按距离增序排列
* int max_nn_chks  - 最大的搜索次数
* int d            - 描述子的大小
*
* 返回值：
* 特征点的kd树
*  
*说明：建立特征点的kd树
********************************************************/
int kdtree_bbf_knn( Kd_Node* kd_root, Key_Point* feat, int k,
				   Key_Point*** nbrs, int max_nn_chks ,int d)
{
	Kd_Node* expl;
	Min_Pq* min_pq;
	Key_Point* tree_feat, ** _nbrs;
	BBF_Data* bbf_data;
	int i, t = 0, n = 0;

	if( !nbrs || !feat || !kd_root )
	{
		/*fprintf( stderr, "Warning: NULL pointer error, %s, line %d\n",
			__FILE__, __LINE__ );*/
		return -1;
	}

	_nbrs =new Key_Point* [k]; 
	min_pq = minpq_init();
	minpq_insert( min_pq, kd_root, 0 );
	//队列里有数据就继续搜，同时控制在t<max_nn_chks 200（即200步内）
	while( min_pq->n > 0  &&  t < max_nn_chks )
	{
		expl = (Kd_Node*)minpq_extract_min( min_pq );//取出最小的，front & pop
		if( ! expl )
		{
			fprintf( stderr, "Warning: PQ unexpectedly empty, %s line %d\n",
				__FILE__, __LINE__ );
			goto fail;
		}

		expl = explore_to_leaf( expl, feat, min_pq ,d);//从该点开始，explore到leaf，路过的“有意义的点”就塞到最小队列min_pq中。
		if( ! expl )
		{
			fprintf( stderr, "Warning: PQ unexpectedly empty, %s line %d\n",
				__FILE__, __LINE__ );
			goto fail;
		}

		for( i = 0; i < expl->n; i++ )
		{
			tree_feat = &expl->features[i];
			bbf_data = new BBF_Data[1];
			if( !bbf_data )
			{
				fprintf( stderr, "Warning: unable to allocate memory,"
					" %s line %d\n", __FILE__, __LINE__ );
				goto fail;
			}
			bbf_data->old_data = tree_feat->feature_data;
			bbf_data->d = descr_dist_sq(feat, tree_feat);//两特征距离
			tree_feat->feature_data = bbf_data;
			 //按从小到大塞到neighbor数组里，到时候取前k个就是KNN， n 每次加1或0，表示目前已有的元素个数
			n += insert_into_nbr_array( tree_feat, _nbrs, n, k );
		}
		t++;
	}

	minpq_release( &min_pq );
	for( i = 0; i < n; i++ )
	{
		bbf_data =(BBF_Data*)_nbrs[i]->feature_data;
		_nbrs[i]->feature_data = bbf_data->old_data;
		delete []bbf_data;
	}
	*nbrs = _nbrs;
	return n;

fail:
	minpq_release( &min_pq );
	for( i = 0; i < n; i++ )
	{
		bbf_data =(BBF_Data*)_nbrs[i]->feature_data;
		_nbrs[i]->feature_data = bbf_data->old_data;
		free( bbf_data );
	}
	free( _nbrs );
	*nbrs = NULL;
	return -1;
}


/****************************************************** 
* 函数名称： 
* kd_node_init()
*  
* 函数参数：
* struct Key_Point* features - 特征点
* int n                      - 特征点的个数
*
* 返回值：
* kd树的一个节点
*  
*说明：初始化kd树的一个节点
********************************************************/
Kd_Node* kd_node_init( Key_Point* features, int n)
{
	Kd_Node* kd_node;

	kd_node =new struct Kd_Node [1];
	memset( kd_node, 0, sizeof( Kd_Node ) );//初始化
	kd_node->ki = -1;
	kd_node->features = features;
	kd_node->n = n;

	return kd_node;
}
/****************************************************** 
* 函数名称： 
* expand_kd_node_subtree()
*  
* 函数参数：
* kd_node* kd_node - kd树节点
* int d            - 向量的维数 128
*
* 返回值：
* 
*  
*说明：以kd_node为根节点扩展子树，建立kd树(递归的过程)
********************************************************/
void expand_kd_node_subtree( Kd_Node* kd_node,int d )
{
	//判断是否为叶子节点
	if( kd_node->n == 1  ||  kd_node->n == 0 )
	{
		kd_node->leaf = 1;
		return;
	}

	assign_part_key( kd_node,d );//找方差最大的维数及最大维数对应的中值
	partition_features( kd_node );//超平面分割

	if( kd_node->kd_left )
		expand_kd_node_subtree( kd_node->kd_left,d );
	if( kd_node->kd_right )
		expand_kd_node_subtree( kd_node->kd_right,d );
}

/****************************************************** 
* 函数名称： 
* partition_features()
*  
* 函数参数：
* kd_node* kd_node - kd树节点
*
*
* 返回值：
* 
*  
*说明：以ki为维度，kv为中值二分特征点
********************************************************/
void partition_features( Kd_Node* kd_node )
{
	Key_Point* features, tmp;
	float kv;
	int n, ki, p, i, j = -1;

	features = kd_node->features;
	n = kd_node->n;
	ki = kd_node->ki;
	kv = kd_node->kv;
	for( i = 0; i < n; i++ )
	{
		if( features[i].descriptor[ki] <= kv )
		{
			tmp = features[++j];
			features[j] = features[i];
			features[i] = tmp;
			if( features[j].descriptor[ki] == kv )
				p = j;//确定根节点所在的序号
		}
	}
	tmp = features[p];
	features[p] = features[j];
	features[j] = tmp;

	//如果所有的点都在同一侧则该节点为叶子节点
	if( j == n - 1 )
	{
		kd_node->leaf = 1;
		return;
	}

	kd_node->kd_left = kd_node_init( features, j + 1 );
	kd_node->kd_right = kd_node_init( features + ( j + 1 ), ( n - j - 1 ) );
}
/****************************************************** 
* 函数名称： 
* assign_part_key()
*  
* 函数参数：
* kd_node* kd_node - kd树节点
* int d  -           描述子的维数 128
*
* 返回值：
* 
*  
*说明：找方差最大的维数（ki）及最大维数对应的中值（kv）
********************************************************/
void assign_part_key( Kd_Node* kd_node,int d )
{
	Key_Point* features;
	float kv, x, mean, var, var_max = 0;
	int n, ki = 0;

	features = kd_node->features;
	n = kd_node->n;  //特征点的个数

	// 找最大方差所在的维度
	for(int j = 0; j < d; j++ )
	{
		mean = var = 0;
		for(int i = 0; i < n; i++ )
			mean += features[i].descriptor[j];
		mean /= n;
		for(int i = 0; i < n; i++ )
		{
			x = features[i].descriptor[j] - mean;
			var += x * x;
		}
		var /= n;

		if( var > var_max )
		{
			ki = j;
			var_max = var;
		}
	}

	// 以方差最大的维度的中值为超平面分割点
	float* tmp = new float [n];
	for(int i = 0; i < n; i++ )
		tmp[i] = features[i].descriptor[ki];
	kv = median_select( tmp, n );
	delete []tmp;

	kd_node->ki = ki;
	kd_node->kv = kv;
}
/****************************************************** 
* 函数名称： 
* median_select()
*  
* 函数参数：
* float* array - n个特征点的ki维数组成的数组
* int n  -       特征点个数
*
* 返回值：
* 中值
*  
*说明：找最大维数对应的中值（kv）
********************************************************/
float median_select( float* array, int n )
{
	return rank_select( array, n, (n - 1) / 2 );
}

/****************************************************** 
* 函数名称： 
* rank_select()
*  
* 函数参数：
* float* array - n个特征点的ki维数组成的数组
* int n  -       特征点个数
* int r  -       数组排序后待查找的值所在的位置
* 返回值：
* 中值
*  
*说明：排序找中值
********************************************************/
float rank_select( float* array, int n, int r )
{
	//只有一个数据
	if( n == 1 )
		return array[0];

	//将数据array每5个数组成一组并排序 groups of 5
	int gr_5 = n / 5;
	int gr_tot = int(ceil( n / 5.0 ));
	int rem_elts = n % 5;
	float *tmp = array;
	for(int i = 0; i < gr_5; i++ )
	{
		insertion_sort( tmp, 5 );
		tmp += 5;
	}
	insertion_sort( tmp, rem_elts );

	// 对以上排好序的gr_tot组，取出每组的中值组成新的数组 
	tmp = new float [gr_tot];
	int i,j;
	for( i = 0, j = 2; i < gr_5; i++, j += 5 )
		tmp[i] = array[j];
	if( rem_elts )
		tmp[i++] = array[n - 1 - rem_elts/2];

	// 对每组的中值组成新的数组，再次进行排序、选中值
	float med;
	med = rank_select( tmp, i, ( i - 1 ) / 2 );
	delete []tmp;

	/* partition around median of medians and recursively select if necessary */
	
	j = partition_array( array, n, med );
	if( r == j )
		return med;
	else if( r < j )
		return rank_select( array, j, r );
	else
	{
		array += j+1;
		return rank_select( array, ( n - j - 1 ), ( r - j - 1 ) );
	}
}

/****************************************************** 
* 函数名称： 
* insertion_sort()
*  
* 函数参数：
* float* array - 待排序的数组
* int n  -       数据个数
* 
* 返回值：
* 
*  
*说明：用差值的方法进行由小到大排序
********************************************************/
void insertion_sort( float* array, int n )
{
	for(int i = 1; i < n; i++ )
	{	
		float tmp = array[i];
		int j = i-1;
		while( j >= 0  &&  array[j] > tmp )
		{
			array[j+1] = array[j];
			j -= 1;
		}
		array[j+1] = tmp;
	}
}

/****************************************************** 
* 函数名称： 
* partition_array()
*  
* 函数参数：
* float* array - 待分割的数组
* int n  -       数据个数
* float pivot -  中值
*
* 返回值：
* 中值在数组中的序号
*  
*说明：以pivot为中值对array进行分割
********************************************************/
int partition_array( float* array, int n, float pivot )
{
	float tmp;
	int p, i, j;

	i = -1;
	for( j = 0; j < n; j++ )
	{
		if( array[j] <= pivot )
		{
			tmp = array[++i];
			array[i] = array[j];
			array[j] = tmp;
			if( array[i] == pivot )
				p = i;
		}
	}
	array[p] = array[i];
	array[i] = pivot;
	return i;
}

/****************************************************** 
* 函数名称： 
* explore_to_leaf()
*  
* 函数参数：
* Kd_Node* kd_node - 待探测子树的根节点
* Key_Point* feat  - 数据个数
* Min_Pq* min_pq   - 极小值的优先队列，用于存放满足条件的树的节点
* int d            - 描述子的大小，128
*
* 返回值：
* 最优队列
*  
*说明：探测并存储最优队列
********************************************************/
Kd_Node* explore_to_leaf( Kd_Node* kd_node, Key_Point* feat, Min_Pq* min_pq ,int d )
{
	Kd_Node* unexpl, * expl = kd_node;
	float kv;
	int ki;

	while( expl  &&  ! expl->leaf )
	{
		ki = expl->ki;
		kv = expl->kv;

		if( ki >= d )
		{
			fprintf( stderr, "Warning: comparing imcompatible descriptors, %s" \
				" line %d\n", __FILE__, __LINE__ );
			return NULL;
		}
		if( feat->descriptor[ki] <= kv )
		{
			unexpl = expl->kd_right;
			expl = expl->kd_left;
		}
		else
		{
			unexpl = expl->kd_left;
			expl = expl->kd_right;
		}

		if( minpq_insert( min_pq, unexpl, abs(int( kv - feat->descriptor[ki] )) ) )
		{
			fprintf( stderr, "Warning: unable to insert into PQ, %s, line %d\n",
				__FILE__, __LINE__ );
			return NULL;
		}
	}

	return expl;
}

/****************************************************** 
* 函数名称： 
* insert_into_nbr_array()
*  
* 函数参数：
* Key_Point* feat  - 待nbr中的特征，它的feature_data指针指向bbf_data
* Key_Point** nbrs - 最近邻数组
* int n            - nbrs中已有的点的个数
* int k            - nbr中最大容纳点的个数
*
* 返回值：
* 成功插入到最近邻数组中返回1否则返回0
*  
*说明：将特征插入到最近邻数组中，并保证数组按照与搜索特征的距离由小到大的排列
********************************************************/
int insert_into_nbr_array( Key_Point* feat, Key_Point** nbrs,
								  int n, int k )
{
	BBF_Data* fdata, * ndata;
	float dn, df;
	int i, ret = 0;

	if( n == 0 )
	{
		nbrs[0] = feat;
		return 1;
	}

	/* check at end of array */
	fdata = (BBF_Data*)feat->feature_data;
	df = fdata->d;
	ndata = (BBF_Data*)nbrs[n-1]->feature_data;
	dn = ndata->d;
	if( df >= dn )
	{
		if( n == k )
		{
			feat->feature_data = fdata->old_data;
			free( fdata );
			return 0;
		}
		nbrs[n] = feat;
		return 1;
	}

	/* find the right place in the array */
	if( n < k )
	{
		nbrs[n] = nbrs[n-1];
		ret = 1;
	}
	else
	{
		nbrs[n-1]->feature_data = ndata->old_data;
		free( ndata );
	}
	i = n-2;
	while( i >= 0 )
	{
		ndata = (BBF_Data*)nbrs[i]->feature_data;
		dn = ndata->d;
		if( dn <= df )
			break;
		nbrs[i+1] = nbrs[i];
		i--;
	}
	i++;
	nbrs[i] = feat;

	return ret;
}


//释放内存空间
void kdtree_release( Kd_Node* kd_root )
{
	if( ! kd_root )
		return;
	kdtree_release( kd_root->kd_left );
	kdtree_release( kd_root->kd_right );
	free( kd_root );
}