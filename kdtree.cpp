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

/******************************* �������� *******************************/

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

/******************************** ������ **********************************/

/****************************************************** 
* �������ƣ� 
* kdtree_build()
*  
* ����������
* Key_Point* features - ������
* int n                      - ������ĸ���
* int d                      - ����������ά��
*
* ����ֵ��
* �������kd��
*  
*˵���������������kd��
********************************************************/
Kd_Node* kdtree_build( Key_Point* features, int n, int d)
{
	Kd_Node* kd_root;//�洢���ڵ�

	if( !features || n<= 0 )
	{
		/*cout<<"û��ƥ��㣡"<<endl;*/
		return NULL;
	}

	//��ʼ�����ڵ�
	kd_root = kd_node_init( features, n);
	//��kd_rootΪ���ڵ㽨��kd��
	expand_kd_node_subtree( kd_root,d );

	return kd_root;
}
/****************************************************** 
* �������ƣ� 
* kdtree_bbf_knn()
*  
* ����������
* kd_node* kd_root - ���ڵ�
* Key_Point* feat  - ������
* int k            - ���������
* Key_Point*** nbrs- �洢������Ľ��ڣ���������������
* int max_nn_chks  - ������������
* int d            - �����ӵĴ�С
*
* ����ֵ��
* �������kd��
*  
*˵���������������kd��
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
	//�����������ݾͼ����ѣ�ͬʱ������t<max_nn_chks 200����200���ڣ�
	while( min_pq->n > 0  &&  t < max_nn_chks )
	{
		expl = (Kd_Node*)minpq_extract_min( min_pq );//ȡ����С�ģ�front & pop
		if( ! expl )
		{
			fprintf( stderr, "Warning: PQ unexpectedly empty, %s line %d\n",
				__FILE__, __LINE__ );
			goto fail;
		}

		expl = explore_to_leaf( expl, feat, min_pq ,d);//�Ӹõ㿪ʼ��explore��leaf��·���ġ�������ĵ㡱��������С����min_pq�С�
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
			bbf_data->d = descr_dist_sq(feat, tree_feat);//����������
			tree_feat->feature_data = bbf_data;
			 //����С��������neighbor�������ʱ��ȡǰk������KNN�� n ÿ�μ�1��0����ʾĿǰ���е�Ԫ�ظ���
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
* �������ƣ� 
* kd_node_init()
*  
* ����������
* struct Key_Point* features - ������
* int n                      - ������ĸ���
*
* ����ֵ��
* kd����һ���ڵ�
*  
*˵������ʼ��kd����һ���ڵ�
********************************************************/
Kd_Node* kd_node_init( Key_Point* features, int n)
{
	Kd_Node* kd_node;

	kd_node =new struct Kd_Node [1];
	memset( kd_node, 0, sizeof( Kd_Node ) );//��ʼ��
	kd_node->ki = -1;
	kd_node->features = features;
	kd_node->n = n;

	return kd_node;
}
/****************************************************** 
* �������ƣ� 
* expand_kd_node_subtree()
*  
* ����������
* kd_node* kd_node - kd���ڵ�
* int d            - ������ά�� 128
*
* ����ֵ��
* 
*  
*˵������kd_nodeΪ���ڵ���չ����������kd��(�ݹ�Ĺ���)
********************************************************/
void expand_kd_node_subtree( Kd_Node* kd_node,int d )
{
	//�ж��Ƿ�ΪҶ�ӽڵ�
	if( kd_node->n == 1  ||  kd_node->n == 0 )
	{
		kd_node->leaf = 1;
		return;
	}

	assign_part_key( kd_node,d );//�ҷ�������ά�������ά����Ӧ����ֵ
	partition_features( kd_node );//��ƽ��ָ�

	if( kd_node->kd_left )
		expand_kd_node_subtree( kd_node->kd_left,d );
	if( kd_node->kd_right )
		expand_kd_node_subtree( kd_node->kd_right,d );
}

/****************************************************** 
* �������ƣ� 
* partition_features()
*  
* ����������
* kd_node* kd_node - kd���ڵ�
*
*
* ����ֵ��
* 
*  
*˵������kiΪά�ȣ�kvΪ��ֵ����������
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
				p = j;//ȷ�����ڵ����ڵ����
		}
	}
	tmp = features[p];
	features[p] = features[j];
	features[j] = tmp;

	//������еĵ㶼��ͬһ����ýڵ�ΪҶ�ӽڵ�
	if( j == n - 1 )
	{
		kd_node->leaf = 1;
		return;
	}

	kd_node->kd_left = kd_node_init( features, j + 1 );
	kd_node->kd_right = kd_node_init( features + ( j + 1 ), ( n - j - 1 ) );
}
/****************************************************** 
* �������ƣ� 
* assign_part_key()
*  
* ����������
* kd_node* kd_node - kd���ڵ�
* int d  -           �����ӵ�ά�� 128
*
* ����ֵ��
* 
*  
*˵�����ҷ�������ά����ki�������ά����Ӧ����ֵ��kv��
********************************************************/
void assign_part_key( Kd_Node* kd_node,int d )
{
	Key_Point* features;
	float kv, x, mean, var, var_max = 0;
	int n, ki = 0;

	features = kd_node->features;
	n = kd_node->n;  //������ĸ���

	// ����󷽲����ڵ�ά��
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

	// �Է�������ά�ȵ���ֵΪ��ƽ��ָ��
	float* tmp = new float [n];
	for(int i = 0; i < n; i++ )
		tmp[i] = features[i].descriptor[ki];
	kv = median_select( tmp, n );
	delete []tmp;

	kd_node->ki = ki;
	kd_node->kv = kv;
}
/****************************************************** 
* �������ƣ� 
* median_select()
*  
* ����������
* float* array - n���������kiά����ɵ�����
* int n  -       ���������
*
* ����ֵ��
* ��ֵ
*  
*˵���������ά����Ӧ����ֵ��kv��
********************************************************/
float median_select( float* array, int n )
{
	return rank_select( array, n, (n - 1) / 2 );
}

/****************************************************** 
* �������ƣ� 
* rank_select()
*  
* ����������
* float* array - n���������kiά����ɵ�����
* int n  -       ���������
* int r  -       �������������ҵ�ֵ���ڵ�λ��
* ����ֵ��
* ��ֵ
*  
*˵������������ֵ
********************************************************/
float rank_select( float* array, int n, int r )
{
	//ֻ��һ������
	if( n == 1 )
		return array[0];

	//������arrayÿ5�������һ�鲢���� groups of 5
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

	// �������ź����gr_tot�飬ȡ��ÿ�����ֵ����µ����� 
	tmp = new float [gr_tot];
	int i,j;
	for( i = 0, j = 2; i < gr_5; i++, j += 5 )
		tmp[i] = array[j];
	if( rem_elts )
		tmp[i++] = array[n - 1 - rem_elts/2];

	// ��ÿ�����ֵ����µ����飬�ٴν�������ѡ��ֵ
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
* �������ƣ� 
* insertion_sort()
*  
* ����������
* float* array - �����������
* int n  -       ���ݸ���
* 
* ����ֵ��
* 
*  
*˵�����ò�ֵ�ķ���������С��������
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
* �������ƣ� 
* partition_array()
*  
* ����������
* float* array - ���ָ������
* int n  -       ���ݸ���
* float pivot -  ��ֵ
*
* ����ֵ��
* ��ֵ�������е����
*  
*˵������pivotΪ��ֵ��array���зָ�
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
* �������ƣ� 
* explore_to_leaf()
*  
* ����������
* Kd_Node* kd_node - ��̽�������ĸ��ڵ�
* Key_Point* feat  - ���ݸ���
* Min_Pq* min_pq   - ��Сֵ�����ȶ��У����ڴ���������������Ľڵ�
* int d            - �����ӵĴ�С��128
*
* ����ֵ��
* ���Ŷ���
*  
*˵����̽�Ⲣ�洢���Ŷ���
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
* �������ƣ� 
* insert_into_nbr_array()
*  
* ����������
* Key_Point* feat  - ��nbr�е�����������feature_dataָ��ָ��bbf_data
* Key_Point** nbrs - ���������
* int n            - nbrs�����еĵ�ĸ���
* int k            - nbr��������ɵ�ĸ���
*
* ����ֵ��
* �ɹ����뵽����������з���1���򷵻�0
*  
*˵�������������뵽����������У�����֤���鰴�������������ľ�����С���������
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


//�ͷ��ڴ�ռ�
void kdtree_release( Kd_Node* kd_root )
{
	if( ! kd_root )
		return;
	kdtree_release( kd_root->kd_left );
	kdtree_release( kd_root->kd_right );
	free( kd_root );
}