#ifndef OUTPUT_H_
#define OUTPUT_H_

#include "siftmatch.h"
#include <vector>

extern int output_point( MatchPoints* matchpoints,ImageFile* imagefile,int image_num,int pointnum,double **trans_h );
extern int readtxt( string abspath,int &num,vector<string> &vect );
extern  int cal_accuracy_error( MatchPoints* matchpoints,double** trans_h,int img_num );
#endif
