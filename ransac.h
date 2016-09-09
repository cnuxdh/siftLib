



#ifndef RANSAC_H_
#define RANSAC_H_

#include "siftmatch.h"


DLL_EXPORT MatchPoints del_gross_error(MatchPoints matchpoints,int imgwidth,int imgheight);
extern double** inline_cal_trans(MatchPoint* pot,int &num);     


#endif



