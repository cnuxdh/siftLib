#ifndef WALLISFILTER_H
#define WALLISFILTER_H

#define GRIDWINDOWSIZE 4
#define FILTERWINDOWSIZE GRIDWINDOWSIZE/**2+1*/
#define VALUEB 0.7f
#define VALUEC 0.8f
#define VALUEMEAN 120.0f
#define VALUESIGMA 60.0f

unsigned char* wallisfilter(unsigned char* grayimgdata,int width,int height);

#endif