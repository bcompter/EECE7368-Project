/*********************************************************************
 * pyramid.h
 *********************************************************************/

#ifndef _PYRAMID_H_
#define _PYRAMID_H_

#include "klt_util.h"

typedef struct  {
  int subsampling;
  int nLevels;
  _KLT_FloatImage *img;
  int *ncols, *nrows;
}  _KLT_PyramidRec, *_KLT_Pyramid;

#endif