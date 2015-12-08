// Design
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>       /* log */
#include <assert.h>
#include <sim.sh>
#include <string.h>

#include "main.sh"

import "i_receiver";

#include "klt_util.h"
#include "klt.h"
#include "convolve.h"
#include "pyramid.h"

#define SWAP3(list, i, j)               \
{register int *pi, *pj, tmp;            \
     pi=list+3*(i); pj=list+3*(j);      \
                                        \
     tmp=*pi;    \
     *pi++=*pj;  \
     *pj++=tmp;  \
                 \
     tmp=*pi;    \
     *pi++=*pj;  \
     *pj++=tmp;  \
                 \
     tmp=*pi;    \
     *pi=*pj;    \
     *pj=tmp;    \
}

behavior Design(i_receiver bytesFromStimulus) 
{    

	typedef enum {SELECTING_ALL, REPLACING_SOME} selectionMode;

 const int mindist = 10;
 const int window_size = 3;
 const int min_eigenvalue = 1;
 const float min_determinant = 0.01f;
 const float min_displacement = 0.1f;
 const int max_iterations = 10;
 const float max_residue = 10.0f;
 const float grad_sigma = 1.0f;
 const float smooth_sigma_fact = 0.1f;
 const float pyramid_sigma_fact = 0.9f;
 const float step_factor = 1.0f;
 const KLT_BOOL sequentialMode = FALSE;
 const KLT_BOOL lighting_insensitive = FALSE;

/* for affine mapping*/
 const int affineConsistencyCheck = -1;
 const int affine_window_size = 15;
 const int affine_max_iterations = 10;
 const float affine_max_residue = 10.0;
 const float affine_min_displacement = 0.02f;
 const float affine_max_displacement_differ = 1.5f;

 const KLT_BOOL smoothBeforeSelecting = TRUE;
 const KLT_BOOL writeInternalImages = FALSE;
 const int search_range = 15;
 const int nSkippedPixels = 0;
 
 /* GLOBAL VARIABLES */
 float sigma_last;
 
 /* Kernels */
 ConvolutionKernel gauss_kernel;
 ConvolutionKernel gaussderiv_kernel;


/*********************************************************************
 * KLTError
 * 
 * Prints an error message and dies.
 * 
 * INPUTS
 * exactly like printf
 */

void KLTError(char *fmt, ...)
{
  printf("*** ERROR...\n");
  exit(1);
}


/*********************************************************************
 * KLTWarning
 * 
 * Prints a warning message.
 * 
 * INPUTS
 * exactly like printf
 */

void KLTWarning(char *fmt, ...)
{
  printf("WARNING...\n");
}

/*********************************************************************/

void _fillFeaturemap(
  int x, int y, 
  uchar *featuremap, 
  int l_mindist, 
  int ncols, 
  int nrows)
{
  int ix, iy;

  for (iy = y - l_mindist ; iy <= y + l_mindist ; iy++)
    for (ix = x - l_mindist ; ix <= x + l_mindist ; ix++)
      if (ix >= 0 && ix < ncols && iy >= 0 && iy < nrows)
        featuremap[iy*ncols+ix] = 1;
}


/*********************************************************************
 * _enforceMinimumDistance
 *
 * Removes features that are within close proximity to better features.
 *
 * INPUTS
 * featurelist:  A list of features.  The nFeatures property
 *               is used.
 *
 * OUTPUTS
 * featurelist:  Is overwritten.  Nearby "redundant" features are removed.
 *               Writes -1's into the remaining elements.
 *
 * RETURNS
 * The number of remaining features.
 */

void _enforceMinimumDistance(
  int *pointlist,              /* featurepoints */
  int npoints,                 /* number of featurepoints */
  KLT_FeatureList featurelist, /* features */
  int ncols, int nrows,        /* size of images */
  int l_mindist,                 /* min. dist b/w features */
  int l_min_eigenvalue,          /* min. eigenvalue */
  KLT_BOOL overwriteAllFeatures)
{
  int indx;          /* Index into features */
  int x, y, val;     /* Location and trackability of pixel under consideration */
  uchar *featuremap; /* Boolean array recording proximity of features */
  int *ptr;
	
  /* Cannot add features with an eigenvalue less than one */
  if (l_min_eigenvalue < 1)  l_min_eigenvalue = 1;

  /* Allocate memory for feature map and clear it */
  featuremap = (uchar *) malloc(ncols * nrows * sizeof(uchar));
  memset(featuremap, 0, ncols*nrows);
	
  /* Necessary because code below works with (mindist-1) */
  l_mindist--;

  /* If we are keeping all old good features, then add them to the featuremap */
  if (!overwriteAllFeatures)
    for (indx = 0 ; indx < featurelist->nFeatures ; indx++)
      if (featurelist->feature[indx]->val >= 0)  {
        x   = (int) featurelist->feature[indx]->x;
        y   = (int) featurelist->feature[indx]->y;
        _fillFeaturemap(x, y, featuremap, l_mindist, ncols, nrows);
      }

  /* For each feature point, in descending order of importance, do ... */
  ptr = pointlist;
  indx = 0;
  while (1)  {

    /* If we can't add all the points, then fill in the rest
       of the featurelist with -1's */
    if (ptr >= pointlist + 3*npoints)  {
      while (indx < featurelist->nFeatures)  {	
        if (overwriteAllFeatures || 
            featurelist->feature[indx]->val < 0) {
          featurelist->feature[indx]->x   = -1;
          featurelist->feature[indx]->y   = -1;
          featurelist->feature[indx]->val = KLT_NOT_FOUND;
	  featurelist->feature[indx]->aff_img = 0;
	  featurelist->feature[indx]->aff_img_gradx = 0;
	  featurelist->feature[indx]->aff_img_grady = 0;
	  featurelist->feature[indx]->aff_x = -1.0;
	  featurelist->feature[indx]->aff_y = -1.0;
	  featurelist->feature[indx]->aff_Axx = 1.0;
	  featurelist->feature[indx]->aff_Ayx = 0.0;
	  featurelist->feature[indx]->aff_Axy = 0.0;
	  featurelist->feature[indx]->aff_Ayy = 1.0;
        }
        indx++;
      }
      break;
    }

    x   = *ptr++;
    y   = *ptr++;
    val = *ptr++;
		
    /* Ensure that feature is in-bounds */
    assert(x >= 0);
    assert(x < ncols);
    assert(y >= 0);
    assert(y < nrows);
	
    while (!overwriteAllFeatures && 
           indx < featurelist->nFeatures &&
           featurelist->feature[indx]->val >= 0)
      indx++;

    if (indx >= featurelist->nFeatures)  break;

    /* If no neighbor has been selected, and if the minimum
       eigenvalue is large enough, then add feature to the current list */
    if (!featuremap[y*ncols+x] && val >= l_min_eigenvalue)  {
      featurelist->feature[indx]->x   = (KLT_locType) x;
      featurelist->feature[indx]->y   = (KLT_locType) y;
      featurelist->feature[indx]->val = (int) val;
      featurelist->feature[indx]->aff_img = 0;
      featurelist->feature[indx]->aff_img_gradx = 0;
      featurelist->feature[indx]->aff_img_grady = 0;
      featurelist->feature[indx]->aff_x = -1.0;
      featurelist->feature[indx]->aff_y = -1.0;
      featurelist->feature[indx]->aff_Axx = 1.0;
      featurelist->feature[indx]->aff_Ayx = 0.0;
      featurelist->feature[indx]->aff_Axy = 0.0;
      featurelist->feature[indx]->aff_Ayy = 1.0;
      indx++;

      /* Fill in surrounding region of feature map, but
         make sure that pixels are in-bounds */
      _fillFeaturemap(x, y, featuremap, l_mindist, ncols, nrows);
    }
  }

  /* Free feature map  */
  free(featuremap);
}

void _quicksort(int *pointlist, int n)
{
  unsigned int i, j, ln, rn;

  while (n > 1)
  {
    SWAP3(pointlist, 0, n/2);
    for (i = 0, j = n; ; )
    {
      do
        --j;
      while (pointlist[3*j+2] < pointlist[2]);
      do
        ++i;
      while (i < j && pointlist[3*i+2] > pointlist[2]);
      if (i >= j)
        break;
      SWAP3(pointlist, i, j);
    }
    SWAP3(pointlist, j, 0);
    ln = j;
    rn = n - ++j;
    if (ln < rn)
    {
      _quicksort(pointlist, ln);
      pointlist += 3*j;
      n = rn;
    }
    else
    {
      _quicksort(pointlist + 3*j, rn);
      n = ln;
    }
  }
}

/*********************************************************************
 * _sortPointList
 */

void _sortPointList(
  int *pointlist,
  int npoints)
{
#ifdef KLT_USE_QSORT
  qsort(pointlist, npoints, 3*sizeof(int), _comparePoints);
#else
  _quicksort(pointlist, npoints);
#endif
}

/*********************************************************************
 * _minEigenvalue
 *
 * Given the three distinct elements of the symmetric 2x2 matrix
 *                     [gxx gxy]
 *                     [gxy gyy],
 * Returns the minimum eigenvalue of the matrix.  
 */

float _minEigenvalue(float gxx, float gxy, float gyy)
{
  return (float) ((gxx + gyy - sqrt((gxx - gyy)*(gxx - gyy) + 4*gxy*gxy))/2.0f);
}

/*********************************************************************
 * _KLTComputeSmoothSigma
 */
float _KLTComputeSmoothSigma(
  KLT_TrackingContext* tc)
{
  return (tc->smooth_sigma_fact * max(tc->window_width, tc->window_height));
}

/*********************************************************************
 * _computeKernels
 */

void _computeKernels(
  float sigma,
  ConvolutionKernel *gauss,
  ConvolutionKernel *gaussderiv)
{
  const float factor = 0.01f;   /* for truncating tail */
  int hw = MAX_KERNEL_WIDTH / 2;
  float max_gauss, max_gaussderiv;
  float den;
  int i;

  assert(MAX_KERNEL_WIDTH % 2 == 1);
  assert(sigma >= 0.0);

  /* Compute kernels, and automatically determine widths */
  {
    
    max_gauss = 1.0f;
    max_gaussderiv = (float) (sigma*exp(-0.5f));
	
    /* Compute gauss and deriv */
    for (i = -hw ; i <= hw ; i++)  {
      gauss->data[i+hw]      = (float) exp(-i*i / (2*sigma*sigma));
      gaussderiv->data[i+hw] = -i * gauss->data[i+hw];
    }

    /* Compute widths */
    gauss->width = MAX_KERNEL_WIDTH;
    for (i = -hw ; fabs(gauss->data[i+hw] / max_gauss) < factor ; 
         i++, gauss->width -= 2);
    gaussderiv->width = MAX_KERNEL_WIDTH;
    for (i = -hw ; fabs(gaussderiv->data[i+hw] / max_gaussderiv) < factor ; 
         i++, gaussderiv->width -= 2);
    if (gauss->width == MAX_KERNEL_WIDTH || 
        gaussderiv->width == MAX_KERNEL_WIDTH)
      KLTError("(_computeKernels) MAX_KERNEL_WIDTH %d is too small for "
               "a sigma of %f", MAX_KERNEL_WIDTH, sigma);
  }

  /* Shift if width less than MAX_KERNEL_WIDTH */
  for (i = 0 ; i < gauss->width ; i++)
    gauss->data[i] = gauss->data[i+(MAX_KERNEL_WIDTH-gauss->width)/2];
  for (i = 0 ; i < gaussderiv->width ; i++)
    gaussderiv->data[i] = gaussderiv->data[i+(MAX_KERNEL_WIDTH-gaussderiv->width)/2];
  /* Normalize gauss and deriv */
  {	
    hw = gaussderiv->width / 2;
    den = 0.0;
    for (i = 0 ; i < gauss->width ; i++)  den += gauss->data[i];
    for (i = 0 ; i < gauss->width ; i++)  gauss->data[i] /= den;
    den = 0.0;
    for (i = -hw ; i <= hw ; i++)  den -= i*gaussderiv->data[i+hw];
    for (i = -hw ; i <= hw ; i++)  gaussderiv->data[i+hw] /= den;
  }

  sigma_last = sigma;
}

/*********************************************************************
 * _KLTFreeFloatImage
 */

void _KLTFreeFloatImage(
  _KLT_FloatImage floatimg)
{
  free(floatimg);
}

/*********************************************************************
 * _KLTCreateFloatImage
 */

_KLT_FloatImage _KLTCreateFloatImage(
  int ncols,
  int nrows)
{
  _KLT_FloatImage floatimg;
  int nbytes;
  nbytes = sizeof(_KLT_FloatImageRec) +
    ncols * nrows * sizeof(float);

  floatimg = (_KLT_FloatImage)  malloc(nbytes);
  if (floatimg == NULL)
    KLTError("(_KLTCreateFloatImage)  Out of memory");
  floatimg->ncols = ncols;
  floatimg->nrows = nrows;
  floatimg->data = (float *)  (floatimg + 1);

  return(floatimg);
}

/*********************************************************************
 * _convolveImageHoriz
 */

void _convolveImageHoriz(
  _KLT_FloatImage imgin,
  ConvolutionKernel kernel,
  _KLT_FloatImage imgout)
{
  float *ptrrow;           /* Points to row's first pixel */
  register float *ptrout, /* Points to next output pixel */
    *ppp;
  register double sum;
  register int radius;
  register int ncols, nrows;
  register int i, j, k;

  ptrrow = imgin->data;           /* Points to row's first pixel */
  ptrout = imgout->data; /* Points to next output pixel */
  radius = kernel.width / 2;
  ncols = imgin->ncols;
  nrows = imgin->nrows;

  /* Kernel width must be odd */
  assert(kernel.width % 2 == 1);

  /* Must read from and write to different images */
  assert(imgin != imgout);

  /* Output image must be large enough to hold result */
  assert(imgout->ncols >= imgin->ncols);
  assert(imgout->nrows >= imgin->nrows);

  /* For each row, do ... */
  for (j = 0 ; j < nrows ; j++)  {

    /* Zero leftmost columns */
    for (i = 0 ; i < radius ; i++)
      *ptrout++ = 0.0;

    /* Convolve middle columns with kernel */
    for ( ; i < ncols - radius ; i++)  {
      ppp = ptrrow + i - radius;
      sum = 0.0;
      for (k = kernel.width-1 ; k >= 0 ; k--)
        sum += *ppp++ * kernel.data[k];
      *ptrout++ = sum;
    }

    /* Zero rightmost columns */
    for ( ; i < ncols ; i++)
      *ptrout++ = 0.0;

    ptrrow += ncols;
  }
}

/*********************************************************************
 * _convolveImageVert
 */

void _convolveImageVert(
  _KLT_FloatImage imgin,
  ConvolutionKernel kernel,
  _KLT_FloatImage imgout)
{
  float *ptrcol;            /* Points to row's first pixel */
  register float *ptrout,  /* Points to next output pixel */
    *ppp;
  register double sum;
  register int radius;
  register int ncols, nrows;
  register int i, j, k;

  ptrcol = imgin->data;            /* Points to row's first pixel */
  ptrout = imgout->data;  /* Points to next output pixel */
  radius = kernel.width / 2;
  ncols = imgin->ncols, nrows = imgin->nrows;

  /* Kernel width must be odd */
  assert(kernel.width % 2 == 1);

  /* Must read from and write to different images */
  assert(imgin != imgout);

  /* Output image must be large enough to hold result */
  assert(imgout->ncols >= imgin->ncols);
  assert(imgout->nrows >= imgin->nrows);

  /* For each column, do ... */
  for (i = 0 ; i < ncols ; i++)  {

    /* Zero topmost rows */
    for (j = 0 ; j < radius ; j++)  {
      *ptrout = 0.0;
      ptrout += ncols;
    }

    /* Convolve middle rows with kernel */
    for ( ; j < nrows - radius ; j++)  {
      ppp = ptrcol + ncols * (j - radius);
      sum = 0.0;
      for (k = kernel.width-1 ; k >= 0 ; k--)  {
        sum += *ppp * kernel.data[k];
        ppp += ncols;
      }
      *ptrout = sum;
      ptrout += ncols;
    }

    /* Zero bottommost rows */
    for ( ; j < nrows ; j++)  {
      *ptrout = 0.0;
      ptrout += ncols;
    }

    ptrcol++;
    ptrout -= nrows * ncols - 1;
  }
}

/*********************************************************************
 * _convolveSeparate
 */

void _convolveSeparate(
  _KLT_FloatImage imgin,
  ConvolutionKernel horiz_kernel,
  ConvolutionKernel vert_kernel,
  _KLT_FloatImage imgout)
{
  /* Create temporary image */
  _KLT_FloatImage tmpimg;
  tmpimg = _KLTCreateFloatImage(imgin->ncols, imgin->nrows);
  
  /* Do convolution */
  _convolveImageHoriz(imgin, horiz_kernel, tmpimg);
  _convolveImageVert(tmpimg, vert_kernel, imgout);

  /* Free memory */
  _KLTFreeFloatImage(tmpimg);
}

	
/*********************************************************************
 * _KLTComputeGradients
 */

void _KLTComputeGradients(
  _KLT_FloatImage img,
  float sigma,
  _KLT_FloatImage gradx,
  _KLT_FloatImage grady)
{
				
  /* Output images must be large enough to hold result */
  assert(gradx->ncols >= img->ncols);
  assert(gradx->nrows >= img->nrows);
  assert(grady->ncols >= img->ncols);
  assert(grady->nrows >= img->nrows);

  /* Compute kernels, if necessary */
  if (fabs(sigma - sigma_last) > 0.05)
    _computeKernels(sigma, &gauss_kernel, &gaussderiv_kernel);
	
  _convolveSeparate(img, gaussderiv_kernel, gauss_kernel, gradx);
  _convolveSeparate(img, gauss_kernel, gaussderiv_kernel, grady);

}
	

/*********************************************************************
 * _KLTComputeSmoothedImage
 */

void _KLTComputeSmoothedImage(
  _KLT_FloatImage img,
  float sigma,
  _KLT_FloatImage smooth)
{
  /* Output image must be large enough to hold result */
  assert(smooth->ncols >= img->ncols);
  assert(smooth->nrows >= img->nrows);

  /* Compute kernel, if necessary; gauss_deriv is not used */
  if (fabs(sigma - sigma_last) > 0.05)
    _computeKernels(sigma, &gauss_kernel, &gaussderiv_kernel);

  _convolveSeparate(img, gauss_kernel, gauss_kernel, smooth);
}

/*********************************************************************
 * _KLTToFloatImage
 *
 * Given a pointer to image data (probably unsigned chars), copy
 * data to a float image.
 */

void _KLTToFloatImage(
  KLT_PixelType *img,
  int ncols, int nrows,
  _KLT_FloatImage floatimg)
{
  KLT_PixelType *ptrend;
  float *ptrout;

  ptrend = img + ncols*nrows;
  ptrout = floatimg->data;

  /* Output image must be large enough to hold result */
  assert(floatimg->ncols >= ncols);
  assert(floatimg->nrows >= nrows);

  floatimg->ncols = ncols;
  floatimg->nrows = nrows;

  while (img < ptrend)  *ptrout++ = (float) *img++;
}



/*********************************************************************/

void _KLTSelectGoodFeatures(
  KLT_TrackingContext *tc,
  KLT_PixelType *img, 
  int ncols, 
  int nrows,
  KLT_FeatureList featurelist,
  selectionMode mode)
{
  _KLT_FloatImage floatimg, gradx, grady;
  int window_hw, window_hh;
  int *pointlist;
  int npoints;
  KLT_BOOL overwriteAllFeatures;
  KLT_BOOL floatimages_created;
  
  npoints = 0;
  overwriteAllFeatures = (mode == SELECTING_ALL) ? TRUE : FALSE;
  floatimages_created = FALSE;

  /* Check window size (and correct if necessary) */
  if (tc->window_width % 2 != 1) {
    tc->window_width = tc->window_width+1;
    KLTWarning("Tracking context's window width must be odd.  "
               "Changing to %d.\n", tc->window_width);
  }
  if (tc->window_height % 2 != 1) {
    tc->window_height = tc->window_height+1;
    KLTWarning("Tracking context's window height must be odd.  "
               "Changing to %d.\n", tc->window_height);
  }
  if (tc->window_width < 3) {
    tc->window_width = 3;
    KLTWarning("Tracking context's window width must be at least three.  \n"
               "Changing to %d.\n", tc->window_width);
  }
  if (tc->window_height < 3) {
    tc->window_height = 3;
    KLTWarning("Tracking context's window height must be at least three.  \n"
               "Changing to %d.\n", tc->window_height);
  }
  window_hw = tc->window_width/2; 
  window_hh = tc->window_height/2;
		
  /* Create pointlist, which is a simplified version of a featurelist, */
  /* for speed.  Contains only integer locations and values. */
  pointlist = (int *) malloc(ncols * nrows * 3 * sizeof(int));

  /* Create temporary images, etc. */
  if (mode == REPLACING_SOME && 
      tc->sequentialMode && tc->pyramid_last != NULL)  {
    floatimg = ((_KLT_Pyramid) tc->pyramid_last)->img[0];
    gradx = ((_KLT_Pyramid) tc->pyramid_last_gradx)->img[0];
    grady = ((_KLT_Pyramid) tc->pyramid_last_grady)->img[0];
    assert(gradx != NULL);
    assert(grady != NULL);
  } else  {
    floatimages_created = TRUE;
    floatimg = _KLTCreateFloatImage(ncols, nrows);
    gradx    = _KLTCreateFloatImage(ncols, nrows);
    grady    = _KLTCreateFloatImage(ncols, nrows);
    if (tc->smoothBeforeSelecting)  {
      _KLT_FloatImage tmpimg;
      tmpimg = _KLTCreateFloatImage(ncols, nrows);
      _KLTToFloatImage(img, ncols, nrows, tmpimg);
      _KLTComputeSmoothedImage(tmpimg, _KLTComputeSmoothSigma(tc), floatimg);
      _KLTFreeFloatImage(tmpimg);
    } else _KLTToFloatImage(img, ncols, nrows, floatimg);
 
    /* Compute gradient of image in x and y direction */
    _KLTComputeGradients(floatimg, tc->grad_sigma, gradx, grady);
  }
	
  /* Write internal images */
  if (tc->writeInternalImages)  {
    //_KLTWriteFloatImageToPGM(floatimg, "kltimg_sgfrlf.pgm");
    //_KLTWriteFloatImageToPGM(gradx, "kltimg_sgfrlf_gx.pgm");
    //_KLTWriteFloatImageToPGM(grady, "kltimg_sgfrlf_gy.pgm");
  }

  /* Compute trackability of each image pixel as the minimum
     of the two eigenvalues of the Z matrix */
  {
    register float gx, gy;
    register float gxx, gxy, gyy;
    register int xx, yy;
    register int *ptr;
    float val;
    unsigned int limit;
    int x, y;
    int i;
    int borderx;	/* Must not touch cols */
    int bordery;	/* lost by convolution */
    
    limit = 1;
    borderx = tc->borderx;	/* Must not touch cols */
    bordery = tc->bordery;	/* lost by convolution */
    
    if (borderx < window_hw)  borderx = window_hw;
    if (bordery < window_hh)  bordery = window_hh;

    /* Find largest value of an int */
    for (i = 0 ; i < sizeof(int) ; i++)  limit *= 256;
    limit = limit/2 - 1;
		
    /* For most of the pixels in the image, do ... */
    ptr = pointlist;
    for (y = bordery ; y < nrows - bordery ; y += tc->nSkippedPixels + 1)
      for (x = borderx ; x < ncols - borderx ; x += tc->nSkippedPixels + 1)  {

        /* Sum the gradients in the surrounding window */
        gxx = 0;  gxy = 0;  gyy = 0;
        for (yy = y-window_hh ; yy <= y+window_hh ; yy++)
          for (xx = x-window_hw ; xx <= x+window_hw ; xx++)  {
            gx = *(gradx->data + ncols*yy+xx);
            gy = *(grady->data + ncols*yy+xx);
            gxx += gx * gx;
            gxy += gx * gy;
            gyy += gy * gy;
          }

        /* Store the trackability of the pixel as the minimum
           of the two eigenvalues */
        *ptr++ = x;
        *ptr++ = y;
        val = _minEigenvalue(gxx, gxy, gyy);
        if (val > limit)  {
          KLTWarning("(_KLTSelectGoodFeatures) minimum eigenvalue %f is "
                     "greater than the capacity of an int; setting "
                     "to maximum value", val);
          val = (float) limit;
        }
        *ptr++ = (int) val;
        npoints++;
      }
  }
			
  /* Sort the features  */
  _sortPointList(pointlist, npoints);

  /* Check tc->mindist */
  if (tc->mindist < 0)  {
    KLTWarning("(_KLTSelectGoodFeatures) Tracking context field tc->mindist "
               "is negative (%d); setting to zero", tc->mindist);
    tc->mindist = 0;
  }

  /* Enforce minimum distance between features */
  _enforceMinimumDistance(
    pointlist,
    npoints,
    featurelist,
    ncols, nrows,
    tc->mindist,
    tc->min_eigenvalue,
    overwriteAllFeatures);

  /* Free memory */
  free(pointlist);
  if (floatimages_created)  {
    _KLTFreeFloatImage(floatimg);
    _KLTFreeFloatImage(gradx);
    _KLTFreeFloatImage(grady);
  }
}

/*********************************************************************
 * KLTSelectGoodFeatures
 *
 * Main routine, visible to the outside.  Finds the good features in
 * an image.  
 * 
 * INPUTS
 * tc:	Contains parameters used in computation (size of image,
 *        size of window, min distance b/w features, sigma to compute
 *        image gradients, # of features desired).
 * img:	Pointer to the data of an image (probably unsigned chars).
 * 
 * OUTPUTS
 * features:	List of features.  The member nFeatures is computed.
 */

void KLTSelectGoodFeatures(
  KLT_TrackingContext * tc,
  KLT_PixelType *img,
  int ncols, 
  int nrows,
  KLT_FeatureList fl
  )
{
  _KLTSelectGoodFeatures(tc, img, ncols, nrows, fl, SELECTING_ALL);
}

/*********************************************************************
 * _createArray2D
 *
 * Creates a two-dimensional array.
 *
 * INPUTS
 * ncols:      no. of columns
 * nrows:      no. of rows
 * nbytes:     no. of bytes per entry
 *
 * RETURNS
 * Pointer to an array.  Must be coerced.
 *
 * EXAMPLE
 * char **ar;
 * ar = (char **) createArray2D(8, 5, sizeof(char));
 */

void** _createArray2D(int ncols, int nrows, int nbytes)
{
  char **tt;
  int i;

  tt = (char **) malloc(nrows * sizeof(void *) +
                        ncols * nrows * nbytes);
  if (tt == NULL)
    KLTError("(createArray2D) Out of memory");

  for (i = 0 ; i < nrows ; i++)
    tt[i] = ((char *) tt) + (nrows * sizeof(void *) +
                             i * ncols * nbytes);

  return((void **) tt);
}

/*********************************************************************
 * KLTCreateFeatureTable
 *
 */

KLT_FeatureTable KLTCreateFeatureTable(
  int nFrames,
  int nFeatures)
{
  KLT_FeatureTable ft;
  KLT_Feature first;
  int nbytes;
  int i, j;
  
  nbytes = sizeof(KLT_FeatureTableRec);
	
  /* Allocate memory for feature history */
  ft = (KLT_FeatureTable)  malloc(nbytes);
	
  /* Set parameters */
  ft->nFrames = nFrames; 
  ft->nFeatures = nFeatures; 
	
  /* Set pointers */
  ft->feature = (KLT_Feature **) 
    _createArray2D(nFrames, nFeatures, sizeof(KLT_Feature));
  first = (KLT_Feature) malloc(nFrames * nFeatures * sizeof(KLT_FeatureRec));
  for (j = 0 ; j < nFeatures ; j++)
    for (i = 0 ; i < nFrames ; i++)
      ft->feature[j][i] = first + j*nFrames + i;

  /* Return feature table */
  return(ft);
}

/*********************************************************************
 * KLTCreateFeatureList
 *
 */

KLT_FeatureList KLTCreateFeatureList(int nFeatures)
{
  KLT_FeatureList fl;
  KLT_Feature first;
  int nbytes, i;
  nbytes = sizeof(KLT_FeatureListRec) +
    nFeatures * sizeof(KLT_Feature) +
    nFeatures * sizeof(KLT_FeatureRec);
	
  /* Allocate memory for feature list */
  fl = (KLT_FeatureList)  malloc(nbytes);
	
  /* Set parameters */
  fl->nFeatures = nFeatures; 

  /* Set pointers */
  fl->feature = (KLT_Feature *) (fl + 1);
  first = (KLT_Feature) (fl->feature + nFeatures);
  for (i = 0 ; i < nFeatures ; i++) {
    fl->feature[i] = first + i;
    fl->feature[i]->aff_img = 0;           /* initialization fixed by Sinisa Segvic */
    fl->feature[i]->aff_img_gradx = 0;
    fl->feature[i]->aff_img_grady = 0;
  }
  /* Return feature list */
  return(fl);
}


/**
 * _pyramidSigma
 */
float _pyramidSigma(
  KLT_TrackingContext* tc)
{
  return (tc->pyramid_sigma_fact * tc->subsampling);
}

/*********************************************************************
 * _KLTGetKernelWidths
 *
 */

void _KLTGetKernelWidths(
  float sigma,
  int *gauss_width,
  int *gaussderiv_width)
{
  _computeKernels(sigma, &gauss_kernel, &gaussderiv_kernel);
  *gauss_width = gauss_kernel.width;
  *gaussderiv_width = gaussderiv_kernel.width;
}

/*********************************************************************
 * Updates border, which is dependent upon 
 * smooth_sigma_fact, pyramid_sigma_fact, window_size, and subsampling
 */

void KLTUpdateTCBorder(
  KLT_TrackingContext * tc)
{
  float val;
  int pyramid_gauss_hw;
  int smooth_gauss_hw;
  int gauss_width, gaussderiv_width;
  int num_levels;
  int n_invalid_pixels;
  int window_hw;
  int ss;
  int ss_power;
  int border;
  int i;
  
  num_levels = tc->nPyramidLevels;
  ss = tc->subsampling;

  /* Check window size (and correct if necessary) */
  if (tc->window_width % 2 != 1) {
    tc->window_width = tc->window_width+1;
    KLTWarning("(KLTUpdateTCBorder) Window width must be odd.  "
               "Changing to %d.\n", tc->window_width);
  }
  if (tc->window_height % 2 != 1) {
    tc->window_height = tc->window_height+1;
    KLTWarning("(KLTUpdateTCBorder) Window height must be odd.  "
               "Changing to %d.\n", tc->window_height);
  }
  if (tc->window_width < 3) {
    tc->window_width = 3;
    KLTWarning("(KLTUpdateTCBorder) Window width must be at least three.  \n"
               "Changing to %d.\n", tc->window_width);
  }
  if (tc->window_height < 3) {
    tc->window_height = 3;
    KLTWarning("(KLTUpdateTCBorder) Window height must be at least three.  \n"
               "Changing to %d.\n", tc->window_height);
  }
  window_hw = max(tc->window_width, tc->window_height)/2;

  /* Find widths of convolution windows */
  _KLTGetKernelWidths(_KLTComputeSmoothSigma(tc),
                      &gauss_width, &gaussderiv_width);
  smooth_gauss_hw = gauss_width/2;
  _KLTGetKernelWidths(_pyramidSigma(tc),
                      &gauss_width, &gaussderiv_width);
  pyramid_gauss_hw = gauss_width/2;

  /* Compute the # of invalid pixels at each level of the pyramid.
     n_invalid_pixels is computed with respect to the ith level   
     of the pyramid.  So, e.g., if n_invalid_pixels = 5 after   
     the first iteration, then there are 5 invalid pixels in   
     level 1, which translated means 5*subsampling invalid pixels   
     in the original level 0. */
  n_invalid_pixels = smooth_gauss_hw;
  for (i = 1 ; i < num_levels ; i++)  {
    val = ((float) n_invalid_pixels + pyramid_gauss_hw) / ss;
    n_invalid_pixels = (int) (val + 0.99);  /* Round up */
  }

  /* ss_power = ss^(num_levels-1) */
  ss_power = 1;
  for (i = 1 ; i < num_levels ; i++)
    ss_power *= ss;

  /* Compute border by translating invalid pixels back into */
  /* original image */
  border = (n_invalid_pixels + window_hw) * ss_power;

  tc->borderx = border;
  tc->bordery = border;
}

/*********************************************************************
 * KLTChangeTCPyramid
 *
 */

void KLTChangeTCPyramid(KLT_TrackingContext * tc)
{
  float window_halfwidth;
  float subsampling;
  float val;

  /* Check window size (and correct if necessary) */
  if (tc->window_width % 2 != 1) {
    tc->window_width = tc->window_width+1;
    KLTWarning("(KLTChangeTCPyramid) Window width must be odd.  "
               "Changing to %d.\n", tc->window_width);
  }
  if (tc->window_height % 2 != 1) {
    tc->window_height = tc->window_height+1;
    KLTWarning("(KLTChangeTCPyramid) Window height must be odd.  "
               "Changing to %d.\n", tc->window_height);
  }
  if (tc->window_width < 3) {
    tc->window_width = 3;
    KLTWarning("(KLTChangeTCPyramid) Window width must be at least three.  \n"
               "Changing to %d.\n", tc->window_width);
  }
  if (tc->window_height < 3) {
    tc->window_height = 3;
    KLTWarning("(KLTChangeTCPyramid) Window height must be at least three.  \n"
               "Changing to %d.\n", tc->window_height);
  }
  window_halfwidth = min(tc->window_width,tc->window_height)/2.0f;

  subsampling = ((float) search_range) / window_halfwidth;

  if (subsampling < 1.0)  {		/* 1.0 = 0+1 */
    tc->nPyramidLevels = 1;
  } else if (subsampling <= 3.0)  {	/* 3.0 = 2+1 */
    tc->nPyramidLevels = 2;
    tc->subsampling = 2;
  } else if (subsampling <= 5.0)  {	/* 5.0 = 4+1 */
    tc->nPyramidLevels = 2;
    tc->subsampling = 4;
  } else if (subsampling <= 9.0)  {	/* 9.0 = 8+1 */
    tc->nPyramidLevels = 2;
    tc->subsampling = 8;
  } else {
    /* The following lines are derived from the formula:
       search_range = 
       window_halfwidth * \sum_{i=0}^{nPyramidLevels-1} 8^i,
       which is the same as:
       search_range = 
       window_halfwidth * (8^nPyramidLevels - 1)/(8 - 1).
       Then, the value is rounded up to the nearest integer. */
    val = (float) (log(7.0*subsampling+1.0)/log(8.0));
    tc->nPyramidLevels = (int) (val + 0.99);
    tc->subsampling = 8;
  }
  if (tc->nPyramidLevels > 2)
	  printf("Error! nPyramidLevels > 2, please change "
			  "the size statically allocated for pyramid!\n");
}


/*********************************************************************
 * KLTCreateTrackingContext
 *
 */
KLT_TrackingContext KLTCreateTrackingContext()
{
  KLT_TrackingContext tc;

  /* Set values to default values */
  tc.mindist = mindist;
  tc.window_width = window_size;
  tc.window_height = window_size;
  tc.sequentialMode = sequentialMode;
  tc.smoothBeforeSelecting = smoothBeforeSelecting;
  tc.writeInternalImages = writeInternalImages;
  tc.lighting_insensitive = lighting_insensitive;
  tc.min_eigenvalue = min_eigenvalue;
  tc.min_determinant = min_determinant;
  tc.max_iterations = max_iterations;
  tc.min_displacement = min_displacement;
  tc.max_residue = max_residue;
  tc.grad_sigma = grad_sigma;
  tc.smooth_sigma_fact = smooth_sigma_fact;
  tc.pyramid_sigma_fact = pyramid_sigma_fact;
  tc.step_factor = step_factor;
  tc.nSkippedPixels = nSkippedPixels;
  tc.pyramid_last = NULL;
  tc.pyramid_last_gradx = NULL;
  tc.pyramid_last_grady = NULL;
  
  /* for affine mapping */
  tc.affineConsistencyCheck = affineConsistencyCheck;
  tc.affine_window_width = affine_window_size;
  tc.affine_window_height = affine_window_size;
  tc.affine_max_iterations = affine_max_iterations;
  tc.affine_max_residue = affine_max_residue;
  tc.affine_min_displacement = affine_min_displacement;
  tc.affine_max_displacement_differ = affine_max_displacement_differ;

  /* Change nPyramidLevels and subsampling */
  KLTChangeTCPyramid(&tc);
	
  /* Update border, which is dependent upon  */
  /* smooth_sigma_fact, pyramid_sigma_fact, window_size, and subsampling */
  KLTUpdateTCBorder(&tc);

  return(tc);
}

    // Launch 
  	void main(void) 
  	{
  		// Initialize data
  		unsigned char img1 [NUM_COLS*NUM_ROWS];
  		unsigned char img2 [NUM_COLS*NUM_ROWS];
  		KLT_TrackingContext tc;
  		KLT_FeatureList fl;
  		KLT_FeatureTable ft;
  		int nFeatures, nFrames;
  		
  		// Loop variables
  		int i, ii;
  		
  		// More initialization
  		tc = KLTCreateTrackingContext();
  		nFeatures = 1024;
  		nFrames = 510;
  		
  		/* DEBUG, print out the Tracking Contents */
  		//KLTPrintTrackingContext(tc);
  		
  		fl = KLTCreateFeatureList(nFeatures);
  		
  		ft = KLTCreateFeatureTable(nFrames, nFeatures);
  		tc.sequentialMode = TRUE;
  		tc.writeInternalImages = FALSE;
  		sigma_last = -10.0;
  		
  		/* set this to 2 to turn on affine consistency check */
  		tc.affineConsistencyCheck = -1;  
  		
  		// Receive and store the first image data
  		for (i = 0; i < NUM_ROWS*NUM_COLS; i++)
  		{
  			bytesFromStimulus.receive(&img1[i], sizeof(char));
  		}
  		
  		// Select Good Features
  		KLTSelectGoodFeatures(&tc, img1, NUM_COLS, NUM_ROWS, fl);
  	
    	
    	{      		
      		// Receive Image behavior
      		for (i = 0; i < NUM_ROWS*NUM_COLS; i++)
  			{
  				bytesFromStimulus.receive(&img2[i], sizeof(char));
  			}
      		
      		// Track features
      		// ...
      		
      		// Store features?  or send to monitor???
      		// ...
      		
    	}  // end par
    
  	}  // end void main void
  
};  // end behavior
