// Select Good Features
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>       /* log */
#include <assert.h>
#include <sim.sh>
#include <string.h>

#include "main.sh"

#include "klt_util.h"
#include "klt.h"
#include "pyramid.h"
#include "convolve.h"



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

/**
 * Select Good Features
 */
behavior Select(KLT_TrackingContext tc, KLT_FeatureList fl, unsigned char img1[NUM_ROWS*NUM_COLS*sizeof(unsigned char)])
{

typedef enum {SELECTING_ALL, REPLACING_SOME} selectionMode;

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
int i;
  KLT_PixelType *ptrend;
  float *ptrout;

  ptrend = img + ncols*nrows;
  ptrout = floatimg->data;

  /* Output image must be large enough to hold result */
  assert(floatimg->ncols >= ncols);
  assert(floatimg->nrows >= nrows);

  floatimg->ncols = ncols;
  floatimg->nrows = nrows;

  i = 0;
  while (img < ptrend)  
  {
  	*ptrout++ = (float) *img++;
  	i++;
  }
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
int i;
				
  /* Output images must be large enough to hold result */
  assert(gradx->ncols >= img->ncols);
  assert(gradx->nrows >= img->nrows);
  assert(grady->ncols >= img->ncols);
  assert(grady->nrows >= img->nrows);

  /* Compute kernels, if necessary */
  if (fabs(sigma - sigma_last) > 0.05)
  {
    _computeKernels(sigma, &gauss_kernel, &gaussderiv_kernel);
  }
	
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

uchar featuremap[NUM_COLS*NUM_ROWS]; /* Boolean array recording proximity of features */
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
  int x, y, val, i;     /* Location and trackability of pixel under consideration */
  
  int *ptr;
	
  /* Cannot add features with an eigenvalue less than one */
  if (l_min_eigenvalue < 1)  l_min_eigenvalue = 1;

  /* Initialize variables */
  for (i = 0; i < NUM_COLS*NUM_ROWS; i++)
  {
  	featuremap[i] = 0;
  }
	
  /* Necessary because code below works with (mindist-1) */
  l_mindist--;

  /* If we are keeping all old good features, then add them to the featuremap */
  if (!overwriteAllFeatures)
    for (indx = 0 ; indx < featurelist->nFeatures ; indx++)
      if (featurelist->feature[indx]->val >= 0)  {
        x   = (int) featurelist->feature[indx]->x;
        y   = (int) featurelist->feature[indx]->y;
        _fillFeaturemap(x, y, &featuremap[0], l_mindist, ncols, nrows);
      }

  /* For each feature point, in descending order of importance, do ... */
  ptr = pointlist;
  indx = 0;
  while (1)  
  {

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
    if (!featuremap[y*ncols+x] && val >= l_min_eigenvalue)  
    {
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
      _fillFeaturemap(x, y, &featuremap[0], l_mindist, ncols, nrows);
    }
    
  }  // end while

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
  int window_hw, window_hh, i;
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
  if (mode == REPLACING_SOME && tc->sequentialMode && tc->pyramid_last != NULL)  
  {
    floatimg = ((_KLT_Pyramid) tc->pyramid_last)->img[0];
    gradx = ((_KLT_Pyramid) tc->pyramid_last_gradx)->img[0];
    grady = ((_KLT_Pyramid) tc->pyramid_last_grady)->img[0];
    assert(gradx != NULL);
    assert(grady != NULL);
  } 
  else  
  {
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


	/**
	 * Main behavior entry
	 */
	void main(void)
  	{
  		sigma_last = -10.0;
  	
  		// Select Good Features
  		printf("SELECT::KLTSelectGoodFeatures\n");
  		KLTSelectGoodFeatures(&tc, img1, NUM_COLS, NUM_ROWS, fl);
  		//KLTSelectGoodFeatures(&tc, img1);
  	}

};  // end behavior