// Design
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
#include "convolve.h"
#include "pyramid.h"

import "i_sender";
import "i_receiver";
import "i_send";
import "i_receive";
import "read.sc";
import "track.sc";
import "c_queue";

behavior Design(
				unsigned char img [NUM_COLS*NUM_ROWS],
				i_receive start,
				i_send ready,
				i_sender imageBytesToMonitor, 
				i_sender featureBytesToMonitor) 
{    

 const int mindist = 10;
 const int window_size = 4;
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

 const KLT_BOOL smoothBeforeSelecting = FALSE;
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
    
    printf("KLTCreateFeatureList:: nbytes = %d\n", nbytes);
	
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
 * KLTPrintTrackingContext
 */

void KLTPrintTrackingContext(
  KLT_TrackingContext *tc)
{
  printf("\n\nTracking context:\n\n");
  printf("\tmindist = %d\n", tc->mindist);
  printf("\twindow_width = %d\n", tc->window_width);
  printf("\twindow_height = %d\n", tc->window_height);
  printf("\tsequentialMode = %s\n",
          tc->sequentialMode ? "TRUE" : "FALSE");
  printf("\tsmoothBeforeSelecting = %s\n",
          tc->smoothBeforeSelecting ? "TRUE" : "FALSE");
  printf("\twriteInternalImages = %s\n",
          tc->writeInternalImages ? "TRUE" : "FALSE");

  printf("\tmin_eigenvalue = %d\n", tc->min_eigenvalue);
  printf("\tmin_determinant = %f\n", tc->min_determinant);
  printf("\tmin_displacement = %f\n", tc->min_displacement);
  printf("\tmax_iterations = %d\n", tc->max_iterations);
  printf("\tmax_residue = %f\n", tc->max_residue);
  printf("\tgrad_sigma = %f\n", tc->grad_sigma);
  printf("\tsmooth_sigma_fact = %f\n", tc->smooth_sigma_fact);
  printf("\tpyramid_sigma_fact = %f\n", tc->pyramid_sigma_fact);
  printf("\tnSkippedPixels = %d\n", tc->nSkippedPixels);
  printf("\tborderx = %d\n", tc->borderx);
  printf("\tbordery = %d\n", tc->bordery);
  printf("\tnPyramidLevels = %d\n", tc->nPyramidLevels);
  printf("\tsubsampling = %d\n", tc->subsampling);
  
  printf("\n\tpyramid_last = %s\n", (tc->pyramid_last!=NULL) ?
          "points to old image" : "NULL");
  printf("\tpyramid_last_gradx = %s\n", 
          (tc->pyramid_last_gradx!=NULL) ?
          "points to old image" : "NULL");
  printf("\tpyramid_last_grady = %s\n",
          (tc->pyramid_last_grady!=NULL) ?
          "points to old image" : "NULL");
  printf("\n\n");
}

/*********************************************************************
 * KLTCreateTrackingContext
 *
 */
KLT_TrackingContext KLTCreateTrackingContext(KLT_TrackingContext * tc)
{
  /* Set values to default values */
  tc->mindist = mindist;
  tc->window_width = window_size;
  tc->window_height = window_size;
  tc->sequentialMode = sequentialMode;
  tc->smoothBeforeSelecting = smoothBeforeSelecting;
  tc->writeInternalImages = writeInternalImages;
  tc->lighting_insensitive = lighting_insensitive;
  tc->min_eigenvalue = min_eigenvalue;
  tc->min_determinant = min_determinant;
  tc->max_iterations = max_iterations;
  tc->min_displacement = min_displacement;
  tc->max_residue = max_residue;
  tc->grad_sigma = grad_sigma;
  tc->smooth_sigma_fact = smooth_sigma_fact;
  tc->pyramid_sigma_fact = pyramid_sigma_fact;
  tc->step_factor = step_factor;
  tc->nSkippedPixels = nSkippedPixels;
  tc->pyramid_last = NULL;
  tc->pyramid_last_gradx = NULL;
  tc->pyramid_last_grady = NULL;
  
  /* for affine mapping */
  tc->affineConsistencyCheck = affineConsistencyCheck;
  tc->affine_window_width = affine_window_size;
  tc->affine_window_height = affine_window_size;
  tc->affine_max_iterations = affine_max_iterations;
  tc->affine_max_residue = affine_max_residue;
  tc->affine_min_displacement = affine_min_displacement;
  tc->affine_max_displacement_differ = affine_max_displacement_differ;

  /* Change nPyramidLevels and subsampling */
  KLTChangeTCPyramid(tc);
	
  /* Update border, which is dependent upon  */
  /* smooth_sigma_fact, pyramid_sigma_fact, window_size, and subsampling */
  KLTUpdateTCBorder(tc);

}

	KLT_TrackingContext tc;
  	KLT_FeatureList fl;

	unsigned char img1 [NUM_COLS*NUM_ROWS];
	unsigned char img2 [NUM_COLS*NUM_ROWS];
	
	// Queues between the design and the monitor
	const unsigned long qSize = 2048;
	c_queue dataToTrack(qSize);

	// Behaviors From Imports
  	Read read(start, ready, img, dataToTrack);
  	Track track(tc, fl, dataToTrack, imageBytesToMonitor, featureBytesToMonitor);

    // Launch main behavior code
  	void main(void) 
  	{
  		
  		// Loop variables
  		int i, ii;

  		printf("Starting design\n");
  		
  		// More initialization
  		printf("Creating tracking context\n");
  		KLTCreateTrackingContext(&tc);
  		
  		/* DEBUG, print out the Tracking Contents */
  		KLTPrintTrackingContext(&tc);
  		
  		printf("DESIGN::Creating feature list\n");
  		fl = KLTCreateFeatureList(NUM_FEATURES);
  		
  		tc.sequentialMode = TRUE;
  		tc.writeInternalImages = FALSE;
  		sigma_last = -10.0;
  		
  		/* set this to 2 to turn on affine consistency check */
  		tc.affineConsistencyCheck = -1;  
  	    
    	par
    	{      		
      		// Receive Image data
      		read;
      		
      		// Track features
      		track;
      		
    	}  // end par
    
  	}  // end void main void
  
};  // end behavior
