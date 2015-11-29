/**********************************************************************
Finds the 100 best features in an image, and tracks these
features to the next image.  Saves the feature
locations (before and after tracking) to text files and to PPM files, 
and prints the features to the screen.
**********************************************************************/

#include <stdlib.h>
#include "pnmio.h"
#include "klt.h"
#include <math.h>
//#include "cv.h"
#include "highgui/highgui_c.h"

int trigger = 0;
void calc_motion_vector(KLT_FeatureTable ft, int idx);

int main()
{
  unsigned char *img1, *img2;
  KLT_TrackingContext tc;
  KLT_FeatureList fl;
  char fnamein[100], fnameout[100];
  //int nFeatures = (512+128), nFrames = 61;
  int nFeatures = 1024, nFrames = 510;
  int ncols, nrows;
  int i;
  // timing variables
  struct timeval start,end;
  struct timeval imgRead_start, imgRead_end;
  double total_time = 0; 
  double imgRead_time = 0;

  tc = KLTCreateTrackingContext();
  KLTPrintTrackingContext(tc);
  fl = KLTCreateFeatureList(nFeatures);
  KLT_FeatureTable ft;
  ft = KLTCreateFeatureTable(nFrames, nFeatures);
  tc->sequentialMode = TRUE;
  tc->writeInternalImages = FALSE;
  /* set this to 2 to turn on affine consistency check */
  tc->affineConsistencyCheck = -1;  

  img1 = pgmReadFile("huntington_1280/huntington_1080p_60fps_1.pgm",
      NULL, &ncols, &nrows);
  img2 = (unsigned char *) malloc(ncols*nrows*sizeof(unsigned char));
	printf("");

  KLTSelectGoodFeatures(tc, img1, ncols, nrows, fl);

  // start timing
  gettimeofday(&start, NULL);
  
  for (i = 1 ; i < nFrames ; i++)  
  {
    gettimeofday(&imgRead_start, NULL);
    sprintf(fnamein, "huntington_1280/huntington_1080p_60fps_%d.pgm", i+1);
    pgmReadFile(fnamein, img2, &ncols, &nrows);

    gettimeofday(&imgRead_end, NULL);
    imgRead_time += (double)(imgRead_end.tv_sec - imgRead_start.tv_sec) 
      + (double)(imgRead_end.tv_usec - imgRead_start.tv_usec)/1000000;

    KLTTrackFeatures(tc, img1, img2, ncols, nrows, fl);
#ifdef REPLACE
    KLTReplaceLostFeatures(tc, img2, ncols, nrows, fl);
#endif
    KLTStoreFeatureList(fl, ft, i);
  }  // end for number of frames...

  // end timing
  gettimeofday(&end, NULL);
  total_time += (double)(end.tv_sec - start.tv_sec) 
    + (double)(end.tv_usec - start.tv_usec)/1000000;

  //printf("OpticalFlow time is %f sec.\n", total_time - imgRead_time);
  // calc the motion vector (frame 0 -> frame 1)
  calc_motion_vector(ft, nFrames);

  KLTWriteFeatureTable(ft, "features.txt", "%5.1f");
  //KLTWriteFeatureTable(ft, "features.ft", NULL);
  //KLTWriteFeatureListToPPM(fl, img1, ncols, nrows, "feat1.ppm");
  //KLTWriteFeatureList(fl, "feat1.txt", "%3d");

  KLTFreeFeatureTable(ft);
  KLTFreeFeatureList(fl);
  KLTFreeTrackingContext(tc);
  free(img1);
  free(img2);
  return 0;
}


// calc motion vector and store it into files
void calc_motion_vector(KLT_FeatureTable ft, int nFrame)
{
//	FILE *fp;
//	fp = fopen("13-window.txt", "w");
//
//	int i, val;
//	float x, y;
//	for (i = 0; i < ft->nFeatures; i ++) {
//		val = ft->feature[i][2]->val;
//		if (val == KLT_SMALL_DET || val == KLT_OOB || val == KLT_LARGE_RESIDUE
//				|| val == KLT_MAX_ITERATIONS) {
//			x = 0;
//			y = 0;
//		} else {
//			x = fabs(ft->feature[i][2]->x - ft->feature[i][1]->x);
//			y = fabs(ft->feature[i][2]->y - ft->feature[i][1]->y);
//		}
//		fprintf(fp, "%f, %f,\n", x, y);
//	}
//
//	fclose(fp);

	/* use OpenCV to write bmp file */
	IplImage * pCVFile;
	char filename[32];
	int i, val, f_idx;
	float x, y;

	pCVFile = cvCreateImage(cvSize(ft->nFeatures, 2*(nFrame-2)), 8, 1);

	for (f_idx = 2; f_idx < nFrame; f_idx ++) {
		for(i = 0; i < ft->nFeatures; i ++)
		{
			val = ft->feature[i][2]->val;
			if (val == KLT_SMALL_DET || val == KLT_OOB || val == KLT_LARGE_RESIDUE
					|| val == KLT_MAX_ITERATIONS) {
				x = 0;
				y = 0;
			} else {
				x = fabs(ft->feature[i][f_idx]->x - ft->feature[i][f_idx - 1]->x) * 100;
				y = fabs(ft->feature[i][f_idx]->y - ft->feature[i][f_idx - 1]->y) * 100;
			}

			if (x > 255) {
				pCVFile->imageData[(f_idx-2)*2*nFrame + 2*i] = (char)(255);
			} else {
				pCVFile->imageData[(f_idx-2)*2*nFrame + 2*i] = (char)x;
			}

			if (y > 255) {
				pCVFile->imageData[(f_idx-2)*2*nFrame + 2*i+1] = (char)(255);
			} else {
				pCVFile->imageData[(f_idx-2)*2*nFrame + 2*i+1] = (char)y;
			}
		}
	}

//	sprintf(filename, "11-window-%04d.bmp", idx);
//	cvSaveImage(filename, pCVFile);
	cvSaveImage("result/3_window.bmp", pCVFile, NULL);
	cvReleaseImage(&pCVFile);
}

