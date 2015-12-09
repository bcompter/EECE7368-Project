// Monitor
// 

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sim.sh>

/*Our includes */
#include main.sh



import "i_receiver";

behavior Monitor(i_receiver bytes) 
{
	
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
	
	
	/*allows table to be written to file*/
	void KLTWriteFeatureTable(
  KLT_FeatureTable ft,
  char *fname, 
  char *fmt)
{
  FILE *fp;
  char format[100];
  char type;
  int i, j;

  if (KLT_verbose >= 1 && fname != NULL)  {
    fprintf(stderr,  
            "(KLT) Writing feature table to %s file: '%s'\n", 
            (fmt == NULL ? "binary" : "text"), fname);
  }

  if (fmt != NULL) {  /* text file or stderr */
    fp = _printSetupTxt(fname, fmt, format, &type);
    _printHeader(fp, format, FEATURE_TABLE, ft->nFrames, ft->nFeatures);

    for (j = 0 ; j < ft->nFeatures ; j++)  {
      fprintf(fp, "%7d | ", j);
      for (i = 0 ; i < ft->nFrames ; i++)
        _printFeatureTxt(fp, ft->feature[j][i], format, type);
      fprintf(fp, "\n");
    }
    _printShutdown(fp);
  } else {  /* binary file */
    fp = _printSetupBin(fname);
    fwrite(binheader_ft, sizeof(char), BINHEADERLENGTH, fp); 
    fwrite(&(ft->nFrames), sizeof(int), 1, fp);
    fwrite(&(ft->nFeatures), sizeof(int), 1, fp);
    for (j = 0 ; j < ft->nFeatures ; j++)  {
      for (i = 0 ; i < ft->nFrames ; i++)  {
        _printFeatureBin(fp, ft->feature[j][i]);
      }
    }
    fclose(fp);
  }
}

	
	
  void main(void) 
  {
	bytes.receive()//todo, need to figure out how to get from bytes.receive to get the ft.
	//can make nFrame static value here but would prefer to pass value to keep code in sync
	calc_motion_vector (ft, nFrame);
	writeFeatures(ft, "features.txt", "%5.1f");	  
	  
  
  
    printf ("Frames comparison written to file!\n");
  }  // end main(void)
  
};  // end Monitor