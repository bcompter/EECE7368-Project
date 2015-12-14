// Monitor
// 

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sim.sh>
#include <string.h>
#include "main.sh"
#include "klt_util.h"
#include "klt.h"

import "i_receiver";

behavior Monitor(i_receiver imageData, i_receiver featureData) 
{

	KLT_FeatureList fl;


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
 * ppmWrite
 */

void ppmWrite(
  FILE *fp,
  unsigned char *redimg,
  unsigned char *greenimg,
  unsigned char *blueimg,
  int ncols, 
  int nrows)
{
  int i, j;

  /* Write header */
  fprintf(fp, "P6\n");
  fprintf(fp, "%d %d\n", ncols, nrows);
  fprintf(fp, "255\n");

  /* Write binary data */
  for (j = 0 ; j < nrows ; j++)  {
    for (i = 0 ; i < ncols ; i++)  {
      fwrite(redimg, 1, 1, fp); 
      fwrite(greenimg, 1, 1, fp);
      fwrite(blueimg, 1, 1, fp);
      redimg++;  greenimg++;  blueimg++;
    }
  }
}


/*********************************************************************
 * ppmWriteFileRGB
 */

void ppmWriteFileRGB(
  char *fname, 
  unsigned char *redimg,
  unsigned char *greenimg,
  unsigned char *blueimg,
  int ncols, 
  int nrows)
{
  FILE *fp;

  /* Open file */
  if ( (fp = fopen(fname, "wb")) == NULL)
  {
    printf("(ppmWriteFileRGB) Can't open file named '%s' for writing\n", fname);
    KLTError("(ppmWriteFileRGB) Can't open file named '%s' for writing\n", fname);
  }
  
  /* Write to file */
  ppmWrite(fp, redimg, greenimg, blueimg, ncols, nrows);

  /* Close file */
  fclose(fp);
}


/*********************************************************************
 * KLTWriteFeatureListToPPM
 */

void KLTWriteFeatureListToPPM(
  KLT_FeatureList featurelist,
  KLT_PixelType *greyimg,
  int ncols,
  int nrows,
  char *filename)
{
  int nbytes;
  uchar *redimg, *grnimg, *bluimg;
  int offset;
  int x, y, xx, yy;
  int i;
  
  nbytes = ncols * nrows * sizeof(char);

  /* Allocate memory for component images */
  redimg = (uchar *)  malloc(nbytes);
  grnimg = (uchar *)  malloc(nbytes);
  bluimg = (uchar *)  malloc(nbytes);
  if (redimg == NULL || grnimg == NULL || bluimg == NULL)
    KLTError("(KLTWriteFeaturesToPPM)  Out of memory\n");

  /* Copy grey image to component images */
  if (sizeof(KLT_PixelType) != 1)
    KLTWarning("(KLTWriteFeaturesToPPM)  KLT_PixelType is not uchar");
  memcpy(redimg, greyimg, nbytes);
  memcpy(grnimg, greyimg, nbytes);
  memcpy(bluimg, greyimg, nbytes);
	
  /* Overlay features in red */
  for (i = 0 ; i < featurelist->nFeatures ; i++)
    if (featurelist->feature[i]->val >= 0)  {
      x = (int) (featurelist->feature[i]->x + 0.5);
      y = (int) (featurelist->feature[i]->y + 0.5);
      for (yy = y - 1 ; yy <= y + 1 ; yy++)
        for (xx = x - 1 ; xx <= x + 1 ; xx++)  
          if (xx >= 0 && yy >= 0 && xx < ncols && yy < nrows)  {
            offset = yy * ncols + xx;
            *(redimg + offset) = 255;
            *(grnimg + offset) = 0;
            *(bluimg + offset) = 0;
          }
    }
	
  /* Write to PPM file */
  ppmWriteFileRGB(filename, redimg, grnimg, bluimg, ncols, nrows);

  /* Free memory */
  free(redimg);
  free(grnimg);
  free(bluimg);
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


unsigned char img [NUM_COLS*NUM_ROWS];
  
  void main(void) 
  {

  		int frame, i;
  		
  		char filename[100];
		
  		// Buffer for feature list transmission
  		unsigned char buffer [sizeof(Feature)];
  		Feature tempFeature;
  		
  		printf("MONITOR::Starting\n");
  		
  		printf("MONITOR::Creating feature list\n");
  		fl = KLTCreateFeatureList(NUM_FEATURES);
  
  		frame = 0;
  		while(frame < NUM_FRAMES)
  		{
  			// Receive image data
  			printf("MONITOR::Receiving Frame %d\n", frame);
  			for (i = 0; i < NUM_ROWS * NUM_COLS; i++)
  			{
  				imageData.receive(&img[i], sizeof(char));
  			}
  			printf("MONITOR::Received Frame %d\n", frame);
  			
  			// Receive feature data
  			printf("MONITOR::Receiving Features %d\n", frame);
  			for (i = 0; i < NUM_FEATURES; i++)
  			{
  				// Receive the data
  				featureData.receive(&tempFeature, sizeof(Feature));
  				
  				// Now cram into our local feature list
  				fl->feature[i]->x = tempFeature.x;
  				fl->feature[i]->y = tempFeature.y;
  				fl->feature[i]->val = tempFeature.value;
  				
  			}  			
  			printf("MONITOR::Received Features %d\n", frame);
  			
  			// Write our feature list to a PPM file
  			sprintf(filename, "results/output_frame_%d.ppm", frame+1);
  			KLTWriteFeatureListToPPM(fl, img, NUM_COLS, NUM_ROWS, filename);
  			
  			frame++;
  		}
  		
  		
  		// All done, exit the program...
  		exit(0);
    
  }  // end main(void)
  
};  // end Monitor