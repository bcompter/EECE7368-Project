// Monitor
// 

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sim.sh>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include "main.sh"
#include "klt_util.h"
#include "klt.h"

import "i_receiver";

typedef enum {FEATURE_LIST, FEATURE_HISTORY, FEATURE_TABLE} structureType;

#define BINHEADERLENGTH	6


behavior Monitor(i_receiver imageData, i_receiver featureData) 
{

	KLT_FeatureList fl;

	char binheader_fl[BINHEADERLENGTH+1] = "KLTFL1";
	char binheader_fh[BINHEADERLENGTH+1] = "KLTFH1";
	char binheader_ft[BINHEADERLENGTH+1] = "KLTFT1";


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


FILE* _printSetupTxt(
  char *fname, 	/* Input: filename, or NULL for stderr */
  char *fmt,	/* Input: format (e.g., %5.1f or %3d) */
  char *format,	/* Output: format (e.g., (%5.1f,%5.1f)=%3d) */
  char *type)	/* Output: either 'f' or 'd', based on input format */
{
  FILE *fp;
  const int val_width = 5;
  int i;

  /* Either open file or use stderr */
  if (fname == NULL)  fp = stderr;
  else  fp = fopen(fname, "wb");
  if (fp == NULL)
    KLTError("(KLTWriteFeatures) "
             "Can't open file '%s' for writing\n", fname);

  /* Parse format */
  if (fmt[0] != '%')
    KLTError("(KLTWriteFeatures) Bad Format: %s\n", fmt);
  i = 0;  while (fmt[i] != '\0') i++;  *type = fmt[i-1];
  if (*type != 'f' && *type != 'd')
    KLTError("(KLTWriteFeatures) Format must end in 'f' or 'd'.");

  /* Construct feature format */
  sprintf(format, "(%s,%s)=%%%dd ", fmt, fmt, val_width);
     
  return fp;
}


FILE* _printSetupBin(
  char *fname) 	/* Input: filename */
{
  FILE *fp;
  if (fname == NULL) 
    KLTError("(KLTWriteFeatures) Can't write binary data to stderr");
  fp = fopen(fname, "wb");
  if (fp == NULL)
    KLTError("(KLTWriteFeatures) "
             "Can't open file '%s' for writing", fname);
  return fp;
}


void _printNhyphens(
  FILE *fp,
  int n)
{
  int i;
  for (i = 0 ; i < n ; i++)
    fprintf(fp, "-");
}

void _printInteger(
  FILE *fp,
  int integer,
  int width)
{
  char fmt[80];
  sprintf(fmt, "%%%dd", width);
  fprintf(fp, fmt, integer);
}


KLT_BOOL _isCharInString(
  char c,
  char *str)
{
  int width;
  int i;
  width = strlen(str);

  for (i = 0 ; i < width ; i++)
    if (c == str[i])  return TRUE;

  return FALSE;
}


/*********************************************************************
 * _findStringWidth
 *
 * Calculates the length of a string after expansion.  E.g., the
 * length of "(%6.1f)" is eight -- six for the floating-point number,
 * and two for the parentheses.
 */

int _findStringWidth(
  char *str)
{
	int width, add, maxi, i;

  width = 0;
  maxi = strlen(str) - 1;
  i = 0;
	
  while (str[i] != '\0')  {
    if (str[i] == '%')  {
      if (isdigit(str[i+1]))  {
        sscanf(str+i+1, "%d", &add);
        width += add;
        i += 2;
        while (!_isCharInString(str[i], "diouxefgn"))  {
          i++;
          if (i > maxi)
            KLTError("(_findStringWidth) Can't determine length "
                     "of string '%s'", str);
        }
        i++;
      } else if (str[i+1] == 'c')  {
        width++;
        i += 2;
      } else 
        KLTError("(_findStringWidth) Can't determine length "
                 "of string '%s'", str);
    } else  {
      i++;
      width++;
    }
  }
	
  return width;
}


void _printHeader(
  FILE *fp,
  char *format,
  structureType id,
  int nFrames,
  int nFeatures)
{
  int i, width;
  width = _findStringWidth(format);
	
  assert(id == FEATURE_LIST || id == FEATURE_HISTORY || id == FEATURE_TABLE);

  if (fp != stderr)  {
    fprintf(fp, "Feel free to place comments here.\n\n\n");
    fprintf(fp, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    fprintf(fp, "WARNING");
    fprintf(fp, "\n");
  }
  fprintf(fp, "------------------------------\n");
  switch (id)  {
     case FEATURE_LIST: fprintf(fp, "KLT Feature List\n");    break;
     case FEATURE_HISTORY: fprintf(fp, "KLT Feature History\n"); break;
     case FEATURE_TABLE: fprintf(fp, "KLT Feature Table\n");   break;
  }

  fprintf(fp, "------------------------------\n\n");
  switch (id)  {
     case FEATURE_LIST: fprintf(fp, "nFeatures = %d\n\n", nFeatures);    break;
     case FEATURE_HISTORY: fprintf(fp, "nFrames = %d\n\n", nFrames); break;
     case FEATURE_TABLE: fprintf(fp, "nFrames = %d, nFeatures = %d\n\n",
                                 nFrames, nFeatures);   break;
  }

  switch (id)  {
     case FEATURE_LIST: fprintf(fp, "feature | (x,y)=val\n");
       fprintf(fp, "--------+-");
       _printNhyphens(fp, width);
       fprintf(fp, "\n");   
       break;
     case FEATURE_HISTORY: fprintf(fp, "frame | (x,y)=val\n");
       fprintf(fp, "------+-");
       _printNhyphens(fp, width);
       fprintf(fp, "\n");   
       break;
     case FEATURE_TABLE: fprintf(fp, "feature |          frame\n");   
       fprintf(fp, "        |");   
       for (i = 0 ; i < nFrames ; i++) _printInteger(fp, i, width);
       fprintf(fp, "\n--------+-");   
       for (i = 0 ; i < nFrames ; i++) _printNhyphens(fp, width);
       fprintf(fp, "\n");   
       break;
  }
}


void _printFeatureTxt(
  FILE *fp,
  KLT_Feature feat,
  char *format,
  char type)
{
	float x, y;

  assert(type == 'f' || type == 'd');

  if (type == 'f')
    fprintf(fp, format, (float) feat->x, (float) feat->y, feat->val);
  else if (type == 'd')  
  {
    /* Round x & y to nearest integer, unless negative */
    x = feat->x;
    y = feat->y;
    if (x >= 0.0) x += 0.5;
    if (y >= 0.0) y += 0.5;
    fprintf(fp, format,(int) x, (int) y, feat->val);
  }
}

void _printFeatureBin(
  FILE *fp,
  KLT_Feature feat)
{
  fwrite(&(feat->x), sizeof(KLT_locType), 1, fp);
  fwrite(&(feat->y), sizeof(KLT_locType), 1, fp);
  fwrite(&(feat->val), sizeof(int), 1, fp);
}


void _printShutdown(
  FILE *fp)
{
  /* Close file, if necessary */
  if (fp != stderr)
    fclose(fp);
}


/*********************************************************************
 * KLTWriteFeatureList()
 * KLTWriteFeatureHistory()
 * KLTWriteFeatureTable()
 * 
 * Writes features to file or to screen.
 *
 * INPUTS
 * fname: name of file to write data; if NULL, then print to stderr
 * fmt:   format for printing (e.g., "%5.1f" or "%3d");
 *        if NULL, and if fname is not NULL, then write to binary file.
 */

void KLTWriteFeatureList(
  KLT_FeatureList fl,
  char *fname, 
  char *fmt)
{
  FILE *fp;
  char format[100];
  char type;
  int i;


  if (fmt != NULL) 
  {  /* text file or stderr */
    fp = _printSetupTxt(fname, fmt, format, &type);
    //_printHeader(fp, format, FEATURE_LIST, 0, fl->nFeatures);
	
    for (i = 0 ; i < fl->nFeatures ; i++)  
    {
      fprintf(fp, "%d,", i);
      fprintf(fp, "%d,%d,%d",(int) fl->feature[i]->x, (int) fl->feature[i]->y, fl->feature[i]->val);
      //_printFeatureTxt(fp, fl->feature[i], format, type);
      fprintf(fp, "\n");
    }
    _printShutdown(fp);
  } 
  else 
  {  /* binary file */
    fp = _printSetupBin(fname);
    fwrite(binheader_fl, sizeof(char), BINHEADERLENGTH, fp); 
    fwrite(&(fl->nFeatures), sizeof(int), 1, fp);
    for (i = 0 ; i < fl->nFeatures ; i++)  
    {
      _printFeatureBin(fp, fl->feature[i]);
    }
    fclose(fp);
  }
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
  			//sprintf(filename, "results/output_frame_%d.ppm", frame+1);
  			//KLTWriteFeatureListToPPM(fl, img, NUM_COLS, NUM_ROWS, filename);
  			
  			// Write out our feature list to file
  			//for (i = 0; i < NUM_FEATURES; i++)
  			{
  				sprintf(filename, "results/output_features_%d.txt", frame+1);
  				KLTWriteFeatureList(fl, filename, "%3d");
  			}
  			
  			frame++;
  		}
  		
  		
  		// All done, exit the program...
  		exit(0);
    
  }  // end main(void)
  
};  // end Monitor