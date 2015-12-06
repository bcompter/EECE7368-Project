//  Stimulus
//

#include <stdio.h>
#include <stdlib.h>
#include <sim.sh>

/* Our common defines */
#include "main.sh"

import "i_sender";

behavior Stimulus(unsigned char imageBuffer[NUM_ROWS*NUM_COLS*sizeof(unsigned char)], 
	i_sender bytesToDesign)
{  

/*********************************************************************
 * pnmio.c
 *
 * Various routines to manipulate PNM files.
 *********************************************************************/


/* Standard includes */
#include <stdio.h>   /* FILE  */
#include <stdlib.h>  /* malloc(), atoi() */

/* Our includes */
//#include "error.h"

#define LENGTH 80

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
  printf("ERROR...\n");
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

void _getNextString(
  FILE *fp,
  char *line)
{
  int i;

  line[0] = '\0';

  while (line[0] == '\0')  {
    fscanf(fp, "%s", line);
    i = -1;
    do  {
      i++;
      if (line[i] == '#')  {
        line[i] = '\0';
        while (fgetc(fp) != '\n') ;
      }
    }  while (line[i] != '\0');
  }
}


/*********************************************************************
 * pnmReadHeader
 */

void pnmReadHeader(
  FILE *fp, 
  int *magic, 
  int *ncols, int *nrows, 
  int *maxval)
{
  char line[LENGTH];
	
  /* Read magic number */
  _getNextString(fp, line);
  if (line[0] != 'P')
    KLTError("(pnmReadHeader) Magic number does not begin with 'P', "
             "but with a '%c'", line[0]);
  sscanf(line, "P%d", magic);
	
  /* Read size, skipping comments */
  _getNextString(fp, line);
  *ncols = atoi(line);
  _getNextString(fp, line);
  *nrows = atoi(line);
  if (*ncols < 0 || *nrows < 0 || *ncols > 10000 || *nrows > 10000)
    KLTError("(pnmReadHeader) The dimensions %d x %d are unacceptable",
             *ncols, *nrows);
	
  /* Read maxval, skipping comments */
  _getNextString(fp, line);
  *maxval = atoi(line);
  fread(line, 1, 1, fp); /* Read newline which follows maxval */
	
  if (*maxval != 255)
    KLTWarning("(pnmReadHeader) Maxval is not 255, but %d", *maxval);
}


/*********************************************************************
 * pgmReadHeader
 */

void pgmReadHeader(
  FILE *fp, 
  int *magic, 
  int *ncols, int *nrows, 
  int *maxval)
{
  pnmReadHeader(fp, magic, ncols, nrows, maxval);
  if (*magic != 5)
    KLTError("(pgmReadHeader) Magic number is not 'P5', but 'P%d'", *magic);
}


/*********************************************************************
 * ppmReadHeader
 */

void ppmReadHeader(
  FILE *fp, 
  int *magic, 
  int *ncols, int *nrows, 
  int *maxval)
{
  pnmReadHeader(fp, magic, ncols, nrows, maxval);
  if (*magic != 6)
    KLTError("(ppmReadHeader) Magic number is not 'P6', but 'P%d'", *magic);
}


/*********************************************************************
 * pgmReadHeaderFile
 */

void pgmReadHeaderFile(
  char *fname, 
  int *magic, 
  int *ncols, int *nrows, 
  int *maxval)
{
  FILE *fp;

  /* Open file */
  if ( (fp = fopen(fname, "rb")) == NULL)
    KLTError("(pgmReadHeaderFile) Can't open file named '%s' for reading\n", fname);

  /* Read header */
  pgmReadHeader(fp, magic, ncols, nrows, maxval);

  /* Close file */
  fclose(fp);
}


/*********************************************************************
 * ppmReadHeaderFile
 */

void ppmReadHeaderFile(
  char *fname, 
  int *magic, 
  int *ncols, int *nrows, 
  int *maxval)
{
  FILE *fp;

  /* Open file */
  if ( (fp = fopen(fname, "rb")) == NULL)
    KLTError("(ppmReadHeaderFile) Can't open file named '%s' for reading\n", fname);

  /* Read header */
  ppmReadHeader(fp, magic, ncols, nrows, maxval);

  /* Close file */
  fclose(fp);
}


/*********************************************************************
 * pgmRead
 *
 * NOTE:  If img is NULL, memory is allocated.
 */

void pgmRead(FILE *fp)
{
  int magic, maxval;
  int i, index;
  int nCols, nRows;

  /* Read header */
  pgmReadHeader(fp, &magic, &nCols, &nRows, &maxval);

  /* Read binary image data */
  for (i = 0 ; i < NUM_ROWS ; i++)  
  {
    fread(&imageBuffer[index], NUM_COLS, 1, fp);
    index += NUM_COLS;
  }

}


/*********************************************************************
 * pgmReadFile
 *
 * NOTE:  If img is NULL, memory is allocated.
 */

void pgmReadFile(char *fname)
{
  FILE *fp;

  /* Open file */
  if ( (fp = fopen(fname, "rb")) == NULL)
    KLTError("(pgmReadFile) Can't open file named '%s' for reading\n", fname);

  /* Read file */
  pgmRead(fp);

  /* Close file */
  fclose(fp);
}


/*********************************************************************
 * pgmWrite
 */

void pgmWrite(
  FILE *fp,
  unsigned char *img, 
  int ncols, 
  int nrows)
{
  int i;

  /* Write header */
  fprintf(fp, "P5\n");
  fprintf(fp, "%d %d\n", ncols, nrows);
  fprintf(fp, "255\n");

  /* Write binary data */
  for (i = 0 ; i < nrows ; i++)  {
    fwrite(img, ncols, 1, fp);
    img += ncols;
  }
}


/*********************************************************************
 * pgmWriteFile
 */

void pgmWriteFile(
  char *fname, 
  unsigned char *img, 
  int ncols, 
  int nrows)
{
  FILE *fp;

  /* Open file */
  if ( (fp = fopen(fname, "wb")) == NULL)
    KLTError("(pgmWriteFile) Can't open file named '%s' for writing\n", fname);

  /* Write to file */
  pgmWrite(fp, img, ncols, nrows);

  /* Close file */
  fclose(fp);
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
    KLTError("(ppmWriteFileRGB) Can't open file named '%s' for writing\n", fname);

  /* Write to file */
  ppmWrite(fp, redimg, greenimg, blueimg, ncols, nrows);

  /* Close file */
  fclose(fp);
}

  /** 
   * Read input files and pass frame data to the design behavior using a queue
   */
  void main(void)
  {
  
  	char fnamein[100], fnameout[100];
  	int i, ii;
  	
  	// Loop over all frames
  	for (i = 0; i < NUM_FRAMES; i++)
  	{
  		// Load the frame data
  		sprintf(fnamein, "huntington_1280/huntington_1080p_60fps_%d.pgm", i+1);
    	pgmReadFile(fnamein);
  		
  		// Send the to queue
  		for (ii = 0; ii < NUM_ROWS*NUM_COLS; ii++)
  		{
  			bytesToDesign.send(&imageBuffer[ii], sizeof(char));
  			
  		}  // end of loading the queue
  		
  	}  // end for each frame
	
    return;
  }

};  // end behavior