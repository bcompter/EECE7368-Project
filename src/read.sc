#include <stdio.h>
#include <stdlib.h>
#include <math.h>       /* log */
#include <assert.h>
#include <sim.sh>
#include <string.h>

#include "main.sh"

import "i_receiver";

behavior Read ( i_receiver bytesFromStimulus,
				unsigned char img [NUM_COLS*NUM_ROWS])
				
{  
	
  void main(void) 
  {
  		int i;
  
    	// Receive and store image data
  		printf("Receiving image data\n");
  		for (i = 0; i < NUM_ROWS*NUM_COLS; i++)
  		{
  			bytesFromStimulus.receive(&img[i], sizeof(char));
  		}
  }  // end void main void
  
};  // end behavior
