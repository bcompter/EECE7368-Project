#include <stdio.h>
#include <stdlib.h>
#include <math.h>       /* log */
#include <assert.h>
#include <sim.sh>
#include <string.h>

#include "main.sh"

import "i_sender";
import "i_receive";
import "i_send";

behavior Read ( i_receive start,
				i_send ready,
				unsigned char img [NUM_COLS*NUM_ROWS],
				i_sender dataToTrack)
				
{  
	
  void main(void) 
  {
  		int i;
  		
  		// Notify Stimulus that we are ready to begin
  		ready.send();
  		
  		while (1)
  		{
  			// Wait for the stimulus to load data for us
  			start.receive();
  			
  			// Send data to the Track behavior
  			for (i = 0; i < NUM_ROWS*NUM_COLS; i++)
  			{
  				dataToTrack.send(&img[i], sizeof(char));
  			}
  			
  			// Ready for the next frame
  			ready.send();
  			
  		}
    	
  }  // end void main void
  
};  // end behavior
