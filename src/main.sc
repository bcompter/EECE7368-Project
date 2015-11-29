// Lucas-Kanade Project
//
// Final Project
// Group Members: 
//   Brian Compter, 001986259
//
//

#include <stdio.h>
#include <sim.sh>
#include "main.sh"

import "stimulus";
import "design";
import "monitor";

import "c_queue";
import "c_handshake";
import "c_double_handshake";

behavior Main 
{  
	// Global storage for frame data
	unsigned char imageBuffer[NUM_ROWS*NUM_COLS*sizeof(unsigned char)];

	// Queue between stimulus and design
	const unsigned long qSize = 1024;
  	c_queue bytesToDesign(qSize);

	// Trigger Read to start
  	c_handshake start;

    // Behaviors
  	Stimulus stimulus(imageBuffer, bytesToDesign);
	Design design(bytesToDesign);
  	Monitor monitor();

  // Main application entry point
  int main(void) 
  {
    par
    {
      	stimulus;
      	design;
      	monitor;
    }

    return 0;
  }  // end int main void
  
};  // end behavior
