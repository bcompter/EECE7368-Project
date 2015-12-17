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

behavior Main 
{  
	// Global storage for frame data
	unsigned char imageBuffer[NUM_ROWS*NUM_COLS*sizeof(unsigned char)];

	// Queues between the design and the monitor
	const unsigned long qSize = 2048;
	c_queue imageBytesToMonitor(qSize);
	c_queue featureBytesToMonitor(qSize);
	
	// Handshakes for synchronization
	c_handshake start;
	c_handshake ready;

    // Behaviors
  	Stimulus stimulus(imageBuffer, start, ready);
	Design design(imageBuffer, start, ready, imageBytesToMonitor, featureBytesToMonitor);
  	Monitor monitor(imageBytesToMonitor, featureBytesToMonitor);

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
