// Common definitions for the Lucas Kanade algorithm
//

#define NUM_COLS 1280
#define NUM_ROWS 720

#define NUM_FRAMES 5

#define NUM_FEATURES 1024

#define FEATURE_BYTES 53256

#define min(a,b) (((a)<(b))?(a):(b))
#define max(a,b) (((a)>(b))?(a):(b))

#ifndef uchar
#define uchar unsigned char
#endif

// Shorthand form of a KLT feature
typedef struct {
	int x;
	int y;
	float value;
} Feature;



