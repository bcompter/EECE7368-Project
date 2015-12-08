#define MAX_KERNEL_WIDTH 	71

typedef struct  {
  int width;
  float data[MAX_KERNEL_WIDTH];
}  ConvolutionKernel;