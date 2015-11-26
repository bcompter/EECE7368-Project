######################################################################
# Choose your favorite C compiler
CC = gcc -O3 #-g

######################################################################
# -DNDEBUG prevents the assert() statements from being included in 
# the code.  If you are having problems running the code, you might 
# want to comment this line to see if an assert() statement fires.
FLAG1 = -DNDEBUG

######################################################################
# -DKLT_USE_QSORT forces the code to use the standard qsort() 
# routine.  Otherwise it will use a quicksort routine that takes
# advantage of our specific data structure to greatly reduce the
# running time on some machines.  Uncomment this line if for some
# reason you are unhappy with the special routine.
# FLAG2 = -DKLT_USE_QSORT

######################################################################
# Add your favorite C flags here.
CFLAGS = $(FLAG1) $(FLAG2) 
INC_PATH=-I/Users/Grad/chulian/OpenCV/include/opencv2 \
	 -I/Users/Grad/chulian/OpenCV/include \
	 -I.

######################################################################
# There should be no need to modify anything below this line (but
# feel free to if you want).

EXAMPLES = sequential.c 
ARCH = convolve.c error.c pnmio.c pyramid.c selectGoodFeatures.c \
       storeFeatures.c trackFeatures.c klt.c klt_util.c writeFeatures.c
LIB = -L/usr/local/lib -L/usr/lib -L/Users/Grad/chulian/OpenCV/lib


.SUFFIXES:  .c .o

all:  lib $(EXAMPLES:.c=)

.c.o:
	$(CC) ${INC_PATH} -c $(CFLAGS) $<

lib: $(ARCH:.c=.o)
	rm -f libklt.a
	ar ruv libklt.a $(ARCH:.c=.o)
	rm -f *.o

sequential: libklt.a
	$(CC) ${INC_PATH} -O3 $(CFLAGS) -o $@ $@.c -L. -lklt $(LIB) -lm -lopencv_core -lopencv_highgui

depend:
	makedepend $(ARCH) $(EXAMPLES)

clean:
	rm -f *.o *.a $(EXAMPLES:.c=) *.tar *.tar.gz libklt.a



