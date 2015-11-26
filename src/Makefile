# src/Makefile 
# 
#    build specification code 
# Author: Gunar Schirner  10/06/09


# name of the design and of the top level file 
DESIGN          = digicam


# list all your *.sc files
#   -- auto determined 
SOURCES := $(wildcard *.sc)

# list all header files *.sh 
#   -- auto determined
HEADERS := $(wildcard *.sh)

SRC             = .


SCC		= scc
SCCDEF		= 
SCCOPT		= -ww -vv -d -g -G -xlx -I$(SPECC)/inc -P$(SPECC)/imp
SCCINC		= -I$(SRC)
SCCIMP		= -P$(SRC)
SCCLIB		=



# get common definitions
include ../Makefile.macros 



# default rule 
#  reqeust the final executable (same as design name)
default: chkspecc $(DESIGN)


# build the design using pattern rules
#   make design depend on all sources and headers to trigger rebuild
$(DESIGN): $(HEADERS) $(SOURCES)


# compare output against golden reference
test: $(DESIGN)
	./$(DESIGN) 
	$(CMP) $(TESTFILE) $(GOLDFILE)

clean:
	-$(RM) *~ *.o
	-$(RM) *.si
	-$(RM) *.sir
	-$(RM) $(DESIGN)
