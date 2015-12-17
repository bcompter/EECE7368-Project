#
# Top level makefile for lucas kanade example
#

# include default definitions
include Makefile.macros 


# paths 
REF_DIR  = ref
SPEC_DIR = src

SUBDIRS  = $(REF_DIR) $(SPEC_DIR)

# default rule 

.PHONY: default all spec ref clean checkspecc

default: spec 

all: ref spec

spec: chkspecc
	make -C src

ref:
	make -C ref

clean:
	set -e ; for d in $(SUBDIRS); do		\
	  $(MAKE) -C $$d $@ ;	\
	done
	-$(RM) *~

