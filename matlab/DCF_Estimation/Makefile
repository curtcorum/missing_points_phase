#
# Make file, should work on unix/linux as is
# otherwise use cmake
#


MEX=mex
MKOCTFILE=mkoctfile
CFLAGS:=$(CFLAGS) -fopenmp -lgomp
LDFLAGS:=$(LDFLAGS) -fopenmp -lgomp
#CFLAGS:=$(CFLAGS) -fopenmp
#LDFLAGS:=$(LDFLAGS) -fopenmp -lgomp

all:
	$(MAKE) matlab 
	$(MAKE) octavemex
	#$(MAKE) octave

matlab: cGrid.mexa64 cInvGridNew.mexa64

octavemex: cGrid.mex cInvGridNew.mex

#eventually make c++ oct code...CAC 160813
#octave: cGrid.oct cInvGridNew.oct

%.mex:%.c
	@if [ -f /`which $(MKOCTFILE)` ] ; then \
		echo "detected octave, building octave mex files for gridding" ; \
		CFLAGS="$(CFLAGS)" LDFLAGS="$(LDFLAGS)" $(MKOCTFILE) --mex -v $< ; \
	else \
		echo "did not detect mkoctfile, NOT building octave mex files for gridding" ; \
	fi

%.mexa64:%.c
	@if [ -f /`which $(MEX)` ] ; then \
		echo "detected matlab, building matlab mex files for gridding" ; \
		$(MEX) -v CFLAGS="\$$CFLAGS $(CFLAGS)" LDFLAGS="\$$LDFLAGS $(LDFLAGS)" $< ; \
	else \
		echo "did not detect matlab, NOT building matlab mex files for gridding" ; \
	fi

clean:
	-rm -f *.mexa64 *.oct *.mex *.o
