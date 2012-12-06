#
# Makefile
#

########################################################################### 
#
# Compiling and Linking Support
#
########################################################################### 

#
# Note: the variables CC, CFLAGS, LDFLAGS, etc., are usually left to be set
# outside of the makefile. Most system administrators set these variables so
# that compilations work "outside the box" with no developer intervention.
# Because we're compiling MPI programs using specific compilers and specific
# MPI libraries, they're redefined here for each system.
#

# --------------------------------
# Gordon Configuration
# --------------------------------

#
# Compilers
#

# Intel Compilers (Default)
CC = mpicc
CXX = mpicxx

CFLAGS = -g -Wall -shared-intel -DMPICH_IGNORE_CXX_SEEK -I /usr/local/include/opencv -L /usr/local/lib -lopencv_highgui -lopencv_core -lopencv_imgproc
CXXFLAGS = $(CFLAGS)

########################################################################### 
#
# File Declarations
#
########################################################################### 

#
# C and Object Files
#
C_FILES = houghTransform.cpp pgm.cpp pprintf.c
O_FILES = houghTransform.o   pgm.o   pprintf.o

#
# Main Targets
#
all: pHough

pHough: $(O_FILES)
	$(CXX) $(CXXFLAGS) -o pHough $(O_FILES) $(LDFLAGS)

#
# Housekeeping Targets
#
.PHONY: clean
clean:		
	/bin/rm -f core $(O_FILES) pHough output.pgm

#
# Dependencies
#

# All of the object files depend on the globals, so rebuild everything if they
# change!
*.o: globals.h

# Nothing really depends on the pprintf prototypes, but just be safe
*.o: pprintf.h

# Conway depends on PGM utilities
houghTransform.o: pgm.h

