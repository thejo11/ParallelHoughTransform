ParallelHoughTransform
======================

A parallel implementation of the Hough Transform algorithm using MPI.

A sample image is included in this project, which is the image this program was tested with.

Usage:

First, in houghTransform.cpp set part_method to SERIAL, BLOCKROW, or CHECKERBOARD, depending on how you want to partition the image.

Next use the Makefile:
make

Then to run the program:
mpiexec -n {processor count} ./pHough