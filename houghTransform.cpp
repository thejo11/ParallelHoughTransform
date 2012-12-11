// Parallel Hough Transform
// Main Executable Program
//
// CSCI 4576/5576 High Performance Scientific Computing
// Jackie Myrose

// System includes
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sstream>
#include <math.h>

//OpenCV includes
#include <cv.h>
#include <highgui.h>

#include "mpi.h"

#include <fstream>
using std::ofstream;

#include <iostream>
using std::endl;

// Include global variables
#define __MAIN
#include "globals.h"
#undef __MAIN

// User includes
#include "pprintf.h"
#include "pgm.h"

//compile-time defines
//#define _DEBUG       //uncomment for debugging print statements
#define SERIAL       1
#define BLOCKROW     2
#define CHECKERBOARD 3


//
// main
//
int main(int argc, char* argv[]) {
	
    // Initialize MPI
    MPI_Init(&argc, &argv);
	
    // Get the communicator and process information
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    MPI_Get_processor_name( my_name, &my_name_len );
	
    // Initialize the pretty printer
    init_pprintf( rank );
    pp_set_banner( "main" );
    if( rank==0 )
        pprintf( "Welcome to Conway's Game of Life!\n" );
	
    //set options for the game
    /* set to SERIAL, BLOCKROW, or CHECKERBOARD */
    const int part_method = BLOCKROW;
    /* total iterations for simulation  */
    const int iterations = 1;
    /* input file name */
    const char input_file_name[] = "test-500x500.pgm";
    
    //
    // Determine the partitioning
    //
    if( part_method == SERIAL ) {
		nrows = 1;
		ncols = 1;
		if( rank == 0 ) pprintf( "Using SERIAL partitioning\n" );
    }
    else if( part_method == BLOCKROW ) {
		nrows = np;
		ncols = 1;
		if( rank == 0 ) pprintf( "Using BLOCKROW partitioning\n" );
    }
    else if( part_method == CHECKERBOARD ) {
		nrows = sqrt( np );
		ncols = sqrt( np );
		if( rank == 0 ) pprintf( "Using CHECKERBOARD partitioning\n" );
    }
    else {
        if( rank==0 )
            pprintf("Error: unknown partitioning method %d\n",
					part_method );
        MPI_Finalize();
        return 1;
    }
    
    //check partitioning
    if( np != nrows * ncols )
    {
        if( rank==0 )
            pprintf("Error: %ix%i partitioning requires %i np (%i provided)\n",
					nrows, ncols, nrows * ncols, np );
        MPI_Finalize();
        return 1;
    }
	
    // Now, calculate neighbors (N, S, E, W, NW, NE, SW, SE)
    MPI_Comm cart_comm;
    int coords[2];
    int test_coords[2];
    int my_cart_rank;
    int dim_sizes[2];
    int wrap_around[2];
    int reorder;
    int N, S, E, W, NW, NE, SW, SE;
    reorder = 0;
    dim_sizes[0] = nrows;
    dim_sizes[1] = ncols;
    wrap_around[0] = 0;
    wrap_around[1] = 0;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dim_sizes, wrap_around,
					reorder, &cart_comm);
    MPI_Comm_rank(cart_comm, &my_cart_rank);
    MPI_Cart_coords(cart_comm, my_cart_rank, 2, coords);
	
#if defined(_DEBUG)
    pprintf("I'm cart rank %d at cart (%d,%d)\n",
			my_cart_rank, coords[0], coords[1]);
#endif
	
    /* find north neighbor */
    MPI_Cart_shift(cart_comm, 1, 1, &W, &E);
    MPI_Cart_shift(cart_comm, 0, 1, &N, &S);
	
    // get NW rank
    test_coords[0] = coords[0] - 1;
    test_coords[1] = coords[1] - 1;
    if(test_coords[0] >= 0 && test_coords[1] >= 0)
		MPI_Cart_rank(cart_comm, test_coords, &NW);
    else
		NW = MPI_PROC_NULL;
    /* get NE rank */
    test_coords[0] = coords[0] - 1;
    test_coords[1] = coords[1] + 1;
    if(test_coords[0] >= 0 && test_coords[1] < dim_sizes[1])
		MPI_Cart_rank(cart_comm, test_coords, &NE);
    else
		NE = MPI_PROC_NULL;
    /* get SW rank */
    test_coords[0] = coords[0] + 1;
    test_coords[1] = coords[1] - 1;
    if(test_coords[0] < dim_sizes[0] && test_coords[1] >= 0)
		MPI_Cart_rank(cart_comm, test_coords, &SW);
    else
		SW = MPI_PROC_NULL;
    /* get SE rank */
    test_coords[0] = coords[0] + 1;
    test_coords[1] = coords[1] + 1;
    if(test_coords[0] < dim_sizes[0] && test_coords[1] < dim_sizes[1])
		MPI_Cart_rank(cart_comm, test_coords, &SE);
    else
		SE = MPI_PROC_NULL;
	
#if defined(_DEBUG)
    pprintf("my neighbor mapping\n%2d %2d %2d\n%2d %2d %2d\n%2d %2d %2d\n",
			NW, N, NE, W, rank, E, SW, S, SE);
#endif
	
    // set local coords (for readpgm())
    my_col = coords[1];
    my_row = coords[0];
	
    // Read the PGM file. The readpgm() routine reads the PGM file and, based
    // on the previously set nrows, ncols, my_row, and my_col variables, loads
    // just the local part of the field onto the current processor. The
    // variables local_width, local_height, field_width, field_height, as well
    // as the fields (field_a, field_b) are allocated and filled.
    if(!readpgm( input_file_name ))
    {
        if( rank==0 )
            pprintf( "An error occured while reading the pgm file\n" );
        MPI_Finalize();
        return 1;
    }
    
    MPI_Datatype datatype;
    MPI_Type_vector(local_height, local_width, field_width, MPI_CHAR, &datatype);
    
    MPI_Datatype filetype;
	
	int gsizes[2], distribs[2], dargs[2], psizes[2];
    
    gsizes[0] = local_height*nrows;
	gsizes[1] = local_width*ncols;
	
	distribs[0] = MPI_DISTRIBUTE_BLOCK;
	distribs[1] = MPI_DISTRIBUTE_BLOCK;
    
	dargs[0] = MPI_DISTRIBUTE_DFLT_DARG;
	dargs[1] = MPI_DISTRIBUTE_DFLT_DARG;
	
	psizes[0] = nrows;
	psizes[1] = ncols;
	
	MPI_Type_create_darray(np, rank, 2, gsizes, distribs, dargs,
						   psizes, MPI_ORDER_C, MPI_CHAR, &filetype);
	MPI_Type_commit(&filetype);
	
	std::string header = ("P5\n");
	std::string cols;
	std::stringstream out1;
	out1 << local_width*ncols;
	cols = out1.str();
	header.append(cols);
	header.append(" ");
	std::string rows;
	std::stringstream out2;
	out2 << local_height*nrows;
	rows = out2.str();
	header.append(rows);
	header.append("\n255\n");
    char *header2 = (char*)header.c_str();
    
    int offset = strlen(header2);
	
    /* assign 2d pointers to field_a */
    unsigned char **board_ptr_a = (unsigned char **) malloc(field_height * sizeof(unsigned char *));
    for (int i=0; i < field_height; i++) {
		board_ptr_a[i] = field_a + (i*field_width);
    }
    /* assign 2d pointers to field_b */
    unsigned char **board_ptr_b = (unsigned char **) malloc(field_height * sizeof(unsigned char *));
    for (int i=0; i < field_height; i++) {
		board_ptr_b[i] = field_b + (i * field_width);
    }
	
	unsigned char *final_data = (unsigned char*)malloc(field_width*field_height*sizeof(unsigned char));
	
    /** column type */
    MPI_Datatype col_type;
    MPI_Type_vector( local_height, 1, field_width, MPI_CHAR, &col_type);
    MPI_Type_commit(&col_type);
	
    /**
     * Main game loop
     */
    int iteration;
    for (iteration=1; iteration<=iterations; iteration++) {
		
		/**
		 * send rows to neighbors, recv ghost rows
		 */
		MPI_Request requests[16];
		//sending my top row to N
		MPI_Isend( &board_ptr_a[ 1 ][ 1 ], local_width, MPI_CHAR,
				  N, 0, cart_comm, requests );
		//recving my top ghost row from N
		MPI_Irecv( &board_ptr_a[ 0 ][ 1 ], local_width, MPI_CHAR,
				  N, 0, cart_comm, requests+1 );
		//sending my bottom row to S
		MPI_Isend( &board_ptr_a[ field_height-2 ][ 1 ], local_width, MPI_CHAR,
				  S, 0, cart_comm, requests+2 );
		//recving my bottom ghost row from S
		MPI_Irecv( &board_ptr_a[ field_height-1 ][ 1 ], local_width, MPI_CHAR,
				  S, 0, cart_comm, requests+3 );
		
		//sending my left col to W
		MPI_Isend( &board_ptr_a[ 1 ][ 1 ], 1, col_type,
				  W, 0, cart_comm, requests+4 );
		//recving my left ghost row from W
		MPI_Irecv( &board_ptr_a[ 1 ][ 0 ], 1, col_type,
				  W, 0, cart_comm, requests+5 );
		//sending my right col to E
		MPI_Isend( &board_ptr_a[ 1 ][ field_width-2 ], 1, col_type,
				  E, 0, cart_comm, requests+6 );
		//recving my right ghost col from E
		MPI_Irecv( &board_ptr_a[ 1 ][ field_width-1 ], 1, col_type,
				  E, 0, cart_comm, requests+7 );
		
		//sending my NW corner to NW
		MPI_Isend( &board_ptr_a[ 1 ][ 1 ], 1, MPI_CHAR,
				  NW, 0, cart_comm, requests+8 );
		//recving my NW ghost corner from NW
		MPI_Irecv( &board_ptr_a[ 0 ][ 0 ], 1, MPI_CHAR,
				  NW, 0, cart_comm, requests+9 );
		//sending my NE corner to NE
		MPI_Isend( &board_ptr_a[ 1 ][ field_width-2 ], 1, MPI_CHAR,
				  NE, 0, cart_comm, requests+10 );
		//recving my NE ghost corner from NE
		MPI_Irecv( &board_ptr_a[ 0 ][ field_width-1 ], 1, MPI_CHAR,
				  NE, 0, cart_comm, requests+11 );
		//sending my SW corner to SW
		MPI_Isend( &board_ptr_a[ field_height-2 ][ 1 ], 1, MPI_CHAR,
				  SW, 0, cart_comm, requests+12 );
		//recving my SW ghost corner from SW
		MPI_Irecv( &board_ptr_a[ field_height-1 ][ 0 ], 1, MPI_CHAR,
				  SW, 0, cart_comm, requests+13 );
		//sending my SE corner to SE
		MPI_Isend( &board_ptr_a[ field_height-2 ][ field_width-2 ], 1, MPI_CHAR,
				  SE, 0, cart_comm, requests+14 );
		//recving my SE ghost corner from SE
		MPI_Irecv( &board_ptr_a[ field_height-1 ][ field_width-1 ], 1, MPI_CHAR,
				  SE, 0, cart_comm, requests+15 );
		
		MPI_Waitall( 16, requests, MPI_STATUSES_IGNORE );
		
		/**
		 * Perform Hough transform
		 */
		
		
		//First write rank's part of image to file
		MPI_File tempOutFile;
		char fileName[13] = "tempOut";
		char intString[3];
		sprintf(intString, "%d", rank);
		strcat(fileName, intString);
		strcat(fileName, ".pgm");
		
		MPI_File_open(MPI_COMM_SELF, fileName,
					  MPI_MODE_RDWR | MPI_MODE_CREATE,
					  MPI_INFO_NULL, &tempOutFile);
		
		std::string tempHeader = ("P5\n");
		std::string tempCols;
		std::stringstream tempOut1;
		tempOut1 << local_width;
		tempCols = tempOut1.str();
		tempHeader.append(tempCols);
		tempHeader.append(" ");
		std::string tempRows;
		std::stringstream tempOut2;
		tempOut2 << local_height;
		tempRows = tempOut2.str();
		tempHeader.append(tempRows);
		tempHeader.append("\n255\n");
		char *tempHeader2 = (char*)tempHeader.c_str();
		
		int tempOffset = strlen(tempHeader2);
		
		MPI_File_write(tempOutFile, tempHeader2, tempOffset, MPI_CHAR, MPI_STATUS_IGNORE);
		
		MPI_File_set_view(tempOutFile, tempOffset, MPI_CHAR, MPI_CHAR, "native", MPI_INFO_NULL);
		
		MPI_File_write(tempOutFile, &field_a[field_width+1], 1, datatype, MPI_STATUS_IGNORE);
		
		MPI_File_close(&tempOutFile);
		
		//Next do the Canney and Hough transforms on the smaller image
		IplImage* src;
		src=cvLoadImage(fileName, 0);
		IplImage* dst = cvCreateImage( cvGetSize(src), 8, 1 );
		IplImage* color_dst = cvCreateImage( cvGetSize(src), 8, 3 );
		IplImage* final_dst = cvCreateImage( cvGetSize(src), 8, 1 );
		CvMemStorage* storage = cvCreateMemStorage(0);
		CvSeq* lines = 0;
		int i;
		cvCanny( src, dst, 50, 200, 3 );
		cvCvtColor( dst, color_dst, CV_GRAY2BGR );
		lines = cvHoughLines2( dst,
							  storage,
							  CV_HOUGH_STANDARD,
							  1,
							  CV_PI/180,
							  100,
							  0,
							  0 );
		
		for( i = 0; i < MIN(lines->total,100); i++ )
        {
			float* line = (float*)cvGetSeqElem(lines,i);
			float rho = line[0];
			float theta = line[1];
			CvPoint pt1, pt2;
			double a = cos(theta), b = sin(theta);
			double x0 = a*rho, y0 = b*rho;
			pt1.x = cvRound(x0 + 1000*(-b));
			pt1.y = cvRound(y0 + 1000*(a));
			pt2.x = cvRound(x0 - 1000*(-b));
			pt2.y = cvRound(y0 - 1000*(a));
			cvLine( color_dst, pt1, pt2, CV_RGB(255,0,0), 2, 8 );
        }
		
		cvCvtColor(color_dst, final_dst, CV_BGR2GRAY);
		
		cvSaveImage(fileName, final_dst);
		
		FILE *fp = fopen( fileName, "r" );
		char header[10];
		int width, height, depth;
		fscanf( fp, "%6s\n%i %i\n%i\n", header, &width, &height, &depth );
		
		//clear array
		for (int y=0; y<field_width*field_height; y++) {
			*(final_data + y) = 0;
		}
		
		int b, ll, lx, ly;
		
		for( int y=0; y<height; y++ )
		{
			for( int x=0; x<width; x++ )
			{
				b = fgetc(fp);
				if( b == EOF )
				{
					pprintf( "Error: Encountered EOF at [%i,%i]\n", y,x );
					return false;
				}
				// Calculate the local pixels (+1 for ghost row,col)
				lx = x + 1;
				ly = y + 1;
				ll = (ly * field_width + lx );
				final_data[ ll ] = b;
				
			} // for x
		} // for y
		
    }
	
    
    
    MPI_File file;
    char fileName[13] = "output";
    char intString[3];
    sprintf(intString, "%d", 1);
    strcat(fileName, intString);
    strcat(fileName, ".pgm");
    
    MPI_File_open(MPI_COMM_WORLD, fileName,
                  MPI_MODE_RDWR | MPI_MODE_CREATE,
                  MPI_INFO_NULL, &file);
    if(rank == 0){
        //First write header to file
        MPI_File_write(file, header2, offset, MPI_CHAR, MPI_STATUS_IGNORE);
    }
    
    MPI_File_set_view(file, offset, MPI_CHAR, filetype, "native", MPI_INFO_NULL);
    
    MPI_File_write_all(file, &final_data[field_width+1], 1,
                       datatype, MPI_STATUS_IGNORE);
	
    
    MPI_File_close(&file);
	
    // free the board pointers
    if( board_ptr_a != NULL ) free( board_ptr_a );
    if( board_ptr_b != NULL ) free( board_ptr_b );
	
    // Free the fields
    if( field_a != NULL ) free( field_a );
    if( field_b != NULL ) free( field_b );
	
    // Finalize MPI and terminate
    if( rank==0 )
        pprintf( "Terminating normally\n" );
    MPI_Finalize();
    return 0;
}
