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

/**
 * prints local game field, including ghost cells
 * (should only be used with small boards)
 */
void print_local_field(unsigned char *local_field) {
	for (int r=0; r<np; r++) {
		if (rank == r) {
			pprintf("local field\n");
			for (int y=0; y<field_height; y++) {
				for (int x=0; x<field_width; x++) {
					printf("%2x ", local_field[ y*field_width + x ]);
				}
				printf("\n");
			}
		}
		fflush(stdout);
		MPI_Barrier(MPI_COMM_WORLD);  //forces serialized printing
	}
}

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
    const int part_method = SERIAL; 
    /* total iterations for simulation  */
    const int iterations = 1;
    /* frequency at which to count bugs */
    const int count_total_freq = 1;
    /* input file name */
    const char input_file_name[] = "vatican-city-900x600.pgm";
    
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
	
	std::string header = ("P2\n");
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
	
    //
    // Count the buggies. Note that we count from [1,1] - [height+1,width+1];
    // we need to ignore the ghost row!
    //
    /*int i=0;
	 for( int y=1; y<local_height+1; y++ )
	 {
	 for( int x=1; x<local_width+1; x++ )
	 {
	 if( field_a[ y * field_width + x ] == ALIVE )
	 {
	 i++;
	 }
	 }
	 }*/
	
    int total;
    //MPI_Allreduce( &i, &total, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
	
    unsigned char cell;
    int num_neighbors = 0;
    int local_left = 0;
	
    /* assign 2d pointers to field_a */
    int **board_ptr_a = (int **) malloc(field_height * sizeof(int *));
    for (int i=0; i < field_height; i++) {
		board_ptr_a[i] = field_a + (i*field_width);
    }
    /* assign 2d pointers to field_b */
    int **board_ptr_b = (int **) malloc(field_height * sizeof(int *));
    for (int i=0; i < field_height; i++) {
		board_ptr_b[i] = field_b + (i * field_width);
    }
	
    /** column type */
    MPI_Datatype col_type;
    MPI_Type_vector( local_height, 1, field_width, MPI_CHAR, &col_type);
    MPI_Type_commit(&col_type);
	
    /**
     * Main game loop
     */
    int iteration;
    for (iteration=1; iteration<=iterations && total>0; iteration++) {
		
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
		 * update stencil section
		 */
		/*local_left = 0; //keep running total
		 for(int y=1; y<local_height+1; y++ ) {
		 for(int x=1; x<local_width+1; x++ ) {*/
		
		/** counting neighbors */
		/*num_neighbors = 0;
		 // top of stencil
		 if( board_ptr_a[ y-1 ][ x-1 ] == ALIVE ) num_neighbors++;
		 if( board_ptr_a[ y-1 ][ x   ] == ALIVE ) num_neighbors++;
		 if( board_ptr_a[ y-1 ][ x+1 ] == ALIVE ) num_neighbors++;
		 // middle of stencil
		 // cell_ptr = cell_ptr + field_width;
		 if( board_ptr_a[ y   ][ x-1 ] == ALIVE ) num_neighbors++;
		 if( board_ptr_a[ y   ][ x+1 ] == ALIVE ) num_neighbors++;
		 // bottom of stencil
		 // cell_ptr = cell_ptr + field_width;
		 if( board_ptr_a[ y+1 ][ x-1 ] == ALIVE ) num_neighbors++;
		 if( board_ptr_a[ y+1 ][ x   ] == ALIVE ) num_neighbors++;
		 if( board_ptr_a[ y+1 ][ x+1 ] == ALIVE ) num_neighbors++;
		 
		 /** actual update section */
		/*cell = board_ptr_a[ y ][ x ];
		 
		 if( num_neighbors < 2 || num_neighbors > 3 ) {
		 board_ptr_b[ y ][ x ] = DEAD;
		 }
		 else if( num_neighbors == 3 ) {
		 board_ptr_b[ y ][ x ] = ALIVE;
		 local_left++;
		 }
		 else {
		 board_ptr_b[ y ][ x ] = cell;
		 if (cell == ALIVE) local_left++;
		 }
		 
		 }
		 } /** end update section */
		
		/** sum and print total bugs every count_total_freq iterations */
		/*if (iteration % count_total_freq == 0) {
#if defined(_DEBUG)
			pprintf( "local number alive after update iteration %d: %d\n",
					iteration, local_left );
#endif
			
			/** get a board count */
			/*MPI_Allreduce(&local_left, &total, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
			
			if (rank == 0) {
				pprintf( "Total number alive after iteration %d: %d\n",
						iteration, total );
			}
		}*/
		
		/**
		 * copy the updated field_b to field_a
		 * (note: can avoid this copy by doing a pointer swap)
		 */
		memcpy ( field_a, field_b, field_width*field_height*sizeof(unsigned char) );
		
    }
	
	char fileName[13] = "output";
	strcat(fileName, ".pgm");
	
	ofstream outdata;
	if(rank == 0){
		
		outdata.open(fileName);
		
		outdata << header;
		
		for(int i = 0; i < local_height*local_width; i++){
			outdata << field_a[i] << endl;
		}
		outdata.close();
	}
	
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