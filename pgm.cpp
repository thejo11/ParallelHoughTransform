
// System includes
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <iostream>
#include <fstream>
#include "mpi.h"

// User includes
#include "globals.h"
#include "pprintf.h"

//
// readpgm
//
bool readpgm( const char *filename )
{
    //
    // This procedure reads a PGM file into the local processor.
    //
    // Input: char *filename, name of file to read
    // Returns: True if file read successfully, False otherwise
    //
    // Preconditions:
    //  * global variables nrows, ncols, my_row, my_col must be set
    //
    // Side effects:
    //  * sets global variables local_width, local_height to local image size
    //  * sets global variables field_width, field_height to local field size
    //  * allocates global variables field_a and field_b
	
    //
    pp_set_banner( "pgm:readpgm" );
	
    // Open the file
    if( rank==0 )
        pprintf( "Opening file %s\n", filename );
    FILE *fp = fopen( filename, "r" );
    if( !fp )
    {
        pprintf( "Error: The file '%s' could not be opened.\n", filename );
        return false;
    }
	
    // Read the PGM header, which looks like this:
    //  |P5                magic version number
    //  |900 900           width height
    //  |255               depth
    char header[10];
    int width, height, depth;
    int rv = fscanf( fp, "%6s\n%i %i\n%i\n", header, &width, &height, &depth );
    if( rv != 4 )
    {
        if(rank==0) 
            pprintf( "Error: The file '%s' did not have a valid PGM header\n", 
					filename );
        return false;
    }
    if( rank==0 )
        pprintf( "%s: %s %i %i %i\n", filename, header, width, height, depth );
	
    // Make sure the header is valid
    if( strcmp( header, "P2") )
    {
        if(rank==0) 
            pprintf( "Error: PGM file is not a valid P2 pixmap.\n" );
        return false;
    }
    if( depth != 255 )
    {
        if(rank==0) 
            pprintf( "Error: PGM file has depth=%i, require depth=255 \n", 
					depth );
        return false;
    }
	
    //
    // Make sure that the width and height are divisible by the number of
    // processors in x and y directions
    //
#if defined(_DEBUG)
    pprintf("nrows %d, ncols %d, my_row %d, my_col %d\n",
			nrows, ncols, my_row, my_col);
#endif
	
    if( width % ncols )
    {
        if( rank==0 )
            pprintf( "Error: %i pixel width cannot be divided into %i cols\n", 
                    width, ncols );
        return false;
    }
    if( height % nrows )
    {
        if( rank==0 )
            pprintf( "Error: %i pixel height cannot be divided into %i rows\n",
                    height, nrows );
        return false;
    }
	
    // Divide the total image among the local processors
    local_width = width / ncols;
    local_height = height / nrows;
	
    // Find out where my starting range it
    int start_x = local_width * my_col;
    int start_y = local_height * my_row;
	
    pprintf( "Hosting data for x:%03i-%03i y:%03i-%03i\n", 
            start_x, start_x + local_width,
            start_y, start_y + local_height );
	
    // Create the array!
    field_width = local_width + 2;
    field_height = local_height + 2;
    field_a = (int *)malloc(field_width * field_height * sizeof(int));
	field_b = (int *)malloc(field_width * field_height * sizeof(int));
	
    // start with completely blank board
    for (int y=0; y<field_width*field_height; y++) {
		*(field_a + y) = 0;
		*(field_b + y) = 0;
    }
	
    //
    // Read the data from the file. Save the local data to the local array.
    //
    int b, ll, lx, ly;
	std::string oldB;
	std::string line;
	std::ifstream myfile(filename);
	
	if(myfile.is_open()){
		getline(myfile, line);
		getline(myfile, line);
		getline(myfile, line);
		for( int y=0; y<height; y++ )
		{
			for( int x=0; x<width; x++ )
			{
				// Read the next number
				getline(myfile, oldB);
				
				b = atoi(oldB.c_str());
				
				// If the character is local, then save it!
				if( x >= start_x && x < start_x + local_width &&
				   y >= start_y && y < start_y + local_height )
				{
					// Calculate the local pixels (+1 for ghost row,col)
					lx = x - start_x;
					ly = y - start_y;
					ll = (ly * local_width + lx );
					field_a[ ll ] = b;
				} // save local point
				
			} // for x
		} // for y
	}
	   
	   fclose( fp );
	   
	   pp_reset_banner();
	   return true;
	   
	   }
