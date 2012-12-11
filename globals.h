// Parallel Hough Transform
// Global variable include file
//
// CSCI 4576/5576 High Performance Scientific Computing
// Matthew Woitaszek
// 12 February 2006
//
// Modified by Jackie Myrose

// How it works:
//  * One .cpp file -- usually the one that contains main(), includes this file
//    within #define __MAIN, like this:
//      #define __MAIN
//      #include globals.h
//      #undef __MAIN
//  * The other files just "#include globals.h"  

#ifdef __MAIN
int rank;
int np;
int my_name_len;
char my_name[255];
#else
extern int rank;
extern int np;
extern int my_name_len;
extern char *my_name;
#endif

//
// Conway globals
//
#ifdef __MAIN
int nrows;          // Number of rows in our partitioning
int ncols;          // Number of columns in our partitioning
int my_row;         // My row number
int my_col;         // My column number

// Local logical game size
int local_width;    // Width and height of game on this processor
int local_height;

// Local physical field size
int field_width;        // Width and height of field on this processor
int field_height;       // (should be local_width+2, local_height+2)
unsigned char *field_a = NULL;      // The local data fields
unsigned char *field_b = NULL;

#else
extern int nrows;   
extern int ncols;   
extern int my_row;  
extern int my_col;  

extern int local_width;    
extern int local_height;

extern int field_width;
extern int field_height;
extern unsigned char *field_a;
extern unsigned char *field_b;

#endif

