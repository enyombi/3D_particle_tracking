/*
Author: Eru K., 5-Dec-2014
modification: 
recieves input .txt file for the dimensions of the container 
holding the particles

generates strings to automatically name povray (*_v.pov, *_p.pov) and voronoi volume (*_voro.txt) files after the name of the packing analyzed  

voro++ libray and function used here were written by:
// Radical Voronoi tessellation example code
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

*/

/*
to execute:
$g++ generateResults.cc -I/usr/local/include/voro++ /usr/local/lib/libvoro++.a -o generateResults.exe
$./generateResults.exe
 */

#include "voro++.hh"
#include <cstring>
#include <cstdio>
#include <cmath> //for pow(,)

using namespace voro;

char *concat(const char *str1, const char *str2){
  char *result = (char *) malloc(strlen(str1)+strlen(str2)+1); //assigns char-type memory to 'result'
  strcpy(result,str1);//copy 'str1' into 'result'
  strcat(result,str2);//add 'str2' to 'result'
  return result;
}

int main(int argc, char **argv) {
  const char *file1 = &(*(*(argv+1)));
  const char *file2 = &(*(*(argv+2)));

  if ( argc != 3 ){
    printf("You input %d files. Please enter 2!\n",argc-1);
    printf("e.g.,\n");
    printf("$ %s <file1> <file2>\n",&(*(*(argv+0))));
    printf("\nwhere,\n");
    printf("<file1> = packing info\n");
    printf("\tparticle# x y z radius\n");
    printf("<file2> = dimensions of container\n");
    printf("\tx_max y_max z_max\n");
    printf("\tx_min y_min z_min\n");
    return 0;
  }

  if( argc == 3 ){
    int i = 0;
    double temp = 0, dimC[6] = {0,0,0,0,0,0}; 

    FILE *fid2 = fopen(file2,"r");
    while( fscanf(fid2,"%lf",&temp) > 0 ){
      dimC[i++] = temp;
    }
    fclose(fid2);

    printf("\n");
    for( int j = 0; j < 6; j++){
      printf("dimC[%d] = %lf\n",j,dimC[j]);
    }

    int numP = 0;
    char tempStr[1000];
    FILE *fid1 = fopen(file1,"r");
    while(fgets(tempStr,sizeof(tempStr),fid1))
      {
	numP++;
	//printf("%d) %s\n",numP++,tempStr);
      }
    fclose(fid1);
    printf("\nnumber of particles = %d\n",numP);

    /*
      ref: http://math.lbl.gov/voro++/doc/refma/n
      Voro++ makes use of internal computational grid of blocks that
      are used to configure the code for maximum efficiency. As 
      discussed on the library website, the best performance is 
      achieved for around 5 particles per block, with anything in the
      range from 3 to 12 giving good performance. Usually the size of
      the grid can be chosen by ensuring that the number of blocks is
      equal to the number of particles divided by 5.
    */ 
    int dimB = roundf((float) pow((numP/5),(1.0/3)));
    printf("particles per computational block(n_x,n_y,n_z) = (number of particles/5)^(1.0/3) = %d\n\n",dimB);

    // Set up constants for the container geometry
    const double x_min=dimC[3],x_max=dimC[0];
    const double y_min=dimC[4],y_max=dimC[1];
    const double z_min=dimC[5],z_max=dimC[2];
    
    // Set up the number of blocks that the container is divided
    // into.
    const int n_x=dimB,n_y=dimB,n_z=dimB;

    const char *suffixV = "_v.pov";
    const char *suffixP = "_p.pov";
    const char *suffixVoro = "_voro.txt";

    char *fileV = concat(file1,suffixV);
    char *fileP = concat(file1,suffixP);
    char *fileVoro = concat(file1,suffixVoro);
    
    printf("The following will be written: \n");
    printf("%s\n",fileV);
    printf("%s\n",fileP);
    printf("%s\n\n",fileVoro);


    // Create a container with the geometry given above, and make it
    // non-periodic in each of the three coordinates. Allocate space for
    // eight particles within each computational block. Import
    // the test packing and output the Voronoi tessellation in POV-Ray formats.
  
    // Create a container for polydisperse particles using the same
    // geometry as above. Import the polydisperse test packing and
    // output the Voronoi radical tessellation in gnuplot and POV-Ray
    // formats.

    //for polydisperse (use container_poly):
    container_poly conp(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,
			false,false,false,8);
    conp.import(argv[1]);
  
    //write the following files:
    conp.draw_cells_pov(fileV);//povray voronoi file
    conp.draw_particles_pov(fileP);//povray particle file
    conp.print_custom("%i %r %v",fileVoro);//text file of voronoi volumes
  }
}
