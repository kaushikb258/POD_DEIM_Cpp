#include <iostream>
#include <math.h>
#include <string>
#include <sstream>
#include <fstream>
#include "consts.h"
#include "classes.h"
#include "funcs.h"
#include "out_vtk.h"
using namespace std;

#define Nx imax
#define Ny jmax


void out_vtk(fluid *fl)
{
 string filename = "output_paraview.vtk";
 string title = "H2O2Flame"; 
 string cell_size_string, node_num_string;
 int node;
 int i, j, k, n;
 stringstream s_node_num, s_cells, s_imax, s_jmax, s_kmax;
 ofstream f;
 int Nz = 1;

 s_node_num << (Nx*Ny);
 s_cells << (Nx-1)*(Ny-1);
 s_imax << Nx;
 s_jmax << Ny;
 s_kmax << Nz;

 // initialize coordinates
 double coord[Nx*Ny*3];
 
 n = 0;
 for (j=0; j<Ny; j++)
 {
 for (i=0; i<Nx; i++)
 {
  coord[n] = ((double) (i) + 0.5)*dx; //x
  coord[n+1] = (double)(Ny)*dy - ((double) (j) + 0.5)*dy; //y
  coord[n+2] = 0.0; //z = 0 as 2D grid
  n += 3;
 }
 }
 

 f.open (filename.c_str(),ios_base::out);
 f<< "# vtk DataFile Version 2.0\n";
 f<< title<<"\n";
 f<< "ASCII\n";
 f<< "DATASET STRUCTURED_GRID\n";
 f<< "DIMENSIONS "<<"\t"<<s_imax.str()<<"\t\t"<<s_jmax.str()<<"\t\t"<<s_kmax.str()<<"\n";
 f<< "POINTS "<<"\t"<<s_node_num.str()<<"\t"<<"double\n";

 n = 0;
 for (j=0; j<Ny; j++)
 {
 for (i=0; i<Nx; i++)
 { 
  f<<coord[n]<<" "<<coord[n+1]<<" "<<coord[n+2]<<"\n";
  n += 3;
 }
 }

 f<< "CELL_DATA "<<"\t"<<s_cells.str()<<"\n";
 f<< "POINT_DATA "<<"\t"<<s_node_num.str()<<"\n";
 
 f<< "SCALARS YH2 double \n";
 f<< "LOOKUP_TABLE default \n";
 for (j=0; j<Ny; j++)
 {
 for (i=0; i<Nx; i++)
 {
  f<< fl->getYH2(i+j*Nx)<<" ";
 } 
 }
 f<<"\n";

 f<< "SCALARS YO2 double \n";
 f<< "LOOKUP_TABLE default \n";
 for (j=0; j<Ny; j++)
 {
 for (i=0; i<Nx; i++)
 {
  f<< fl->getYO2(i+j*Nx)<<" ";
 }
 }
 f<<"\n";
 
 f<< "SCALARS H2O double \n";
 f<< "LOOKUP_TABLE default \n";
 for (j=0; j<Ny; j++)
 {
 for (i=0; i<Nx; i++)
 {
  f<< fl->getYH2O(i+j*Nx)<<" ";
 }
 }
 f<<"\n"; 

 f<< "SCALARS T double \n";
 f<< "LOOKUP_TABLE default \n";
 for (j=0; j<Ny; j++)
 {
 for (i=0; i<Nx; i++)
 {
  f<< fl->getT(i+j*Nx)*Tref<<" ";
 }
 }
 f<<"\n";

 f.close();

} 
