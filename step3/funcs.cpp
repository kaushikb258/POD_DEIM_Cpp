#include <iostream>
#include <math.h>
#include <string>
#include <fstream>
#include <cstring>
#include <iomanip>
#include "consts.h"
#include "funcs.h"

using namespace std;

//----------------------------------------------

void read_qm(double *qm)
{
 int i, j, k, l;
 ifstream f;
 string filename; 
 
 filename = "../step2/pod_basis";
 f.open (filename.c_str(),ios_base::in);  

 k = -1;
 for (i=0; i<ndof; i++)
 {
  for (j=0; j<nmodes; j++)
  {
   k++;
   f>>qm[k];
  } 
 }  

 f.close();
}

//----------------------------------------------

void read_qtil0(double *qtil0)
{
 int i, l;
 ifstream f;
 string filename; 
 
 filename = "../step2/qtil0";
 f.open (filename.c_str(),ios_base::in);  

 for (i=0; i<nmodes; i++)
 {
  f>>qtil0[i]; 
 }  

 f.close();
}

//----------------------------------------------

void read_B(double *B_deim)
{
 int i, j, k, l;
 ifstream f;
 string filename; 
 
 filename = "../step2/B_DEIM";
 f.open (filename.c_str(),ios_base::in);  

 k = -1;
 for (i=0; i<nmodes; i++)
 {
  for (j=0; j<ndeim; j++)
  {
   k++;
   f>>B_deim[k];
  } 
 }  

 f.close();
}

//----------------------------------------------

void read_P(double *P_deim)
{
 int i, j, k, l;
 ifstream f;
 string filename; 
 
 filename = "../step2/P_DEIM";
 f.open (filename.c_str(),ios_base::in);  

 k = -1;
 for (i=0; i<ndof; i++)
 {
  for (j=0; j<ndeim; j++)
  {
   k++;
   f>>P_deim[k];
  } 
 }  

 f.close();
}

//----------------------------------------------

void read_pvector(int *pvector)
{
 int i, l;
 ifstream f;
 string filename; 
 
 filename = "../step2/pvector";
 f.open (filename.c_str(),ios_base::in);  

 for (i=0; i<ndeim; i++)
 {
  f>>pvector[i]; 
 }  

 f.close();
}

//----------------------------------------------

void compute_ijloc(int *pvector, int *iloc, int *jloc, int *eqnno)
{
 int i, j, k, l, kk;
 int nentries; 

 for (k=0; k<ndeim; k++)
 {
  iloc[k] = 0;
  jloc[k] = 0;
  eqnno[k] = 0;
 } 

 nentries = 0;

 for (k=0; k<ndeim; k++)
 {
  kk = 0;
  for (j=0; j<jmax; j++)
  {
   for (i=0; i<imax; i++)
   {
    for (l=0; l<4; l++)
    {
     kk++;
     if(kk == pvector[k])
     {
      iloc[k] = i;
      jloc[k] = j;
      eqnno[k] = l;
      nentries++;
     } 
    }
   }
  }
 }
 
 if(nentries != ndeim) 
 {
  cout<<"nentries = "<<nentries<<"; ndeim = "<<ndeim<<"\n";
 } 

}

//----------------------------------------------

void computeTranspose(double *A, double *AT, int rows, int cols)
{
 
 int i, j, k, l, ll;

 setToZero(AT,rows*cols);
 
 k = -1;
 ll = 0;
 l = ll; 
 for (i=0; i<rows; i++)
 {
  for (j=0; j<cols; j++)
  {
   k++;
   AT[l] = A[k];
   l += rows;    
   if(l > ll + (cols-1)*rows)
   {
    ll ++;
    l = ll;
   }
  
  } 
 } 

}

//----------------------------------------------

void matrixMult(double *A, double *B, double *C, int n1, int n2, int n3)
{
 int i, j, k;
 
 // set C = 0
  setToZero(C,n1*n3);

 for(i=0; i<n1; ++i)
 for(j=0; j<n3; ++j)
 for(k=0; k<n2; ++k)
 {
  C[j+i*n3]+=A[k+i*n2]*B[k*n3+j];
 }

}

//----------------------------------------------

void setToZero(double *A, int nA)
{
 for (int i=0; i<nA; i++)
 {
  A[i] = 0.0;
 }
}

//----------------------------------------------

void makeGTZero(double *A, int nA)
{
 for (int i=0; i<nA; i++)
 {
  A[i] = A[i]>0.0? A[i] : 0.0;
 }
}

//----------------------------------------------

void update(double *qtil, double *VTLin, double *BPTF)
{

 for (int i=0; i<nmodes; i++)
 {
  qtil[i] += dt*(VTLin[i] + BPTF[i]);
 }

}

//----------------------------------------------



