#ifndef CLASSES_H
#define CLASSES_H

#include <iostream>
#include <math.h>
#include "consts.h"

class fluid
{

 private:

  double YH2[imax*jmax];
  double YO2[imax*jmax];
  double YH2O[imax*jmax];
  double T[imax*jmax]; 

 public:

  // constructor 
  fluid()
  {
   init();
  } 

  // destructor
  ~fluid()
  {
  } 

  // copy constructor 
  fluid(fluid &f)
  {
   int i, j, k;
   k = -1;
   for (j=0; j<jmax; j++)
   {
    for (i=0; i<imax; i++)
    {
     k++;
     this->YH2[k] = f.YH2[k];
     this->YO2[k] = f.YO2[k];
     this->YH2O[k] = f.YH2O[k];
     this->T[k] = f.T[k]; 
    } 
   }
  }

  // init()
  void init()
  {
   int i, j, k;
   k = -1;
   for (j=0; j<jmax; j++)
   {
    for (i=0; i<imax; i++)
    {
     k++;
     YH2[k] = 0.0;
     YO2[k] = 0.0;
     YH2O[k] = 0.0;
     T[k] = 300.0/Tref; 
    } 
   }
  }

  // output values

  double getYH2(int k)
  {
   return YH2[k];
  }

  double getYO2(int k)
  {
   return YO2[k];
  }
  
  double getYH2O(int k)
  {
   return YH2O[k];
  }

  double getT(int k)
  {
   return T[k];
  }

  
  // set values

  void setYH2(int k, double val)
  {
   YH2[k] = val;
  }

  void setYO2(int k, double val)
  {
   YO2[k] = val;
  }
  
  void setYH2O(int k, double val)
  {
   YH2O[k] = val;
  }

  void setT(int k, double val)
  {
   T[k] = val;
  }

};

#endif
