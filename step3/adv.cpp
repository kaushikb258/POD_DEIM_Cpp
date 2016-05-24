#include <iostream>
#include <math.h>
#include "consts.h"
#include "funcs.h"

using namespace std;

//----------------------------------------------

void advance(double *q, double *Lin, double *PTF, int *pvector, int *iloc, int *jloc, int *eqnno)
{

 int i, j, k, l;  
 double *YH2 = new double[imax*jmax]; 
 double *YO2 = new double[imax*jmax];
 double *YH2O = new double[imax*jmax];
 double *T = new double[imax*jmax];
 
 double vi, vip1, vim1, vjp1, vjm1;
 double d2dx2, d2dy2, ddx;
 double CH2, CO2, TT, RR;

//-------------------

  setToZero(Lin,ndof);  
  setToZero(PTF,ndeim);

//-------------------

 k = -1;
 l = -1;
 for (j=0; j<jmax; j++)
 {
  for (i=0; i<imax; i++)
  {
   k++;
   l++;
   YH2[k] = q[l];
   l++;
   YO2[k] = q[l];
   l++;
   YH2O[k] = q[l];
   l++;
   T[k] = q[l];
  }
 } 

//-------------------

  k = -1;
  l = -1;
  for (j=0; j<jmax; j++)
  {
   for (i=0; i<imax; i++)
   {
    k++;
//-------------------
    // H2

    vi = YH2[k]; 

    if(i < imax-1)
    {
     vip1 = YH2[k+1];
    }
    else
    {
     vip1 = 2.0*YH2[k] - YH2[k-1];
    }

    if(i > 0)
    {
     vim1 = YH2[k-1];
    }    
    else
    {
     if(j < jmax/3 || j > 2*jmax/3)
     {
      vim1 = yh2_w;
     }
     else
     {
      vim1 = yh2_i;
     }  
    }

    if(j < jmax-1)
    {
     vjp1 = YH2[k+imax];
    }  
    else
    {
     vjp1 = 2.0*YH2[k] - YH2[k-imax];
    }

    if(j > 0)
    {
     vjm1 = YH2[k-imax];
    }
    else
    {
     vjm1 = 2.0*YH2[k] - YH2[k+imax];
    }

    d2dx2 = (vip1 - 2.0*vi + vim1)/(dx*dx);
    d2dy2 = (vjp1 - 2.0*vi + vjm1)/(dy*dy);
    ddx = (vip1 - vim1)/(2.0*dx);
  
    l++;
    Lin[l] = kappa*(d2dx2 + d2dy2) - w*ddx;

//-------------------
    // O2

    vi = YO2[k]; 

    if(i < imax-1)
    {
     vip1 = YO2[k+1];
    }
    else
    {
     vip1 = 2.0*YO2[k] - YO2[k-1];
    }

    if(i > 0)
    {
     vim1 = YO2[k-1];
    }    
    else
    {
     if(j < jmax/3 || j > 2*jmax/3)
     {
      vim1 = yo2_w;
     }
     else
     {
      vim1 = yo2_i;
     }  
    }

    if(j < jmax-1)
    {
     vjp1 = YO2[k+imax];
    }  
    else
    {
     vjp1 = 2.0*YO2[k] - YO2[k-imax];
    }

    if(j > 0)
    {
     vjm1 = YO2[k-imax];
    }
    else
    {
     vjm1 = 2.0*YO2[k] - YO2[k+imax];
    }

    d2dx2 = (vip1 - 2.0*vi + vim1)/(dx*dx);
    d2dy2 = (vjp1 - 2.0*vi + vjm1)/(dy*dy);
    ddx = (vip1 - vim1)/(2.0*dx);
   
    l++;
    Lin[l] = kappa*(d2dx2 + d2dy2) - w*ddx;

//-------------------
    // H2O

    vi = YH2O[k]; 

    if(i < imax-1)
    {
     vip1 = YH2O[k+1];
    }
    else
    {
     vip1 = 2.0*YH2O[k] - YH2O[k-1];
    }

    if(i > 0)
    {
     vim1 = YH2O[k-1];
    }    
    else
    {
     if(j < jmax/3 || j > 2*jmax/3)
     {
      vim1 = yh2o_w;
     }
     else
     {
      vim1 = yh2o_i;
     }  
    }

    if(j < jmax-1)
    {
     vjp1 = YH2O[k+imax];
    }  
    else
    {
     vjp1 = 2.0*YH2O[k] - YH2O[k-imax];
    }

    if(j > 0)
    {
     vjm1 = YH2O[k-imax];
    }
    else
    {
     vjm1 = 2.0*YH2O[k] - YH2O[k+imax];
    }

    d2dx2 = (vip1 - 2.0*vi + vim1)/(dx*dx);
    d2dy2 = (vjp1 - 2.0*vi + vjm1)/(dy*dy);
    ddx = (vip1 - vim1)/(2.0*dx);
    
    l++;
    Lin[l] = kappa*(d2dx2 + d2dy2) - w*ddx;
    
//-------------------
    // T

    vi = T[k]; 

    if(i < imax-1)
    {
     vip1 = T[k+1];
    }
    else
    {
     vip1 = 2.0*T[k] - T[k-1];
    }

    if(i > 0)
    {
     vim1 = T[k-1];
    }    
    else
    {
     if(j < jmax/3 || j > 2*jmax/3)
     {
      vim1 = T_w;
     }
     else
     {
      vim1 = T_i;
     }  
    }

    if(j < jmax-1)
    {
     vjp1 = T[k+imax];
    }  
    else
    {
     vjp1 = 2.0*T[k] - T[k-imax];
    }

    if(j > 0)
    {
     vjm1 = T[k-imax];
    }
    else
    {
     vjm1 = 2.0*T[k] - T[k+imax];
    }

    d2dx2 = (vip1 - 2.0*vi + vim1)/(dx*dx);
    d2dy2 = (vjp1 - 2.0*vi + vjm1)/(dy*dy);
    ddx = (vip1 - vim1)/(2.0*dx);

    l++;
    Lin[l] = kappa*(d2dx2 + d2dy2) - w*ddx;
    
//---------------------

   }
  }

//----------------------

// Chemistry - DEIM
  
   for (l=0; l<ndeim; l++)
   {
    k = iloc[l] + jloc[l]*imax;

    CH2 = rho*YH2[k]/2.0;
    CO2 = rho*YO2[k]/32.0;
    TT = T[k]*Tref;
    RR = mu1*exp(-mu2/Runiv/TT)*(CH2*CH2)*CO2;
 

    if(RR < 0.0)
    {
     cout<<"RR < 0 "<<RR<<"\n";
    } 

    if(eqnno[l]==0) { PTF[l] = - 2.0*RR*(2.0/rho);}
    else if(eqnno[l]==1) {PTF[l] = - RR*(32.0/rho);}
    else if(eqnno[l]==2) {PTF[l] = 2.0*RR*(18.0/rho);}
    else if(eqnno[l]==3) {PTF[l] = Q*(2.0*RR*18.0/rho);}   
   }
 
//---------------------

 delete [] YH2;
 delete [] YO2;
 delete [] YH2O;
 delete [] T;

}

//----------------------------------------------
