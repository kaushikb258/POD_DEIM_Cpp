#include <iostream>
#include <math.h>
#include <string>
#include <fstream>
#include <iomanip>
#include "consts.h"
#include "classes.h"
#include "funcs.h"

using namespace std;

void advance(fluid *fl)
{
 int n, i, j, k;
 double d2dx2, d2dy2, ddx;
 double vip1, vi, vim1, vjp1, vjm1;
 double conv_diff[4], source[4]; 
 double CH2, CO2, TT, RR;
 double rhs1[imax*jmax], rhs2[imax*jmax], rhs3[imax*jmax], rhs4[imax*jmax];
 double s1[imax*jmax], s2[imax*jmax], s3[imax*jmax], s4[imax*jmax]; 
  
 double oldval, newval; 
 string filename1, filename2;
 ofstream f1, f2;

 filename1 = "param1/snapshot_u";
 filename2 = "param1/snapshot_rhs";
 f1.open (filename1.c_str(),ios_base::out | ios_base::app);  
 f2.open (filename2.c_str(),ios_base::out | ios_base::app); 

 for (n=0; n<nmax; n++)
 {
    
  k = -1;
  for (j=0; j<jmax; j++)
  {
   for (i=0; i<imax; i++)
   {
    k++;

//-------------------
    // H2
    vi = fl->getYH2(k); 

    if(i < imax-1)
    {
     vip1 = fl->getYH2(k+1);
    }
    else
    {
     vip1 = 2.0*fl->getYH2(k) - fl->getYH2(k-1);
    }

    if(i > 0)
    {
     vim1 = fl->getYH2(k-1);
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
     vjp1 = fl->getYH2(k+imax);
    }  
    else
    {
     vjp1 = 2.0*fl->getYH2(k) - fl->getYH2(k-imax);
    }

    if(j > 0)
    {
     vjm1 = fl->getYH2(k-imax);
    }
    else
    {
     vjm1 = 2.0*fl->getYH2(k) - fl->getYH2(k+imax);
    }

    d2dx2 = (vip1 - 2.0*vi + vim1)/(dx*dx);
    d2dy2 = (vjp1 - 2.0*vi + vjm1)/(dy*dy);
    ddx = (vip1 - vim1)/(2.0*dx);
    conv_diff[0] = kappa*(d2dx2 + d2dy2) - w*ddx;
//-------------------
    // O2
    vi = fl->getYO2(k); 

    if(i < imax-1)
    {
     vip1 = fl->getYO2(k+1);
    }
    else
    {
     vip1 = 2.0*fl->getYO2(k) - fl->getYO2(k-1);
    }

    if(i > 0)
    {
     vim1 = fl->getYO2(k-1);
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
     vjp1 = fl->getYO2(k+imax);
    }  
    else
    {
     vjp1 = 2.0*fl->getYO2(k) - fl->getYO2(k-imax);
    }

    if(j > 0)
    {
     vjm1 = fl->getYO2(k-imax);
    }
    else
    {
     vjm1 = 2.0*fl->getYO2(k) - fl->getYO2(k+imax);
    }

    d2dx2 = (vip1 - 2.0*vi + vim1)/(dx*dx);
    d2dy2 = (vjp1 - 2.0*vi + vjm1)/(dy*dy);
    ddx = (vip1 - vim1)/(2.0*dx);
    conv_diff[1] = kappa*(d2dx2 + d2dy2) - w*ddx;
//------------------- 
    // H2O
    vi = fl->getYH2O(k); 

    if(i < imax-1)
    {
     vip1 = fl->getYH2O(k+1);
    }
    else
    {
     vip1 = 2.0*fl->getYH2O(k) - fl->getYH2O(k-1);
    }

    if(i > 0)
    {
     vim1 = fl->getYH2O(k-1);
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
     vjp1 = fl->getYH2O(k+imax);
    }  
    else
    {
     vjp1 = 2.0*fl->getYH2O(k) - fl->getYH2O(k-imax);
    }

    if(j > 0)
    {
     vjm1 = fl->getYH2O(k-imax);
    }
    else
    {
     vjm1 = 2.0*fl->getYH2O(k) - fl->getYH2O(k+imax);
    }

    d2dx2 = (vip1 - 2.0*vi + vim1)/(dx*dx);
    d2dy2 = (vjp1 - 2.0*vi + vjm1)/(dy*dy);
    ddx = (vip1 - vim1)/(2.0*dx);
    conv_diff[2] = kappa*(d2dx2 + d2dy2) - w*ddx;
//------------------- 
    // T
    vi = fl->getT(k); 

    if(i < imax-1)
    {
     vip1 = fl->getT(k+1);
    }
    else
    {
     vip1 = 2.0*fl->getT(k) - fl->getT(k-1);
    }

    if(i > 0)
    {
     vim1 = fl->getT(k-1);
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
     vjp1 = fl->getT(k+imax);
    }  
    else
    {
     vjp1 = 2.0*fl->getT(k) - fl->getT(k-imax);
    }

    if(j > 0)
    {
     vjm1 = fl->getT(k-imax);
    }
    else
    {
     vjm1 = 2.0*fl->getT(k) - fl->getT(k+imax);
    }

    d2dx2 = (vip1 - 2.0*vi + vim1)/(dx*dx);
    d2dy2 = (vjp1 - 2.0*vi + vjm1)/(dy*dy);
    ddx = (vip1 - vim1)/(2.0*dx);
    conv_diff[3] = kappa*(d2dx2 + d2dy2) - w*ddx;
//------------------- 

    // combustion

    CH2 = rho*fl->getYH2(k)/2.0;
    CO2 = rho*fl->getYO2(k)/32.0;
    TT = fl->getT(k)*Tref;
    RR = mu1*exp(-mu2/Runiv/TT)*(CH2*CH2)*CO2;
 
    source[0] = - 2.0*RR*(2.0/rho);
    source[1] = - RR*(32.0/rho);
    source[2] = 2.0*RR*(18.0/rho);
    source[3] = Q*(2.0*RR*18.0/rho);

//-------------------
    
    rhs1[k] = conv_diff[0] + source[0];
    rhs2[k] = conv_diff[1] + source[1]; 
    rhs3[k] = conv_diff[2] + source[2];
    rhs4[k] = conv_diff[3] + source[3];    

    s1[k] = source[0];
    s2[k] = source[1]; 
    s3[k] = source[2];  
    s4[k] = source[3];

//-------------------

   } 
  }

// --------------------
// update

   k = -1;
   for (j=0; j<jmax; j++)
   {
    for (i=0; i<imax; i++)
    {
     k++;

     // H2
     oldval = fl->getYH2(k);
     newval = oldval + rhs1[k]*dt;
     fl->setYH2(k,newval);

     // O2
     oldval = fl->getYO2(k);
     newval = oldval + rhs2[k]*dt;
     fl->setYO2(k,newval);
  
     // H2O
     oldval = fl->getYH2O(k);
     newval = oldval + rhs3[k]*dt;
     fl->setYH2O(k,newval);

     // T
     oldval = fl->getT(k);
     newval = oldval + rhs4[k]*dt;
     fl->setT(k,newval);

     // write snapshots
     if(n%2 == 0)
     {
      f1<<std::fixed<<std::setprecision(16)<<fl->getYH2(k)<<endl;
      f1<<std::fixed<<std::setprecision(16)<<fl->getYO2(k)<<endl;
      f1<<std::fixed<<std::setprecision(16)<<fl->getYH2O(k)<<endl;
      f1<<std::fixed<<std::setprecision(16)<<fl->getT(k)<<endl;

      f2<<std::fixed<<std::setprecision(16)<<s1[k]<<endl;
      f2<<std::fixed<<std::setprecision(16)<<s2[k]<<endl;
      f2<<std::fixed<<std::setprecision(16)<<s3[k]<<endl;
      f2<<std::fixed<<std::setprecision(16)<<s4[k]<<endl;
     }
  
    } 
   } 

//---------------------
  
 }

  f1.close();
  f2.close();

}

//----------------------------------------------

void writeData(string filename, fluid *fl)
{
 int i, j, k;
 ofstream f;
 f.open (filename.c_str(),ios_base::out);  

 k = -1;
 for (j=0; j<jmax; j++)
 {
  for (i=0; i<imax; i++)
  {
   k++;
   f<<std::fixed<<std::setprecision(16)<<fl->getYH2(k)<<endl;
   f<<fl->getYO2(k)<<endl;
   f<<fl->getYH2O(k)<<endl;
   f<<fl->getT(k)<<endl;
  } 
 }  

 f.close();
}

//----------------------------------------------

