#include <iostream>
#include <math.h>
#include <string>
#include "consts.h"
#include "classes.h"
#include "funcs.h"
#include "out_vtk.h"
using namespace std;


int main()
{
 string filename;
 fluid *fl = new fluid;

 filename = "param1/init"; 
 writeData(filename,fl); 

 advance(fl); 
 out_vtk(fl);

 delete fl;

}
