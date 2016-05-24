#ifndef FUNCS_H
#define FUNCS_H

#include <string>

void read_qm(double*);
void read_qtil0(double*);
void read_B(double*);
void read_P(double*);
void read_pvector(int*);
void compute_ijloc(int*,int*,int*,int*);
void computeTranspose(double*,double*,int,int);
void matrixMult(double*,double*,double*,int,int,int);
void setToZero(double*,int);
void makeGTZero(double*,int);
void compute_PTF(double*,int*,double*);
void update(double*,double*,double*);

#endif
