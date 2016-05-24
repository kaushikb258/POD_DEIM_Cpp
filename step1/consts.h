#ifndef CONSTS_H
#define CONSTS_H

const int imax = 80;
const int jmax = 40;
const int nmax = 600;

const double mu1 = 2.3375e12; 
const double mu2 = 5.6255e3;
//const double mu1 = 6.5e12; 
//const double mu2 = 9.0e3;

const double Runiv = 8.314e0;
const double rho = 1.39e-3;
const double Tref = 2000.0e0;
const double Q = 9800.0e0/Tref;

const double xlen = 1.8e0;
const double ylen = 0.9e0;
const double dx = xlen/((double) imax);
const double dy = ylen/((double) jmax);
const double dt = 5.0e-5;

const double kappa = 2.0e0;
const double w = 50.0e0;

const double yh2_w = 0.0;
const double yo2_w = 0.0;
const double yh2o_w = 0.0;
const double T_w = 300.0/Tref;
const double yh2_i = 0.0282;
const double yo2_i = 0.2259;
const double yh2o_i = 0.0;
const double T_i = 950.0/Tref;

#endif

