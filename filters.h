#ifndef FILTERS_H    // To make sure you don't declare the function more than once by including the header multiple times.
#define FILTERS_H

double* construct_xkernel( double alpha, double beta, double a, double b,\
			   double c, double d, int nx, int ny);
double* construct_ykernel( double alpha, double beta, double a, double b,\
			   double c, double d, int nx, int ny);
void xfilter(double *x, double* k3, int nx, int ny);
void yfilter(double *x, double * k3, int nx, int ny);
void source_substep2(double *vxo, double *vyo,double *dno, double cv,   \
                     double dx,double dy, double dt, int nx, int ny);
#endif
