#ifndef VARIABLE_TRANSFORMATION_H    // To make sure you don't declare the function more than once by including the header multiple times.
#define VARIABLE_TRANSFORMATION_H
void u_to_s(double * s, double * u, double * rho, int nx, int ny);
void s_to_u(double *u, double * s, double * rho, int nx, int ny);

#endif
