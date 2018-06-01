#ifndef EOS_H    // To make sure you don't declare the function more than once by including the header multiple times.
#define EOS_H

void EOS(double* p, double * cs, double * E, double * rho,\
         double * ux, double * uy, int nx,int ny, double gamma, int EOStype);
void ideal(double * p , double * cs, double * E, double * rho,          \
           double * ux, double * uy, int nx, int ny, double gamma);
void isothermal(double * p, double * cs, double* rho, int nx, int ny);

#endif
