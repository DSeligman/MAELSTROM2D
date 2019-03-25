#ifndef GRAVITY_H    // To make sure you don't declare the function more than once by including the header multiple times.
#define GRAVITY_H
//2d gravity function
void gravity2D(double *u,  int nx, int ny, double dt, double g);
void gravity(double * u, int nx, int ny, int nz, double dt, double g);
#endif
