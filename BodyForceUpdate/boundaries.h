#ifndef BOUNDARIES_H    // To make sure you don't declare the function more than once by including the header multiple times.
#define BOUNDARIES_H

void BC(double *q,double *qi, int nx, int ny,			\
        int xbtype, int ybtype,  int vartype,int nghosts);
void boundary(double *q, double *qi,int nx, int ny,int dim=1, int side =1, int nghosts=2,int btype=1, \
	      int vartype=0);
void periodic(double *q, int nx, int ny, int dim, int side, int nghosts);
void reflect(double *q, int nx, int ny, int nz, int dim, int side, int nghosts);
void reflectux(double *q, int nx, int ny, int nz, int dim, int side, int nghosts);
void reflectuy(double *q, int nx, int ny, int nz, int dim, int side, int nghosts);
void reflectuz(double *q, int nx, int ny, int nz, int dim, int side, int nghosts);
void analytic(double *q,double*qi, int nx, int ny, int dim, int side, int nghosts);
#endif
