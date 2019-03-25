//MAELSTROM3D
//Darryl Seligman
//August 4th 2017 v1.0
//convective.cxx file
#include "convective.h"  // player.h must be in the current directory. or use relative or absolute path to it. e.g #include "include/player.h"
#include <iostream>
#include "operators.h"
using namespace std;
//Here we will define the functions to compute the convective residuals at any dual time step
//note to compile MAELSTROM3D compiler must be linked to convective.cxx
/*these functions will require arrays already allocated in the memory
  to fill with the convective residuals
  q1 -- to compute rho
  q2 -- to compute s1
  q3 -- to compute s2
  q4 -- to compute s3*/
void convective_rho(double * q1, double * sx, double * sy,	\
		    double dx, double dy,			\
		    int nx, int ny ){

  for(int i=1;i<nx-1;i++){
    for(int j=1;j<ny-1;j++){
      //for(int k=1;k<nz-1;k++){
	q1[i*ny+(j)]=1./dx*d1mu2mu1mu2(sx,i,j,nx,ny)\
	  +1./dy*d2mu1mu1mu2(sy,i,j,nx,ny);

    }}
}
void convective_s1(double * q2, double * rhoux2_plus_p, double * rhouxuy,\
		   double *rhoax, double dx, double dy,			\
		   int nx, int ny){
  for(int i=1;i<nx-1;i++){
    for(int j=1;j<ny-1;j++){
      
      q2[i*ny+(j)]=1./dx*d1mu2mu1mu2(rhoux2_plus_p,i,j,nx,ny) \
	+1./dy*d2mu1mu1mu2(rhouxuy,i,j,nx,ny)\
	+1.*mu1mu2mu1mu2(rhoax,i,j,nx,ny);
    }}
}
void convective_s2(double * q3,double * rhouxuy,double * rhouy2_plus_p,\
		   double *rhoay,double dx, double dy,		       \
                   int nx, int ny){
  for(int i=1;i<nx-1;i++){
    for(int j=1;j<ny-1;j++){
      
      q3[i*ny+(j)]=1./dx*d1mu2mu1mu2(rhouxuy,i,j,nx,ny)			\
	+1./dy*d2mu1mu1mu2(rhouy2_plus_p,i,j,nx,ny)			\
	+1.*mu1mu2mu1mu2(rhoay,i,j,nx,ny);
    }}
}
void convective_E(double * q4,double * rhouxE_plus_pux, double *rhouyE_plus_puy, \
		  double dx, double dy,					\
		  int nx, int ny){
  for(int i=1;i<nx-1;i++){
    for(int j=1;j<ny-1;j++){
      
      q4[i*ny+(j)]=1./dx*d1mu2mu1mu2(rhouxE_plus_pux,i,j,nx,ny)		\
	+1./dy*d2mu1mu1mu2(rhouyE_plus_puy,i,j,nx,ny);
    }}
}


