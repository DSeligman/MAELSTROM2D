//MAELSTROM3D
//Darryl Seligman
//August 4th 2017 v1.0
//dissipative.cxx file
#include "dissipative.h"  // player.h must be in the current directory. or use relative or absolute path to it. e.g #include "include/player.h"
#include <iostream>
#include "operators.h"
using namespace std;
//Here we will define the functions to compute the dissipative residuals at any dual time step
//note to compile MAELSTROM3D compiler must be linked to dissipative.cxx
/*these functions will require arrays already allocated in the memory
  to fill with the dissipative residuals
  d1 -- to compute rho
  d2 -- to compute s1
  d3 -- to compute s2
  d4 -- to compute s3*/
/* note that for these residuals, we need stored values of primitive variables
   at the previous two timesteps; therefore denote
   q_np1 --- q^{n+1}
   q_n   --- q^{n}
   q_nm1 --- q^{n-1}
*/
double h1D1r(double *w_np1,double *w_n, double * w_nm1,	\
	     double *f1_n, double *f2_n,		\
	     int i, int j,\
	     int nx, int ny,\
	     double dx, double dy, double dt){
  double f;
  f= 3./(2.*dt)*d1mu2mu1mu2(w_np1,i,j,nx,ny)			 \
    -2./dt*d1mu2mu1mu2(w_n,i,j,nx,ny)			 \
    +1./(2.*dt)*d1mu2mu1mu2(w_nm1,i,j,nx,ny)		 \
    +1./dx*d1mu2d1mu2(f1_n,i,j,nx,ny)			 \
    +1./dy*d1mu2d2mu1(f2_n,i,j,nx,ny);//		 \
    //+1.*d1mu2mu1mu2(w_np1,i,j,nx,ny) 
  
  return(f);
}
double h2D2r(double *w_np1,double *w_n, double * w_nm1, \
             double *f1_n, double *f2_n,\
             int i, int j,\
             int nx, int ny,\
             double dx, double dy, double dt){
  double f;
  f= 3./(2.*dt)*d2mu1mu1mu2(w_np1,i,j,nx,ny)             \
    -2./dt*d2mu1mu1mu2(w_n,i,j,nx,ny)                    \
    +1./(2.*dt)*d2mu1mu1mu2(w_nm1,i,j,nx,ny)             \
    +1./dx*d2mu1d1mu2(f1_n,i,j,nx,ny)                    \
    +1./dy*d2mu1d2mu1(f2_n,i,j,nx,ny);			 
  
  return(f);
}

void dissipative_rho(double * d1,				\
		     double *ux, double *uy,				\
		     double *rho_np1, double * rho_n, double *rho_nm1,	\
		     double *s1_np1, double * s1_n, double * s1_nm1,\
		     double *s2_np1, double * s2_n, double * s2_nm1,	\
		     double * rhouxuy,					\
		     double *rhoux2_plus_p, double *rhouy2_plus_p,\
		     double * rhoax, double * rhoay, double dx, double dy,double dt, \
		     int nx, int ny){
  
  for(int i=1;i<nx-1;i++){
    for(int j=1;j<ny-1;j++){
      
      d1[i*ny+(j)]=						\
	.5*ux[i*ny+(j)]*						\
	h1D1r(rho_np1,rho_n,rho_nm1,s1_n,s2_n,				\
	      i,j,nx,ny,dx,dy,dt)					\
	+.5*h1D1r(s1_np1,s1_n,s1_nm1,rhoux2_plus_p,rhouxuy,	\
		  i,j,nx,ny,dx,dy,dt)				\
	+.5*uy[i*ny+(j)]*						\
	h2D2r(rho_np1,rho_n,rho_nm1,s1_n,s2_n,				\
	      i,j,nx,ny,dx,dy,dt)					\
	+.5*h2D2r(s2_np1,s2_n,s2_nm1,rhouxuy,rhouy2_plus_p,		\
		  i,j,nx,ny,dx,dy,dt)					\
	+.5*1.*d1mu2mu1mu2(rhoax,i,j,nx,ny) \
	+.5*1.*d2mu1mu1mu2(rhoay,i,j,nx,ny);
      
      
      
    }
  }
}
void dissipative_s1(double *d2,				\
		    double *ux, double *uy,				\
		    double *rho_np1, double * rho_n, double *rho_nm1,	\
		    double *s1_np1, double * s1_n, double * s1_nm1,\
		    double *s2_np1, double * s2_n, double * s2_nm1,    \
		    double * rhouxuy,				       \
		    double *rhoux2_plus_p, double *rhouy2_plus_p,      \
		    double *cs,					       \
		    double * rhoax, double * rhoay,double dx, double dy,double dt, \
		    int nx, int ny){
  
  for(int i=1;i<nx-1;i++){
    for(int j=1;j<ny-1;j++){
      
      d2[i*ny+(j)]=						\
	.5*(ux[i*ny+(j)]*ux[i*ny+(j)]+cs[i*ny+(j)]*cs[i*ny+(j)])*	\
	h1D1r(rho_np1,rho_n,rho_nm1,s1_n,s2_n,				\
	      i,j,nx,ny,dx,dy,dt)+					\
	.5*2.*ux[i*ny+(j)]*						\
	h1D1r(s1_np1,s1_n,s1_nm1,rhoux2_plus_p,rhouxuy,			\
	      i,j,nx,ny,dx,dy,dt)+					\
	.5*ux[i*ny+(j)]*uy[i*ny+(j)]*					\
	h2D2r(rho_np1,rho_n,rho_nm1,s1_n,s2_n,			\
	      i,j,nx,ny,dx,dy,dt)+				\
	.5*uy[i*ny+(j)]*						\
	h2D2r(s1_np1,s1_n,s1_nm1,rhoux2_plus_p,rhouxuy,			\
	      i,j,nx,ny,dx,dy,dt)+					\
	.5*ux[i*ny+(j)]*						\
	h2D2r(s2_np1,s2_n,s2_nm1,rhouxuy,rhouy2_plus_p,			\
	      i,j,nx,ny,dx,dy,dt)\
	+.5*2.*ux[i*ny+(j)]*1.*d1mu2mu1mu2(rhoax,i,j,nx,ny)\
	+.5*uy[i*ny+(j)]*1.*d2mu1mu1mu2(rhoax,i,j,nx,ny)\
	+.5*ux[i*ny+(j)]*1.*d2mu1mu1mu2(rhoay,i,j,nx,ny);

      
    }
  }
}

void dissipative_s2(double *d3,				\
                    double *ux, double *uy,				\
                    double *rho_np1, double * rho_n, double *rho_nm1,	\
                    double *s1_np1, double * s1_n, double * s1_nm1,	\
                    double *s2_np1, double * s2_n, double * s2_nm1,    \
                    double * rhouxuy,				       \
                    double *rhoux2_plus_p, double *rhouy2_plus_p,      \
                    double *cs,					       \
                    double * rhoax, double * rhoay,double dx, double dy, double dt, \
                    int nx, int ny){	  
  for(int i=1;i<nx-1;i++){
    for(int j=1;j<ny-1;j++){
      
      d3[i*ny+(j)]=					\
	.5*ux[i*ny+(j)]*uy[i*ny+(j)]*			\
	h1D1r(rho_np1,rho_n,rho_nm1,s1_n,s2_n,		\
	      i,j,nx,ny,dx,dy,dt)+			\
	.5*uy[i*ny+(j)]*						\
	h1D1r(s1_np1,s1_n,s1_nm1,rhoux2_plus_p,rhouxuy,			\
	      i,j,nx,ny,dx,dy,dt)+					\
	.5*ux[i*ny+(j)]*						\
	h1D1r(s2_np1,s2_n,s2_nm1,rhouxuy,rhouy2_plus_p,			\
	      i,j,nx,ny,dx,dy,dt)+					\
	.5*(uy[i*ny+(j)]*uy[i*ny+(j)]+cs[i*ny+(j)]*cs[i*ny+(j)])*	\
	h2D2r(rho_np1,rho_n,rho_nm1,s1_n,s2_n,				\
	      i,j,nx,ny,dx,dy,dt)+					\
	.5*2*uy[i*ny+(j)]*						\
	h2D2r(s2_np1,s2_n,s2_nm1,rhouxuy,rhouy2_plus_p,			\
	      i,j,nx,ny,dx,dy,dt)\
	+.5*uy[i*ny+(j)]*1.*d1mu2mu1mu2(rhoax,i,j,nx,ny)	\
        +.5*ux[i*ny+(j)]*1.*d1mu2mu1mu2(rhoay,i,j,nx,ny)	\
        +.5*2.*uy[i*ny+(j)]*1.*d2mu1mu1mu2(rhoay,i,j,nx,ny);
      
    }
  }
}
				  
void dissipative_E(double *d4,                         \
		   double *ux, double *uy,				\
		   double *E,\
		   double *rho_np1, double * rho_n, double *rho_nm1,	\
		   double *s1_np1, double * s1_n, double * s1_nm1,	\
		   double *s2_np1, double * s2_n, double * s2_nm1,	\
		   double *rhoE_np1, double * rhoE_n, double * rhoE_nm1,\
		   double * rhouxuy,					\
		   double *rhoux2_plus_p, double *rhouy2_plus_p,	\
		   double * cs,double *rhouxE_plus_pux, double *rhouyE_plus_puy, \
		   double dx, double dy,double dt,		\
		   int nx, int ny){				  
						   

  for(int i=1;i<nx-1;i++){
    for(int j=1;j<ny-1;j++){
      
      d4[i*ny+(j)]=					\
	.5*(ux[i*ny+j]*E[i*ny+j]+ux[i*ny+j]*cs[i*ny+j])*\
	h1D1r(rho_np1,rho_n,rho_nm1,s1_n,s2_n,          \
              i,j,nx,ny,dx,dy,dt)+			\
	.5*(E[i*ny+j]+cs[i*ny+j])*\
	h1D1r(s1_np1,s1_n,s1_nm1,rhoux2_plus_p,rhouxuy,                 \
              i,j,nx,ny,dx,dy,dt)+\
	.5*ux[i*ny+j]*\
	h1D1r(rhoE_np1,rhoE_n,rhoE_nm1,rhouxE_plus_pux,rhouyE_plus_puy,                 \
              i,j,nx,ny,dx,dy,dt)+\
	.5*(uy[i*ny+j]*E[i*ny+j]+uy[i*ny+j]*cs[i*ny+j])*\
        h2D2r(rho_np1,rho_n,rho_nm1,s1_n,s2_n,          \
              i,j,nx,ny,dx,dy,dt)+                      \
	.5*(E[i*ny+j]+cs[i*ny+j])*\
	h2D2r(s2_np1,s2_n,s2_nm1,rhouxuy,rhouy2_plus_p,	\
              i,j,nx,ny,dx,dy,dt)+\
	.5*uy[i*ny+j]*							\
	h1D1r(rhoE_np1,rhoE_n,rhoE_nm1,rhouxE_plus_pux,rhouyE_plus_puy,                 \
              i,j,nx,ny,dx,dy,dt);
      
      
	
    }
  }
}
