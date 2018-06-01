//MAELSTROM3D
//Darryl Seligman
//August 4th 2017 v1.0
//dualtime.cxx file
#include "dualtime.h"  // player.h must be in the current directory. or use relative or absolute path to it. e.g #include "include/player.h"
#include "dissipative.h"
#include "convective.h"
#include <iostream>
using namespace std;
//Here we will define the dual time steppig functions for the JAMESON iteration
//note to compile MAELSTROM3D compiler must be linked to dualtime.cxx
void update_primitive(double *q_np1, double *q_n, double *q_nm1, double *qc, double *qd, \
		      double dt, int nx, int ny){
  for(int i=1;i<nx-1;i++){
    for(int j=1;j<ny-1;j++){
      
      q_np1[i*ny+(j)]=4./3.*q_n[i*ny+(j)]-1./3.*q_nm1[i*ny+(j)]	\
	-2./3.*dt*(qc[i*ny+(j)]-qd[i*ny+(j)]);
      }
    }
  }

		      
void dualtimestep_1(double * rho_np1, double *s1_np1, double *s2_np1,double * rhoE_np1,	\
		    double * d1, double * d2, double * d3, double * d4, \
		    double * rho_n, double *s1_n, double *s2_n, double * rhoE_n, \
		    double * rho_nm1, double *s1_nm1, double *s2_nm1, double * rhoE_nm1, \
		    double * rhouxuy,\
		    double * rhoux2_plus_p, double * rhouy2_plus_p,\
		    double * rhouxE_plus_pux, double * rhouyE_plus_puy,	\
		    double * ux, double *uy, double *E,		      \
		    double * cs,					\
		    int nx, int ny,  double dx, double dy, \
		    double dt, double alpha){
  /*the first iteration of the Jameson scheme will populate the arrays required
    for the next iteration; namely:
    rho_np1 --- density field at the next time step
    s1_np1  --- x momentum at the next time step
    s2_np1  --- y momentum at the next time step
    s3_np1  --- z momentum at the next time step
    d1      --- dissipative residual used at every following iteration
    d2      --- dissipative residual used at every following iteration 
    d3      --- dissipative residual used at every following iteration
    d4      --- dissipative residual used at every following iteration
  */
  //find the appropriate fictitous time step
  double dts = dt*alpha;
  //instantiate arrays to hold the convective residuals
  double * q1=NULL;
  q1 = new double[nx*ny];
  
  double * q2=NULL;
  q2 = new double[nx*ny];
  
  double * q3=NULL;
  q3 = new double[nx*ny];

  double * q4=NULL;
  q4 = new double[nx*ny];

  //call functions from convective.h to build the convective residuals
  convective_rho(q1,s1_n,s2_n,dx,dy,nx,ny);
  convective_s1(q2,rhoux2_plus_p,rhouxuy,dx,dy,nx,ny);
  convective_s2(q3,rhouxuy,rhouy2_plus_p,dx,dy,nx,ny);
  convective_E(q4,rhouxE_plus_pux,rhouyE_plus_puy,dx,dy,nx,ny);
  
  //call functions from dissipative.h to build dissipative residuals
  
  dissipative_rho(d1,ux,uy,rho_np1,rho_n,rho_nm1,s1_np1,s1_n,s1_nm1,	\
		  s2_np1,s2_n,s2_nm1,rhouxuy,rhoux2_plus_p,rhouy2_plus_p, \
		  dx,dy,dts,nx,ny);
  
  dissipative_s1(d2,ux,uy,rho_np1,rho_n,rho_nm1,s1_np1,s1_n,s1_nm1,    \
		 s2_np1,s2_n,s2_nm1,rhouxuy,rhoux2_plus_p,rhouy2_plus_p, \
		 cs,dx,dy,dts,nx,ny);
  
  dissipative_s2(d3,ux,uy,rho_np1,rho_n,rho_nm1,s1_np1,s1_n,s1_nm1,    \
                 s2_np1,s2_n,s2_nm1,rhouxuy,rhoux2_plus_p,rhouy2_plus_p, \
		 cs,dx,dy,dts,nx,ny);
  dissipative_E(d4,ux,uy,E,rho_np1,rho_n,rho_nm1,s1_np1,s1_n,s1_nm1,    \
		s2_np1,s2_n,s2_nm1,rhoE_np1,rhoE_n,rhoE_nm1,\
		rhouxuy,rhoux2_plus_p,rhouy2_plus_p,cs,rhouxE_plus_pux,rhouyE_plus_puy,\
		dx,dy,dts,nx,ny);
  
  //update the density and momentum
  update_primitive(rho_np1,rho_n,rho_nm1,q1,d1,dts,nx,ny);
  update_primitive(s1_np1,s1_n,s1_nm1,q2,d2,dts,nx,ny);
  update_primitive(s2_np1,s2_n,s2_nm1,q3,d3,dts,nx,ny);
  update_primitive(rhoE_np1,rhoE_n,rhoE_nm1,q4,d4,dts,nx,ny);
  //delete arrays containing convective residuals
  delete q1;
  delete q2;
  delete q3;
  delete q4;
  return;
}
void dualtimestep_2(double * rho_np1, double *s1_np1, double *s2_np1,double * rhoE_np1, \
                    double * d1, double * d2, double * d3, double * d4, \
                    double * rho_n, double *s1_n, double *s2_n, double * rhoE_n, \
                    double * rho_nm1, double *s1_nm1, double *s2_nm1, double * rhoE_nm1, \
                    double * rhouxuy,\
                    double * rhoux2_plus_p, double * rhouy2_plus_p,\
                    double * rhouxE_plus_pux, double * rhouyE_plus_puy, \
                    int nx, int ny,  double dx, double dy,		\
                    double dt, double alpha){

  /*the following iterations of the Jameson scheme will populate the arrays required
    for the next iteration; namely:                                       
    rho_np1 --- density field at the next time step                       
    s1_np1  --- x momentum at the next time step                          
    s2_np1  --- y momentum at the next time step                          
    s3_np1  --- z momentum at the next time step                          
    they require as input the dissipative residual calculated in the first fictitous time
    step, namely:
    d1      --- dissipative residual used at every following iteration    
    d2      --- dissipative residual used at every following iteration    
    d3      --- dissipative residual used at every following iteration    
    d4      --- dissipative residual used at every following iteration    
  */
  //find the appropriate fictitous time step
  double dts = dt*alpha;
  //instantiate arrays to hold the convective residuals  
  double * q1=NULL;
  q1 = new double[nx*ny];

  double * q2=NULL;
  q2 = new double[nx*ny];

  double * q3=NULL;
  q3 = new double[nx*ny];

  double * q4=NULL;
  q4 = new double[nx*ny];

  //call functions from convective.h to build the convective residuals
  
  convective_rho(q1,s1_n,s2_n,dx,dy,nx,ny);
  convective_s1(q2,rhoux2_plus_p,rhouxuy,dx,dy,nx,ny);
  convective_s2(q3,rhouxuy,rhouy2_plus_p,dx,dy,nx,ny);
  convective_E(q4,rhouxE_plus_pux,rhouyE_plus_puy,dx,dy,nx,ny);
  
  
  //update the density and momentum
  update_primitive(rho_np1,rho_n,rho_nm1,q1,d1,dts,nx,ny);
  update_primitive(s1_np1,s1_n,s1_nm1,q2,d2,dts,nx,ny);
  update_primitive(s2_np1,s2_n,s2_nm1,q3,d3,dts,nx,ny);
  update_primitive(rhoE_np1,rhoE_n,rhoE_nm1,q4,d4,dts,nx,ny);
  //delete arrays containing convective residuals
  delete q1;
  delete q2;
  delete q3;
  delete q4;
}
