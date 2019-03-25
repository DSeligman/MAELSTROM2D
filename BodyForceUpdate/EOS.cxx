//Maelstrom2d
//Darryl Seligman
//October 7th 2017 v1.0
//EOS.cxx file
#include "EOS.h"  // player.h must be in the current directory. or use relative or absolute path to it. e.g #include "include/player.h"
#include <cmath>
//here we will define the EOS functions 
//note to compile MAELSTROM3D compiler must be linked to operators.cxx
void EOS(double * p, double * cs, double * E, double * rho,\
	 double * ux, double * uy, int nx,int ny, double gamma, int EOStype){
  //wrapper function for all equations of state
  /*
    EOStype -- tells the code which Equation of State to use
    EOS = 0 --- ideal equation of state P = (\gamma -1) \epsilon
        where \epsilon is the internal energy density
	\epsilon = E - 1/2 rho |v|^2
     
    EOS = 1 -- isothermal P = cs^2 rho
    
   */
  switch(EOStype){                                                                           
  case 0 :                                                                                   
    ideal(p,cs, E, rho,ux,uy,nx,ny,gamma);
    break;                                                                                   
  case 1 :                                                                                   
    isothermal(p,cs,rho,nx,ny);
    break;                                                                                   
  }
}

void ideal(double * p , double * cs, double * E, double * rho,		\
	   double * ux, double * uy, int nx, int ny, double gamma){
  double epsilon;
  for(int i=0;i<nx;i++){
    for(int j=0;j<ny;j++){
      epsilon = E[i*ny+j]-.5*rho[i*ny+j]*(ux[i*ny+j]*ux[i*ny+j]+uy[i*ny+j]*uy[i*ny+j]);
      p[i*ny+j]=(gamma-1.)*epsilon;
      cs[i*ny+j]=gamma*(gamma-1.)*epsilon;
    }
  }
}
void isothermal(double * p, double * cs, double * rho, int nx, int ny){
  for(int i=0;i<nx;i++){
    for(int j=0;j<ny;j++){
      p[i*ny+j]=cs[i*ny+j]*cs[i*ny+j]*rho[i*ny+j];
    }
  }
}
