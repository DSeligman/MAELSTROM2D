//MAELSTROM3D
//Darryl Seligman
//August 4th 2017 v1.0
//variable_transformation.cxx file
#include "variable_transformation.h"  // player.h must be in the current directory. or use relative or absolute path to it. e.g #include "include/player.h"
#include <iostream>
using namespace std;
//Here we will define the primitive variable transformation functions
//These are three dimensional boundaries
//note to compile MAELSTROM3D compiler must be linked to variable_transformation.cxx
void u_to_s(double * s, double * u, double * rho, int nx, int ny){
  for(int i=0;i<nx;i++){
    for(int j=0;j<ny;j++){
      
      s[i*ny+(j)]=rho[i*ny+(j)]*u[i*ny+(j)];}}
}
void s_to_u(double *u, double * s, double * rho, int nx, int ny){
  for(int i=0;i<nx;i++){
    for(int j=0;j<ny;j++){
      
      //if(rho[i+ny*(j+nz*k)]!=0.){
      u[i*ny+(j)]=s[i*ny+(j)]/rho[i*ny+(j)];
	  //else{
	  //u[i+ny*(j+nz*k)]=u[i+ny*(j+nz*k)];}}}}
      }}
}
