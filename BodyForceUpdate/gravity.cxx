//MAELSTROM3D                                                                             
//Darryl Seligman                                                                  
//August 4th 2017 v1.0                                                               
//operators.cxx file                                                                    
#include "gravity.h"
//here we will include functions to update the velocities 
//in the presence of a gravitational potential
//note to compile MAELSTROM3D compiler must be linked to gravity.cxc

//gravity2D will take just the 2dimensional paramaters
// u is an arbitrary velocity field, and denotes the direction the gravity acts in
// generally for 2D - -- pass uy as double * u -- gravity acts in the y direction
// denote we will have plus g - so sign of the gravity should be negative for downwards

void gravity(double * u, int nx, int ny, int nz, double dt, double g){
  for(int i=0;i<nx;i++){
    for(int j=0;j<ny;j++){
      for(int k=0;k<nz;k++){
	u[i+ny*(j+nz*k)]=u[i+ny*(j+nz*k)]+g*dt;
      }}}
}

void gravity2D(double * u,  int nx, int ny, double dt, double g){
  //
  for(int i=0;i<nx;i++){
    for(int j=0;j<ny;j++){
        u[i*ny+j]=u[i*ny+j]+g*dt;
    }}
}
