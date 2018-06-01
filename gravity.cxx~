//MAELSTROM3D                                                                             
//Darryl Seligman                                                                  
//August 4th 2017 v1.0                                                               
//operators.cxx file                                                                    
#include "gravity.h"
//here we will include functions to update the velocities 
//in the presence of a gravitational potential
//note to compile MAELSTROM3D compiler must be linked to gravity.cxc

void gravity(double * u, int nx, int ny, int nz, double dt, double g){
  for(int i=0;i<nx;i++){
    for(int j=0;j<ny;j++){
      for(int k=0;k<nz;k++){
	u[i+ny*(j+nz*k)]=u[i+ny*(j+nz*k)]+g*dt;
      }}}
}
