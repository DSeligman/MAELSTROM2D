//Maelstrom2d
//Darryl Seligman
//October5 2017 v1.0
//boundaries.cxx file
#include "boundaries.h"  // player.h must be in the current directory. or use relative or absolute path to it. e.g #include "include/player.h"
#include <iostream>
using namespace std;
//Here we will define the boundary condition functions
//These are three dimensional boundaries
//note to compile MAELSTROM3D compiler must be linked to boundaries.cxx
void BC(double *q, int nx, int ny, \
	int xbtype, int ybtype,  int vartype){
  /* the parameter btype tells the code what type of boundary condition to apply          
     btype = 1 ---- Periodic WRAP boundary condition                                       
     btype = 2 ---- Reflecting Boundary condition  */
  
  boundary(q,nx,ny,1,1,2,xbtype,vartype);
  boundary(q,nx,ny,1,2,2,xbtype,vartype);
  boundary(q,nx,ny,2,1,2,ybtype,vartype);
  boundary(q,nx,ny,2,2,2,ybtype,vartype);


    
}



void boundary(double *q, int nx, int ny,int dim, int side, int nghosts, int btype, \
	      int vartype){
  //wrapper function for all  boundary condition
  /*optional flags and default settings are:
    int dim = 1
    int side =1
    int nghosts = 2
    int vartype = 0
    dim tells the code which dimension to perform the wrap BC on
    dim = 1 is x dimension
    dim = 2 is y dimension
    dim = 3 is z dimension
    similarly side tells the code which edge of the 
    side = 1 means 0,1 left side
    side = 2 means nx-1,nx-2 right side (similarly for y z)
    nghosts tells the code how many ghost zones to fill
    the default for this is 2*/
  /* the parameter btype tells the code what type of boundary condition to apply
     btype = 1 ---- Periodic WRAP boundary condition
     btype = 2 ---- Reflecting Boundary condition

     vartype tells the code what type of variable it is applying the boundary condition to
     vartype = 0 ---- any primitive variable besides velocity
     vartype = 1 ---- vx
     vartype = 2 ---- vy
     vartype = 3 ---- vz
  */
  
  //error flags
  if((dim < 1) ||(dim >3)){
    cout<<"wrong number of dimensions in call to periodic boundary condition"<<endl;
    return;
  }
  if((side < 1) ||(side >2)){
    cout<<"wrong number of sides in call to periodic boundary condition"<<endl;
    return;
  }
  if((nghosts < 1) ||(nghosts >2)){
    cout<<"wrong number of ghost cells in call to periodic boundary condition"<<endl;
    return;
  }
  //periodic boundary
  if(btype==1){
    periodic(q,nx,ny,dim,side,nghosts);
  }
  /*if(btype==2){
    switch(vartype){
    case 0 :
      reflect(q,nx,ny,dim,side,nghosts);
      break;
    case 1 :
      reflectux(q,nx,ny,dim,side,nghosts);
      break;
    case 2:
      reflectuy(q,nx,ny,dim,side,nghosts);
      break;
    case 3:
      reflectuz(q,nx,ny,dim,side,nghosts);
      break;
      }}*/

}

void periodic(double *q, int nx, int ny, int dim, int side, int nghosts){
  //apply periodic boundary condition to array q
  //x boundary
  if(dim==1){
    if(side == 1){
      for(int j=0;j<ny;j++){
	//for(int k=0;k<nz;k++){
	q[0*ny+(j)]=q[(nx-4)*ny+(j)];
	q[1*ny+(j)]=q[(nx-3)*ny+(j)];}}
    if(side ==2){
      for(int j=0;j<ny;j++){
        //for(int k=0;k<nz;k++){
          q[(nx-1)*ny+(j)]=q[3*ny+(j)];
          q[(nx-2)*ny+(j)]=q[2*ny+(j)];}}
  }
  
  if(dim==2){
    if(side == 1){
      for(int i=0;i<nx;i++){
        //for(int k=0;k<nz;k++){
	q[i*ny+(0)]=q[i*ny+(ny-4)];
	q[i*ny+(1)]=q[i*ny+(ny-3)];}}
    if(side ==2){
      for(int i=0;i<nx;i++){
        //for(int k=0;k<nz;k++){
          q[i*ny+(ny-1)]=q[i*ny+(3)];
          q[i*ny+(ny-2)]=q[i*ny+(2)];}}
  }



	    
}
void reflect(double *q, int nx, int ny, int nz, int dim, int side, int nghosts){
  // here we implement the REFLECTING BOUNDARY CONDITIONS
  // reflect function reflect all primitive variables besides velocities
  // as described in Stone & Norman 92
  /*all zone centered variables and tangential components of velocity in the
    ghost zones are set equal to the corresponding values for their images among the active zones
    Normal component of velocity is set to zeron on the boundary
    and reflected for the second ghost zone
   */

  if(dim==1){
    if(side == 1){
      for(int j=0;j<ny;j++){
        for(int k=0;k<nz;k++){
          q[0+ny*(j+nz*k)]=q[3+ny*(j+nz*k)];
          q[1+ny*(j+nz*k)]=q[2+ny*(j+nz*k)];}}}
    if(side ==2){
      for(int j=0;j<ny;j++){
        for(int k=0;k<nz;k++){
          q[nx-1+ny*(j+nz*k)]=q[nx-4+ny*(j+nz*k)];
          q[nx-2+ny*(j+nz*k)]=q[nx-3+ny*(j+nz*k)];}}}
  }

  if(dim==2){
    if(side == 1){
      for(int i=0;i<nx;i++){
        for(int k=0;k<nz;k++){
          q[i+ny*(0+nz*k)]=q[i+ny*(3+nz*k)];
          q[i+ny*(1+nz*k)]=q[i+ny*(2+nz*k)];}}}
    if(side ==2){
      for(int i=0;i<nx;i++){
        for(int k=0;k<nz;k++){
          q[i+ny*(ny-1+nz*k)]=q[i+ny*(ny-4+nz*k)];
          q[i+ny*(ny-2+nz*k)]=q[i+ny*(ny-3+nz*k)];}}}
  }

  if(dim==3){
    if(side == 1){
      for(int i=0;i<nx;i++){
        for(int j=0;j<ny;j++){
          q[i+ny*(j+nz*(0))]=q[i+ny*(j+nz*(3))];
          q[i+ny*(j+nz*(1))]=q[i+ny*(j+nz*(2))];}}}
    if(side ==2){
      for(int i=0;i<nx;i++){
        for(int j=0;j<ny;j++){
          q[i+ny*(j+nz*(nz-1))]=q[i+ny*(j+nz*(nz-4))];
          q[i+ny*(j+nz*(nz-2))]=q[i+ny*(j+nz*(nz-3))];}}}
  }

}
void reflectux(double *q, int nx, int ny, int nz, int dim, int side, int nghosts){
  //reflecting BC for x velocity
  if(dim==1){
    if(side == 1){
      for(int j=0;j<ny;j++){
        for(int k=0;k<nz;k++){
          q[0+ny*(j+nz*k)]=0.;
          q[1+ny*(j+nz*k)]=-1.*q[2+ny*(j+nz*k)];}}}
    if(side ==2){
      for(int j=0;j<ny;j++){
        for(int k=0;k<nz;k++){
          q[nx-1+ny*(j+nz*k)]=0.;
          q[nx-2+ny*(j+nz*k)]=-1.*q[nx-3+ny*(j+nz*k)];}}}
  }
  if(dim==2){
    if(side == 1){
      for(int i=0;i<nx;i++){
        for(int k=0;k<nz;k++){
          q[i+ny*(0+nz*k)]=q[i+ny*(3+nz*k)];
          q[i+ny*(1+nz*k)]=q[i+ny*(2+nz*k)];}}}
    if(side ==2){
      for(int i=0;i<nx;i++){
        for(int k=0;k<nz;k++){
          q[i+ny*(ny-1+nz*k)]=q[i+ny*(ny-4+nz*k)];
          q[i+ny*(ny-2+nz*k)]=q[i+ny*(ny-3+nz*k)];}}}
  }

  if(dim==3){
    if(side == 1){
      for(int i=0;i<nx;i++){
        for(int j=0;j<ny;j++){
          q[i+ny*(j+nz*(0))]=q[i+ny*(j+nz*(3))];
          q[i+ny*(j+nz*(1))]=q[i+ny*(j+nz*(2))];}}}
    if(side ==2){
      for(int i=0;i<nx;i++){
        for(int j=0;j<ny;j++){
          q[i+ny*(j+nz*(nz-1))]=q[i+ny*(j+nz*(nz-4))];
          q[i+ny*(j+nz*(nz-2))]=q[i+ny*(j+nz*(nz-3))];}}}
  }

}
void reflectuy(double *q, int nx, int ny, int nz, int dim, int side, int nghosts){

  if(dim==1){
    if(side == 1){
      for(int j=0;j<ny;j++){
        for(int k=0;k<nz;k++){
          q[0+ny*(j+nz*k)]=q[3+ny*(j+nz*k)];
          q[1+ny*(j+nz*k)]=q[2+ny*(j+nz*k)];}}}
    if(side ==2){
      for(int j=0;j<ny;j++){
        for(int k=0;k<nz;k++){
          q[nx-1+ny*(j+nz*k)]=q[nx-4+ny*(j+nz*k)];
          q[nx-2+ny*(j+nz*k)]=q[nx-3+ny*(j+nz*k)];}}}
  }

  if(dim==2){
    if(side == 1){
      for(int i=0;i<nx;i++){
        for(int k=0;k<nz;k++){
          q[i+ny*(0+nz*k)]=0.;
          q[i+ny*(1+nz*k)]=-1.*q[i+ny*(2+nz*k)];}}}
    if(side ==2){
      for(int i=0;i<nx;i++){
        for(int k=0;k<nz;k++){
          q[i+ny*(ny-1+nz*k)]=0.;
          q[i+ny*(ny-2+nz*k)]=-1.*q[i+ny*(ny-3+nz*k)];}}}
  }

  if(dim==3){
    if(side == 1){
      for(int i=0;i<nx;i++){
        for(int j=0;j<ny;j++){
          q[i+ny*(j+nz*(0))]=q[i+ny*(j+nz*(3))];
          q[i+ny*(j+nz*(1))]=q[i+ny*(j+nz*(2))];}}}
    if(side ==2){
      for(int i=0;i<nx;i++){
        for(int j=0;j<ny;j++){
          q[i+ny*(j+nz*(nz-1))]=q[i+ny*(j+nz*(nz-4))];
          q[i+ny*(j+nz*(nz-2))]=q[i+ny*(j+nz*(nz-3))];}}}
  }

}
void reflectuz(double *q, int nx, int ny, int nz, int dim, int side, int nghosts){
  if(dim==1){
    if(side == 1){
      for(int j=0;j<ny;j++){
	for(int k=0;k<nz;k++){
	  q[0+ny*(j+nz*k)]=q[3+ny*(j+nz*k)];
	  q[1+ny*(j+nz*k)]=q[2+ny*(j+nz*k)];}}}
    if(side ==2){
      for(int j=0;j<ny;j++){
	for(int k=0;k<nz;k++){
	  q[nx-1+ny*(j+nz*k)]=q[nx-4+ny*(j+nz*k)];
	  q[nx-2+ny*(j+nz*k)]=q[nx-3+ny*(j+nz*k)];}}}
  }
  
  if(dim==2){
    if(side == 1){
      for(int i=0;i<nx;i++){
	for(int k=0;k<nz;k++){
	  q[i+ny*(0+nz*k)]=q[i+ny*(3+nz*k)];
	  q[i+ny*(1+nz*k)]=q[i+ny*(2+nz*k)];}}}
    if(side ==2){
      for(int i=0;i<nx;i++){
	for(int k=0;k<nz;k++){
	  q[i+ny*(ny-1+nz*k)]=q[i+ny*(ny-4+nz*k)];
	  q[i+ny*(ny-2+nz*k)]=q[i+ny*(ny-3+nz*k)];}}}
  }
  
  if(dim==3){
    if(side == 1){
      for(int i=0;i<nx;i++){
	for(int j=0;j<ny;j++){
	  q[i+ny*(j+nz*(0))]=0.;
	  q[i+ny*(j+nz*(1))]=-1.*q[i+ny*(j+nz*(2))];}}}
    if(side ==2){
      for(int i=0;i<nx;i++){
	for(int j=0;j<ny;j++){
	  q[i+ny*(j+nz*(nz-1))]=0.;
	  q[i+ny*(j+nz*(nz-2))]=-1.*q[i+ny*(j+nz*(nz-3))];}}}
  }
  
}
