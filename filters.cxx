//Maelstrom2d
//Darryl Seligman
//October5 2017 v1.0
//filters.cxx file
#include "filters.h"  // player.h must be in the current directory. or use relative or absolute path to it. e.g #include "include/player.h"
#include <iostream>
#include <cmath> 
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
using namespace std;
//Here we will set up filters for adding artificial viscosity
//Using Lele Filters
//note to compile MAELSTROM3D compiler must be linked to boundaries.cxx

double* construct_xkernel( double alpha, double beta, double a, double b, double c, double d, int nx, int ny){


  //create the two kernels
  double * k1=NULL;
  k1 = new double[nx*ny];
  double * k2=NULL;
  k2 = new double[nx*ny];
  for(int j=0;j<nx*ny;j++){
    k1[j]=0.;
    k2[j]=0.;}
  //fill first three rows
  k2[0*ny+0]=15./16.;
  k2[0*ny+1]=4./16.;
  k2[0*ny+2]=-6./16.;
  k2[0*ny+3]=4./16.;
  k2[0*ny+4]=-1./16.;
  
  k2[1*ny+0]=1./16.;
  k2[1*ny+1]=3./4.;
  k2[1*ny+2]=6./16.;
  k2[1*ny+3]=-4./16.;
  k2[1*ny+4]=1./16.;
  
  k2[2*ny+0]=-1./16.;
  k2[2*ny+1]=4./16.;
  k2[2*ny+2]=5./8.;
  k2[2*ny+3]=4./16.;
  k2[2*ny+4]=-1./16.;

  k1[0*ny+0]=1.;
  k1[1*ny+1]=1.;
  k1[2*ny+2]=1.;
  //fill the last three rows

  k2[(nx-1)*ny+(ny-1)]=15./16.;
  k2[(nx-1)*ny+(ny-2)]=4./16.;
  k2[(nx-1)*ny+(ny-3)]=-6./16.;
  k2[(nx-1)*ny+(ny-4)]=4./16.;
  k2[(nx-1)*ny+(ny-5)]=-1./16.;
  
  k2[(nx-2)*ny+(ny-1)]=1./16.;
  k2[(nx-2)*ny+(ny-2)]=3./4.;
  k2[(nx-2)*ny+(ny-3)]=6./16.;
  k2[(nx-2)*ny+(ny-4)]=-4./16.;
  k2[(nx-2)*ny+(ny-5)]=1./16.;

  k2[(nx-3)*ny+(ny-1)]=-1./16.;
  k2[(nx-3)*ny+(ny-2)]=4./16.;
  k2[(nx-3)*ny+(ny-3)]=5./8.;
  k2[(nx-3)*ny+(ny-4)]=4./16.;
  k2[(nx-3)*ny+(ny-5)]=-1./16.;

  k1[(nx-1)*ny+(ny-1)]=1.;
  k1[(nx-2)*ny+(ny-2)]=1.;
  k1[(nx-3)*ny+(ny-3)]=1.;

  //now fill in the rest

  for(int i=3;i<(nx-3);i++){
    k2[i*ny+i]=a;
    k2[i*ny+i-1]=b/2.;
    k2[i*ny+i+1]=b/2.;

    k2[i*ny+i-2]=c/2.;
    k2[i*ny+i+2]=c/2.;

    k2[i*ny+i-3]=d/2.;
    k2[i*ny+i+3]=d/2.;

    k1[i*ny+i]=1.;
    
    k1[i*ny+i-1]=alpha;
    k1[i*ny+i+1]=alpha;

    k1[i*ny+i-2]=beta;
    k1[i*ny+i+2]=beta;
  }
  //  double invk1[nx*ny];
  double * invk1=NULL;
  invk1 = new double[nx*ny];
  int s, i, j;
  cout<<"constructing matrices"<<endl;
  gsl_matrix_view m = gsl_matrix_view_array(k1, nx, ny);
  gsl_matrix_view mk2= gsl_matrix_view_array(k2, nx, ny);
  gsl_matrix_view inv= gsl_matrix_view_array(invk1,nx,ny);
  gsl_permutation * p = gsl_permutation_alloc (nx);
  cout<<"about to invert"<<endl;
  gsl_linalg_LU_decomp (&m.matrix, p, &s);
  gsl_linalg_LU_invert (&m.matrix, p, &inv.matrix);
  cout<<"inversion done"<<endl;
  //  gsl_matrix_mul_elements(&inv.matrix,&mk2.matrix);
  //multiply the two matrixes
  double * kernel=NULL;
  kernel = new double[nx*ny];
  for(int i=0;i<nx*ny;i++){
    kernel[i]=0.;}
  
  gsl_matrix_view Kernel = gsl_matrix_view_array(kernel,nx,ny);
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
                  1.0, &inv.matrix, &mk2.matrix,
                  0.0, &Kernel.matrix);
  
  gsl_permutation_free (p);
  return(kernel);
}
double* construct_ykernel( double alpha, double beta, double a, double b, double c, double d, int nx, int ny){
  double * k1=NULL;
  k1 = new double[nx*ny];
  double * k2=NULL;
  k2 = new double[nx*ny];
  for(int j=0;j<nx*ny;j++){
    k1[j]=0.;
    k2[j]=0.;}
  //fill first three columns
  
  k2[0*ny+0]=15./16.;
  k2[1*ny+0]=4./16.;
  k2[2*ny+0]=-6./16.;
  k2[3*ny+0]=4./16.;
  k2[4*ny+0]=-1./16.;

  k2[0*ny+1]=1./16.;
  k2[1*ny+1]=3./4.;
  k2[2*ny+1]=6./16.;
  k2[3*ny+1]=-4./16.;
  k2[4*ny+1]=1./16.;

  k2[0*ny+2]=-1./16.;
  k2[1*ny+2]=4./16.;
  k2[2*ny+2]=5./8.;
  k2[3*ny+2]=4./16.;
  k2[4*ny+2]=-1./16.;

  k1[0*ny+0]=1.;
  k1[1*ny+1]=1.;
  k1[2*ny+2]=1.;
  //fill the last three rows                                                                       

  k2[(nx-1)*ny+(ny-1)]=15./16.;
  k2[(nx-2)*ny+(ny-1)]=4./16.;
  k2[(nx-3)*ny+(ny-1)]=-6./16.;
  k2[(nx-4)*ny+(ny-1)]=4./16.;
  k2[(nx-5)*ny+(ny-1)]=-1./16.;

  k2[(nx-1)*ny+(ny-2)]=1./16.;
  k2[(nx-2)*ny+(ny-2)]=3./4.;
  k2[(nx-3)*ny+(ny-2)]=6./16.;
  k2[(nx-4)*ny+(ny-2)]=-4./16.;
  k2[(nx-5)*ny+(ny-2)]=1./16.;

  k2[(nx-1)*ny+(ny-3)]=-1./16.;
  k2[(nx-2)*ny+(ny-3)]=4./16.;
  k2[(nx-3)*ny+(ny-3)]=5./8.;
  k2[(nx-4)*ny+(ny-3)]=4./16.;
  k2[(nx-5)*ny+(ny-3)]=-1./16.;

  k1[(nx-1)*ny+(ny-1)]=1.;
  k1[(nx-2)*ny+(ny-2)]=1.;
  k1[(nx-3)*ny+(ny-3)]=1.;
  //now fill in the rest
  for(int i=3;i<(nx-3);i++){
    k2[i*ny+i]=a;
    k2[(i-1)*ny+i]=b/2.;
    k2[(i+1)*ny+i]=b/2.;

    k2[(i-2)*ny+i]=c/2.;
    k2[(i+2)*ny+i]=c/2.;

    k2[(i-3)*ny+i]=d/2.;
    k2[(i+3)*ny+i]=d/2.;

    k1[i*ny+i]=1.;

    k1[(i-1)*ny+i]=alpha;
    k1[(i+1)*ny+i]=alpha;

    k1[(i-2)*ny+i]=beta;
    k1[(i+2)*ny+i]=beta;
  }

  
  int s;
  //double invk1[nx*ny];
  double * invk1=NULL;
  invk1 = new double[nx*ny];
  gsl_matrix_view m = gsl_matrix_view_array(k1, nx, ny);
  gsl_matrix_view mk2= gsl_matrix_view_array(k2, nx, ny);
  gsl_matrix_view inv= gsl_matrix_view_array(invk1,nx,ny);
  gsl_permutation * p = gsl_permutation_alloc (nx);
  gsl_linalg_LU_decomp (&m.matrix, p, &s);
  gsl_linalg_LU_invert (&m.matrix, p, &inv.matrix);
  double * kernel=NULL;
  kernel = new double[nx*ny];
  for(int i=0;i<nx*ny;i++){
    kernel[i]=0.;}
  
  gsl_matrix_view Kernel = gsl_matrix_view_array(kernel,nx,ny);
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
                  1.0, &mk2.matrix, &inv.matrix,
                  0.0, &Kernel.matrix);
  
  gsl_permutation_free (p);
  return(kernel);
  
}
void xfilter(double *x, double* k3, int nx, int ny){

  gsl_matrix_view xmat = gsl_matrix_view_array(x, nx, ny);
  gsl_matrix_view kern = gsl_matrix_view_array(k3, nx, ny);
  double * xfilt=NULL;
  xfilt = new double[nx*ny];
  for(int i=0;i<nx*ny;i++){
    xfilt[i]=0.;}

  gsl_matrix_view Xfilt = gsl_matrix_view_array(xfilt,nx,ny);
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
                  1.0, &kern.matrix, &xmat.matrix,
                  0.0, &Xfilt.matrix);


  for(int i=0;i<nx;i++){
    for(int j=0;j<ny;j++){

      x[i*ny+j]=xfilt[i*ny+j];
    }}
  delete xfilt;
}
void yfilter(double *x, double * k3, int nx, int ny){
  /*mat xmat(nx,ny);
  for(int i=0;i<nx;i++){
    for(int j=0;j<nx;j++){
      xmat(i,j)=x[i*ny+j];}}

  mat xmatf(nx,ny);


  xmatf=xmat*k3;

  for(int i=0;i<nx;i++){
    for(int j=0;j<nx;j++){

      x[i*ny+j]=xmatf(i,j);
      }}*/
  
  gsl_matrix_view xmat = gsl_matrix_view_array(x, nx, ny);
  gsl_matrix_view kern = gsl_matrix_view_array(k3, nx, ny);
  double * xfilt=NULL;
  xfilt = new double[nx*ny];
  for(int i=0;i<nx*ny;i++){
    xfilt[i]=0.;}

  gsl_matrix_view Xfilt = gsl_matrix_view_array(xfilt,nx,ny);
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
                  1.0, &xmat.matrix, &kern.matrix,
                  0.0, &Xfilt.matrix);


  for(int i=0;i<nx;i++){
    for(int j=0;j<ny;j++){

      x[i*ny+j]=xfilt[i*ny+j];
    }}
  delete xfilt;

}

void source_substep2(double *vxo, double *vyo,double *dno, double cv,   \
                     double dx,double dy, double dt, int nx, int ny){

  double * q1=NULL;
  q1 = new double[nx*ny];
  double * q2=NULL;
  q2 = new double[nx*ny];
  for(int i=0;i<nx*ny;i++){
    q1[i]=0.;
    q2[i]=0.;}
  for(int i=1;i<nx-1;i++){
    for(int j=1;j<ny-1;j++){
      if((vxo[(i+1)*ny+j]-vxo[(i)*ny+j])<0.){
        q1[i*ny+j]=cv*dno[i*ny+j]*pow((vxo[(i+1)*ny+j]-vxo[i*ny+j]),2);}
      else{
        q1[i*ny+j]=0.;}
      if((vyo[(i)*ny+j+1]-vyo[(i)*ny+j])<0.){
        q2[i*ny+j]=cv*dno[i*ny+j]*pow((vyo[(i)*ny+j+1]-vyo[i*ny+j]),2);}
      else{
        q2[i*ny+j]=0.;}
    }
  }

  for(int i=2;i<nx-2;i++){
    for(int j=2;j<ny-2;j++){
      vxo[i*ny+j]=vxo[i*ny+j]-dt*((q1[i*ny+j]-q1[(i-1)*ny+j])/  \
                                  (dx*(dno[i*ny+j]+dno[(i-1)*ny+j])/2.));
      vyo[i*ny+j]=vyo[i*ny+j]-dt*((q2[i*ny+j]-q2[(i)*ny+j-1])/        \
                                  (dy*(dno[i*ny+j]+dno[(i)*ny+j-1])/2.));
    }
  }
  delete q1;
  delete q2;
}
