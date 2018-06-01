//MAELSTROM3D
//Darryl Seligman
//August 4th 2017 v1.0
//operators.cxx file
#include "operators.h"  // player.h must be in the current directory. or use relative or absolute path to it. e.g #include "include/player.h"
//Here we will define the combined difference and average operators
//note that these operators are each 6 1d operators
//note to compile MAELSTROM3D compiler must be linked to operators.cxx
double d1mu2mu1mu2(double *q,int j,int k,int nx, int ny){
  double f=1./8.*(q[((j+1)*ny)+k+1]-q[((j-1)*ny)+k+1]+q[((j+1)*ny)+k-1]\
                  -q[((j-1)*ny)+k-1]+2.*q[((j+1)*ny)+k]-2.*q[((j-1)*ny)+k]);
  return (f);
}
double d2mu1mu1mu2(double *q,int j,int k,int nx, int ny){
  double f=1./8.*(q[((j+1)*ny)+k+1]+2.*q[((j)*ny)+k+1]-q[((j+1)*ny)+k-1]\
                  -2.*q[((j)*ny)+k-1]+q[((j-1)*ny)+k+1]-q[((j-1)*ny)+k-1]);
  return (f);
}
double d1mu2d1mu2(double *q,int j,int k,int nx, int ny){
  double f=1./4.*(q[((j+1)*ny)+k+1]-2.*q[((j)*ny)+k+1]+2.*q[((j+1)*ny)+k]\
                  -4.*q[((j)*ny)+k]+q[((j-1)*ny)+k+1]+q[((j+1)*ny)+k-1]\
                  +2.*q[((j-1)*ny)+k]-2.*q[((j)*ny)+k-1]+q[((j-1)*ny)+k-1]);
  return (f);
}
double d1mu2d2mu1(double *q,int j,int k,int nx, int ny){
  double f=1./4.*(q[((j+1)*ny)+k+1]-q[((j+1)*ny)+k-1]-q[((j-1)*ny)+k+1]\
                  +q[((j-1)*ny)+k-1]);
  return (f);
}
double d2mu1d1mu2(double *q,int j,int k,int nx, int ny){
  double f=1./4.*(q[((j+1)*ny)+k+1]-q[((j-1)*ny)+k+1]-q[((j+1)*ny)+k-1]\
                  +q[((j-1)*ny)+k-1]);
  return (f);
}
double d2mu1d2mu1(double *q,int j,int k,int nx, int ny){
  double f=1./4.*(q[((j+1)*ny)+k+1]-2.*q[((j+1)*ny)+k]+2.*q[((j)*ny)+k+1]\
                  +q[((j+1)*ny)+k-1]+2.*q[((j)*ny)+k-1]-4.*q[((j)*ny)+k]\
                  +q[((j-1)*ny)+k+1]-2.*q[((j-1)*ny)+k]+q[((j-1)*ny)+k-1]);
  return (f);
}
