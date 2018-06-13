//Version to computer a skewed shear flow - skewed by 45 degrees as in Karim Sharif notes from May 1st 2018

//MAELSTROM 2D v1.1
//Vorticity Conserving Scheme for the Euler Equations in 3 Dimensions
//v1.1 06/01/2018
//The set up here is a 1d Shear Wave
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <sstream>
#include <time.h>
#include <cstdlib> //for drand48 
#include "operators.h"
#include "boundaries.h"
#include "convective.h"
#include "dissipative.h"
#include "dualtime.h"
#include "variable_transformation.h"
#include "gravity.h"
#include "EOS.h"
#include "filters.h"
using namespace std;

/*Note on indexing
  we use flat arrays to keep contiguous memory allocation
  therefore
  Flat[x + WIDTH * (y + DEPTH * z)] = Original[x, y, z]
  therefore for example, with nx ny and nz and indices i j k
  for arrays
  double * dno[nx*ny*nz]
  double DENSITY[nx][ny][nz]
  dno[i+ny*(j+nz*k)]=DENSITY[i][j][k]
*/
/* MAELSTROM3D requires the input file
   params.txt
 */
int main()
{
  //read in from the input file
  clock_t t;
  t = clock();
  string line;
  ifstream myfile ("params_SkewShear.txt");
  double * inputs=NULL;
  inputs = new double[20];
  int nx, ny, nsteps,xbtype,ybtype,EOStype,nghosts,nxmghost,nymghost,viscosity;
  
  double dx, dy,g;
  double extentx,extenty,xstart,ystart,cadence;
  int i=0;

  if (myfile.is_open())
    {
      while ( getline (myfile,line) )
        {
          istringstream ss(line);
          string name;
          ss >> name >> inputs[i] ;
	  i++;
        }
      myfile.close();
    }
  
  else cout << "Unable to open file";
  nxmghost  = inputs[0];
  nymghost  = inputs[1];
  extentx  = inputs[2];
  extenty  = inputs[3];
  //extentx = extentx*2.;
  //extenty=extenty*2.;
    
  nsteps = inputs[4];
  xbtype = inputs[10];
  ybtype = inputs[11];
  g = inputs[12];
  xstart = inputs[13];
  ystart = inputs[14];
  cadence = inputs[15];
  EOStype = inputs[16];
  nghosts = inputs[17];
  viscosity = inputs[18];
  double gamma = 5./3.;
  //int nsteps = 10;

  //instantiate the grid
  nx = nxmghost+2*nghosts;
  ny=nymghost+2*nghosts;
  //nx=nxmghost;
  //ny=nymghost;
  bool have_E=false;
  if(EOStype == 0){
    have_E = true;}

  //dx=extentx/((double)nx-(double)1);
  //dy=extenty/((double)ny-(double)1);
  dx=extentx/((double)nxmghost);
  dy=extenty/((double)nymghost);
  
  double * xg=NULL;
  xg = new double[nx];


  double * yg=NULL;
  yg = new double[ny];


  
  for(int i=0;i<nx;i++){
    //xg[i]=xstart+dx*i;
    xg[i]=xstart+dx*(i+.5-nghosts);
    //xg[i]=xstart+dx*(i+.5);
  }


  for(int j=0;j<ny;j++){
    //yg[j]=ystart+dy*j;
    yg[j]=ystart+dy*(j+.5-nghosts);
    //yg[j]=ystart+dy*(j+.5);
  }

  cout<<"Grid Constructed"<<endl;
  cout<<"Nx with ghosts = "<<nx<<" on domain "<< xg[0]<<"  up   to   "<<xg[nx-1]<<endl;
  cout<<"Ny with ghosts = "<<ny<<" on domain "<< yg[0]<<"  up   to   "<<yg[ny-1]<<endl;
  
  cout<<"Nx without ghosts = "<<nxmghost<<" on domain "<< xg[nghosts]<<"  up   to   "<<xg[nx-1-nghosts]<<endl;
  cout<<"Ny without ghosts= "<<nymghost<<" on domain "<< yg[nghosts]<<"  up   to   "<<yg[ny-1-nghosts]<<endl;
  double dt=10.;
  //double cs = 1.; //sound speed
  //instantiate density and velocities and pressure
  double * rho=NULL;
  rho = new double[nx*ny];
  
  double * ux=NULL;
  ux = new double[nx*ny];
  
  double * uy=NULL;
  uy = new double[nx*ny];
  
  double * rhoi=NULL;
  rhoi = new double[nx*ny];

  double * uxi=NULL;
  uxi = new double[nx*ny];

  double * uyi=NULL;
  uyi = new double[nx*ny];
  
  
  double * p=NULL;
  p = new double[nx*ny];
  
  //internal energy 
  double * E = NULL;
  E= new double[nx*ny];
  
  double * Ei = NULL;
  Ei= new double[nx*ny];
  
  double * cs = NULL;
  cs=new double[nx*ny];
  //initial conditions

  double pi = 3.14159265358979323846;
  double M=2.*pi;
  
  for(int i=0;i<nx;i++){
    for(int j=0;j<ny;j++){

      rho[i*ny+j]=1.0;
      //ux[i*ny+j]=.1*cos(-1.*xg[i])*cos(yg[j])/sqrt(2.);
      ux[i*ny+j]=.1*cos(-1.*xg[i]+yg[j])/sqrt(2.);
      uy[i*ny+j]=.1*cos(-1.*xg[i]+yg[j])/sqrt(2.);
      //uy[i*ny+j]=.1*cos(-1.*xg[i])*cos(yg[j])/sqrt(2.);
      E[i*ny+j]=1.0/(gamma*(gamma-1.));
      p[i*ny+j]=0.;
      cs[i*ny+j]=1.;
    }}
  if(have_E==false){
    for(int i=0;i<nx*ny;i++){
      cs[i]=1.0;}}
  
  //EOS call
  EOS(p,cs,E,rho,ux,uy,nx,ny,gamma,EOStype);
  

  for(int i=0;i<nx*ny;i++){
    Ei[i]=E[i];
    uxi[i]=ux[i];
    uyi[i]=uy[i];
    rhoi[i]=rho[i];
  }

  //CFL
  double speed=0.;
  for(int i=0;i<nx;i++){
    for(int j=0;j<ny;j++){

	double speed1=sqrt(ux[i*ny+j]*ux[i*ny+j]+\
			   uy[i*ny+j]*uy[i*ny+j]+\
			   cs[i*ny+j]*cs[i*ny+j]);
	if(speed1>=speed){
	  speed=speed1;}
      }
    }

  double dtTest=dx/speed;
  cout<<"speedmax"<<speed<<endl;
  if(dtTest<dt){
    dt=dtTest/3.;}
  cout<<"initial timestep:"<<dt<<endl;

  double totalT=0;
  //set readout cadence based on input
  double cad = cadence/(dt);
  int modo = (int) cad;
  cout<<"output files generated every  "<<modo<<" timesteps"<<endl;
  cout<<"corresponding to a readout cadence of "<<modo*dt<<" time units"<<endl;

  ofstream outputFile("output_SkewShear/rho.txt"); // output filestream                           
  ofstream outputFile2("output_SkewShear/ux.txt");    //read out the initial velocity flow           
  ofstream outputFile3("output_SkewShear/uy.txt"); ////////                                          
  ofstream outputFile4("output_SkewShear/E.txt");

  for(int i=0;i<nx*ny;i++){

    outputFile <<rho[i]<<",";
    outputFile2 <<ux[i]<<",";
    outputFile3 <<uy[i]<<",";
    if(have_E==false){
    outputFile4 <<E[i]<<",";
    }
  }
  outputFile<<endl;
  outputFile2<<endl;
  outputFile3<<endl;
  if(have_E==false){
    outputFile4<<endl;}
  //alpha values for JAMESON Scheme
  double alpha1=0.6;
  double alpha2=0.6;
  double alpha3=1.;

  //implement lele filter
  double *kernelx;
  double *kernely; 
  if(viscosity == 1){
    //constants for kernels for filtering                                                               
    double alpha = 0.6522474;
    double beta = 0.1702929;
    double a=0.9891856;
    double b=1.321180;
    double c=0.3333548;
    double d=0.001359850;
    //double *kernelx;
    kernelx=construct_xkernel(  alpha,  beta,  a,  b,  c, d,  nx,  ny);
    //double *kernely;
    kernely=construct_ykernel(  alpha,  beta,  a,  b,  c, d,  nx,  ny);
    cout<<"kernels created"<<endl;}


  //instantiate arrays for density and momentum at the next time step and previous
  double * rho_np1=NULL;
  rho_np1 = new double[nx*ny];
  
  double * rho_nm1=NULL;
  rho_nm1 = new double[nx*ny];
  
  double * s1_np1=NULL;
  s1_np1 = new double[nx*ny];
  
  double * s1=NULL;
  s1 = new double[nx*ny];
  
  double * s1_nm1=NULL;
  s1_nm1 = new double[nx*ny];

  double * s2_np1=NULL;
  s2_np1 = new double[nx*ny];
  
  double * s2=NULL;
  s2 = new double[nx*ny];
  
  double * s2_nm1=NULL;
  s2_nm1 = new double[nx*ny];
  
  double * rhoE_np1=NULL;
  rhoE_np1 = new double[nx*ny];

  double * rhoE=NULL;
  rhoE = new double[nx*ny];
  
  double * rhoE_nm1=NULL;
  rhoE_nm1 = new double[nx*ny];
  

  
  u_to_s(s1_np1,ux,rho,nx,ny);
  u_to_s(s2_np1,uy,rho,nx,ny);
  u_to_s(rhoE_np1,E,rho,nx,ny);
  
  u_to_s(s1_nm1,ux,rho,nx,ny);
  u_to_s(s2_nm1,uy,rho,nx,ny);
  u_to_s(rhoE_nm1,E,rho,nx,ny);

  for(int i=0;i<nx;i++){
    for(int j=0;j<ny;j++){

        rho_nm1[i*ny+(j)]=rho[i*ny+(j)];
	rho_np1[i*ny+(j)]=rho[i*ny+(j)];
	
      }}

  //instantiate arrays for the dissipative residuals
  double * d1=NULL;
  d1 = new double[nx*ny];

  double * d2=NULL;
  d2 = new double[nx*ny];

  double * d3=NULL;
  d3 = new double[nx*ny];
  
  double * d4=NULL;
  d4 = new double[nx*ny];

  //other needed arrays
  double *rhouxuy =NULL;
  rhouxuy = new double[nx*ny];
  

  double *rhoux2_plus_p =NULL;
  rhoux2_plus_p = new double[nx*ny];

  double *rhouy2_plus_p =NULL;
  rhouy2_plus_p = new double[nx*ny];

  double *rhouxE_plus_pux =NULL;
  rhouxE_plus_pux = new double[nx*ny];

  double *rhouyE_plus_puy =NULL;
  rhouyE_plus_puy = new double[nx*ny];


  //time loop for the simulation
  for(int k=0;k<nsteps;k++){
    //update time
    totalT=totalT+dt;

    //update velocity by gravity
    //gravity(uz, nx, ny, nz, dt, g);
      //boundary conditions
    BC(ux,uxi,nx,ny,xbtype,ybtype,1,nghosts);
    BC(uy,uyi,nx,ny,xbtype,ybtype,1,nghosts);

    BC(rho,rhoi,nx,ny,xbtype,ybtype,0,nghosts);
    BC(E,Ei,nx,ny,xbtype,ybtype,0,nghosts);

    u_to_s(s1,ux,rho,nx,ny);
    u_to_s(s2,uy,rho,nx,ny);
    u_to_s(rhoE,E,rho,nx,ny);
    
    for(int i=0;i<nx;i++){
      for(int j=0;j<ny;j++){

	rhouxuy[i*ny+j]=rho[i*ny+j]*ux[i*ny+j]*uy[i*ny+j];
	rhoux2_plus_p[i*ny+j]=rho[i*ny+j]*ux[i*ny+j]*ux[i*ny+j]+p[i*ny+j];
	
	rhouy2_plus_p[i*ny+j]=rho[i*ny+j]*uy[i*ny+j]*uy[i*ny+j]+p[i*ny+j];
	rhouxE_plus_pux[i*ny+j]=rho[i*ny+j]*E[i*ny+j]*ux[i*ny+j]+p[i*ny+j]*ux[i*ny+j];
	rhouyE_plus_puy[i*ny+j]=rho[i*ny+j]*uy[i*ny+j]*E[i*ny+j]+p[i*ny+j]*uy[i*ny+j];

	}}

    dualtimestep_1(rho_np1,s1_np1,s2_np1,rhoE_np1,d1,d2,d3,d4,	\
		   rho,s1,s2,rhoE,rho_nm1,s1_nm1,s2_nm1,rhoE_nm1,\
		   rhouxuy,rhoux2_plus_p,rhouy2_plus_p,rhouxE_plus_pux,rhouyE_plus_puy,\
		   ux,uy,E,cs,nx,ny,dx,dy,dt,alpha1);
    
    

    s_to_u(ux,s1_np1,rho_np1,nx,ny);
    s_to_u(uy,s2_np1,rho_np1,nx,ny);
    s_to_u(E,rhoE_np1,rho_np1,nx,ny);

    BC(ux,uxi,nx,ny,xbtype,ybtype,1,nghosts);
    BC(uy,uyi,nx,ny,xbtype,ybtype,2,nghosts);
    BC(E,Ei,nx,ny,xbtype,ybtype,3,nghosts);
    BC(rho_np1,rhoi,nx,ny,xbtype,ybtype,0,nghosts);
    //changed made here rho -> rho_np1
    u_to_s(s1_np1,ux,rho_np1,nx,ny);
    u_to_s(s2_np1,uy,rho_np1,nx,ny);
    u_to_s(rhoE_np1,E,rho_np1,nx,ny);
    //EOS call                                                                                                    
    EOS(p,cs,E,rho_np1,ux,uy,nx,ny,gamma,EOStype);
    
    for(int i=0;i<nx;i++){
      for(int j=0;j<ny;j++){

	  rhouxuy[i*ny+j]=rho_np1[i*ny+j]*ux[i*ny+j]*uy[i*ny+j];
	  rhoux2_plus_p[i*ny+j]=rho_np1[i*ny+j]*ux[i*ny+j]*ux[i*ny+j]+p[i*ny+j];

	  rhouy2_plus_p[i*ny+j]=rho_np1[i*ny+j]*uy[i*ny+j]*uy[i*ny+j]+p[i*ny+j];
	  rhouxE_plus_pux[i*ny+j]=rho_np1[i*ny+j]*E[i*ny+j]*ux[i*ny+j]+p[i*ny+j]*ux[i*ny+j];
	  rhouyE_plus_puy[i*ny+j]=rho_np1[i*ny+j]*uy[i*ny+j]*E[i*ny+j]+p[i*ny+j]*uy[i*ny+j];
	  
      }}

		   
    dualtimestep_2(rho_np1,s1_np1,s2_np1,rhoE_np1,d1,d2,d3,d4,\
		   rho,s1,s2,rhoE,rho_nm1,s1_nm1,s2_nm1,rhoE_nm1,		\
		   rhouxuy,rhoux2_plus_p,rhouy2_plus_p,rhouxE_plus_pux,rhouyE_plus_puy,\
		   nx,ny,dx,dy,dt,alpha2);


    
    s_to_u(ux,s1_np1,rho_np1,nx,ny);
    s_to_u(uy,s2_np1,rho_np1,nx,ny);
    s_to_u(E,rhoE_np1,rho_np1,nx,ny);
    BC(ux,uxi,nx,ny,xbtype,ybtype,1,nghosts);
    BC(uy,uyi,nx,ny,xbtype,ybtype,2,nghosts);
    BC(E,Ei,nx,ny,xbtype,ybtype,0,nghosts);
    BC(rho_np1,rhoi,nx,ny,xbtype,ybtype,0,nghosts);

    u_to_s(s1_np1,ux,rho_np1,nx,ny);
    u_to_s(s2_np1,uy,rho_np1,nx,ny);
    u_to_s(rhoE_np1,E,rho_np1,nx,ny);

    EOS(p,cs,E,rho_np1,ux,uy,nx,ny,gamma,EOStype);


    for(int i=0;i<nx;i++){
      for(int j=0;j<ny;j++){

	  rhouxuy[i*ny+j]=rho_np1[i*ny+j]*ux[i*ny+j]*uy[i*ny+j];
          rhoux2_plus_p[i*ny+j]=rho_np1[i*ny+j]*ux[i*ny+j]*ux[i*ny+j]+p[i*ny+j];

          rhouy2_plus_p[i*ny+j]=rho_np1[i*ny+j]*uy[i*ny+j]*uy[i*ny+j]+p[i*ny+j];
          rhouxE_plus_pux[i*ny+j]=rho_np1[i*ny+j]*E[i*ny+j]*ux[i*ny+j]+p[i*ny+j]*ux[i*ny+j];
          rhouyE_plus_puy[i*ny+j]=rho_np1[i*ny+j]*uy[i*ny+j]*E[i*ny+j]+p[i*ny+j]*uy[i*ny+j];
	  
	
      }}

    dualtimestep_2(rho_np1,s1_np1,s2_np1,rhoE_np1,d1,d2,d3,d4,\
                   rho,s1,s2,rhoE,rho_nm1,s1_nm1,s2_nm1,rhoE_nm1,               \
                   rhouxuy,rhoux2_plus_p,rhouy2_plus_p,rhouxE_plus_pux,rhouyE_plus_puy,\
                   nx,ny,dx,dy,dt,alpha3);

    s_to_u(ux,s1_np1,rho_np1,nx,ny);
    s_to_u(uy,s2_np1,rho_np1,nx,ny);
    s_to_u(E,rhoE_np1,rho_np1,nx,ny);

    if(viscosity == 1){
      // now filter each primitive variable                                                             
      xfilter(ux, kernelx,  nx,  ny);
      yfilter(ux, kernely,  nx,  ny);
      xfilter(uy, kernelx,  nx,  ny);
      yfilter(uy, kernely,  nx,  ny);
      xfilter(rho_np1, kernelx,  nx,  ny);
      yfilter(rho_np1, kernely,  nx,  ny);
      xfilter(E, kernelx,  nx,  ny);
      yfilter(E, kernely,  nx,  ny);}
    if(viscosity == 2){
      double cv = 6.;
      source_substep2(ux,uy,rho,cv,dx,dy,dt,nx,ny);}
    BC(ux,uxi,nx,ny,xbtype,ybtype,1,nghosts);
    BC(uy,uyi,nx,ny,xbtype,ybtype,2,nghosts);
    BC(E,Ei,nx,ny,xbtype,ybtype,0,nghosts);
    BC(rho_np1,rhoi,nx,ny,xbtype,ybtype,0,nghosts);

    u_to_s(s1_np1,ux,rho_np1,nx,ny);
    u_to_s(s2_np1,uy,rho_np1,nx,ny);
    u_to_s(rhoE_np1,E,rho_np1,nx,ny);

    EOS(p,cs,E,rho_np1,ux,uy,nx,ny,gamma,EOStype);
    ////////
    


    //refill the arrays
    for(int i=0;i<nx;i++){
      for(int j=0;j<ny;j++){

	  
	  rho_nm1[i*ny+j]=rho[i*ny+j];
	  s1_nm1[i*ny+j]=s1[i*ny+j];
	  s2_nm1[i*ny+j]=s2[i*ny+j];
	  rhoE_nm1[i*ny+j]=rhoE[i*ny+j];
	  
      }}

    for(int i=0;i<nx;i++){
      for(int j=0;j<ny;j++){
        rho[i*ny+j]=rho_np1[i*ny+j];
	s1[i*ny+j]=s1_np1[i*ny+j];
	s2[i*ny+j]=s2_np1[i*ny+j];
	rhoE[i*ny+j]=rhoE_np1[i*ny+j];
	
      }}
    EOS(p,cs,E,rho,ux,uy,nx,ny,gamma,EOStype);
    //that concludes all updates
    //readout
    double tstopper = fmod(totalT,cadence);
    if(abs(tstopper) < dt){
      //if(k%modo==0){
      //if(k%10==0){
      cout<<"reading out to file..."<<endl;
      cout<<"timestep number: "<<k<<"/"<<nsteps<<endl;
      cout<<"simulation time: "<<totalT<<endl;
      t = clock();
      cout<<"physical time: "<<((float)t)/CLOCKS_PER_SEC<<endl;
      for(int i=0;i<nx*ny;i++){
	
	outputFile <<rho[i]<<",";
	outputFile2 <<ux[i]<<",";
	outputFile3 <<uy[i]<<",";
	outputFile4 <<E[i]<<",";
	
      }
      outputFile<<endl;
      outputFile2<<endl;
      outputFile3<<endl;
      outputFile4<<endl;
    }
  }
  return 0;
}
