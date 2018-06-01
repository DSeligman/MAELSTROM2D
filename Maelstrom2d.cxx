//MAELSTROM 3D v1.0
//Vorticity Conserving Scheme for the Euler Equations in 3 Dimensions
//v1.0 08/03/2017
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
  ifstream myfile ("params.txt");
  double * inputs=NULL;
  inputs = new double[20];
  int nx, ny, nsteps,xbtype,ybtype,EOStype;
  
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
  nx  = inputs[0];
  ny  = inputs[1];
  extentx  = inputs[2];
  extenty  = inputs[3];
  nsteps = inputs[4];
  xbtype = inputs[10];
  ybtype = inputs[11];
  g = inputs[12];
  xstart = inputs[13];
  ystart = inputs[14];
  cadence = inputs[15];
  EOStype = inputs[16];
  double gamma = 5./3.;
  //int nsteps = 10;
  
  //instantiate the grid
  
  dx=extentx/((double)nx-(double)1);
  dy=extenty/((double)ny-(double)1);

  
  double * xg=NULL;
  xg = new double[nx];


  double * yg=NULL;
  yg = new double[ny];

  
  for(int i=0;i<nx;i++){
    xg[i]=xstart+dx*i;
  }


  for(int j=0;j<ny;j++){
    yg[j]=ystart+dy*j;
  }
  

  
  double dt=10.;
  //double cs = 1.; //sound speed
  //instantiate density and velocities and pressure
  double * rho=NULL;
  rho = new double[nx*ny];
  
  double * ux=NULL;
  ux = new double[nx*ny];
  
  double * uy=NULL;
  uy = new double[nx*ny];
  
  
  
  double * p=NULL;
  p = new double[nx*ny];
  
  //internal energy 
  double * E = NULL;
  E= new double[nx*ny];
  
  double * cs = NULL;
  cs=new double[nx*ny];
  double M = 2.*3.14;
  //initial conditions
  for(int i=0;i<nx;i++){
    for(int j=0;j<ny;j++){
            
      rho[i*ny+j]=1.0;
      ux[i*ny+j]=.1*sin(M *yg[j]) ;//.1*drand48();
      uy[i*ny+j]=0.;//.1*sin(M *yg[j]) ;
      E[i*ny+j]=1.0/(gamma*(gamma-1.));
      p[i*ny+j]=0.;
      cs[i*ny+j]=1.;
    }}
  //EOS call
  EOS(p,cs,E,rho,ux,uy,nx,ny,gamma,EOStype);


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

  ofstream outputFile("rho.txt"); // output filestream                           
  ofstream outputFile2("ux.txt");    //read out the initial velocity flow           
  ofstream outputFile3("uy.txt"); ////////                                          
  ofstream outputFile4("E.txt");

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
  //alpha values for JAMESON Scheme
  double alpha1=0.6;
  double alpha2=0.6;
  double alpha3=1.;

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
    BC(ux,nx,ny,xbtype,ybtype,1);
    BC(uy,nx,ny,xbtype,ybtype,1);

    BC(rho,nx,ny,xbtype,ybtype,0);
    BC(E,nx,ny,xbtype,ybtype,0);

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

    BC(ux,nx,ny,xbtype,ybtype,1);
    BC(uy,nx,ny,xbtype,ybtype,2);
    BC(E,nx,ny,xbtype,ybtype,3);
    BC(rho_np1,nx,ny,xbtype,ybtype,0);
    
    u_to_s(s1_np1,ux,rho,nx,ny);
    u_to_s(s2_np1,uy,rho,nx,ny);
    u_to_s(rhoE_np1,E,rho,nx,ny);
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
    BC(ux,nx,ny,xbtype,ybtype,1);
    BC(uy,nx,ny,xbtype,ybtype,2);
    BC(E,nx,ny,xbtype,ybtype,0);
    BC(rho_np1,nx,ny,xbtype,ybtype,0);

    u_to_s(s1_np1,ux,rho,nx,ny);
    u_to_s(s2_np1,uy,rho,nx,ny);
    u_to_s(rhoE_np1,E,rho,nx,ny);

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
    BC(ux,nx,ny,xbtype,ybtype,1);
    BC(uy,nx,ny,xbtype,ybtype,2);
    BC(E,nx,ny,xbtype,ybtype,0);
    BC(rho_np1,nx,ny,xbtype,ybtype,0);

    u_to_s(s1_np1,ux,rho,nx,ny);
    u_to_s(s2_np1,uy,rho,nx,ny);
    u_to_s(rhoE_np1,E,rho,nx,ny);

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
    if(k%modo==0){
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
