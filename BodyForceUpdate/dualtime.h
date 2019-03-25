#ifndef DUALTIME_H    // To make sure you don't declare the function more than once by including the header multiple times.
#define DUALTIME_H
void update_primitive(double *q_np1, double *q_n, double *q_nm1, double *qc, double *qd, \
                      double dt, int nx, int ny);

void dualtimestep_1(double * rho_np1, double *s1_np1, double *s2_np1,double * rhoE_np1, \
                    double * d1, double * d2, double * d3, double * d4, \
                    double * rho_n, double *s1_n, double *s2_n, double * rhoE_n, \
                    double * rho_nm1, double *s1_nm1, double *s2_nm1, double * rhoE_nm1, \
                    double * rhouxuy,\
                    double * rhoux2_plus_p, double * rhouy2_plus_p,\
                    double * rhouxE_plus_pux, double * rhouyE_plus_puy, \
                    double * ux, double *uy, double *E,               \
                    double * cs,                   double *rhoax, double *rhoay, \
                    int nx, int ny,  double dx, double dy, \
                    double dt, double alpha);
void dualtimestep_2(double * rho_np1, double *s1_np1, double *s2_np1,double * rhoE_np1, \
                    double * d1, double * d2, double * d3, double * d4, \
                    double * rho_n, double *s1_n, double *s2_n, double * rhoE_n, \
                    double * rho_nm1, double *s1_nm1, double *s2_nm1, double * rhoE_nm1, \
                    double * rhouxuy,\
                    double * rhoux2_plus_p, double * rhouy2_plus_p,\
                    double * rhouxE_plus_pux, double * rhouyE_plus_puy, \
                    double *rhoax, double *rhoay,  int nx, int ny,  double dx, double dy, \
                    double dt, double alpha);

#endif
