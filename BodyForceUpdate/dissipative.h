#ifndef DISSIPATIVE_H    // To make sure you don't declare the function more than once by including the header multiple times.
#define DISSIPATIVE_H

double h1D1r(double *w_np1,double *w_n, double * w_nm1, \
             double *f1_n, double *f2_n,		\
             int i, int j,\
             int nx, int ny,\
             double dx, double dy, double dt);
double h2D2r(double *w_np1,double *w_n, double * w_nm1, \
             double *f1_n, double *f2_n,\
             int i, int j,\
             int nx, int ny,\
             double dx, double dy, double dt);
void dissipative_rho(double * d1,                               \
                     double *ux, double *uy,				\
                     double *rho_np1, double * rho_n, double *rho_nm1,  \
                     double *s1_np1, double * s1_n, double * s1_nm1,\
                     double *s2_np1, double * s2_n, double * s2_nm1,    \
                     double * rhouxuy,                                  \
                     double *rhoux2_plus_p, double *rhouy2_plus_p,\
                     double * rhoax, double * rhoay,double dx, double dy,double dt, \
                     int nx, int ny);
void dissipative_s1(double *d2,                         \
                    double *ux, double *uy,                             \
                    double *rho_np1, double * rho_n, double *rho_nm1,   \
                    double *s1_np1, double * s1_n, double * s1_nm1,\
                    double *s2_np1, double * s2_n, double * s2_nm1,    \
                    double * rhouxuy,                                  \
                    double *rhoux2_plus_p, double *rhouy2_plus_p,      \
                    double *cs,                                        \
                    double * rhoax, double * rhoay,double dx, double dy,double dt, \
                    int nx, int ny);
void dissipative_s2(double *d3,						\
                    double *ux, double *uy,                             \
                    double *rho_np1, double * rho_n, double *rho_nm1,   \
                    double *s1_np1, double * s1_n, double * s1_nm1,     \
                    double *s2_np1, double * s2_n, double * s2_nm1,    \
                    double * rhouxuy,                                  \
                    double *rhoux2_plus_p, double *rhouy2_plus_p,      \
                    double *cs,                                        \
                    double * rhoax, double * rhoay,double dx, double dy, double dt, \
                    int nx, int ny);
void dissipative_E(double *d4,                         \
                   double *ux, double *uy,                              \
                   double *E,\
                   double *rho_np1, double * rho_n, double *rho_nm1,    \
                   double *s1_np1, double * s1_n, double * s1_nm1,      \
                   double *s2_np1, double * s2_n, double * s2_nm1,      \
                   double *rhoE_np1, double * rhoE_n, double * rhoE_nm1,\
                   double * rhouxuy,                                    \
                   double *rhoux2_plus_p, double *rhouy2_plus_p,        \
                   double * cs,double *rhouxE_plus_pux, double *rhouyE_plus_puy, \
                   double dx, double dy,double dt,              \
                   int nx, int ny);

#endif
