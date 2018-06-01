#ifndef CONVECTIVE_H    // To make sure you don't declare the function more than once by including the header multiple times.
#define CONVECTIVE_H
void convective_rho(double * q1, double * sx, double * sy,      \
                    double dx, double dy,                       \
                    int nx, int ny );
void convective_s1(double * q2, double * rhoux2_plus_p, double * rhouxuy,\
                   double dx, double dy,                                \
                   int nx, int ny);
void convective_s2(double * q3,double * rhouxuy,double * rhouy2_plus_p,\
                   double dx, double dy,                               \
                   int nx, int ny);
void convective_E(double * q4,double * rhouxE_plus_pux, double *rhouyE_plus_puy, \
                  double dx, double dy,                                 \
                  int nx, int ny);

#endif
