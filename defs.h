#include <math.h>

#define MINIMAL_FOR_COMPARE  1.e-12

#define VISC      0.01
#define PP        10
#define MAX_ITER  1000

#define SQR(A) ((A) * (A))
#ifndef PI
#define PI      3.141592653589793
#endif

double debug_rho (double t, double x, double y)
{
  return exp (t + x - y);
}

double debug_v1 (double t, double x, double y)
{
  return exp (-t) * x * (x - 0.4) * (x - 0.6) * (1 - x) * y * (1 - y) * (y - 0.25);
}

double debug_v2 (double t, double x, double y)
{
  return log (1 + t) * sin (PI * x) * sin (PI * y) * (x - 0.4) * (x - 0.6) * (y - 0.25);
}

double init_rho (double x, double y)
{
  return debug_rho (0, x, y);
}

double init_v1 (double x, double y)
{
  return debug_v1 (0, x, y);
}

double init_v2 (double x, double y)
{
  return debug_v2 (0, x, y);
}

double calc_f1 (double t, double x1, double x2)
{
  return
  exp(-t - x1 + x2)*(exp(t + x1 - x2)*(PP - SQR(((1 - x1))*SQR((x1 - 0.6))*SQR((x1 - 0.4))*SQR(x1)*SQR((1 - x2))*SQR((x2 - 0.25))*SQR(x2))/(double)(3*exp(2*t))) +
  (1/(double)3)*(exp(-t + x1 - x2)*SQR((1 - x1))*SQR((x1 - 0.6))*SQR((x1 - 0.4))*SQR(x1)*SQR((1 - x2))*SQR((x2 - 0.25))*SQR(x2) - 2*exp(-t + x1 - x2)*(1 - x1)*SQR((x1 - 0.6))*SQR((x1 - 0.4))*SQR(x1)*
  SQR((1 - x2))*SQR((x2 - 0.25))*SQR(x2) + 2*exp(-t + x1 - x2)*SQR((1 - x1))*(x1 - 0.6)*SQR((x1 - 0.4))*SQR(x1)*SQR((1 - x2))*SQR((x2 - 0.25))*SQR(x2) +
  2*exp(-t + x1 - x2)*SQR((1 - x1))*SQR((x1 - 0.6))*(x1 - 0.4)*SQR(x1)*SQR((1 - x2))*SQR((x2 - 0.25))*SQR(x2) + 2*exp(-t + x1 - x2)*SQR((1 - x1))*SQR((x1 - 0.6))*SQR((x1 - 0.4))*x1*
  SQR((1 - x2))*SQR((x2 - 0.25))*SQR(x2) + exp(x1 - x2)*(((1 - x1)*(x1 - 0.6)*(x1 - 0.4)*(1 - x2)*(x2 - 0.25)*x2)/(double)exp(t) + ((1 - x1)*(x1 - 0.6)*x1*(1 - x2)*(x2 - 0.25)*x2)/(double)
  exp(t) + ((1 - x1)*(x1 - 0.4)*x1*(1 - x2)*(x2 - 0.25)*x2)/(double)exp(t) - ((x1 - 0.6)*(x1 - 0.4)*x1*(1 - x2)*(x2 - 0.25)*x2)/(double)exp(t))*(1 - x1)*(x1 - 0.6)*(x1 - 0.4)*x1*
  (1 - x2)*(x2 - 0.25)*x2) - exp(x1 - x2)*(1 - x1)*(x1 - 0.6)*(x1 - 0.4)*x1*(1 - x2)*(x2 - 0.25)*x2 -
  VISC*((4/(double)3)*((2*(1 - x1)*(x1 - 0.6)*(1 - x2)*(x2 - 0.25)*x2)/(double)exp(t) + (2*(1 - x1)*(x1 - 0.4)*(1 - x2)*(x2 - 0.25)*x2)/(double)exp(t) -
  (2*(x1 - 0.6)*(x1 - 0.4)*(1 - x2)*(x2 - 0.25)*x2)/(double)exp(t) + (2*(1 - x1)*x1*(1 - x2)*(x2 - 0.25)*x2)/(double)exp(t) - (2*(x1 - 0.6)*x1*(1 - x2)*(x2 - 0.25)*x2)/(double)exp(t) -
  (2*(x1 - 0.4)*x1*(1 - x2)*(x2 - 0.25)*x2)/(double)exp(t)) + (2*(1 - x1)*(x1 - 0.6)*(x1 - 0.4)*x1*(1 - x2))/(double)exp(t) - (2*(1 - x1)*(x1 - 0.6)*(x1 - 0.4)*x1*(x2 - 0.25))/(double)exp(t) -
  (2*(1 - x1)*(x1 - 0.6)*(x1 - 0.4)*x1*x2)/(double)exp(t) + (1/(double)3)*SQR((PI)*(x1 - 0.6)*(x1 - 0.4)*(x2 - 0.25)*cos(PI*x1)*cos(PI*x2)*log(t + 1) +
  PI*(x1 - 0.6)*(x2 - 0.25)*cos(PI*x2)*sin(PI*x1)*log(t + 1) + PI*(x1 - 0.4)*(x2 - 0.25)*cos(PI*x2)*sin(PI*x1)*log(t + 1) +
  PI*(x1 - 0.6)*(x1 - 0.4)*cos(PI*x1)*sin(PI*x2)*log(t + 1) + (x1 - 0.6)*sin(PI*x1)*sin(PI*x2)*log(t + 1) + (x1 - 0.4)*sin(PI*x1)*sin(PI*x2)*log(t + 1))) +
  (1/(double)2)*(exp(x1 - x2)*PI*(1 - x1)*SQR((x1 - 0.6))*x1*(1 - x2)*SQR((x2 - 0.25))*x2*cos(PI*x2)*log(t + 1)*sin(PI*x1)*SQR((x1 - 0.4)) +
  exp(x1 - x2)*(1 - x1)*SQR((x1 - 0.6))*x1*(1 - x2)*SQR((x2 - 0.25))*log(t + 1)*sin(PI*x1)*sin(PI*x2)*SQR((x1 - 0.4)) -
  exp(x1 - x2)*(1 - x1)*SQR((x1 - 0.6))*x1*SQR((x2 - 0.25))*x2*log(t + 1)*sin(PI*x1)*sin(PI*x2)*SQR((x1 - 0.4)) + 2*exp(x1 - x2)*(1 - x1)*SQR((x1 - 0.6))*x1*(1 - x2)*
  (x2 - 0.25)*x2*log(t + 1)*sin(PI*x1)*sin(PI*x2)*SQR((x1 - 0.4)) + exp(t + x1 - x2)*(((1 - x1)*(x1 - 0.6)*(x1 - 0.4)*x1*(1 - x2)*(x2 - 0.25))/(double)exp(t) +
  ((1 - x1)*(x1 - 0.6)*(x1 - 0.4)*x1*(1 - x2)*x2)/(double)exp(t) - ((1 - x1)*(x1 - 0.6)*(x1 - 0.4)*x1*(x2 - 0.25)*x2)/(double)exp(t))*(x1 - 0.6)*(x2 - 0.25)*log(t + 1)*sin(PI*x1)*
  sin(PI*x2)*(x1 - 0.4) - exp(x1 - x2)*(1 - x1)*(x1 - 0.6)*x1*(1 - x2)*(x2 - 0.25)*x2*(PI*(x1 - 0.6)*(x1 - 0.4)*(x2 - 0.25)*cos(PI*x2)*log(t + 1)*sin(PI*x1) +
  (x1 - 0.6)*(x1 - 0.4)*log(t + 1)*sin(PI*x2)*sin(PI*x1))*(x1 - 0.4)))
  ;
}

double calc_f2 (double t, double x1, double x2)
{
  return
  exp(-t - x1 + x2)*((-exp(t + x1 - x2))*(PP - ((1/(double)3)*SQR((1 - x1))*SQR((-0.6 + x1))*SQR((-0.4 + x1))*SQR(x1)*SQR((1 - x2))*SQR((-0.25 + x2))*SQR(x2))/(double)exp(2*t)) +
  (1/(double)3)*(2*exp(t + x1 - x2)*PI*SQR((-0.6 + x1))*SQR((-0.4 + x1))*SQR((-0.25 + x2))*cos(PI*x2)*SQR(log(1 + t))*SQR(sin(PI*x1))*sin(PI*x2) +
  2*exp(t + x1 - x2)*SQR((-0.6 + x1))*SQR((-0.4 + x1))*(-0.25 + x2)*SQR(log(1 + t))*SQR(sin(PI*x1))*SQR(sin(PI*x2)) - exp(t + x1 - x2)*SQR((-0.6 + x1))*SQR((-0.4 + x1))*SQR((-0.25 + x2))*
  SQR(log(1 + t))*SQR(sin(PI*x1))*SQR(sin(PI*x2)) + exp(t + x1 - x2)*(-0.6 + x1)*(-0.4 + x1)*(-0.25 + x2)*log(1 + t)*sin(PI*x1)*sin(PI*x2)*
  (PI*(-0.6 + x1)*(-0.4 + x1)*(-0.25 + x2)*cos(PI*x2)*log(1 + t)*sin(PI*x1) + (-0.6 + x1)*(-0.4 + x1)*log(1 + t)*sin(PI*x1)*sin(PI*x2))) +
  (1/(double)2)*(exp(x1 - x2)*PI*(1 - x1)*SQR((-0.6 + x1))*SQR((-0.4 + x1))*x1*(1 - x2)*SQR((-0.25 + x2))*x2*cos(PI*x1)*log(1 + t)*sin(PI*x2) +
  exp(x1 - x2)*(1 - x1)*SQR((-0.6 + x1))*SQR((-0.4 + x1))*(1 - x2)*SQR((-0.25 + x2))*x2*log(1 + t)*sin(PI*x1)*sin(PI*x2) +
  2*exp(x1 - x2)*(1 - x1)*SQR((-0.6 + x1))*(-0.4 + x1)*x1*(1 - x2)*SQR((-0.25 + x2))*x2*log(1 + t)*sin(PI*x1)*sin(PI*x2) +
  2*exp(x1 - x2)*(1 - x1)*(-0.6 + x1)*SQR((-0.4 + x1))*x1*(1 - x2)*SQR((-0.25 + x2))*x2*log(1 + t)*sin(PI*x1)*sin(PI*x2) -
  exp(x1 - x2)*SQR((-0.6 + x1))*SQR((-0.4 + x1))*x1*(1 - x2)*SQR((-0.25 + x2))*x2*log(1 + t)*sin(PI*x1)*sin(PI*x2) - exp(t + x1 - x2)*(-0.6 + x1)*(-0.4 + x1)*(-0.25 + x2)*
  (((1 - x1)*(-0.6 + x1)*(-0.4 + x1)*x1*(1 - x2)*(-0.25 + x2))/(double)exp(t) + ((1 - x1)*(-0.6 + x1)*(-0.4 + x1)*x1*(1 - x2)*x2)/(double)exp(t) -
  ((1 - x1)*(-0.6 + x1)*(-0.4 + x1)*x1*(-0.25 + x2)*x2)/(double)exp(t))*log(1 + t)*sin(PI*x1)*sin(PI*x2) + exp(x1 - x2)*(1 - x1)*(-0.6 + x1)*(-0.4 + x1)*x1*(1 - x2)*
  (-0.25 + x2)*x2*(PI*(-0.6 + x1)*(-0.4 + x1)*(-0.25 + x2)*cos(PI*x1)*log(1 + t)*sin(PI*x2) + (-0.6 + x1)*(-0.25 + x2)*log(1 + t)*sin(PI*x1)*sin(PI*x2) +
  (-0.4 + x1)*(-0.25 + x2)*log(1 + t)*sin(PI*x1)*sin(PI*x2))) -
  VISC*((1/(double)3)*(((1 - x1)*(-0.6 + x1)*(-0.4 + x1)*(1 - x2)*(-0.25 + x2))/(double)exp(t) + ((1 - x1)*(-0.6 + x1)*x1*(1 - x2)*(-0.25 + x2))/(double)exp(t) +
  ((1 - x1)*(-0.4 + x1)*x1*(1 - x2)*(-0.25 + x2))/(double)exp(t) - ((-0.6 + x1)*(-0.4 + x1)*x1*(1 - x2)*(-0.25 + x2))/(double)exp(t) +
  ((1 - x1)*(-0.6 + x1)*(-0.4 + x1)*(1 - x2)*x2)/(double)exp(t) + ((1 - x1)*(-0.6 + x1)*x1*(1 - x2)*x2)/(double)exp(t) + ((1 - x1)*(-0.4 + x1)*x1*(1 - x2)*x2)/(double)exp(t) -
  ((-0.6 + x1)*(-0.4 + x1)*x1*(1 - x2)*x2)/(double)exp(t) - ((1 - x1)*(-0.6 + x1)*(-0.4 + x1)*(-0.25 + x2)*x2)/(double)exp(t) - ((1 - x1)*(-0.6 + x1)*x1*(-0.25 + x2)*x2)/(double)exp(t) -
  ((1 - x1)*(-0.4 + x1)*x1*(-0.25 + x2)*x2)/(double)exp(t) + ((-0.6 + x1)*(-0.4 + x1)*x1*(-0.25 + x2)*x2)/(double)exp(t)) + 2*PI*(-0.6 + x1)*(-0.25 + x2)*cos(PI*x1)*log(1 + t)*
  sin(PI*x2) + 2*PI*(-0.4 + x1)*(-0.25 + x2)*cos(PI*x1)*log(1 + t)*sin(PI*x2) + 2*(-0.25 + x2)*log(1 + t)*sin(PI*x1)*sin(PI*x2) -
  SQR(PI)*(-0.6 + x1)*(-0.4 + x1)*(-0.25 + x2)*log(1 + t)*sin(PI*x1)*sin(PI*x2) + (4/(double)3)*(2*PI*(-0.6 + x1)*(-0.4 + x1)*cos(PI*x2)*log(1 + t)*sin(PI*x1) -
  SQR(PI)*(-0.6 + x1)*(-0.4 + x1)*(-0.25 + x2)*log(1 + t)*sin(PI*x1)*sin(PI*x2))) + exp(t + x1 - x2)*(-0.6 + x1)*(-0.4 + x1)*(-0.25 + x2)*sin(PI*x1)*sin(PI*x2)/(1 + t))
  ;
}
