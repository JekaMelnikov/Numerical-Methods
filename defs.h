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
  return exp (t) * (cos (2 * PI * x) * cos (2 * PI * y + PI / 4) + 1.5);
}

double debug_v1 (double t, double x, double y)
{
  return exp (x - y) * cos (2 * PI * t) * sin (25 * PI * x * x) * sin (32 * PI * y * y);
}

double debug_v2 (double t, double x, double y)
{
  return exp (y - x) * sin (2 * PI * t) * sin (50 * PI * x * x) * sin (16 * PI * y * y);
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

double calc_f1 (double t, double x1, double x2, double x3)
{
  return
  (-2*exp(t + x1 - x2)*PI*(cos(2*PI*x1)*
  cos(2*PI*x2 + PI/(double)4) + 1.5)*sin(2*PI*t)*
  sin(25*PI*SQR(x1))*sin(32*PI*SQR(x2)) -
  2*exp(t)*PI*cos(2*PI*x2 + PI/(double)4)*sin(2*PI*x1)*
  (PP - (1/(double)3)*exp(2*x1 - 2*x2)*SQR(cos(2*PI*t))*
  SQR(sin(25*PI*SQR(x1)))*SQR(sin(32*PI*SQR(x2)))) +
  (1/(double)3)*(2*exp(t + 2*x1 - 2*x2)*SQR(cos(2*PI*t))*
  (cos(2*PI*x1)*cos(2*PI*x2 + PI/(double)4) + 1.5)*
  SQR(sin(25*PI*SQR(x1)))*SQR(sin(32*PI*SQR(x2))) -
  2*exp(t + 2*x1 - 2*x2)*PI*SQR(cos(2*PI*t))*
  cos(2*PI*x2 + PI/(double)4)*sin(2*PI*x1)*
  SQR(sin(25*PI*SQR(x1)))*SQR(sin(32*PI*SQR(x2))) +
  100*exp(t + 2*x1 - 2*x2)*PI*x1*SQR(cos(2*PI*t))*
  cos(25*PI*SQR(x1))*(cos(2*PI*x1)*
  cos(2*PI*x2 + PI/(double)4) + 1.5)*sin(25*PI*SQR(x1))*
  SQR(sin(32*PI*SQR(x2))) + exp(t + x1 - x2)*
  cos(2*PI*t)*(cos(2*PI*x1)*
  cos(2*PI*x2 + PI/(double)4) + 1.5)*sin(25*PI*SQR(x1))*
  (50*exp(x1 - x2)*PI*x1*cos(2*PI*t)*
  cos(25*PI*SQR(x1))*sin(32*PI*SQR(x2)) +
  exp(x1 - x2)*cos(2*PI*t)*sin(25*PI*SQR(x1))*
  sin(32*PI*SQR(x2)))*sin(32*PI*SQR(x2))) +
  (1/(double)2)*(64*exp(t)*PI*x2*cos(2*PI*t)*cos(32*PI*SQR(x2))*
  (cos(2*PI*x1)*cos(2*PI*x2 + PI/(double)4) + 1.5)*
  sin(2*PI*t)*sin(25*PI*SQR(x1))*sin(50*PI*SQR(x1))*
  sin(16*PI*SQR(x2)) + exp(t - x1 + x2)*
  (cos(2*PI*x1)*cos(2*PI*x2 + PI/(double)4) + 1.5)*
  sin(2*PI*t)*sin(50*PI*SQR(x1))*
  (64*exp(x1 - x2)*PI*x2*cos(2*PI*t)*
  cos(32*PI*SQR(x2))*sin(25*PI*SQR(x1)) -
  exp(x1 - x2)*cos(2*PI*t)*sin(25*PI*SQR(x1))*
  sin(32*PI*SQR(x2)))*sin(16*PI*SQR(x2)) +
  32*exp(t)*PI*x2*cos(2*PI*t)*cos(16*PI*SQR(x2))*
  (cos(2*PI*x1)*cos(2*PI*x2 + PI/(double)4) + 1.5)*
  sin(2*PI*t)*sin(25*PI*SQR(x1))*sin(50*PI*SQR(x1))*
  sin(32*PI*SQR(x2)) - exp(t + x1 - x2)*cos(2*PI*t)*
  (cos(2*PI*x1)*cos(2*PI*x2 + PI/(double)4) + 1.5)*
  sin(25*PI*SQR(x1))*(32*exp(x2 - x1)*PI*x2*
  cos(16*PI*SQR(x2))*sin(2*PI*t)*
  sin(50*PI*SQR(x1)) + exp(x2 - x1)*sin(2*PI*t)*
  sin(16*PI*SQR(x2))*sin(50*PI*SQR(x1)))*
  sin(32*PI*SQR(x2))) -
  VISC*(-4096*exp(x1 - x2)*SQR(PI)*cos(2*PI*t)*
  sin(25*PI*SQR(x1))*sin(32*PI*SQR(x2))*SQR(x2) -
  128*exp(x1 - x2)*PI*cos(2*PI*t)*cos(32*PI*SQR(x2))*
  sin(25*PI*SQR(x1))*x2 + 64*exp(x1 - x2)*PI*
  cos(2*PI*t)*cos(32*PI*SQR(x2))*sin(25*PI*SQR(x1)) +
  (1/(double)3)*(3200*exp(x2 - x1)*SQR(PI)*x1*x2*
  cos(50*PI*SQR(x1))*cos(16*PI*SQR(x2))*
  sin(2*PI*t) - 32*exp(x2 - x1)*PI*x2*
  cos(16*PI*SQR(x2))*sin(50*PI*SQR(x1))*
  sin(2*PI*t) + 100*exp(x2 - x1)*PI*x1*
  cos(50*PI*SQR(x1))*sin(16*PI*SQR(x2))*
  sin(2*PI*t) - exp(x2 - x1)*sin(50*PI*SQR(x1))*
  sin(16*PI*SQR(x2))*sin(2*PI*t)) +
  exp(x1 - x2)*cos(2*PI*t)*sin(25*PI*SQR(x1))*
  sin(32*PI*SQR(x2)) + (4/(double)3)*
  (-2500*exp(x1 - x2)*SQR(PI)*cos(2*PI*t)*
  sin(25*PI*SQR(x1))*sin(32*PI*SQR(x2))*SQR(x1) +
  100*exp(x1 - x2)*PI*cos(2*PI*t)*
  cos(25*PI*SQR(x1))*sin(32*PI*SQR(x2))*x1 +
  50*exp(x1 - x2)*PI*cos(2*PI*t)*
  cos(25*PI*SQR(x1))*sin(32*PI*SQR(x2)) +
  exp(x1 - x2)*cos(2*PI*t)*sin(25*PI*SQR(x1))*
  sin(32*PI*SQR(x2)))))/(double)exp(t)/(double)
  (cos(2*PI*x1)*cos(2*PI*x2 + PI/(double)4) + 1.5)
  ;
}

double calc_f2 (double t, double x1, double x2, double x3)
{
  return
  (2*exp(t - x1 + x2)*PI*cos(2*PI*t)*
  (cos(2*PI*x1)*cos(2*PI*x2 + PI/(double)4) + 1.5)*
  sin(50*PI*SQR(x1))*sin(16*PI*SQR(x2)) +
  (1/(double)2)*(100*exp(t)*PI*x1*cos(2*PI*t)*cos(50*PI*SQR(x1))*
  (cos(2*PI*x1)*cos(2*PI*x2 + PI/(double)4) + 1.5)*
  sin(2*PI*t)*sin(25*PI*SQR(x1))*sin(16*PI*SQR(x2))*
  sin(32*PI*SQR(x2)) + 50*exp(t)*PI*x1*cos(2*PI*t)*
  cos(25*PI*SQR(x1))*(cos(2*PI*x1)*
  cos(2*PI*x2 + PI/(double)4) + 1.5)*sin(2*PI*t)*
  sin(50*PI*SQR(x1))*sin(16*PI*SQR(x2))*
  sin(32*PI*SQR(x2)) + exp(t + x1 - x2)*cos(2*PI*t)*
  (cos(2*PI*x1)*cos(2*PI*x2 + PI/(double)4) + 1.5)*
  sin(25*PI*SQR(x1))*(100*exp(x2 - x1)*PI*x1*
  cos(50*PI*SQR(x1))*sin(2*PI*t)*
  sin(16*PI*SQR(x2)) - exp(x2 - x1)*sin(2*PI*t)*
  sin(50*PI*SQR(x1))*sin(16*PI*SQR(x2)))*
  sin(32*PI*SQR(x2)) - exp(t - x1 + x2)*
  (cos(2*PI*x1)*cos(2*PI*x2 + PI/(double)4) + 1.5)*
  sin(2*PI*t)*sin(50*PI*SQR(x1))*sin(16*PI*SQR(x2))*
  (64*exp(x1 - x2)*PI*x2*cos(2*PI*t)*
  cos(32*PI*SQR(x2))*sin(25*PI*SQR(x1)) -
  exp(x1 - x2)*cos(2*PI*t)*sin(25*PI*SQR(x1))*
  sin(32*PI*SQR(x2)))) -
  VISC*(-10000*exp(x2 - x1)*SQR(PI)*sin(2*PI*t)*
  sin(50*PI*SQR(x1))*sin(16*PI*SQR(x2))*SQR(x1) -
  200*exp(x2 - x1)*PI*cos(50*PI*SQR(x1))*sin(2*PI*t)*
  sin(16*PI*SQR(x2))*x1 + 100*exp(x2 - x1)*PI*
  cos(50*PI*SQR(x1))*sin(2*PI*t)*sin(16*PI*SQR(x2)) +
  exp(x2 - x1)*sin(2*PI*t)*sin(50*PI*SQR(x1))*
  sin(16*PI*SQR(x2)) +
  (4/(double)3)*(-1024*exp(x2 - x1)*SQR(PI)*sin(2*PI*t)*
  sin(50*PI*SQR(x1))*sin(16*PI*SQR(x2))*SQR(x2) +
  64*exp(x2 - x1)*PI*cos(16*PI*SQR(x2))*sin(2*PI*t)*
  sin(50*PI*SQR(x1))*x2 + 32*exp(x2 - x1)*PI*
  cos(16*PI*SQR(x2))*sin(2*PI*t)*
  sin(50*PI*SQR(x1)) + exp(x2 - x1)*sin(2*PI*t)*
  sin(50*PI*SQR(x1))*sin(16*PI*SQR(x2))) +
  (1/(double)3)*(3200*exp(x1 - x2)*SQR(PI)*x1*x2*cos(2*PI*t)*
  cos(25*PI*SQR(x1))*cos(32*PI*SQR(x2)) +
  64*exp(x1 - x2)*PI*x2*cos(2*PI*t)*
  sin(25*PI*SQR(x1))*cos(32*PI*SQR(x2)) -
  50*exp(x1 - x2)*PI*x1*cos(2*PI*t)*
  cos(25*PI*SQR(x1))*sin(32*PI*SQR(x2)) -
  exp(x1 - x2)*cos(2*PI*t)*sin(25*PI*SQR(x1))*
  sin(32*PI*SQR(x2)))) - 2*exp(t)*PI*cos(2*PI*x1)*
  (PP - (1/(double)3)*exp(2*x1 - 2*x2)*SQR(cos(2*PI*t))*
  SQR(sin(25*PI*SQR(x1)))*SQR(sin(32*PI*SQR(x2))))*
  sin(2*PI*x2 + PI/(double)4) +
  (1/(double)3)*(2*exp(t - 2*x1 + 2*x2)*
  (cos(2*PI*x1)*cos(2*PI*x2 + PI/(double)4) + 1.5)*
  SQR(sin(2*PI*t))*SQR(sin(16*PI*SQR(x2)))*
  SQR(sin(50*PI*SQR(x1))) + 64*exp(t - 2*x1 + 2*x2)*PI*
  x2*cos(16*PI*SQR(x2))*(cos(2*PI*x1)*
  cos(2*PI*x2 + PI/(double)4) + 1.5)*SQR(sin(2*PI*t))*
  sin(16*PI*SQR(x2))*SQR(sin(50*PI*SQR(x1))) -
  2*exp(t - 2*x1 + 2*x2)*PI*cos(2*PI*x1)*
  SQR(sin(2*PI*t))*SQR(sin(16*PI*SQR(x2)))*
  sin(2*PI*x2 + PI/(double)4)*SQR(sin(50*PI*SQR(x1))) +
  exp(t - x1 + x2)*(cos(2*PI*x1)*
  cos(2*PI*x2 + PI/(double)4) + 1.5)*sin(2*PI*t)*
  sin(16*PI*SQR(x2))*(32*exp(x2 - x1)*PI*x2*
  cos(16*PI*SQR(x2))*sin(2*PI*t)*
  sin(50*PI*SQR(x1)) + exp(x2 - x1)*sin(2*PI*t)*
  sin(16*PI*SQR(x2))*sin(50*PI*SQR(x1)))*
  sin(50*PI*SQR(x1))))/(double)
  (exp(t)*(cos(2*PI*x1)*cos(2*PI*x2 + PI/(double)4) + 1.5))
  ;
}
