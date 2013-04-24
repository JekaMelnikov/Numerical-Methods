#include <math.h>

#define MINIMAL_FOR_COMPARE  1.e-12

#define VISC      0.01
#define PP        10
#define MAX_ITER  1000

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
