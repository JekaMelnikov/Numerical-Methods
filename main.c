#include <stdio.h>

#include <laspack/vector.h>
#include <laspack/qmatrix.h>
#include <laspack/itersolv.h>

#include "spec_mesh.h"
#include "defs.h"

#define DEBUG_PRINT 0

void init (Vector *u, spec_mesh mesh)
{
  int cell;
  for (cell = 0; cell < size (mesh); cell++)
    if (get_up (mesh, cell) == -1 || get_down (mesh, cell) == -1 || 
        get_right (mesh, cell) == -1 || get_left (mesh, cell) == -1)
      {
        V_SetCmp (u, 3 * cell + 0, init_rho (get_x (mesh, cell), get_y (mesh, cell))); /// density
        V_SetCmp (u, 3 * cell + 1, 0.0);                                               /// x velocity
        V_SetCmp (u, 3 * cell + 2, 0.0);                                               /// y velocity
        
        if (DEBUG_PRINT)
          printf ("(%.3lf %.3lf) cell=%d\tu=%d\td=%d\tr=%d\tl=%d\n", get_x (mesh, cell), get_y (mesh, cell), cell, get_up (mesh, cell), get_down (mesh, cell), get_right (mesh, cell), get_left (mesh, cell));
      }
    else
      {
        V_SetCmp (u, 3 * cell + 0, init_rho (get_x (mesh, cell), get_y (mesh, cell))); /// density
        V_SetCmp (u, 3 * cell + 1, init_v1  (get_x (mesh, cell), get_y (mesh, cell))); /// x velocity
        V_SetCmp (u, 3 * cell + 2, init_v2  (get_x (mesh, cell), get_y (mesh, cell))); /// y velocity
        
        if (DEBUG_PRINT)
          printf ("(%.3lf %.3lf) cell=%d\tu=%d\td=%d\tr=%d\tl=%d\n", get_x (mesh, cell), get_y (mesh, cell), cell, get_up (mesh, cell), get_down (mesh, cell), get_right (mesh, cell), get_left (mesh, cell));
      }
}

void create_system (spec_mesh mesh, Vector u, QMatrix *A, Vector *b)
{
  return;
}

/// calculate rezidual in norm C and print report
void rezidual (spec_mesh mesh, Vector* u, double cur_t)
{
  int cell;
  double x, y;
  double rez_rho_c, rez_v1_c, rez_v2_c, cur_rho_c, cur_v1_c, cur_v2_c;
  
  rez_rho_c = rez_v1_c = rez_v2_c = 0.;
  
  for (cell = 0; cell < size (mesh); cell++)
    {
      x = get_x (mesh, cell);
      y = get_y (mesh, cell);
      
      cur_rho_c = fabs (V_GetCmp (u, 3 * cell + 0) - debug_rho (cur_t, x, y));
      cur_v1_c  = fabs (V_GetCmp (u, 3 * cell + 1) - debug_v1  (cur_t, x, y));
      cur_v2_c  = fabs (V_GetCmp (u, 3 * cell + 2) - debug_v2  (cur_t, x, y));
      
      if (cur_rho_c - rez_rho_c > MINIMAL_FOR_COMPARE)
        rez_rho_c = cur_rho_c;
      if (cur_v1_c - rez_v1_c > MINIMAL_FOR_COMPARE)
        rez_v1_c = cur_v1_c;
      if (cur_v2_c - rez_v2_c > MINIMAL_FOR_COMPARE)
        rez_v2_c = cur_v2_c;
    }
    
  /// printf report
  printf ("Rezidual in C:\n\tDensity\t\t%le\n\tVelosity by X\t%le\n\tVelosity by Y\t%le\n", rez_rho_c, rez_v1_c, rez_v2_c);
}

int main (int argc, char* argv[])
{
  if (argc != 2)
    {
      printf ("Usage %s <input_file>\n", argv[0]);
      return -1;
    }
  
  FILE *file = fopen (argv[1], "r");
  if (!file)
    {
      printf ("Cant open input file %s\n", argv[1]);
      return -2;
    }
  
  int M1, M2, T;
  double tau;
  fscanf (file, "%d%d%d%lf", &M1, &M2, &T, &tau);
  fclose (file);
  
  spec_mesh mesh;
  create (&mesh, M1, M2);
  
  /// printf start report
  print_properties (mesh);
  printf ("tau = %.3lf\n", tau);
  
  Vector u;
  V_Constr (&u, "u", 3 * size (mesh), Normal, True);
  init (&u, mesh);
  
  QMatrix A;
  Q_Constr(&A, "A", 3 * size (mesh), False, Rowws, Normal, True);
  Vector b;
  V_Constr (&b, "b", 3 * size (mesh), Normal, True);
  
  double cur_T = 0.0;
  while (T - cur_T > MINIMAL_FOR_COMPARE) /// main loop
    {
      cur_T += tau;
      
      /// fill matrix A and vector b
      create_system (mesh, u, &A, &b);
      
      /// solve Au = b
      BiCGSTABIter (&A, &u, &b, MAX_ITER, NULL, 1.2);
      
      printf ("\nT = %.3lf\n", cur_T);
      /// calculate and print rezidual
      rezidual (mesh, &u, cur_T);
    }
  
  Q_Destr (&A);
  V_Destr (&u);
  V_Destr (&b);
  return 0;
}
