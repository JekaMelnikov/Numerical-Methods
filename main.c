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

/// calculate rezidual in norm C and L2 and print report
void rezidual (spec_mesh mesh, Vector* u, double cur_t)
{
  int cell, up, down, right, left;
  double x, y, x1, y1, x2, y2;
  double rez_rho_c, rez_v1_c, rez_v2_c, rez_rho_l2, rez_v1_l2, rez_v2_l2, cur_rho_c, cur_v1_c, cur_v2_c;
  double rho_1, rho_2, v1_1, v1_2, v2_1, v2_2;
  
  rez_rho_l2 = rez_v1_l2 = rez_v2_l2 = rez_rho_c = rez_v1_c = rez_v2_c = 0.;
  
  for (cell = 0; cell < size (mesh); cell++)
    {
      x = get_x (mesh, cell);
      y = get_y (mesh, cell);
      
      /// calculate rezidual in norm C
      cur_rho_c = fabs (V_GetCmp (u, 3 * cell + 0) - debug_rho (cur_t, x, y));
      cur_v1_c  = fabs (V_GetCmp (u, 3 * cell + 1) - debug_v1  (cur_t, x, y));
      cur_v2_c  = fabs (V_GetCmp (u, 3 * cell + 2) - debug_v2  (cur_t, x, y));
      
      if (cur_rho_c - rez_rho_c > MINIMAL_FOR_COMPARE)
        rez_rho_c = cur_rho_c;
      if (cur_v1_c - rez_v1_c > MINIMAL_FOR_COMPARE)
        rez_v1_c = cur_v1_c;
      if (cur_v2_c - rez_v2_c > MINIMAL_FOR_COMPARE)
        rez_v2_c = cur_v2_c;
      
      /// calculate rezidual in norm L2
      up    = get_up    (mesh, cell);
      down  = get_down  (mesh, cell);
      right = get_right (mesh, cell);
      left  = get_left  (mesh, cell);
      
      if (up != -1 && left != -1)
        {
          x1 = get_x (mesh, up);
          y1 = get_y (mesh, up);
          
          x2 = get_x (mesh, left);
          y2 = get_y (mesh, left);
          
          rho_1 = V_GetCmp (u, 3 * up + 0) - debug_rho (cur_t, x1, y1);
          v1_1  = V_GetCmp (u, 3 * up + 1) - debug_v1  (cur_t, x1, y1);
          v2_1  = V_GetCmp (u, 3 * up + 2) - debug_v2  (cur_t, x1, y1);
          
          rho_1 *= rho_1;
          v1_1  *= v1_1;
          v2_1  *= v2_1;
          
          rho_2 = V_GetCmp (u, 3 * left + 0) - debug_rho (cur_t, x2, y2);
          v1_2  = V_GetCmp (u, 3 * left + 1) - debug_v1  (cur_t, x2, y2);
          v2_2  = V_GetCmp (u, 3 * left + 2) - debug_v2  (cur_t, x2, y2);
          
          rho_2 *= rho_2;
          v1_2  *= v1_2;
          v2_2  *= v2_2;
          
          rez_rho_l2 += sqrt ((2 * get_h1 (mesh) * get_h2 (mesh) * (rho_1 + rho_2 + cur_rho_c * cur_rho_c)) / 3.);
          rez_v1_l2  += sqrt ((2 * get_h1 (mesh) * get_h2 (mesh) * (v1_1  + v1_2  + cur_v1_c  * cur_v1_c )) / 3.);
          rez_v2_l2  += sqrt ((2 * get_h1 (mesh) * get_h2 (mesh) * (v2_1  + v2_2  + cur_v2_c  * cur_v2_c )) / 3.);
        }
      if (down != -1 && right != -1)
        {
          x1 = get_x (mesh, down);
          y1 = get_y (mesh, down);
          
          x2 = get_x (mesh, right);
          y2 = get_y (mesh, right);
          
          rho_1 = V_GetCmp (u, 3 * down + 0) - debug_rho (cur_t, x1, y1);
          v1_1  = V_GetCmp (u, 3 * down + 1) - debug_v1  (cur_t, x1, y1);
          v2_1  = V_GetCmp (u, 3 * down + 2) - debug_v2  (cur_t, x1, y1);
          
          rho_1 *= rho_1;
          v1_1  *= v1_1;
          v2_1  *= v2_1;
          
          rho_2 = V_GetCmp (u, 3 * right + 0) - debug_rho (cur_t, x2, y2);
          v1_2  = V_GetCmp (u, 3 * right + 1) - debug_v1  (cur_t, x2, y2);
          v2_2  = V_GetCmp (u, 3 * right + 2) - debug_v2  (cur_t, x2, y2);
          
          rho_2 *= rho_2;
          v1_2  *= v1_2;
          v2_2  *= v2_2;
          
          rez_rho_l2 += sqrt ((2 * get_h1 (mesh) * get_h2 (mesh) * (rho_1 + rho_2 + cur_rho_c * cur_rho_c)) / 3.);
          rez_v1_l2  += sqrt ((2 * get_h1 (mesh) * get_h2 (mesh) * (v1_1  + v1_2  + cur_v1_c  * cur_v1_c )) / 3.);
          rez_v2_l2  += sqrt ((2 * get_h1 (mesh) * get_h2 (mesh) * (v2_1  + v2_2  + cur_v2_c  * cur_v2_c )) / 3.);
        }
    }
    
  /// printf report
  printf ("Rezidual in C:\n\tDensity\t\t%le\n\tVelosity by X\t%le\n\tVelosity by Y\t%le\n", rez_rho_c, rez_v1_c, rez_v2_c);
  printf ("Rezidual in L2:\n\tDensity\t\t%le\n\tVelosity by X\t%le\n\tVelosity by Y\t%le\n", rez_rho_l2, rez_v1_l2, rez_v2_l2);
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
