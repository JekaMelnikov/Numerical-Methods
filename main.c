#include <stdio.h>

#include <laspack/vector.h>
#include <laspack/qmatrix.h>
#include <laspack/itersolv.h>
#include <laspack/operats.h>

#include "spec_mesh.h"
#include "defs.h"

#define DEBUG_PRINT 300
#define FULL_OUTPUT 1
#define SOLVER      2

void print_vector (Vector *v);

void init (Vector *u, spec_mesh mesh);

/// build matrix A and left part vector b for system Au = b
void create_system (spec_mesh mesh, Vector *u, QMatrix *A, Vector *b, double cur_t, double tau);

/// calculate rezidual in norm C and L2 and print report
void rezidual (spec_mesh mesh, Vector* u, double cur_t);

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
  
  int M1, M2;
  double T, tau;
  fscanf (file, "%d%d%lf%lf", &M1, &M2, &T, &tau);
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
      if (FULL_OUTPUT)
        {
          printf ("\nT = %.3lf\n", cur_T);
          /// calculate and print rezidual
          rezidual (mesh, &u, cur_T);
        }
      
      
      if (DEBUG_PRINT == 100)
        {
          printf ("b = ");
          print_vector (&b);
          printf ("\n");
        }
      
      cur_T += tau;
      
      /// fill matrix A and vector b
      create_system (mesh, &u, &A, &b, cur_T, tau);

      if (DEBUG_PRINT == 100)
        {
          printf ("b = ");
          print_vector (&b);
          printf ("\n");
          
          printf ("u = ");
          print_vector (&u);
          printf ("\n");
        }
      
      /// solve Au = b
      SetRTCAccuracy (1e-16);
      if (SOLVER == 1)
        BiCGSTABIter (&A, &u, &b, MAX_ITER, NULL, 1.2);
      if (SOLVER == 2)
        BiCGIter (&A, &u, &b, MAX_ITER, NULL, 1.2);
      
      if (DEBUG_PRINT == 300)
        {
          Vector *b1;
          b1 = Mul_QV (&A, &u);
          int it;
          double norm_c, norm_l1, cur;
          norm_c = norm_l1 = 0.;
          
          for (it = 1; it <= 3 * size (mesh); it++)
            {
              cur = fabs (V__GetCmp (b1, it) - V__GetCmp (&b, it));
              if (cur > norm_c)
                norm_c = cur;
              norm_l1 += cur;
            }
          printf ("err_c = %le\nerr_l1 = %le\n", norm_c, norm_l1);
        }      
      
      if (DEBUG_PRINT == 100)
        {
          printf ("u = ");
          print_vector (&u);
          printf ("\n");
        }
    }
  
  printf ("Finaly rezult\n");
  rezidual (mesh, &u, cur_T);
  
  Q_Destr (&A);
  V_Destr (&u);
  V_Destr (&b);
  return 0;
}

void print_vector (Vector *v)
{
  int it;
  printf ("(");
  for (it = 1; it < V_GetDim (v); it++)
    printf ("%.4lf; ", V__GetCmp (v, it));
  printf ("%.4lf)", V__GetCmp (v, it));
}

void init (Vector *u, spec_mesh mesh)
{
  int cell;
  for (cell = 0; cell < size (mesh); cell++)
    if (get_up (mesh, cell) == -1 || get_down (mesh, cell) == -1 || 
        get_right (mesh, cell) == -1 || get_left (mesh, cell) == -1)
      {
        V__SetCmp (u, 3 * cell + 1, init_rho (get_x (mesh, cell), get_y (mesh, cell))); /// density
        V__SetCmp (u, 3 * cell + 2, 0.0);                                               /// x velocity
        V__SetCmp (u, 3 * cell + 3, 0.0);                                               /// y velocity
        
        if (DEBUG_PRINT == 100)
          printf ("(%.3lf %.3lf) cell=%d\tu=%d\td=%d\tr=%d\tl=%d\n", get_x (mesh, cell), get_y (mesh, cell), cell, get_up (mesh, cell), get_down (mesh, cell), get_right (mesh, cell), get_left (mesh, cell));
      }
    else
      {
        V__SetCmp (u, 3 * cell + 1, init_rho (get_x (mesh, cell), get_y (mesh, cell))); /// density
        V__SetCmp (u, 3 * cell + 2, init_v1  (get_x (mesh, cell), get_y (mesh, cell))); /// x velocity
        V__SetCmp (u, 3 * cell + 3, init_v2  (get_x (mesh, cell), get_y (mesh, cell))); /// y velocity
        
        if (DEBUG_PRINT == 100)
          printf ("(%.3lf %.3lf) cell=%d\tu=%d\td=%d\tr=%d\tl=%d\n", get_x (mesh, cell), get_y (mesh, cell), cell, get_up (mesh, cell), get_down (mesh, cell), get_right (mesh, cell), get_left (mesh, cell));
      }
}

/// build matrix A and left part vector b for system Au = b
void create_system (spec_mesh mesh, Vector *u, QMatrix *A, Vector *b, double cur_t, double tau)
{
  int cell;
  int neigh1, neigh2, neigh3, neigh4, neigh13, neigh14, neigh23, neigh24;
  double coef;
  
  for (cell = 0; cell < size (mesh); cell++)
    if (get_up (mesh, cell) == -1) /// w3.2
      {
        neigh1 = get_down (mesh, cell);
        neigh2 = get_down (mesh, neigh1);
        
        /// w3.2
        /// fill matrix A
        Q_SetLen (A, 3 * cell + 1, 4);
        /// rho
        coef = get_h1 (mesh) * get_h2 (mesh) / tau + 0.5 * V__GetCmp (u, 3 * cell + 3) * get_h1 (mesh);
        Q__SetEntry (A, 3 * cell + 1, 0, 3 * cell + 1, coef);
        /// w
        coef = V__GetCmp (u, 3 * cell + 1) * get_h1 (mesh);
        Q__SetEntry (A, 3 * cell + 1, 1, 3 * cell + 3, coef);
        /// rho_m-1
        coef = 0.5 * V__GetCmp (u, 3 * neigh1 + 3) * get_h1 (mesh);
        Q__SetEntry (A, 3 * cell + 1, 2, 3 * neigh1 + 1, -coef);
        /// w_m-1
        coef = V__GetCmp (u, 3 * cell + 1) * get_h1 (mesh);
        Q__SetEntry (A, 3 * cell + 1, 3, 3 * neigh1 + 3, -coef);
        
        /// fill vector b
        coef = V__GetCmp (u, 3 * cell + 1) * get_h1 (mesh) * get_h2 (mesh) / tau;
        coef += 0.5 * V__GetCmp (u, 3 * cell + 1) * (V__GetCmp (u, 3 * cell + 3) - V__GetCmp (u, 3 * neigh1 + 3)) * get_h1 (mesh);
        coef += 0.25 * (V__GetCmp (u, 3 * neigh2 + 1) * V__GetCmp (u, 3 * neigh2 + 3) - 2 * V__GetCmp (u, 3 * neigh1 + 1) * V__GetCmp (u, 3 * neigh1 + 3) + V__GetCmp (u, 3 * cell + 1) * V__GetCmp (u, 3 * cell + 3)) * get_h1 (mesh);
        coef += 0.25 * V__GetCmp (u, 3 * cell + 1) * (V__GetCmp (u, 3 * neigh2 + 3) - 2 * V__GetCmp (u, 3 * neigh1 + 3) + V__GetCmp (u, 3 * cell + 3)) * get_h1 (mesh);
        V__SetCmp (b, 3 * cell + 1, coef);
        
        /// boundary conditions
        Q_SetLen (A, 3 * cell + 2, 1);
        Q__SetEntry (A, 3 * cell + 2, 0, 3 * cell + 2, 1.);
        V__SetCmp (b, 3 * cell + 2, 0.);
        
        Q_SetLen (A, 3 * cell + 3, 1);
        Q__SetEntry (A, 3 * cell + 3, 0, 3 * cell + 3, 1.);
        V__SetCmp (b, 3 * cell + 3, 0.);
      }
    else if (get_down (mesh, cell) == -1 || get_down (mesh, get_left (mesh, cell)) == -1 || get_down (mesh, get_right (mesh, cell)) == -1) /// w2.2
      {
        neigh1 = get_up (mesh, cell);
        neigh2 = get_up (mesh, neigh1);
        
        /// w2.2
        /// fill matrix A
        Q_SetLen (A, 3 * cell + 1, 4);
        /// rho_m+1
        coef = 0.5 * V__GetCmp (u, 3 * neigh1 + 3) * get_h1 (mesh);
        Q__SetEntry (A, 3 * cell + 1, 0, 3 * neigh1 + 1, coef);
        /// v2_m+1
        coef = V__GetCmp (u, 3 * cell + 1) * get_h1 (mesh);
        Q__SetEntry (A, 3 * cell + 1, 1, 3 * neigh1 + 3, coef);
        /// rho
        coef = get_h1 (mesh) * get_h2 (mesh) / tau - 0.5 * V__GetCmp (u, 3 * cell + 3) * get_h1 (mesh);
        Q__SetEntry (A, 3 * cell + 1, 2, 3 * cell + 1, coef);
        /// w
        coef = V__GetCmp (u, 3 * cell + 1) * get_h1 (mesh);
        Q__SetEntry (A, 3 * cell + 1, 3, 3 * cell + 3, -coef);
        
        /// fill vector b
        coef = V__GetCmp (u, 3 * cell + 1) * get_h1 (mesh) * get_h2 (mesh) / tau;
        coef += 0.5 * V__GetCmp (u, 3 * cell + 1) * (V__GetCmp (u, 3 * neigh1 + 3) - V__GetCmp (u, 3 * cell + 3)) * get_h1 (mesh);
        coef += 0.25 * (V__GetCmp (u, 3 * neigh2 + 1) * V__GetCmp (u, 3 * neigh2 + 3) - 2. * V__GetCmp (u, 3 * neigh1 + 1) * V__GetCmp (u, 3 * neigh1 + 3) + V__GetCmp (u, 3 * cell + 1) * V__GetCmp (u, 3 * cell + 3)) * get_h1 (mesh);
        coef += 0.25 * V__GetCmp (u, 3 * cell + 1) * (V__GetCmp (u, 3 * neigh2 + 3) - 2. * V__GetCmp (u, 3 * neigh1 + 3) + V__GetCmp (u, 3 * cell + 3)) * get_h1 (mesh);
        V__SetCmp (b, 3 * cell + 1, coef);
        
        /// boundary conditions
        Q_SetLen (A, 3 * cell + 2, 1);
        Q__SetEntry (A, 3 * cell + 2, 0, 3 * cell + 2, 1.);
        V__SetCmp (b, 3 * cell + 2, 0.);
        
        Q_SetLen (A, 3 * cell + 3, 1);
        Q__SetEntry (A, 3 * cell + 3, 0, 3 * cell + 3, 1.);
        V__SetCmp (b, 3 * cell + 3, 0.);
      }
    else if (get_left (mesh, cell) == -1) /// w2.1
      {
        neigh1 = get_right (mesh, cell);
        neigh2 = get_right (mesh, neigh1);
        
        /// w2.1
        /// fill matrix A
        Q_SetLen (A, 3 * cell + 1, 4);
        /// rho
        coef = get_h1 (mesh) * get_h2 (mesh) / tau - 0.5 * V__GetCmp (u, 3 * cell + 2) * get_h2 (mesh);
        Q__SetEntry (A, 3 * cell + 1, 0, 3 * cell + 1, coef);
        /// v1
        coef = V__GetCmp (u, 3 * cell + 1) * get_h2 (mesh);
        Q__SetEntry (A, 3 * cell + 1, 1, 3 * cell + 2, -coef);
        /// rho_n+1
        coef = 0.5 * V__GetCmp (u, 3 * neigh1 + 2) * get_h2 (mesh);
        Q__SetEntry (A, 3 * cell + 1, 2, 3 * neigh1 + 1, coef);
        /// v1_n+1
        coef = V__GetCmp (u, 3 * cell + 1) * get_h2 (mesh);
        Q__SetEntry (A, 3 * cell + 1, 3, 3 * neigh1 + 2, coef);
        
        /// fill vector b
        coef = get_h1 (mesh) * get_h1 (mesh) * V__GetCmp (u, 3 * cell + 1) / tau;
        coef += 0.5 * V__GetCmp (u, 3 * cell + 1) * (V__GetCmp (u, 3 * neigh1 + 2) - V__GetCmp (u, 3 * cell + 2)) * get_h2 (mesh);
        coef += 0.25 * (V__GetCmp (u, 3 * neigh2 + 1) * V__GetCmp (u, 3 * neigh2 + 2) - 2 * V__GetCmp (u, 3 * neigh1 + 1) * V__GetCmp (u, 3 * neigh1 + 2) + V__GetCmp (u, 3 * cell + 1) * V__GetCmp (u, 3 * cell + 2)) * get_h2 (mesh);
        coef += 0.25 * V__GetCmp (u, 3 * cell + 1) * (V__GetCmp (u, 3 * neigh2 + 2) - 2 * V__GetCmp (u, 3 * neigh1 + 2) + V__GetCmp (u, 3 * cell + 2)) * get_h2 (mesh);
        V__SetCmp (b, 3 * cell + 1, coef);
        
        /// boundary conditions
        Q_SetLen (A, 3 * cell + 2, 1);
        Q__SetEntry (A, 3 * cell + 2, 0, 3 * cell + 2, 1.);
        V__SetCmp (b, 3 * cell + 2, 0.);
        
        Q_SetLen (A, 3 * cell + 3, 1);
        Q__SetEntry (A, 3 * cell + 3, 0, 3 * cell + 3, 1.);
        V__SetCmp (b, 3 * cell + 3, 0.);
      }
    else if (get_right (mesh, cell) == -1) /// w3.1
      {
        neigh1 = get_left (mesh, cell);
        neigh2 = get_left (mesh, neigh1);
        
        /// w3.1
        /// fill matrix A
        Q_SetLen (A, 3 * cell + 1, 4);
        /// rho_n-1
        coef = 0.5 * V__GetCmp (u, 3 * neigh1 + 2) * get_h2 (mesh);
        Q__SetEntry (A, 3 * cell + 1, 0, 3 * neigh1 + 1, -coef);
        /// v1_n-1
        coef = V__GetCmp (u, 3 * cell + 1) * get_h2 (mesh);
        Q__SetEntry (A, 3 * cell + 1, 1, 3 * neigh1 + 2, -coef);
        /// rho
        coef = get_h1 (mesh) * get_h2 (mesh) / tau + 0.5 * V__GetCmp (u, 3 * cell + 2) * get_h2 (mesh);
        Q__SetEntry (A, 3 * cell + 1, 2, 3 * cell + 1, coef);
        /// v1
        coef = V__GetCmp (u, 3 * cell + 1) * get_h2 (mesh);
        Q__SetEntry (A, 3 * cell + 1, 3, 3 * cell + 2, coef);
        
        /// fill vector b
        coef = get_h1 (mesh) * get_h2 (mesh) * V__GetCmp (u, 3 * cell + 1) / tau;
        coef += 0.5 * V__GetCmp (u, 3 * cell + 1) * (V__GetCmp (u, 3 * cell + 2) - V__GetCmp (u, 3 * neigh1 + 2)) * get_h2 (mesh);
        coef += 0.25 * (V__GetCmp (u, 3 * neigh2 + 1) * V__GetCmp (u, 3 * neigh2 + 2) - 2 * V__GetCmp (u, 3 * neigh1 + 1) * V__GetCmp (u, 3 * neigh1 + 2) + V__GetCmp (u, 3 * cell + 1) * V__GetCmp (u, 3 * cell + 2)) * get_h2 (mesh);
        coef += 0.25 * V__GetCmp (u, 3 * cell + 1) * (V__GetCmp (u, 3 * neigh2 + 2) - 2 * V__GetCmp (u, 3 * neigh1 + 2) + V__GetCmp (u, 3 * cell + 2)) * get_h2 (mesh);
        V__SetCmp (b, 3 * cell + 1, coef);
        
        /// boundary conditions
        Q_SetLen (A, 3 * cell + 2, 1);
        Q__SetEntry (A, 3 * cell + 2, 0, 3 * cell + 2, 1.);
        V__SetCmp (b, 3 * cell + 2, 0.);
        
        Q_SetLen (A, 3 * cell + 3, 1);
        Q__SetEntry (A, 3 * cell + 3, 0, 3 * cell + 3, 1.);
        V__SetCmp (b, 3 * cell + 3, 0.);
      }
    else /// w1 w4 w5
      {
        neigh1 = get_up    (mesh, cell);
        neigh2 = get_down  (mesh, cell);
        neigh3 = get_left  (mesh, cell);
        neigh4 = get_right (mesh, cell);
        
        neigh13 = get_left  (mesh, neigh1);
        neigh14 = get_right (mesh, neigh1);
        neigh23 = get_left  (mesh, neigh2);
        neigh24 = get_right (mesh, neigh2);
        
        /// w4
        /// fill matrix A
        Q_SetLen (A, 3 * cell + 1, 7);
        /// v
        coef = V__GetCmp (u, 3 * cell + 1) * get_h1 (mesh) * get_h2 (mesh) / tau + (8. * VISC * get_h2 (mesh)) / (3. * get_h1 (mesh)) + (2. * VISC * get_h1 (mesh)) / (get_h2 (mesh));
        Q__SetEntry (A, 3 * cell + 1, 0, 3 * cell + 2, coef);
        /// v_n+1
        coef = (get_h2 (mesh) / 3.) * (0.5 * V__GetCmp (u, 3 * cell + 1) * V__GetCmp (u, 3 * cell + 2) + 0.5 * V__GetCmp (u, 3 * neigh4 + 1) * V__GetCmp (u, 3 * neigh4 + 2) - 4. * VISC / get_h1 (mesh));
        Q__SetEntry (A, 3 * cell + 1, 1, 3 * neigh4 + 2, coef);
        /// v_n-1
        coef = (get_h2 (mesh) / 3.) * (0.5 * V__GetCmp (u, 3 * cell + 1) * V__GetCmp (u, 3 * cell + 2) + 0.5 * V__GetCmp (u, 3 * neigh3 + 1) * V__GetCmp (u, 3 * neigh3 + 2) + 4. * VISC / get_h1 (mesh));
        Q__SetEntry (A, 3 * cell + 1, 2, 3 * neigh3 + 2, -coef);
        /// v_m+1
        coef = get_h1 (mesh) * (0.25 * V__GetCmp (u, 3 * cell + 1) * V__GetCmp (u, 3 * cell + 3) + 0.25 * V__GetCmp (u, 3 * neigh1 + 1) * V__GetCmp (u, 3 * neigh1 + 3) - VISC / get_h2 (mesh));
        Q__SetEntry (A, 3 * cell + 1, 3, 3 * neigh1 + 2, coef);
        /// v_m-1
        coef = get_h1 (mesh) * (0.25 * V__GetCmp (u, 3 * cell + 1) * V__GetCmp (u, 3 * cell + 3) + 0.25 * V__GetCmp (u, 3 * neigh2 + 1) * V__GetCmp (u, 3 * neigh2 + 3) + VISC / get_h2 (mesh));
        Q__SetEntry (A, 3 * cell + 1, 4, 3 * neigh2 + 2, -coef);
        /// p_n+1
        coef = 0.5 * PP * get_h2 (mesh);
        Q__SetEntry (A, 3 * cell + 1, 5, 3 * neigh4 + 1, coef);
        /// p_n-1
        coef = 0.5 * PP * get_h2 (mesh);
        Q__SetEntry (A, 3 * cell + 1, 6, 3 * neigh3 + 1, -coef);
        
        /// fill vector b
        coef = V__GetCmp (u, 3 * cell + 1) * V__GetCmp (u, 3 * cell + 2) * get_h1 (mesh) * get_h2 (mesh) / tau;
        coef += 0.25 * V__GetCmp (u, 3 * cell + 2) * V__GetCmp (u, 3 * cell + 3) * (V__GetCmp (u, 3 * neigh1 + 1) - V__GetCmp (u, 3 * neigh2 + 1)) * get_h1 (mesh);
        coef += 0.25 * V__GetCmp (u, 3 * cell + 1) * V__GetCmp (u, 3 * cell + 2) * (V__GetCmp (u, 3 * neigh1 + 3) - V__GetCmp (u, 3 * neigh2 + 3)) * get_h1 (mesh);
        coef += (0.5 * V__GetCmp (u, 3 * cell + 2) * V__GetCmp (u, 3 * cell + 2) * (V__GetCmp (u, 3 * neigh4 + 1) - V__GetCmp (u, 3 * neigh3 + 1)) * get_h2 (mesh)) / 3.;
        coef += 0.25 * VISC * (V__GetCmp (u, 3 * neigh14 + 3) - V__GetCmp (u, 3 * neigh13 + 3) - V__GetCmp (u, 3 * neigh24 + 3) + V__GetCmp (u, 3 * neigh23 + 3));
        coef += V__GetCmp (u, 3 * cell + 1) * calc_f1 (cur_t, get_x (mesh, cell), get_y (mesh, cell)) * get_h1 (mesh) * get_h2 (mesh);
        V__SetCmp (b, 3 * cell + 1, coef);
        
        /// w5
        /// fill matrix A
        Q_SetLen (A, 3 * cell + 2, 7);
        /// w
        coef = V__GetCmp (u, 3 * cell + 1) * get_h1 (mesh) * get_h2 (mesh) / tau + (8. * VISC * get_h1 (mesh)) / (3. * get_h2 (mesh)) + (2. * VISC * get_h2 (mesh)) / (get_h1 (mesh));
        Q__SetEntry (A, 3 * cell + 2, 0, 3 * cell + 3, coef);
        /// w_m+1
        coef = (get_h1 (mesh) / 3.) * (0.5 * V__GetCmp (u, 3 * cell + 1) * V__GetCmp (u, 3 * cell + 3) + 0.5 * V__GetCmp (u, 3 * neigh1 + 1) * V__GetCmp (u, 3 * neigh1 + 3) - 4. * VISC / get_h2 (mesh));
        Q__SetEntry (A, 3 * cell + 2, 1, 3 * neigh1 + 3, coef);
        /// w_m-1
        coef = (get_h1 (mesh) / 3.) * (0.5 * V__GetCmp (u, 3 * cell + 1) * V__GetCmp (u, 3 * cell + 3) + 0.5 * V__GetCmp (u, 3 * neigh2 + 1) * V__GetCmp (u, 3 * neigh2 + 3) + 4. * VISC / get_h2 (mesh));
        Q__SetEntry (A, 3 * cell + 2, 2, 3 * neigh2 + 3, -coef);
        /// w_n+1
        coef = get_h2 (mesh) * (0.25 * V__GetCmp (u, 3 * cell + 1) * V__GetCmp (u, 3 * cell + 2) + 0.25 * V__GetCmp (u, 3 * neigh4 + 1) * V__GetCmp (u, 3 * neigh4 + 2) - VISC / get_h1 (mesh));
        Q__SetEntry (A, 3 * cell + 2, 3, 3 * neigh4 + 3, coef);
        /// w_n-1
        coef = get_h2 (mesh) * (0.25 * V__GetCmp (u, 3 * cell + 1) * V__GetCmp (u, 3 * cell + 2) + 0.25 * V__GetCmp (u, 3 * neigh3 + 1) * V__GetCmp (u, 3 * neigh3 + 2) + VISC / get_h1 (mesh));
        Q__SetEntry (A, 3 * cell + 2, 4, 3 * neigh3 + 3, -coef);
        /// p_m+1
        coef = 0.5 * PP * get_h1 (mesh);
        Q__SetEntry (A, 3 * cell + 2, 5, 3 * neigh1 + 1, coef);
        /// p_m-1
        coef = 0.5 * PP * get_h1 (mesh);
        Q__SetEntry (A, 3 * cell + 2, 6, 3 * neigh2 + 1, -coef);
        
        /// fill vector b
        coef = V__GetCmp (u, 3 * cell + 1) * V__GetCmp (u, 3 * cell + 3) * get_h1 (mesh) * get_h2 (mesh) / tau;
        coef += 0.25 * V__GetCmp (u, 3 * cell + 2) * V__GetCmp (u, 3 * cell + 3) * (V__GetCmp (u, 3 * neigh4 + 1) - V__GetCmp (u, 3 * neigh3 + 1)) * get_h2 (mesh);
        coef += 0.25 * V__GetCmp (u, 3 * cell + 1) * V__GetCmp (u, 3 * cell + 3) * (V__GetCmp (u, 3 * neigh4 + 2) - V__GetCmp (u, 3 * neigh3 + 2)) * get_h2 (mesh);
        coef += (0.5 * V__GetCmp (u, 3 * cell + 3) * V__GetCmp (u, 3 * cell + 3) * (V__GetCmp (u, 3 * neigh1 + 1) - V__GetCmp (u, 3 * neigh2 + 1)) * get_h1 (mesh)) / 3.;
        coef += 0.25 * VISC * (V__GetCmp (u, 3 * neigh14 + 2) - V__GetCmp (u, 3 * neigh13 + 2) - V__GetCmp (u, 3 * neigh24 + 2) + V__GetCmp (u, 3 * neigh23 + 2));
        coef += V__GetCmp (u, 3 * cell + 1) * calc_f2 (cur_t, get_x (mesh, cell), get_y (mesh, cell)) * get_h1 (mesh) * get_h2 (mesh);
        V__SetCmp (b, 3 * cell + 2, coef);
        
        /// w1
        /// fill matrix A
        Q_SetLen (A, 3 * cell + 3, 9);
        /// p
        coef = get_h1 (mesh) * get_h2 (mesh) / tau;
        Q__SetEntry (A, 3 * cell + 3, 0, 3 * cell + 1, coef);
        /// p_n+1
        coef = 0.25 * get_h2 (mesh) * (V__GetCmp (u, 3 * cell + 2) + V__GetCmp (u, 3 * neigh4 + 2));
        Q__SetEntry (A, 3 * cell + 3, 1, 3 * neigh4 + 1, coef);
        /// p_n-1
        coef = 0.25 * get_h2 (mesh) * (V__GetCmp (u, 3 * cell + 2) + V__GetCmp (u, 3 * neigh3 + 2));
        Q__SetEntry (A, 3 * cell + 3, 2, 3 * neigh3 + 1, -coef);
        /// p_m+1
        coef = 0.25 * get_h1 (mesh) * (V__GetCmp (u, 3 * cell + 3) + V__GetCmp (u, 3 * neigh1 + 3));
        Q__SetEntry (A, 3 * cell + 3, 3, 3 * neigh1 + 1, coef);
        /// p_m-1
        coef = 0.25 * get_h1 (mesh) * (V__GetCmp (u, 3 * cell + 3) + V__GetCmp (u, 3 * neigh2 + 3));
        Q__SetEntry (A, 3 * cell + 3, 4, 3 * neigh2 + 1, -coef);
        /// v_n+1
        coef = 0.5 * get_h2 (mesh) * V__GetCmp (u, 3 * cell + 1);
        Q__SetEntry (A, 3 * cell + 3, 5, 3 * neigh4 + 2, coef);
        /// v_n-1
        coef = 0.5 * get_h2 (mesh) * V__GetCmp (u, 3 * cell + 1);
        Q__SetEntry (A, 3 * cell + 3, 6, 3 * neigh3 + 2, -coef);
        /// w_m+1
        coef = 0.5 * get_h1 (mesh) * V__GetCmp (u, 3 * cell + 1);
        Q__SetEntry (A, 3 * cell + 3, 7, 3 * neigh1 + 3, coef);
        /// w_m-1
        coef = 0.5 * get_h1 (mesh) * V__GetCmp (u, 3 * cell + 1);
        Q__SetEntry (A, 3 * cell + 3, 8, 3 * neigh2 + 3, -coef);
        
        /// fill vector b
        coef = V__GetCmp (u, 3 * cell + 1) * get_h1 (mesh) * get_h2 (mesh) / tau;
        coef += 0.25 * get_h2 (mesh) * V__GetCmp (u, 3 * cell + 1) * (V__GetCmp (u, 3 * neigh4 + 2) - V__GetCmp (u, 3 * neigh3 + 2));
        coef += 0.25 * get_h1 (mesh) * V__GetCmp (u, 3 * cell + 1) * (V__GetCmp (u, 3 * neigh1 + 3) - V__GetCmp (u, 3 * neigh2 + 3));
        V__SetCmp (b, 3 * cell + 3, coef);
      }
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
      cur_rho_c = fabs (V__GetCmp (u, 3 * cell + 1) - debug_rho (cur_t, x, y));
      cur_v1_c  = fabs (V__GetCmp (u, 3 * cell + 2) - debug_v1  (cur_t, x, y));
      cur_v2_c  = fabs (V__GetCmp (u, 3 * cell + 3) - debug_v2  (cur_t, x, y));
      
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
          
          rho_1 = V__GetCmp (u, 3 * up + 1) - debug_rho (cur_t, x1, y1);
          v1_1  = V__GetCmp (u, 3 * up + 2) - debug_v1  (cur_t, x1, y1);
          v2_1  = V__GetCmp (u, 3 * up + 3) - debug_v2  (cur_t, x1, y1);
          
          rho_1 *= rho_1;
          v1_1  *= v1_1;
          v2_1  *= v2_1;
          
          rho_2 = V__GetCmp (u, 3 * left + 1) - debug_rho (cur_t, x2, y2);
          v1_2  = V__GetCmp (u, 3 * left + 2) - debug_v1  (cur_t, x2, y2);
          v2_2  = V__GetCmp (u, 3 * left + 3) - debug_v2  (cur_t, x2, y2);
          
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
          
          rho_1 = V__GetCmp (u, 3 * down + 1) - debug_rho (cur_t, x1, y1);
          v1_1  = V__GetCmp (u, 3 * down + 2) - debug_v1  (cur_t, x1, y1);
          v2_1  = V__GetCmp (u, 3 * down + 3) - debug_v2  (cur_t, x1, y1);
          
          rho_1 *= rho_1;
          v1_1  *= v1_1;
          v2_1  *= v2_1;
          
          rho_2 = V__GetCmp (u, 3 * right + 1) - debug_rho (cur_t, x2, y2);
          v1_2  = V__GetCmp (u, 3 * right + 2) - debug_v1  (cur_t, x2, y2);
          v2_2  = V__GetCmp (u, 3 * right + 3) - debug_v2  (cur_t, x2, y2);
          
          rho_2 *= rho_2;
          v1_2  *= v1_2;
          v2_2  *= v2_2;
          
          rez_rho_l2 += sqrt ((get_h1 (mesh) * get_h2 (mesh) * (rho_1 + rho_2 + cur_rho_c * cur_rho_c)) / 6.);
          rez_v1_l2  += sqrt ((get_h1 (mesh) * get_h2 (mesh) * (v1_1  + v1_2  + cur_v1_c  * cur_v1_c )) / 6.);
          rez_v2_l2  += sqrt ((get_h1 (mesh) * get_h2 (mesh) * (v2_1  + v2_2  + cur_v2_c  * cur_v2_c )) / 6.);
        }
    }
    
  /// printf report
  printf ("Rezidual in C:\n\tDensity\t\t%le\n\tVelosity by X\t%le\n\tVelosity by Y\t%le\n", rez_rho_c, rez_v1_c, rez_v2_c);
  printf ("Rezidual in L2:\n\tDensity\t\t%le\n\tVelosity by X\t%le\n\tVelosity by Y\t%le\n", rez_rho_l2, rez_v1_l2, rez_v2_l2);
}
