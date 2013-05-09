#include <stdio.h>

typedef
struct _spec_mesh
{
  int M1;
  int M2;
} spec_mesh;

void create (spec_mesh *mesh, int M1, int M2);
      
int size (spec_mesh mesh);
      
int boundary_size (spec_mesh mesh);

void print_properties (spec_mesh mesh);
    
int get_right (spec_mesh mesh, int ind);
int get_left  (spec_mesh mesh, int ind);
int get_up    (spec_mesh mesh, int ind);
int get_down  (spec_mesh mesh, int ind);
    
double get_x (spec_mesh mesh, int ind);
double get_y (spec_mesh mesh, int ind);

double get_h1 (spec_mesh mesh);
double get_h2 (spec_mesh mesh);
