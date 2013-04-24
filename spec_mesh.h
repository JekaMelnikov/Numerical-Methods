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
    
int right (spec_mesh mesh, int ind);
int left  (spec_mesh mesh, int ind);
int up    (spec_mesh mesh, int ind);
int down  (spec_mesh mesh, int ind);
    
double x (spec_mesh mesh, int ind);
double y (spec_mesh mesh, int ind);
