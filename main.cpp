#include <stdio.h>

#include "spec_mesh.h"
#include "defs.h"

int main (int argc, char* argv[])
{
  FILE *file = fopen (argv[1], "r");
  if (!file)
    {
      printf ("Cant open initialization file.\n");
      return -1;
    }
  
  int M1, M2, T;
  double tau;
  fscanf (file, "%d%d%d%lf", &M1, &M2, &T, &tau);
  fclose (file);
  
  spec_mesh mesh (M1, M2);
  mesh.print_properties ();
  
  return 0;
}
