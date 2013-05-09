#include "spec_mesh.h"

void create (spec_mesh *mesh, int M1, int M2)
{
  mesh->M1 = M1 * 5;
  mesh->M2 = M2 * 4;
}

int size (spec_mesh mesh)
{
  return (mesh.M1 + 1) * (mesh.M2 + 1) - (mesh.M1 / 5 - 1) * (mesh.M2 / 4);
}

int boundary_size (spec_mesh mesh)
{
  return 2 * (mesh.M1 + mesh.M2) + mesh.M2 / 2;
}

void print_properties (spec_mesh mesh)
{
  printf ("Special mesh %d x %d\nh1 = %.4lf\nh2 = %.4lf\nTotal number of cells %d\n",
          mesh.M1, mesh.M2, (double)1 / mesh.M1, (double)1 / mesh.M2, size (mesh));
}

/// return number of cell neighbor to the top or -1 if this neighbor don't exist
int get_up (spec_mesh mesh, int ind)
{
  int row, col;
  row = ind / (mesh.M1 + 1);
  
  if (row > mesh.M2 - mesh.M2/4)
    {
      row = (ind - (mesh.M1 + 1) * (mesh.M2 + 1 - mesh.M2 / 4)) / (mesh.M1 + 1 - (mesh.M1 / 5 - 1));
      col = (ind - (mesh.M1 + 1) * (mesh.M2 + 1 - mesh.M2 / 4)) % (mesh.M1 + 1 - (mesh.M1 / 5 - 1));
      if (row == 0)
        {
          if (col > (2 * mesh.M1) / 5)
            return (mesh.M1 + 1) * (mesh.M2 - mesh.M2 / 4) + col + (mesh.M1 / 5 - 1);
          else
            return (mesh.M1 + 1) * (mesh.M2 - mesh.M2 / 4) + col;
        }
      else
        return ind - (mesh.M1 + 1 - (mesh.M1 / 5 - 1));
    }
  else
    {
      if (row == 0)
        return -1;
      else
        return ind - (mesh.M1 + 1);
    }
}

/// return number of cell neighbor to the top or -1 if this neighbor don't exist
int get_down (spec_mesh mesh, int ind)
{
  int row, col;
  row = ind / (mesh.M1 + 1);
  
  if (row > mesh.M2 - mesh.M2 / 4)
    {
      row = (ind - (mesh.M1 + 1) * (mesh.M2 + 1 - mesh.M2 / 4)) / (mesh.M1 + 1 - (mesh.M1 / 5 - 1));
      col = (ind - (mesh.M1 + 1) * (mesh.M2 + 1 - mesh.M2 / 4)) % (mesh.M1 + 1 - (mesh.M1 / 5 - 1));
      if (row == mesh.M2 / 4 - 1)
        return -1;
      else
        return ind + mesh.M1 + 1 - (mesh.M1 / 5 - 1);
    }
  else
    {
      col = ind % (mesh.M1 + 1);
      if (row == mesh.M2 - mesh.M2/4)
        {
          if (col <= (2 * mesh.M1) / 5)
            return (mesh.M1 + 1) * (mesh.M2 + 1 - mesh.M2/4) + col;
          else if (col >= (3 * mesh.M1) / 5)
            return (mesh.M1 + 1) * (mesh.M2 + 1 - mesh.M2/4) + col - (mesh.M1 / 5 - 1);
          else
            return -1; 
        }
      else
        return ind + (mesh.M1 + 1);
    }
}

/// return number of cell neighbor to the left or -1 if this neighbor don't exist
int get_left (spec_mesh mesh, int ind)
{
  int row, col;
  row = ind / (mesh.M1 + 1);
  
  if (row > mesh.M2 - mesh.M2 / 4)
    {
      col = (ind - (mesh.M1 + 1) * (mesh.M2 - mesh.M2 / 4 + 1)) % (mesh.M1 + 1 - (mesh.M1 / 5 - 1));
      if (col == 0 || col == (3 * mesh.M1) / 5 - (mesh.M1 / 5 - 1))
        return -1;
      else
        return ind - 1;
    }
  else
    {
      col = ind % (mesh.M1 + 1);
      if (col == 0)
        return -1;
      else
        return ind - 1;
    }
}

/// return number of cell neighbor to the right or -1 if this neighbor don't exist
int get_right (spec_mesh mesh, int ind)
{
  int row, col;
  row = ind / (mesh.M1 + 1);
  
  if (row > mesh.M2 - mesh.M2 / 4)
    {
      col = (ind - (mesh.M1 + 1) * (mesh.M2 - mesh.M2 / 4 + 1)) % (mesh.M1 + 1 - (mesh.M1 / 5 - 1));
      if (col == mesh.M1 - (mesh.M1 / 5 - 1) || col == 2 * (mesh.M1 / 5))
        return -1;
      else
        return ind + 1;
    }
  else
    {
      col = ind % (mesh.M1 + 1);
      if (col == mesh.M1)
        return -1;
      else
        return ind + 1;
    }
}

/// return x coordinate of cell
double get_x (spec_mesh mesh, int ind)
{
  int row, col;
  row = ind / (mesh.M1 + 1);
  
  if (row > mesh.M2 - mesh.M2 / 4)
    {
      col = (ind - (mesh.M1 + 1) * (mesh.M2 - mesh.M2 / 4 + 1)) % (mesh.M1 + 1 - (mesh.M1 / 5 - 1));
      if (col > (2 * mesh.M1) / 5)
        col += mesh.M1 / 5 - 1;
    }
  else
    col = ind % (mesh.M1 + 1);
  
  return (double)(1.0 / mesh.M1) * col;
}

/// return y coordinate of cell
double get_y (spec_mesh mesh, int ind)
{
  int row;
  row = ind / (mesh.M1 + 1);
  
  if (row > mesh.M2 - mesh.M2 / 4)
    row = (ind - (mesh.M1 + 1) * (mesh.M2 - mesh.M2 / 4 + 1)) / (mesh.M1 + 1 - (mesh.M1 / 5 - 1)) + mesh.M2 - mesh.M2 / 4 + 1;
  
  return (double)(1.0 / mesh.M2) * (mesh.M2 - row);
}
