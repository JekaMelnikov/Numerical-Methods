#include "spec_mesh.h"

void spec_mesh::print_properties ()
{
  printf ("Special mesh %d x %d\nh1 = %.4lf\nh2 = %.4lf\nTotal number of cells %d\n",
          m_M1, m_M2, (double)1 / m_M1, (double)1 / m_M2, size ());
}

/// return number of cell neighbor to the top or -1 if this neighbor don't exist
int spec_mesh::up (int ind)
{
  int row, col;
  row = ind / (m_M1 + 1);
  
  if (row > m_M2 - m_M2/4)
    {
      row = (ind - (m_M1 + 1) * (m_M2 + 1 - m_M2 / 4)) / (m_M1 + 1 - (m_M1 / 5 - 1));
      col = (ind - (m_M1 + 1) * (m_M2 + 1 - m_M2 / 4)) % (m_M1 + 1 - (m_M1 / 5 - 1));
      if (row == 0)
        {
          if (col > (m_M1 + 1 - (m_M1 / 5 - 1)) / 2)
            return (m_M1 + 1) * (m_M2 - m_M2 / 4) + col + (m_M1 / 5 - 1);
          else
            return (m_M1 + 1) * (m_M2 - m_M2 / 4) + col;
        }
      else
        return ind - (m_M1 + 1 - (m_M1 / 5 - 1));
    }
  else
    {
      if (row == 0)
        return -1;
      else
        return ind - (m_M1 + 1);
    }
}

/// return number of cell neighbor to the top or -1 if this neighbor don't exist
int spec_mesh::down (int ind)
{
  int row, col;
  row = ind / (m_M1 + 1);
  
  if (row > m_M2 - m_M2 / 4)
    {
      row = (ind - (m_M1 + 1) * (m_M2 + 1 - m_M2 / 4)) / (m_M1 + 1 - (m_M1 / 5 - 1));
      col = (ind - (m_M1 + 1) * (m_M2 + 1 - m_M2 / 4)) % (m_M1 + 1 - (m_M1 / 5 - 1));
      if (row == m_M2 / 4)
        return -1;
      else
        return ind + m_M1 + 1 - (m_M1 / 5 - 1);
    }
  else
    {
      col = ind % (m_M1 + 1);
      if (row == m_M2 - m_M2/4)
        {
          if (col <= (2 * m_M1) / 5)
            return (m_M1 + 1) * (m_M2 + 1 - m_M2/4) + col;
          else if (col >= (3 * m_M1) / 5)
            return (m_M1 + 1) * (m_M2 + 1 - m_M2/4) + col - (m_M1 / 5 - 1);
          else
            return -1; 
        }
      else
        return ind + (m_M1 + 1);
    }
}

/// return number of cell neighbor to the left or -1 if this neighbor don't exist
int spec_mesh::left (int ind)
{
  int row, col;
  row = ind / (m_M1 + 1);
  
  if (row > m_M2 - m_M2 / 4)
    {
      col = (ind - (m_M1 + 1) * (m_M2 - m_M2 / 4 + 1)) % (m_M1 + 1 - (m_M1 / 5 - 1));
      if (col == 0 || col == (3 * m_M2) / 5 - (m_M1 / 5 - 1))
        return -1;
      else
        return ind - 1;
    }
  else
    {
      col = ind % (m_M1 + 1);
      if (col == 0)
        return -1;
      else
        return ind - 1;
    }
}

/// return number of cell neighbor to the right or -1 if this neighbor don't exist
int spec_mesh::right (int ind)
{
  int row, col;
  row = ind / (m_M1 + 1);
  
  if (row > m_M2 - m_M2 / 4)
    {
      col = (ind - (m_M1 + 1) * (m_M2 - m_M2 / 4 + 1)) % (m_M1 + 1 - (m_M1 / 5 - 1));
      if (col == m_M1 + 1 - (m_M1 / 5 - 1) || col == (2 * m_M2) / 5)
        return -1;
      else
        return ind + 1;
    }
  else
    {
      col = ind % (m_M1 + 1);
      if (col == m_M1)
        return -1;
      else
        return ind + 1;
    }
}
