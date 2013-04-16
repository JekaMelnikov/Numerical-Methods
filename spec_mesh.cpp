#include "spec_mesh.h"

void spec_mesh::print_properties ()
{
  printf ("Special mesh %d x %d\nh1 = %.4lf\nh2 = %.4lf\nTotal number of cells %d\n",
          m_M1, m_M2, (double)1 / m_M1, (double)1 / m_M2, size ());
}
