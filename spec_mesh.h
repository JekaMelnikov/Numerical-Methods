#include <stdio.h>

class spec_mesh
{
  private:
    int m_M1;
    int m_M2;
  
  public:
    spec_mesh ()
      {
        m_M1 = 25;
        m_M2 = 20;
      }      
      
    spec_mesh (int M1, int M2)
      {
        m_M1 = M1 * 25;
        m_M2 = M2 * 20;
      }
    
    int size ()
      {
       return (m_M1 + 1) * (m_M2 + 1) - (m_M1 / 5 - 1) * (m_M2 / 4);
      }
    
    int boundary_size ()
      {
        return 2 * (m_M1 + m_M2) + m_M2 / 2;
      }
      
    void print_properties ();
};