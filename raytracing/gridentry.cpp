#include "gridentry.h"
namespace GOAT
{
 namespace raytracing
 { 
  std::ostream & operator << (std::ostream &os, const stepEntry &se)
        {
          os << se.l << "\t" << se.matIndex << std::endl;
          return os;
        }
   std::ostream & operator << (std::ostream &os, const gridEntry &ge)
  { 
   for (stepEntry se : ge.step)
     os << se << "\t";
    os << ge.E;
   return os;  
  }
 }
}

