// #include <unistd.h>
#include <fstream>
#include <stdio.h>
#include <string.h>
#include "mem.h"


namespace GOAT
{
   namespace raytracing
   {
MemInfo memstat ()
{
 MemInfo mem; 
 mem.VmData = 0;
 return mem; 
}

      std::ostream &operator << (std::ostream &os, MemInfo mem)
      {
      os << "Speicherbedarf :" << std::endl;
      os << "VmSize   :" << mem.VmSize << "kB" << std::endl;
      os << "VmLck    :" << mem.VmLck  << "kB" << std::endl; 
      os << "VmRSS    :" << mem.VmRSS  << "kB" << std::endl;
      os << "VmData   :" << mem.VmData << "kB" << std::endl; 
      os << "VmStk    :" << mem.VmStk  << "kB" << std::endl; 
      os << "VmExe    :" << mem.VmExe  << "kB" << std::endl;
      os << "VmLib    :" << mem.VmLib  << "kB" << std::endl; 
      return os;
      }

      SysMemInfo sysmem ()
      {
      SysMemInfo SysMem;
      char dummy[255];
      std::ifstream is;
      is.open ("/proc/meminfo");
      is >> dummy;
      if (strcmp(dummy,"total:")==NULL)
      {
      for (int i=1; i<=6; i++)
      is >> dummy;
      is >> SysMem.total >> SysMem.used >> SysMem.free;
      is >> SysMem.shared >> SysMem.buffers >> SysMem.cached;
      is >> dummy >> SysMem.swapTotal >> SysMem.swapUsed >> SysMem.swapFree;
      SysMem.total/=1024;
      SysMem.used/=1024;
      SysMem.free/=1024;
      SysMem.shared/=1024;
      SysMem.buffers/=1024;
      SysMem.cached/=1024;
      SysMem.swapTotal/=1024;
      SysMem.swapUsed/=1024;
      SysMem.swapFree/=1024;
      }
      else 
      {
      is >> SysMem.total >> dummy;
      is >> dummy >> SysMem.free >> dummy;
      is >> dummy >> SysMem.buffers >> dummy;
      is >> dummy >> SysMem.cached >> dummy;
      for (int i=1; i<=7; i++)
         is >> dummy >> dummy >> dummy;
      is >> dummy >> SysMem.swapTotal >> dummy;
      is >> dummy >> SysMem.swapFree >> dummy;
      SysMem.swapUsed=SysMem.swapTotal-SysMem.swapFree;
      }
      is.close();
      return SysMem;
      }

      std::ostream& operator << (std::ostream &os, SysMemInfo smi)
      {
      os << "Speicher: " << std::endl;
      os << "Total      : " << smi.total     << "        Used: " << smi.used     << "        Free : " << smi.free << std::endl;
      os << "Swap(total): " << smi.swapTotal << "  Swap(Used): " << smi.swapUsed << "   swap(Free): " << smi.swapFree << std::endl;
      os << "Buffers    : " << smi.buffers << "Cached     : " << smi.cached    << std::endl; 
      return os;
      }
   }
}