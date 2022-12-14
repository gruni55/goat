#ifndef __MEM__
#define __MEM__
#include <iostream>
namespace GOAT
{
    namespace raytracing
    {
        typedef struct 
        {
        long int VmSize;
        long int VmLck;
        long int VmRSS;
        long int VmData;
        long int VmStk;
        long int VmExe;
        long int VmLib;
        } MemInfo;

        typedef struct
        {
        long int total;
        long int used;
        long int free;
        long int shared;
        long int buffers;
        long int cached;
        long int swapTotal;
        long int swapUsed;
        long int swapFree;
        } SysMemInfo;

        MemInfo memstat ();
        SysMemInfo sysmem ();
        std::ostream& operator << (std::ostream &os, SysMemInfo smi);
    }
}
#endif
