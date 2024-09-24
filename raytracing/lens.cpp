#include "lens.h"
namespace GOAT
{
    namespace raytracing
    {
        void binWrite(lensSide ls, std::ofstream& os)
        {
            ls.P.binWrite(os);
            os.write((char*)&ls.R, (char)sizeof(ls.R));
            os.write((char*)&ls.shift, (char)sizeof(ls.shift));
            os.write((char*)&ls.curvature, (char)sizeof(ls.curvature));
        }

        void binWrite(lensParms lp, std::ofstream& os)
        {
            binWrite(lp.left, os);
            binWrite(lp.right, os);
            os.write((char*)&lp.offset, (char)sizeof(lp.offset));
            os.write((char*)&lp.radius, (char)sizeof(lp.radius));
        }

        void binRead(lensSide ls, std::ifstream& is)
        {
            ls.P.binRead(is);
            is.read((char*)&ls.R, (char)sizeof(ls.R));
            is.read((char*)&ls.shift, (char)sizeof(ls.shift));
            is.read((char*)&ls.curvature, (char)sizeof(ls.curvature));
        }

        void binRead(lensParms lp, std::ifstream& is)
        {
            binRead(lp.left, is);
            binRead(lp.right, is);
            is.read((char*)&lp.offset, (char)sizeof(lp.offset));
            is.read((char*)&lp.radius, (char)sizeof(lp.radius));
        }
    }
}