
/***************************************************************************
                          misc.h  -  description
                             -------------------
    begin                : Mit Mär 12 2003
    copyright            : (C) 2003 by Thomas Weigel
    email                : weigel@lat.ruhr-uni-bochum.de
 ***************************************************************************/

#include "surface.h"
#include "ellipsoid.h"
#include "vector.h"
#include "tubedray.h"
#include "box.h"
#include "cone.h"
#include "constants.h"
#include "cylinder.h"
#include "vortex.h"

#ifndef MISC_H
#define MISC_H
namespace GOAT
{
    namespace raytracing
    {



        void initInc(ObjectShape* E);
        void setR0(ObjectShape* E, double r0);
        void copyFormList(ObjectShape**& d, ObjectShape** s, int anz);
        void binWriteIncList(std::ofstream& os, ObjectShape** E, int anz);
        void binWriteInc(std::ofstream& os, ObjectShape* E);
        void binReadIncList(std::ifstream& is, ObjectShape**& E, int anz);
        void binReadInc(std::ifstream& is, ObjectShape*& E, bool isNew);
        void copyInc(ObjectShape*& d, ObjectShape* s);
        void deleteInc(ObjectShape* E);
        maths::Vector<double> force(maths::Vector<double> norm, tubedRay Se, tubedRay Sr, tubedRay St, double df);
        double gaussw(double z, double wvl, double w0);
        std::complex<double> gaussphase(maths::Vector<double> P, maths::Vector<double> F, maths::Vector<double> k, double w0, double k0);
        double  NA2w0(double lambda, double NA, std::complex<double> n);
         float readLE_float32(std::istream& is); 
        int readLE_int32(std::istream& is);
    }
}
#endif
