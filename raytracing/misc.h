
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

#ifndef MISC_H
#define MISC_H
namespace GOAT
{
    namespace raytracing
    {
#ifndef c_light
#define c_light 299792458.0   // Vakuumlichtgeschwindigkeit in m/s
#endif

#define Z0 376.730313461    // Impedanz des Vakuums Z0=sqrt (mu_0/epsilon_0)=mu_0*c_light in Ohm 
#define mu0 4.0*M_PI*1E-7   // µ0 in N/A^2
#define eps0 8.854187817E-12 // Epsilon_0
#define Planck_h 6.62606896E-34 // Plancksches Wirkungsquantum in Js
#define Planck_hquer 1.054571628E-34 // h/2PI in Js 


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
