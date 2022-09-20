#ifndef FRESNEL_CC
#define FRESNEL_CC

#include "fresnel.h"
#include "vector.h" 
#include <fstream>
#include <complex>
#define _USE_MATH_DEFINES
#include <math.h>
 
#ifndef double_complex
#define double_complex complex<double>
#endif
#include "goodies.h"

using namespace std;
namespace GOAT
{
    namespace maths
    { 
    complex<double> Fresnel_reflect(int pol, Vector<double> k, Vector<double> n, double_complex n1, double_complex n2)

        /*  Berechnet den (komplexen) Fresnel-Koeffizienten fuer Reflexion
            ( auf die Amplitude bezogen !!!!)
            S : Strahl
            n : Normale auf die reflektierende Flaeche
            n1 : Brechungsindex innen
            n2 : Brechungsindex aussen

        */

    {
        double alpha;
        double_complex w, n12, cosa, sin2a;
        complex<double> wurzel, Erg;
        alpha = acos(k * n);
        if (alpha > M_PI / 2.0) alpha = M_PI - alpha;
        n12 = n2 / n1;
        n12 *= n12;
        cosa = cos(alpha);
        sin2a = sin(alpha);
        sin2a *= sin2a;
        w = n12 - sin2a;
        wurzel = sqrt(w);
        if (pol == SENKRECHT) Erg = (cosa - wurzel) / (cosa + wurzel);
        else Erg = (n12 * cosa - wurzel) / (n12 * cosa + wurzel);
        return Erg;
    }


    complex<double> Fresnel_trans(int pol, Vector<double> k, Vector<double> n,
        double_complex n1, double_complex n2)

        //  Berechnet den (komplexen) Fresnel-Koeffizienten fuer Transmission 
        //  S : Strahl
        //  n : Normale auf die reflektierende Stelle
        //  n1 : Brechungsindex innen 
        //  n2 : Brechungsindex aussen 
    {
        double_complex n12;
        double alpha;
        complex<double> beta, Erg;
        alpha = acos(k * n);
        if (alpha > M_PI / 2.0) alpha = M_PI - alpha;
        n12 = n2 / n1;
        beta = asin((double_complex)sin(alpha) / n12);
        if (pol == SENKRECHT) Erg = 2.0 * sin(beta) * cos(alpha) / sin(alpha + beta);
        else Erg = 2.0 * sin(beta) * cos(alpha) / (sin(alpha + beta) * cos(alpha - beta));
        return Erg;
    }

    double_complex freflect(int pol, double alpha, double_complex n1, double_complex n2)
    {
        double_complex r;
        double_complex n12;
        double_complex beta;
        n12 = n2 / n1;
        beta = asin((double_complex)sin(alpha) / n12);
        if (pol == PARALLEL)
            r = tan(alpha - beta) / tan(alpha + beta);
        else r = -sin(alpha - beta) / sin(alpha + beta);
        return r;
    }

    double_complex ftrans(int pol, double alpha, double_complex n1, double_complex n2)
    {
        double_complex beta, t;
        if (pol == SENKRECHT)
            t = 2.0 * sin(beta) * cos(alpha) / sin(alpha + beta);
        else
            t = 2.0 * sin(beta) * cos(alpha) / (sin(alpha + beta) * cos(alpha - beta));
        return t;
    }
    }
}
#endif
