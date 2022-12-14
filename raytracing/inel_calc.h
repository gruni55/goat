#include <complex>
#include <math.h>
#include "vector.h"
#ifdef WITH_SUPERGITTER
#include "superarray.h"
#endif
#include "grid.h"
#include "plane.h"
namespace GOAT
{
    namespace raytracing
    {
        // Berechnet den Punkt ix,iy eines Einstrahlgrids, welches sich im Abstand r
        // unter einem Winkel theta,phi vor dem Partikel befindet
        maths::Vector<double> startpunkt(int ix, int iy,  grid& git, 
                                double r, double theta, double phi);

        // Berechnet den Punkt ix,iy eines Einstrahlgrids, welches sich
        // unter einem Winkel theta,phi direkt vor dem Partikel befindet
        maths::Vector<double> startpunkt(int ix, int iy, grid& git,
                                double theta, double phi);

        // Dreht den Vektor vein um theta und um phi
        maths::Vector<double> drehvektor(maths::Vector<double> vein, double theta, double phi);

        // Dreht den Vektor vein um phi
        maths::Vector<double> drehphivektor(maths::Vector<double> vein, double phi);

        // Gibt den Ursprung eines Einstrahlgrids wieder, dessen Mittelpunkt um  
        // Winkel theta,phi gedreht wurde und das urspr�nglich um dr �ber der Ebene 
        // z = git.zmax lag. Wird in der Routine startpunkt benoetigt
        maths::Vector<double> ursprungrot(double dr, double theta, double phi, 
                                maths::Vector<double> &exrot, maths::Vector<double> &eyrot, 
                                grid &git);

        // Berechnet die eigentliche Strahlverfolgung im grid git bei Einstrahlung  
        // an dem Startpunkt p0  in Richtung k0


        #ifdef WITH_SUPERGITTER
        SuperArray verfolgung(maths::Vector<double> p0, maths::Vector<double> k0, SuperArray &git);


        maths::Vector<double> pnext(maths::Vector<double> p0, maths::Vector<double> k0, SuperArray &git);

        maths::Vector<double> pnext(maths::Vector<double> p0, maths::Vector<double> k0, SuperArray &git, double eps);
        #endif
        maths::Vector<double> pnext(maths::Vector<double> p0, maths::Vector<double> k0,
        maths::Vector<double> d, double eps);
        maths::Vector<double> pnext(Plane E, maths::Vector<double> p0s, maths::Vector<double> k0s,
        maths::Vector<double> d, double eps);
    }
}
