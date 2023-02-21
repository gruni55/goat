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
#define SGN(x) (x<0) ? -1 : (x>0)
        SuperArray<maths::Vector<std::complex<double> > > verfolgung(maths::Vector<double> p0, maths::Vector<double> k0, SuperArray<maths::Vector<std::complex<double> > >&git);


       template<class T>  maths::Vector<double> pnext(maths::Vector<double> p0, maths::Vector<double> k0, SuperArray<T> &git);

        template<class T> maths::Vector<double> pnext(maths::Vector<double> p0, maths::Vector<double> k0, SuperArray<T> &git, double eps);
        #endif
        maths::Vector<double> pnext(maths::Vector<double> p0, maths::Vector<double> k0, maths::Vector<double> d, double eps=1E-50);
        maths::Vector<double> pnext(Plane E, maths::Vector<double> p0s, maths::Vector<double> k0s, maths::Vector<double> d, double eps=1E-50);


        /* ------------------------------ IMPLEMENTATION ---------------------------- */

#define NUM_EPS 1E-15

        template <class T> maths::Vector<double> pnext(maths::Vector<double> p0, maths::Vector<double> k0, SuperArray<T>& git)
        {

            double lambdax, lambday, lambdaz, lambda;
            double signx, signy, signz;
            double sx, sy, sz;

            signx = SGN(k0[0]);
            signy = SGN(k0[1]);
            signz = SGN(k0[2]);

            int nx = floor((p0[0] + git.r0 + NUM_EPS) / git.d[0]) + signx;
            int ny = floor((p0[1] + git.r0 + NUM_EPS) / git.d[1]) + signy;
            int nz = floor((p0[2] + git.r0 + NUM_EPS) / git.d[2]) + signz;

            lambdax = (nx  * git.d[0] - p0[0] - git.r0) / k0[0];
            lambday = (ny  * git.d[1] - p0[1] - git.r0) / k0[1];
            lambdaz = (nz  * git.d[2] - p0[2] - git.r0) / k0[2];

            // fabs avoids question after k0[i]=0 which can result into -inf !
            if (fabs(lambdax) < fabs(lambday))
                lambda = lambdax;
            else
                lambda = lambday;
            if (fabs(lambdaz) < lambda)
                lambda = lambdaz;
            return  p0 + lambda * k0;
        }

        template <class T> maths::Vector<double> pnext(maths::Vector<double> p0, maths::Vector<double> k0, SuperArray<T>& git, double eps)
        {
            double lambdax, lambday, lambdaz, lambda;
            double signx, signy, signz;
            double sx, sy, sz;

            signx = SGN(k0[0]);
            signy = SGN(k0[1]);
            signz = SGN(k0[2]);

            int nx = floor((p0[0] + git.r0 + NUM_EPS) / git.d[0]) + signx;
            int ny = floor((p0[1] + git.r0 + NUM_EPS) / git.d[1]) + signy;
            int nz = floor((p0[2] + git.r0 + NUM_EPS) / git.d[2]) + signz;

            lambdax = (nx * git.d[0] - p0[0] - git.r0) / k0[0];
            lambday = (ny * git.d[1] - p0[1] - git.r0) / k0[1];
            lambdaz = (nz * git.d[2] - p0[2] - git.r0) / k0[2];

            // fabs avoids question after k0[i]=0 which can result into -inf !
            if (fabs(lambdax) < fabs(lambday))
                lambda = lambdax;
            else
                lambda = lambday;
            if (fabs(lambdaz) < lambda)
                lambda = lambdaz;
            return  p0 + lambda * k0;
        }

    }
}
