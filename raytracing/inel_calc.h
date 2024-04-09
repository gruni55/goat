#include <complex>
#include <math.h>
#include "vector.h"
#ifdef WITH_SUPERGITTER
#include "superarray.h"
#endif
#include "grid.h"
#include "plane.h"
#include <cfloat>
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

        GOAT::maths::Vector<INDEX_TYPE> currentIndex(-1, -1, -1);
        

       template<class T>  maths::Vector<double> pnext(maths::Vector<double> p0, maths::Vector<double> k0, SuperArray<T> &git);

        template<class T> maths::Vector<double> pnext(maths::Vector<double> p0, maths::Vector<double> k0, SuperArray<T> &git, double eps);
        #endif
        maths::Vector<double> pnext(maths::Vector<double> p0, maths::Vector<double> k0, maths::Vector<double> d, double eps=1E-50);
        maths::Vector<double> pnext(Plane E, maths::Vector<double> p0s, maths::Vector<double> k0s, maths::Vector<double> d, double eps=1E-50);
        /**
        * This function describes the next crossing point of a ray. The ray is described by a point P on the ray  and a direction vector k. 
        * (P and k are already transformed into the local coordinate system to make the calculation faster
        * The dimensions of a grid cell is described by the vector d. 
        */
     //   maths::Vector<double> pnext(maths::Vector<double> p0, maths::Vector<double> k0, maths::Vector<double> d);


        /* ------------------------------ IMPLEMENTATION ---------------------------- */

#define NUM_EPS 1E-10

        template <class T> maths::Vector<double> pnext(maths::Vector<double> p0, maths::Vector<double> k0, SuperArray<T>& git)
        {
            double lambdax, lambday, lambdaz, lambda;
            // double signx, signy, signz;
            double sx, sy, sz;
            GOAT::maths::Vector<INDEX_TYPE> index;

            if (currentIndex[0] == -1)
            {

                INDEX_TYPE nx = (p0[0] + git.r0) / git.d[0];
                INDEX_TYPE ny = (p0[1] + git.r0) / git.d[1];
                INDEX_TYPE nz = (p0[2] + git.r0) / git.d[2];

                currentIndex = GOAT::maths::Vector<int>(nx, ny, nz);
            }

            int signx = SGN(k0[0]);
            int signy = SGN(k0[1]);
            int signz = SGN(k0[2]);
            index[0] = currentIndex[0] + signx;
            index[1] = currentIndex[1] + signy;
            index[2] = currentIndex[2] + signz;

            lambdax = (index[0] * git.d[0] - p0[0] - git.r0) / k0[0];
            lambday = (index[1] * git.d[1] - p0[1] - git.r0) / k0[1];
            lambdaz = (index[2] * git.d[2] - p0[2] - git.r0) / k0[2];

            // fabs avoids question after k0[i]=0 which can result into -inf !
            int i = 0;
            if (fabs(lambdax) < fabs(lambday))
                lambda = lambdax;
            else
            {
                lambda = lambday;
                i = 1;
            }
            if (fabs(lambdaz) < lambda)
            {
                lambda = lambdaz;
                i = 2;

            }
            currentIndex[i] = index[i];
            return  p0 + lambda * k0;
            
        }

        template <class T> maths::Vector<double> pnext(maths::Vector<double> p0, maths::Vector<double> k0, SuperArray<T>& git, double eps)
        {
            double lambdax = DBL_MAX; 
            double lambday = DBL_MAX; 
            double lambdaz = DBL_MAX, lambda;
           // double signx, signy, signz;
            double sx, sy, sz;
            GOAT::maths::Vector<INDEX_TYPE> index;

           // if (currentIndex[0] == -1)
            {             

                INDEX_TYPE nx = (p0[0] + git.r0) / git.d[0];
                INDEX_TYPE ny = (p0[1] + git.r0) / git.d[1];
                INDEX_TYPE nz = (p0[2] + git.r0) / git.d[2];
                
                currentIndex = GOAT::maths::Vector<INDEX_TYPE>(nx, ny, nz);
            }
            
                int signx = SGN(k0[0]);
                int signy = SGN(k0[1]);
                int signz = SGN(k0[2]);
                index[0] = currentIndex[0] + signx;
                index[1] = currentIndex[1] + signy;
                index[2] = currentIndex[2] + signz;

                if (k0[0] != 0) lambdax = ((double)index[0] * git.d[0] - p0[0] - git.r0) / k0[0];
                if (k0[1] != 0) lambday = ((double)index[1] * git.d[1] - p0[1] - git.r0) / k0[1];
                if (k0[2] != 0) lambdaz = ((double)index[2] * git.d[2] - p0[2] - git.r0) / k0[2];
            
           //     std::cout << "pnext: " << p0<< ":" << k0 << "/" << lambdax << "," << lambday << "," << lambdaz << std::endl;
            int i = 0;
            if ((lambdax < lambday) && (lambdax > 10.0*DBL_MIN))
                lambda = lambdax;
            else
            {
                lambda = lambday;
                i = 1;
            }
            if ((lambdaz < lambda) && (lambdaz > 10.0 * DBL_MIN) || lambda < 10.0 * DBL_MIN)
            {
                lambda = lambdaz;
                i = 2;
               
            }
             currentIndex[i] = index[i];    
             //std::cout << "currentIndex=" << currentIndex << std::endl;
            return  p0 + lambda * k0;
        }
        /*
        maths::Vector<double> pnext(maths::Vector<double> p0, maths::Vector<double> k0, maths::Vector<double> d)
        {
            double lambdax, lambday, lambdaz, lambda;
            double signx, signy, signz;
            double sx, sy, sz;

            signx = SGN(k0[0]);
            signy = SGN(k0[1]);
            signz = SGN(k0[2]);

            int nx = floor(p0[0] / d[0]) + signx;
            int ny = floor(p0[1] / d[1]) + signy;
            int nz = floor(p0[2] / d[2]) + signz;

            lambdax = (nx * d[0] - p0[0] ) / k0[0];
            lambday = (ny * d[1] - p0[1] ) / k0[1];
            lambdaz = (nz * d[2] - p0[2] ) / k0[2];

            // fabs avoids question after k0[i]=0 which can result into -inf !
            if (fabs(lambdax) < fabs(lambday))
                lambda = lambdax;
            else
                lambda = lambday;
            if (fabs(lambdaz) < lambda)
                lambda = lambdaz;
            return  p0 + lambda * k0;
        }*/

    }
}
