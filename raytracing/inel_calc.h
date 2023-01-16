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
        SuperArray<maths::Vector<std::complex<double> > > verfolgung(maths::Vector<double> p0, maths::Vector<double> k0, SuperArray<maths::Vector<std::complex<double> > >&git);


       template<class T>  maths::Vector<double> pnext(maths::Vector<double> p0, maths::Vector<double> k0, SuperArray<T> &git);

        template<class T> maths::Vector<double> pnext(maths::Vector<double> p0, maths::Vector<double> k0, SuperArray<T> &git, double eps);
        #endif
        maths::Vector<double> pnext(maths::Vector<double> p0, maths::Vector<double> k0,
        maths::Vector<double> d, double eps);
        maths::Vector<double> pnext(Plane E, maths::Vector<double> p0s, maths::Vector<double> k0s,
        maths::Vector<double> d, double eps);


        /* ------------------------------ IMPLEMENTATION ---------------------------- */
        template <class T> maths::Vector<double> pnext(maths::Vector<double> p0, maths::Vector<double> k0, SuperArray<T>& git)
        {

            double lambdax, lambday, lambdaz, lambda;
            double signx, signy, signz;

            signx = copysign(1, k0[0]); signy = copysign(1, k0[1]); signz = copysign(1, k0[2]);

            lambdax = ((floor(p0[0] / git.d[0]) + signx + (signx < 0) * (fmod(p0[0], git.d[0]) != 0)) * git.d[0] - p0[0]) / k0[0];
            lambday = ((floor(p0[1] / git.d[1]) + signy + (signy < 0) * (fmod(p0[1], git.d[1]) != 0)) * git.d[1] - p0[1]) / k0[1];
            lambdaz = ((floor(p0[2] / git.d[2]) + signz + (signz < 0) * (fmod(p0[2], git.d[2]) != 0)) * git.d[2] - p0[2]) / k0[2];

            if (lambdax < lambday)
                lambda = lambdax;
            else
                lambda = lambday;

            if (lambda > lambdaz)
                lambda = lambdaz;

            return  p0 + lambda * k0;

        }

        template <class T> maths::Vector<double> pnext(maths::Vector<double> p0, maths::Vector<double> k0, SuperArray<T>& git, double eps)
        {
            double lambdax, lambday, lambdaz, lambda;
            double signx, signy, signz;
            double sx, sy, sz;

            /***********************************************************************************************
             * Algorithmus zur Bestimmung der Schnittpunkte des Strahls mit den gridebenen               *
             *                                                                                             *
             * Prinzip:                                                                                    *
             * Ausgehend vom Punkt p_0 werden die Schnittpunkte mit Ebenen x=const.,                       *
             * y=const. und z=const. gebildet. Hierzu werden zunaechst die Stuetzpunkte s_x, s_y und s_z,  *
             * welche den Abstand dieser Ebenen vom Ursprung des grids angeben, berechnet. Der           *
             * Schnittpunkt mit dem kleinsten Abstand ist dann der gesuchte Punkt.                         *
             *                                                                                             *
             * Die Stuetzpunkte s_i bestimmen sich folgendermassen:                                        *
             *                                                                                             *
             *                                                                                             *
             * Die i-Komponente des Punktes p_0 ist von der Form                                           *
             *                                                                                             *
             *  p_{0,i} = m*di + a,  m ganz, a >=0                                                         *
             *                                                                                             *
             * wobei di die Schrittweite des grids in i-Richtung ist.                                    *
             *                                                                                             *
             * Ist k_i > 0 und a = 0, so ist der Stuetzpunkt s_i gegeben durch (m+1)*di=(m+sign(k_i)+0)*di *
             * Ist k_i > 0 und a > 0, so ist der Stuetzpunkt s_i gegeben durch (m+1)*di=(m+sign(k_i)+0)*di *
             * Ist k_i < 0 und a = 0, so ist der Stuetzpunkt s_i gegeben durch (m-1)*di=(m+sign(k_i)+0)*di *
             * Ist k_i < 0 und a > 0, so ist der Stuetzpunkt s_i gegeben durch   m*di  =(m+sign(k_i)+1)*di *
             *                                                                                             *
             * Der Schnittpunkt mit dem kleinsten Abstand ist dann der gesuchte Punkt.                     *
             ***********************************************************************************************/

            signx = copysign(1.0, k0[0]); signy = copysign(1.0, k0[1]); signz = copysign(1.0, k0[2]);

            sx = (floor((p0[0] + signx * 2 * eps) / git.d[0]) + signx + (signx < 0) * (fmod(p0[0] + signx * 2 * eps, git.d[0]) != 0)) * git.d[0];
            sy = (floor((p0[1] + signy * 2 * eps) / git.d[1]) + signy + (signy < 0) * (fmod(p0[1] + signy * 2 * eps, git.d[1]) != 0)) * git.d[1];
            sz = (floor((p0[2] + signz * 2 * eps) / git.d[2]) + signz + (signz < 0) * (fmod(p0[2] + signz * 2 * eps, git.d[2]) != 0)) * git.d[2];

            lambdax = (sx - p0[0]) / k0[0];
            lambday = (sy - p0[1]) / k0[1];
            lambdaz = (sz - p0[2]) / k0[2];

            if (lambdax <= lambday)
                lambda = lambdax;
            else
                lambda = lambday;

            if (lambda >= lambdaz)  lambda = lambdaz;


            return  p0 + lambda * k0;

        }

    }
}
