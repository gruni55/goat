#include "ray_pow.h"

namespace GOAT
{
    namespace raytracing
    {

        Ray_pow::Ray_pow() : IRay()
        {
        }

        Ray_pow::Ray_pow(double pow, const maths::Vector<double>& p,
            const maths::Vector<std::complex<double> >& Pol, const maths::Vector<double>& K,
            std::complex<double>  n0, double r0, double k0,
            const int numObjs = 0, ObjectShape** Objs = NULL) : IRay(p, Pol, K, n0, r0, k0, numObjs, Objs)
        {
            this->Pow = pow;
            E1 = Pol / abs(Pol);
            E2 = E1 * sqrt(Pow);
        }

        Ray_pow Ray_pow::reflect(maths::Vector<double> n, std::complex<double> n1, std::complex<double> n2)
            /** Strahl wird an einer Oberflaeche reflektiert. Wird an einer Einschlussoberflaeche
              * reflektiert (objIndex >-1 ), dann wird der transmittierte Strahl zurueckgegeben */
        {
            Ray_pow Erg;
            maths::Matrix <double> D, H, R;
            maths::Matrix<std::complex<double> > FR, FT;
            maths::Vector <double> Ph, e0, e1, e2;
            double alpha, gamma;

            // n=-n;

             /* Erst mal die Fresnelmatrix für Reflexion berechnen */

            double nk = n * k / (abs(n) * abs(k));
            double h;
            if (nk < 0) { n = -n; nk = -nk; }
            getKSystem(n, k, e0, e1, e2);
            if (nk > 1.0) nk = 1.0;
            alpha = std::acos(nk);
            gamma = M_PI - 2.0 * alpha;
            if (alpha > M_PI / 2.0) { gamma = -gamma; alpha = M_PI - alpha; }
            if (objIndex > -1) e0 = -e0;
            trafo(e0, e1, e2, H, R);
            FR = Fresnel_reflect(alpha, n1, n2);

            /* Fresnelmatrix für Transmission */
            FT(0, 0) = sqrt(abs(1 - maths::abs2(FR(0, 0))));
            FT(1, 1) = sqrt(abs(1 - maths::abs2(FR(1, 1))));
            FT(2, 2) = sqrt(abs(1 - maths::abs2(FR(2, 2))));
            D = rotMatrixA(e2, k, gamma);


            Erg.E1 = this->E1;
            Erg.E2 = this->E2;
            Erg.Pow = Pow;


            if (inObject)
            {
                // Strahl will aus dem Einschluss raus
                Erg = *this;
                Erg.OK = maths::dzero;
                Erg.inObject = false;
                Erg.n = n2;
                Erg.r0 = r0;
                Erg.objIndex = -1;
                Erg.refract(FT, n, n1, n2);
                Erg.Pow = abs2(Erg.E2);
            } // if inObject
            else
            {
                if (objIndex != -1)
                {
                    // Einschluss wurde getroffen
                    Erg = *this;
                    Erg.inObject = true;
                    // Erg.n=Ein[objIndex]->n;
                    Erg.n = n2;
                    Erg.OK = Obj[objIndex]->P; // ??
                    Erg.objIndex = objIndex;
                    Erg.refract(FT, n, n1, n2);

                    Erg.Pow = abs2(Erg.E2);
                    objIndex = -1;
                } // if objIndex!=-1
                else // Reflexion an der Partikeloberfläche
                {
                    Erg = *this;
                    Erg.inObject = false;
                    Erg.objIndex = -1;
                    Erg.n = n2;
                    Erg.OK = maths::dzero;
                    Erg.P = P;
                    Erg.E1 = E1;
                    Erg.E2 = E2;
                    Erg.k = k;
                    Erg.refract(FT, n, n1, n2);
                    Erg.Pow = abs2(Erg.E2);

                }
                /* E1 ist "nur" der Polarisationsvektor und zeigt die Richtung des E-Felds an */

            } // else if inObject

            k = D * k;
            iR++;
            E1 = D * E1;
            E2 = H * E2;
            E2 = FR * E2;
            E2 = R * E2;
            E2 = D * E2;

            Pow = abs2(E2);

            return Erg;
        }

        void Ray_pow::reflectRay(RayBase*& tray, maths::Vector<double> n, std::complex<double> n1, std::complex<double> n2)
        {
            Ray_pow t;
            t = reflect(n, n1, n2);
            *(Ray_pow*)tray = Ray_pow(t);
        }

        Ray_pow::~Ray_pow(void)
        {
        }


        void Ray_pow::refract(maths::Matrix<std::complex<double> > T, maths::Vector<double> N, std::complex<double> n1, std::complex<double> n2)
        {
            double s;
            maths::Matrix<double> H, R;
            maths::Matrix<double> D;
            maths::Vector <double> n, e0, e1, e2;
            double alpha, gamma;
            std::complex<double>  beta;

            n = N / abs(N);
            double nk = (n * k) / abs(k);
            // if (nk > 0) { n = -n;  nk = -nk; }
            getKSystem(n, k, e0, e1, e2);

            if (nk >= 1.0) alpha = 0.0;
            else  alpha = std::acos(nk);
            //  cout << "n=" << n << "   nk=" << nk << "   alpha=" << alpha << endl;
            if (alpha > M_PI / 2.0) { alpha = M_PI - alpha; e2 = -e2; }
            beta = asin((std::complex<double>) std::real(n1) / real(n2) * sin(alpha));
            gamma = real(beta) - alpha;
            s = 1.0;
            trafo(e0, e1, e2, H, R);
            D = rotMatrixA(e2, k, gamma);

            k = D * k;
            E1 = H * E1;
            E1 = T * E1;
            E1 = R * E1;
            E1 = D * E1;
            E1 = E1 / abs(E1);   // E1 ist nur der Polarisationsvektor und gibt somit nur die Polarisationsrichtung vor, d.h. |E1|=1

            E2 = H * E2;
            E2 = T * E2;
            E2 = R * E2;
            E2 = D * E2;
        }

        std::ostream& operator << (std::ostream& os, Ray_pow S)
        {
            os << "P=" << S.P << "   k=" << S.k << std::endl;
            os << "E1=" << S.E1 << std::endl;
            os << "E2=" << S.E2 << std::endl;
            os << "Pow=" << S.Pow << std::endl;
            os << "Anzahl Einschlüsse=" << S.numObj << std::endl;
            os << "inObject=" << S.inObject << "    objIndex=" <<
                S.objIndex << std::endl;
            os << "OK=" << S.OK << std::endl;
            return os;
        }
    }
}