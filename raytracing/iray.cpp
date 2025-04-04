/***************************************************************************
                          strahl.cpp  -  description
                             -------------------
    begin                : Fri Oct 15 1999
    copyright            : (C) 1999 by Thomas Weigel
    email                : weigel@lat.ruhr-uni-bochum.de
 ***************************************************************************/

#include "iray.h"
#include "matrix.h"

namespace GOAT
{
    namespace raytracing
    {   
        extern maths::Matrix<double> Dall;
        IRay::IRay()
        {
            iR = 0;
            objIndex = -1;
            OK = maths::dzero;
            isValid = true;
        }

        void  IRay::reflectRay(RayBase*& tray, maths::Vector<double> n, std::complex<double> n1, std::complex<double> n2)
        {
            IRay t;
            t = reflect(n, n1, n2);
            *(IRay*)tray = IRay(t);
        }

        IRay::~IRay() {
        }

        bool IRay::checkObjectIntersection(int& Index, maths::Vector<double>& Pmin)
        {
            maths::Vector<double> P1, n;
            double amin, a;
            bool found;
            bool isInside = (objIndex > -1);
            a = -1;
            amin = 100.0 * r0;
            Pmin = P;
            Index = -1;
            for (int i = 0; i < numObj; i++)
            {
                found = Obj[i]->next(P, k, P1);
                if (found)
                {
                    a = abs(P1 - P);
                    if (a < amin)
                    {
                        amin = a;
                        Pmin = P1;
                        /*if (i!=objIndex)*/ Index = i;
                    }
                }
            }
            found = (a > 0);
            return found;           
        }




        /*Vector<double> IRay::checkEinschluss(const Vector<double>& P,int& index)
        {
         Vector<double> P1,Pmin;
         double amin,a;
         bool found;

         amin=1E+99;
         index=-1;
         Pmin=P;
         for (int i=0; i<numObj; i++)
         {
         //  P1=nextP(P, k, Ein[i].P*r0,Ein[i].a*r0,found);
           found=Ein[i]->next(P,k,P1);
           if (found)
           {
            a=abs(P1-P);
            if (a<amin)
            {
              amin=a;
              Pmin=P1;
              if (i!=objIndex) index=i;
            }
           }
          }
         return Pmin;
        }*/

        IRay::IRay(const maths::Vector<double>& p,
            const maths::Vector<std::complex<double> >& Pol, const maths::Vector<double>& K,
            std::complex<double>  n0, double r0, double k0,
            const int numObj, std::vector<ObjectShape*> obj)
        {
            double l;
            E1 = maths::czero;
            E2 = maths::czero;
            P = p;
            k = K;
            Obj = obj;
            this->numObj = numObj;
            n = n0;
            this->r0 = r0;
            this->k0 = k0;
            l = k0 / (2.0 * M_PI);
            // KORR=1E-10;
             //init_Efeld(Pol);
            inObject = false;
            objIndex = -1;
            OK = maths::dzero;
            isValid = true;
            suppress_phase_progress = false;
        }


        bool IRay::next()
        {
            int Index;
            maths::Vector<double> p, R;
            bool found = true;
            bool inObject = true && objIndex > -1;
            if (!inObject) // Strahl ist ausserhalb eines Einschlusses
            {
                if (numObj == 0) Index = -1;
                else found = checkObjectIntersection(Index, R);  // Schnittpunkt mit Einschluss suchen
                if (Index == -1 || !found)    // Es wurde kein Einschluss gefunden !
                {
                    R = nextP(this->P, k, maths::dzero, r0, found); // Schnittpunkt mit der Außenkugel
                    if (!found) { R = this->P; Index = -2; } // Wenn nicht gefunden => Da stimmt was nicht 
                }

                objIndex = Index;
                if (!suppress_phase_progress)
                {
                    E1 = E1 * exp(I * k0 * n * abs(R - P));
                    E2 = E2 * exp(I * k0 * n * abs(R - P));
                }
                //  this->P=R;
            }
            else
            {
                found = Obj[objIndex]->next(P, k, R);
                if (!suppress_phase_progress)
                {
                    E1 = E1 * exp(I * k0 * Obj[objIndex]->n * abs(R - P));
                    E2 = E2 * exp(I * k0 * Obj[objIndex]->n * abs(R - P));
                }
            }

            this->P = R;
            /*if (abs(this->P-R)>1E-15) this->P=R;
            else objIndex=-3;*/
            return found;
        }

        
        void IRay::refract(maths::Vector<double> N, std::complex<double>  n1, std::complex<double>  n2)
        {
            double s;
            maths::Matrix<double> H, R;
            maths::Matrix<std::complex<double> > T;
            maths::Matrix<double> D;
            maths::Vector <double> n, e0, e1, e2;
            double alpha, gamma;
            std::complex<double>  beta;

            n = N / abs(N);
            getKSystem(n, k, e0, e1, e2);
            double nk = (n * k) / abs(k);
            alpha = std::acos(nk);
            if (alpha > M_PI / 2.0) { alpha = M_PI - alpha; e2 = -e2; }
            beta = asin((std::complex<double>) real(n1) / real(n2) * sin(alpha));
            if (imag(beta) > 1E-10) status = RAYBASE_STATUS_TIR;
            gamma = real(beta) - alpha;
            //if (real(n2)<real(n1)) { e2=-e2;}
            s = 1.0;
            trafo(e0, e1, e2, H, R);
            D = rotMatrixA(e2, k, gamma);
            if (beta != 0.0) T = Fresnel_trans(alpha, beta, n1, n2);
            else
            {
                T(0, 0) = (2.0 * n1) / (n2 + n1);
                T(1, 1) = T(0, 0);
                T(2, 2) = T(0, 0);
            }

            k = D * k;
            E1 = H * E1;
            E1 = T * E1;
            E1 = R * E1;
            E1 = D * E1;

            E2 = H * E2;
            E2 = T * E2;
            E2 = R * E2;
            E2 = D * E2;
        }




        IRay IRay::reflect(maths::Vector<double> n, std::complex<double>  n1, std::complex<double>  n2)
            /* Strahl wird an einer Oberflaeche reflektiert. Wird an einer Einschlussoberflaeche
               reflektiert (objIndex >-1 ), dann wird der transmittierte Strahl zurueckgegeben */
        {
            IRay Erg;
            maths::Matrix <double> D, H, R;
            maths::Matrix<std::complex<double> > F;
            maths::Vector <double> Ph, e0, e1, e2;
            double alpha, gamma;

            Erg.E1 = this->E1;
            Erg.E2 = this->E2;

            if (inObject)
            {
                // Strahl will aus dem Einschluss raus
                Erg = *this;
                Erg.OK = maths::dzero;
                Erg.inObject = false;
                Erg.n = n2;
                Erg.r0 = r0;
                Erg.objIndex = -1;
                Erg.refract(n, n1, n2);
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
                    Erg.OK = Obj[objIndex]->P;
                    Erg.objIndex = objIndex;
                    Erg.refract(n, n1, n2);
                    Erg.suppress_phase_progress = suppress_phase_progress;
                    objIndex = -1;
                } // if objIndex!=-1
                else // Reflexion an der Partikeloberfläche
                {
                    Erg = *this;
                    Erg.suppress_phase_progress = suppress_phase_progress;
                    Erg.inObject = false;
                    Erg.objIndex = -1;
                    Erg.n = n2;
                    Erg.OK = maths::dzero;
                    Erg.P = P;
                    Erg.E1 = E1;
                    Erg.E2 = E2;
                    Erg.k = k;
                    Erg.refract(n, n1, n2);
                }
            } // else if inObject

            n = -n;
            double nk = n * k / (abs(n) * abs(k));
            getKSystem(n, k, e0, e1, e2);
            alpha = std::acos(nk);
            gamma = M_PI - 2.0 * alpha;
            if (alpha > M_PI / 2.0) { gamma = -gamma; alpha = M_PI - alpha; }
            if (objIndex > -1) e0 = -e0;
            trafo(e0, e1, e2, H, R);
            F = Fresnel_reflect(alpha, n1, n2);
            D = rotMatrixA(e2, k, gamma);

            k = D * k;

            E1 = H * E1;
            E1 = F * E1;
            E1 = R * E1;
            E1 = D * E1;

            E2 = H * E2;
            E2 = F * E2;
            E2 = R * E2;
            E2 = D * E2;

            iR++;
            return Erg;
        }

        maths::Matrix<std::complex<double> > IRay::Fresnel_reflect(double alpha, std::complex<double>  n1, std::complex<double>  n2)
        {
            std::complex<double>  rs, rp;
            maths::Matrix<std::complex<double> > R;
            double  n12;
            double  beta;

            if (alpha != 0.0)
            {
                n12 = real(n2) / real(n1);
                beta = real(asin((std::complex<double>) sin(alpha) / n12));
                rp = tan(alpha - beta) / tan(alpha + beta);
                rs = -sin(alpha - beta) / sin(alpha + beta);
            }
            else
            {
                rp = (real(n2 - n1)) / real(n1 + n2);
                rs = -rp;
            }

            R(0, 0) = rp;
            R(1, 1) = rp;
            R(2, 2) = rs;
            return R;
        }

        maths::Matrix<std::complex<double> > IRay::Fresnel_trans(double alpha, std::complex<double>  beta, std::complex<double>  n1, std::complex<double>  n2)
        {
            std::complex<double>  ts, tp;
            maths::Matrix<std::complex<double> > T;

            if (alpha == 0.0)
            {
                ts = 2.0 * real(n1) / real(n1 + n2);
                tp = ts;
            }
            else
            {
                ts = 2.0 * sin(beta) * cos(alpha) / sin(alpha + beta);
                tp = 2.0 * sin(beta) * cos(alpha) / (sin(alpha + beta) * cos(alpha - beta));
            }

            T(0, 0) = tp;
            T(1, 1) = tp;
            T(2, 2) = ts;

            return T;
        }

        void IRay::initElectricField(const Plane& Eb, const maths::Vector<std::complex<double> >& Pol, const int AnzRays)
        {
            double x1, x2, xn, p;
            maths::Vector<double> r;
            k /= abs(k);
            x1 = P * Eb.e1;
            x2 = P * Eb.e2;
            xn = P * Eb.n;
            r = P;
            p = x1 * x1 + x2 * x2;
            P = x1 * Eb.e1 + x2 * Eb.e2 - sqrt(fabs(r0 * r0 - p)) * Eb.n;
            E1 = Pol * exp(I * abs(r - P) * k0);
            E2 = Pol * exp(I * abs(r - P) * k0);
        }

        void IRay::initElectricField(const Plane& Eb, const maths::Vector<std::complex<double> >& Pol1, const maths::Vector<std::complex<double> >& Pol2, const int AnzRays)
        {
            double x1, x2, xn, p;
            maths::Vector<double> r;
            k /= abs(k);
            x1 = P * Eb.e1;
            x2 = P * Eb.e2;
            xn = P * Eb.n;
            r = P;
            p = x1 * x1 + x2 * x2;
            P = x1 * Eb.e1 + x2 * Eb.e2 - sqrt(fabs(r0 * r0 - p)) * Eb.n;
            E1 = Pol1 * exp(I * abs(r - P) * k0);
            E2 = Pol2 * exp(I * abs(r - P) * k0);
        }



        void IRay::initElectricField(
            const maths::Vector<std::complex<double> >& PolS,
            const maths::Vector<std::complex<double> >& PolP,
            const int AnzRays)
        {
            maths::Vector<double> hr, h;
            hr = P;
            h = P - (P * k) * k;
            P = h - sqrt(r0 * r0 - abs(h) * abs(h)) * k;
            E1 = PolS * exp(I * abs(hr - P) * k0);
            E2 = PolP * exp(I * abs(hr - P) * k0);
        }

        void IRay::initElectricFieldGauss(const Plane& Eb,
            const maths::Vector<std::complex<double> >& PolS,
            const maths::Vector<std::complex<double> >& PolP,
            Gauss g)
        {
            k = g.F - P;
            k = k / abs(k);
            double theta;
            double z = abs(P - g.F);
            double r;
            double R;
            double l0 = 2.0 * M_PI / std::real(k0);
            maths::Vector<double> hr, h;
            std::complex<double> E0;
            double w, G;

            double z0 = M_PI * g.w0 * g.w0 / l0;
            theta = atan2(g.w0, z0);
            r = sin(theta / 2.0) * z;
            w = g.w0 * sqrt(1.0 + z * z / (z0 * z0));
            R = z * sqrt(1.0 + z * z / (z0 * z0));
            G = atan2(z, z0);
            E0 = g.w0 / w * exp((-r * r / w * w) - I * k0 * r * r / (2 * R) - I * k0 * z + I * G);
            double l;
            double Pk = P * k;
            l = -Pk - sqrt(Pk * Pk - P * P + r0 * r0);
            P = P + l * k;
            E1 = PolS * E0 * exp(I * k0 * l);
            E2 = PolP * E0 * exp(I * k0 * l);
        }
        void IRay::initElectricFieldGauss(maths::Vector<std::complex<double> >& Pol, Gauss g)
        {
            double w2, w02, r, zi0, h2, alpha, phi, det, l, R;

            maths::Matrix<double> D;
            maths::Vector<double> n, Pkt, Ps, d, hk;

            hk = k;
            d = P - g.F * r0;
            r = sqrt(d[0] * d[0] + d[1] * d[1]);
            w02 = g.w0 * g.w0 * r0 * r0;
            zi0 = std::real(k0) * w02 / 2.0;
            h2 = std::real(k0) * w02 / (2 * d[2]);
            R = d[2] * (1 + h2 * h2);

            // Normale zur Phasenfläche
            n[0] = std::real(k0) * d[0] / R;
            n[1] = std::real(k0) * d[1] / R;
            n[2] = -1 / (zi0 * (1 + d[2] * d[2])) - std::real(k0) + std::real(k0) * r * r / (R * R) * (1 - std::real(k0) * std::real(k0) * w02 * w02 / (2 * d[2] * d[2]));
            k = n / maths::abs(n);

            // Nächsten Punkt auf der Kugel berechnen
            det = (P * k) * (P * k) - (abs2(P) - r0 * r0);
            if (det > 0) // Schnittpunkt gefunden
            {
                l = (-P * k - std::sqrt(det));
                P = P + l * k;
                h2 = 2 * d[2] / (std::real(k0) * w02);
                w2 = w02 * (1 + h2 * h2);
                alpha = abs(maths::ez * k);
                phi = atan(d[2] / zi0) - std::real(k0) * (d[2] + r * r / (2 * R));
                E1 = Pol * g.w0 * r0 / w2 * exp(-r * r / w2) * exp(I * phi);
                D = maths::rotMatrix(k % hk, alpha);
                E1 = D * E1;
                E2 = E1;
            }
        }



        double IRay::cross(const maths::Vector<double> P10, const maths::Vector<double> P11,
            const maths::Vector<double> P20, const maths::Vector<double> P21)
        {
            maths::Vector<double> a0, a1;
            double k1, k2, k;

            a0 = P20 - P10;
            a1 = P21 - P11;
            k1 = a0[1] / a0[0] * (P11[0] - P10[0]) + P10[1] - P11[1];
            k2 = a1[1] - a0[1] / a0[0] * a1[0];
            k = k1 / k2;
            return k;
        }

        maths::Vector<double> IRay::crossPlane(const maths::Vector<double> Pe, const maths::Vector<double> n)
            /* Berechnet den Schnittpunkt des zentralen Strahls mit einer Ebene, die durch den Aufpunkt Pe und die Normale n beschrieben wird
             * Die Routine liefert inf-Vektor, falls kein Schnittpunkt */
        {
            double l = k * n;
            double h = (Pe - P) * n;
            maths::Vector<double> S;

            if ((l == 0) && (h == 0)) return P; // Strahl befindet sich in der Ebene und bewegt sich in der Ebene
            if ((l == 0) && (h != 0)) return maths::Vector<double>(-1, -1, -1);
            l = h / l;
            S = P + l * k;
            return S;
        }


        std::ostream& operator << (std::ostream& os, IRay S)
        {
            os << "P=" << S.P << "   k=" << S.k << std::endl;
            os << "E1=" << S.E1 << std::endl;
            os << "E2=" << S.E2 << std::endl;
            os << "Anzahl Einschlüsse=" << S.numObj << std::endl;
            os << "inObject=" << S.inObject << "    objIndex=" <<
                S.objIndex << std::endl;
            os << "OK=" << S.OK << std::endl;
            return os;
        }



        void IRay::initElectricFieldGauss(double sigma2, maths::Vector<double> focuspos, maths::Vector<std::complex<double> > Pol)
        {
            k = (focuspos - P);
            k = k / abs(k);

            double r = (focuspos - P) * k;
            double h = r - abs(focuspos - P);

            E1 = Pol * exp(-I * k0 * h) * exp(-(r * r - h * h) / sigma2);
            E2 = E1;
        }

        maths::Vector<double> IRay::intersectRect(const maths::Vector<double> Pe, const maths::Vector<double> e1, const maths::Vector<double> e2)
        {
            maths::Matrix<double> M(e1, e2, -k);
            bool inv;
            maths::Matrix<double>  Mi = invert(M, inv);
            if (inv)
            {
                maths::Vector<double> Ps = P - Pe;
                maths::Vector<double> L = Mi * Ps;
                if ((L[0] < 0.0) || (L[0] > 1.0) || (L[1] < 0.0) || (L[1] > 1.0)) return maths::nanV("");
                return P + L[2] * k;
            }
        }
    }
}