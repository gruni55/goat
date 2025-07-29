#include "asphericLens.h"
namespace GOAT
{
    namespace raytracing
    {

        asphericLens::asphericLens(const maths::Vector<double>& P,
            std::complex<double>  n, asphericLensParms lensParms, 
            const maths::Matrix<std::complex<double> > alpha,
            const maths::Vector<double>& Ex,
            const maths::Vector<double>& Ey,
            const maths::Vector<double>& Ez) : ObjectShape(P,n,alpha,Ex,Ey,Ez,OBJECTSHAPE_ASPHERIC_LENS)
        {
            

            // first, let's calculate the shift for the left side

            if (lensParms.left.isPlano)
            {                
                lensParms.left.shift = -lensParms.offset / 2.0;
            }
            else
            {
                if (lensParms.left.R > 0) // is convex
                {
                    lensParms.left.shift = z(lensParms.radius,lensParms.left)-lensParms.offset/2.0;
                }
                else // side is concave
                {
                    lensParms.left.shift = lensParms.offset / 2.0 - z(0, lensParms.left);
                }
            }

            // now, let's calculate the shift for the right side

            if (lensParms.right.isPlano)
            {
                lensParms.right.shift = lensParms.offset / 2.0;
            }
            else
            {
                if (lensParms.right.R > 0) // is convex
                {
                    lensParms.right.shift = z(lensParms.radius,lensParms.right) - lensParms.offset / 2.0;
                }
                else // side is concave
                {
                    lensParms.right.shift = lensParms.offset / 2.0 - z(lensParms.radius, lensParms.right);
                }
            }    
            this->lensParms = lensParms;
        }

        void asphericLens::binWrite() // todo
        {
        }

        void asphericLens::binRead() // todo
        {
        }

        void asphericLens::scale(double sf) // todo
        {
        }



        inline double asphericLens::zleft(maths::Vector<double> P, double lambda)
        {
            double r2 = P[0] * P[0] + P[1] * P[1] + 2.0 * lambda * (P[0] * k[0] + P[1] * k[1]) + lambda * lambda * (k[0] * k[0] + k[1] * k[1]);
            double r = sqrt(r2);
            double rn = 1.0;
            double R = lensParms.left.R;
            double kap = lensParms.left.k;
            double z = r2 / (R * (1 + sqrt(1 - (1 + kap) * r2 / (R * R))));
            for (int l = 0; l < lensParms.left.A.size(); l++)
            {
                rn *= r;
                z += rn * lensParms.left.A[l];
            }

            return z - lensParms.left.shift;
        }

        inline double asphericLens::zright(maths::Vector<double> P, double lambda)
        {
            double r2 = P[0] * P[0] + P[1] * P[1] + 2.0 * lambda * (P[0] * k[0] + P[1] * k[1]) + lambda * lambda * (k[0] * k[0] + k[1] * k[1]);
            double r = sqrt(r2);
            double rn = 1.0;
            double R = lensParms.right.R;
            double kap = lensParms.right.k;
            double z = r2 / (R * (1 + sqrt(1 - (1 + kap) * r2 / (R * R))));
            for (int l = 0; l < lensParms.right.A.size(); l++)
            {
                rn *= r;
                z += rn * lensParms.right.A[l];
            }

            return z + lensParms.left.shift;
        }

        double asphericLens::z(double r, asphericLensSide side)
        {
            double R = lensParms.left.R;
            double erg = r * r / (R * (1 + sqrt(1 - (1 + side.k) * r * r / (R * R))));
            double rn = r;
            for (int i = 0; i < side.A.size(); i++)
            {
                erg += side.A[i] * rn;
                rn *= r;
            }
            return erg;
        }

        double asphericLens::dzdrleft(maths::Vector<double> P, double lambda)
        {
            double r2 = P[0] * P[0] + P[1] * P[1] + 2.0 * lambda * (P[0] * k[0] + P[1] * k[1]) + lambda * lambda * (k[0] * k[0] + k[1] * k[1]);
            double r = sqrt(r2);
            double R = lensParms.left.R;
            double D1 = (1 + lensParms.left.k) * r2 / (R * R);
            double D2 = 1 + sqrt(1 - D1);
            double dz = r / (R * D2) * (2 + D1 / D2);
            double rn = 1.0;
            for (int l = 0; l < lensParms.left.A.size(); l++)
            {                
                dz += rn * l * lensParms.left.A[l];
                rn *= r;
            }
            return dz;
        }

        double asphericLens::dzdrright(maths::Vector<double> P, double lambda)
        {
            double r2 = P[0] * P[0] + P[1] * P[1] + 2.0 * lambda * (P[0] * k[0] + P[1] * k[1]) + lambda * lambda * (k[0] * k[0] + k[1] * k[1]);
            double r = sqrt(r2);
            double R = lensParms.right.R;
            double D1 = (1 + lensParms.right.k) * r2 / (R * R);
            double D2 = 1 + sqrt(1 - D1);
            double dz = r / (R * D2) * (2 + D1 / D2);
            double rn = 1.0;
            for (int l = 0; l < lensParms.right.A.size(); l++)
            {
                dz += rn * l * lensParms.right.A[l];
                rn *= r;
            }
            return dz;
        }


        /*
        inline double asphericLens::dzright(double lambda)
        {
             return dzdrright(lambda)*drdlambda;
        }

        double asphericLens::dzleft(double lambda)
        {
            return dzdrleft(lambda) * drdlambda;
        }
        
        */
        bool asphericLens::next(const maths::Vector<double>& ps, const maths::Vector<double>& ks, maths::Vector<double>& pout)
        {            
            maths::Vector<double> er;
            Side side;
            P = H * ps;
            k = H * ks;
            double lambdaLeft, lambdaRight, lambdaLateral, lambda = 0;
            double dr1 = P[0] * k[0] + P[1] * k[1];
            double dr2 = k[0] * k[0] + k[1] * k[1];
            drdlambda = (dr1 + lambda * dr2) / sqrt(P[0] * P[0] + P[1] * P[1] + 2 * lambda * dr1 + lambda * lambda * dr2);

            // test intersection with left side
            if (lensParms.left.isPlano)
            {                
                lambdaLeft = (ps[2] - lensParms.offset / 2.0) / k[2];
            }
            else
            {
                double lold, lnew;
                lold = 0.0;
                double dl;
                do
                {                    
                    lnew = lold - (zleft(P, lold)  - lold * k[2]) / (dzdrleft(P, lold) - k[2]);
                    dl = fabs(lold - lnew);
                    lold = lnew;
                } while (dl > 1E-100);
                lambdaLeft = lnew;
            }
            
            // test intersection with right side
            
            if (lensParms.right.isPlano)
            {
                lambdaRight = (ps[2] + lensParms.offset / 2.0) / k[2];
            }
            else
            {
                double lold, lnew;
                lold = 0.0;
                double dl;
                do
                {
                    lnew = lold - (zright(P, lold) - lold * k[2]) / (dzdrright(P, lold) - k[2]);
                    dl = fabs(lold - lnew);
                    lold = lnew;
                } while (dl > 1E-100);
                lambdaRight = lnew;
            }

            // test intersection with lateral surface
            double a, b, c, D;
            a = dr2;
            if (a == 0) lambdaLateral = -1;
            else
            {
                b = 2.0 * dr1;
                c = P[0] * P[0] + P[1] * P[1] - lensParms.radius;
                D = b * b - 4.0 * a * c;
                if (D < 0) lambdaLateral = -1;
                else
                {             
                   er=maths::Vector<double>(P[0], P[1], 0);
                    er = er / abs(er);
                    double sqrtD = sqrt(D);
                    lambdaLateral = (-b - sqrtD) / (2.0 * a);
                    if (lambdaLateral<0) lambdaLateral= (-b + sqrtD) / (2.0 * a);
                    if (lambdaLateral > 0)
                    {
                        pout = P + lambdaLateral * k;
                        if ((pout[2] < -lensParms.offset / 2.0) || (pout[2] > lensParms.offset / 2.0))
                            lambdaLateral = -1;
                    }
                }
            }

            if ((lambdaLeft < lambdaRight) && (lambdaLeft > 0)) {
                side = left;  lambda = lambdaLeft;
            }
            else {
                side = right;  lambda = lambdaRight;
            }
            if ((lambda < 0) || ((lambdaLateral < lambda) && (lambdaLateral > 0)))
            {
                lambda = lambdaLateral;
                side = lateral;
            }

            if (lambda < 0) return false;
            pout = P + lambda * k;
            double dzdr;

            switch (side)
            {
              case left: dzdr = dzdrleft(P,lambdaLeft); currentnorm = maths::ez + dzdr * er; break;
              case right: dzdr = dzdrright(P,lambdaRight); currentnorm = maths::ez + dzdr * er; break;
              case lateral: currentnorm = er;
            }
            return true;
        }
        maths::Vector<double> asphericLens::norm(const maths::Vector<double>& P)
        {
            return currentnorm;
        }
        bool asphericLens::isInside(const maths::Vector<double>& p) // todo
        {
            return false;
        }
        double asphericLens::volume() // todo
        {
            return 0.0;
        }
    }
}
