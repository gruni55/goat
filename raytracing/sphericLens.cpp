#include "sphericLens.h"

namespace GOAT
{
    namespace raytracing
    {
        sphericLens::sphericLens(const maths::Vector<double>& P, std::complex<double> n, lensParms lp, GOAT::maths::Matrix<std::complex<double>> alpha, const maths::Vector<double>& Ex, const maths::Vector<double>& Ey, const maths::Vector<double>& Ez, const int type)
            : ObjectShape(P, n, alpha, Ex, Ey, Ez, type)
        {
            this->lp = lp;

            
            init();
        }

        void sphericLens::binWrite(std::ofstream& os)
        {
            P.binWrite(os);
            H.binWrite(os);
            R.binWrite(os);
            os.write((char*)&n, (char)sizeof(n));
            alpha.binWrite(os);
            pul.binWrite(os);
            por.binWrite(os);
            for (int i = 0; i < 3; i++)
                e[i].binWrite(os);
            os.write((char*)&Ealpha, (char)sizeof(Ealpha));
            os.write((char*)&Ebeta, (char)sizeof(Ebeta));
            os.write((char*)&Egamma, (char)sizeof(Egamma));
            os.write((char*)&sf, (char)sizeof(sf));
            os.write((char*)&r0, (char)sizeof(r0));

            GOAT::raytracing::binWrite(lp, os);
        }

        void sphericLens::binRead(std::ifstream& is)
        {
            type = OBJECTSHAPE_SPHERIC_LENS;
            P.binRead(is);
            H.binRead(is);
            R.binRead(is);
            is.read((char*)&n, (char)sizeof(n));
            alpha.binRead(is);
            pul.binRead(is);
            por.binRead(is);
            for (int i = 0; i < 3; i++)
                e[i].binRead(is);
            is.read((char*)&Ealpha, (char)sizeof(Ealpha));
            is.read((char*)&Ebeta, (char)sizeof(Ebeta));
            is.read((char*)&Egamma, (char)sizeof(Egamma));
            is.read((char*)&sf, (char)sizeof(sf));
            is.read((char*)&r0, (char)sizeof(r0));

            GOAT::raytracing::binRead(lp, is);
        }
        void sphericLens::scale(double sf)
        {
            this->sf = sf;
            init();
            initQuad();
        }

        void sphericLens::initQuad()
        {
            double b, l, d;
            if (lp.left.curvature == concave)
            {
                l = sqrt(lp.left.R * lp.left.R - lp.radius * lp.radius);
                d = lp.left.R - l;
                b = lp.offset / 2.0 + d;
            }

            if (lp.left.curvature == convex)
                b = lp.left.R - lp.left.P[2];
            
            if (lp.left.curvature == flat)
                b = lp.offset / 2.0;
            
            pul = maths::Vector<double>(-lp.radius, -lp.radius, -b);


            if (lp.right.curvature == concave)
            {
                l = sqrt(lp.right.R * lp.right.R - lp.radius * lp.radius);
                d = lp.right.R - l;
                b = lp.offset / 2.0 + d;
            }

            if (lp.right.curvature == convex)
                b= lp.left.R - lp.right.P[2];

            por = maths::Vector<double>(lp.radius, lp.radius, b);
        }

        void sphericLens::setr0(double r0)
        {
            this->r0 = r0;
            initQuad();
        }

        bool sphericLens::next(const maths::Vector<double>& Ps, const maths::Vector<double>& K, maths::Vector<double>& pout)
        {
            double lambda;
            double lleft, lright;
            double D, l1, l2;
            
            GOAT::maths::Vector<double> n, k, p, ps;
            GOAT::maths::Vector<double> P1, P2;

            k = H * K;

            // left side
            ps = Ps - P;
            ps = H * ps;
            p = ps - lp.left.P;
            lleft = -1;
            if (lp.left.curvature == flat) 
                lleft = (k[2]==0) ? -1 : - p[2] / k[2];
            else
            {
                double pk = p * k;                
                D = pk * pk - abs2(p) + lp.left.R * lp.left.R;
                if (D < 0) lleft = -1; // no intersection with left side
                else
                {
                    D = sqrt(D);
                    l1 = -pk + D;
                    l2 = -pk - D;
                    P1 = p + l1 * k;
                    P2 = p + l2 * k;
                    if (lp.left.curvature == convex)
                        lleft = P1[2] < 0 ? l1 : l2;
                    else
                        lleft = P1[2] > 0 ? l1 : l2;
                }
            }


            p = ps - lp.right.P;            
            lright = -1;
            if (lp.right.curvature == flat)
               lright = (k[2]==0) ? -1 : - p[2] / k[2];
            
            else
            {
                double pk = p * k;
                D = pk * pk - abs2(p) + lp.right.R * lp.right.R;
                if (D < 0) lright= -1; // no intersection with left side
                else
                {
                    D = sqrt(D);
                    l1 = -pk + D;
                    l2 = -pk - D;
                    P1 = p + l1 * k;
                    P2 = p + l2 * k;
                    if (lp.right.curvature == concave)
                        lright = P1[2] < 0 ? l1 : l2;
                    else 
                        lright = P1[2] > 0 ? l1 : l2;
                }
            }

            // lateral surface
            
            p = ps;
            double llateral;
            double h1 = k[0] * p[0] + k[1] * p[1];
            double h2 = (k[0] * k[0] + k[1] * k[1]);
            D = 2 * k[0] * k[1] * p[0] * p[1] - k[0] * k[0] * p[1] * p[1] - k[1] * k[1] * p[0] * p[0] + lp.radius * lp.radius * h2;
           
            if (D <= 0) llateral = -1; // no intersection with lateral surface
            else
            {                
                D = sqrt(D);
                l1 = (-h1 + D) / h2;
                l2 = (-h1 - D) / h2;
                llateral = ((l2 < 1E-10) || (l1 < l2)) ? l1 : l2;
            }

            // finally everthing is prepared, let's have a look, which intersection point is nearer 
            // first, compare the side faces
            if (((lleft > 0) && (lleft < lright))|| (lright<=1E-10))
            {
                lambda = lleft;
                pout = Ps + lleft * K;
                if (lp.left.curvature == flat) currentnorm = -R * maths::ez;
                else
                {
                    currentnorm = R * (p + lleft * k - lp.left.P);
                    currentnorm = currentnorm / abs(currentnorm);
                }
            }
            else
            {               
                lambda = lright;
                pout = Ps + lambda * K;
                if (lp.right.curvature == flat) currentnorm = maths::ez;
                else
                {                    
                    currentnorm = R * (p + lambda * k  - lp.right.P);
                    currentnorm = currentnorm / abs(currentnorm);
                }
            }
           
            if (lambda < 0)
            {
                if (llateral < 0) { pout = Ps; return false; }
                else 
                {                     
                    pout = Ps + llateral * K; 
                    if (abs(Ps - P) > lp.offset / 2.0) { pout = Ps; return false; }                    
                    currentnorm = p + llateral * k;
                    currentnorm[2] = 0;
                    currentnorm = R * currentnorm / abs(currentnorm);                    
                    return true; 
                }
            }
            else
            {
                if ((llateral < lambda) && (llateral > 0))
                {
                    pout = p + llateral * k;
                    currentnorm = pout;
                    currentnorm[2] = 0;
                    currentnorm = R * currentnorm / abs(currentnorm);
                    if (fabs(pout[2]) < lp.offset/2.0) { pout = Ps + llateral * K; return true; }
                }

                pout = p + lambda * k; // check if the ray hits the lens
                if (pout[0] * pout[0] + pout[1] * pout[1] > lp.radius * lp.radius) { pout = Ps; return false; }
                pout = Ps + lambda * K;
                return true;                
            }
            pout = Ps;
            return false;
        }

        double sphericLens::volume()
        {
            double Vleft, Vlateral, Vright;

            // left part          
            double h, l;
            double R2 = lp.radius * lp.radius;

            l = sqrt(lp.left.R * lp.left.R - lp.radius * lp.radius);
            h = lp.left.R - l;
            if (lp.left.curvature == convex)
                Vleft = M_PI / 3.0 * h * h * (3.0 * lp.left.R - h);
            else
            {
                if (lp.left.curvature == flat)
                    Vleft = 0.0;
                else
                    Vleft = M_PI * R2 * h - M_PI / 3.0 * h * h * (3.0 * lp.left.R - h);
            }

            // right part          
            
            l = sqrt(lp.right.R * lp.right.R - lp.radius * lp.radius);
            h = lp.right.R - l;
            if (lp.right.curvature == convex)
                Vright = M_PI / 3.0 * h * h * (3.0 * lp.right.R - h);
            else
            {
                if (lp.right.curvature == flat)
                    Vright = 0.0;
                else
                    Vright = M_PI * R2 * h - M_PI / 3.0 * h * h * (3.0 * lp.right.R - h);
            }

            // lateral surface
            Vlateral = M_PI * R2 * lp.offset;
            return Vleft + Vlateral + Vright;
        }

        void sphericLens::setPos(maths::Vector<double> r)
        {
            P = r;
        }

        void sphericLens::setPos(double x, double y, double z)
        {
            P = maths::Vector < double>(x, y, z);
        }

        void sphericLens::init()
        {                        
            // Set the centers of the spherical boundary surfaces         
            // ---- left side ----
            double zl;
            // the radius of the lens must be smaller or equal to the smallest curvature radius
            if ((lp.radius > lp.left.R) && (lp.left.curvature!=flat)) lp.radius = lp.left.R;
            if ((lp.radius > lp.right.R) && (lp.right.curvature!=flat)) lp.radius = lp.right.R;
            if (lp.left.curvature == convex) zl = sqrt(lp.left.R * lp.left.R - lp.radius * lp.radius) - lp.offset / 2.0;
            if (lp.left.curvature == concave) zl = -(lp.left.R + lp.offset / 2.0);
            if (lp.left.curvature == flat) zl = -lp.offset / 2.0;

            lp.left.P = maths::Vector<double>(0, 0, zl);

            // ---- right side ----
            double zr;
            if (lp.right.curvature == convex) zr = lp.offset / 2.0 - sqrt(lp.right.R * lp.right.R - lp.radius * lp.radius);
            if (lp.right.curvature == concave) zr = lp.offset / 2.0 + lp.right.R;
            if (lp.right.curvature == flat) zr = lp.offset / 2.0;

            lp.right.P = maths::Vector<double>(0, 0, zr);
        }
    }
}