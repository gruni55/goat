#include "ellipsoid.h"
#include "matrix.h"

namespace GOAT
{
    namespace raytracing
    {
        Ellipsoid::Ellipsoid()
            : ObjectShape()
        {
            type = OBJECTSHAPE_ELLIPSOID;
        }

        Ellipsoid::Ellipsoid(const ObjectShape& F) : ObjectShape(F)
        {
            type = OBJECTSHAPE_ELLIPSOID;
        }

        Ellipsoid::Ellipsoid(const Ellipsoid& E) :ObjectShape(E)
        {
            r = E.r;
            P = E.P;
            r_2 = E.r_2;
            type = OBJECTSHAPE_ELLIPSOID;
        }


        Ellipsoid::Ellipsoid(
            const GOAT::maths::Vector<double>& P,
            const GOAT::maths::Vector<double>& r,
            std::complex<double>  n,
            double r0,
            const GOAT::maths::Matrix<std::complex<double> > alpha,
            const GOAT::maths::Vector<double>& Ex,
            const GOAT::maths::Vector<double>& Ey,
            const GOAT::maths::Vector<double>& Ez)
            : ObjectShape(P, n, alpha, Ex, Ey, Ez, OBJECTSHAPE_ELLIPSOID)
        {
            GOAT::maths::Vector<double> h[8];
            trafo(Ex, Ey, Ez, H, R);
            this->r0 = r0;
            this->r = r;
            this->P = P;
            type = OBJECTSHAPE_ELLIPSOID;
            P2 = GOAT::maths::Vector<double>(P[0] * P[0], P[1] * P[1], P[2] * P[2]);
            r_2 = GOAT::maths::Vector<double>(1.0 / (r[0] * r[0]), 1.0 / (r[1] * r[1]), 1.0 / (r[2] * r[2]));
            initQuad();
        }

        void Ellipsoid::scale(double sf)
        {
            double sfold = this->sf;
            this->sf = sf;
            r = r * sf / sfold;
            r_2 = GOAT::maths::Vector<double>(1.0 / (r[0] * r[0]), 1.0 / (r[1] * r[1]), 1.0 / (r[2] * r[2]));
            initQuad();
        }

        bool Ellipsoid::next(const GOAT::maths::Vector<double>& Ps, const GOAT::maths::Vector<double>& K,
            GOAT::maths::Vector<double>& pout)
        {            
            double A, B, C;
            double l1, l2, l;
            GOAT::maths::Vector<double> n, k, p = Ps - P;
            p = H * p;
            k = H * K;
            pout = Ps;
            n = norm(Ps);
            /* cout.precision (10);
             cout << "P=" << P << "   Ps=" << Ps << "  k=" << K << "   n=" << n << "  |p|=" << abs(p)  << "  phi=" << acos (fabs(K*n))/M_PI*180.0 << endl;*/

            GOAT::maths::Vector<double> k2 = GOAT::maths::Vector<double>(k[0] * k[0], k[1] * k[1], k[2] * k[2]);
            GOAT::maths::Vector<double> p2 = GOAT::maths::Vector<double>(p[0] * p[0], p[1] * p[1], p[2] * p[2]);
            A = k2 * r_2;
            B = 2.0 * emult(p, k) * r_2;
            C = p2 * r_2 - 1.0;

            l1 = B * B - 4.0 * A * C;
            if (l1 <= 0.0) { /*cout << "l1=" << l1 << "  kein Schnittpunkt" << endl;*/ return false; }
            l2 = (-B + sqrt(l1)) / (2.0 * A);
            l1 = (-B - sqrt(l1)) / (2.0 * A);
            //  cout << "l1=" << l1 << "   l2=" << l2;
            if (l1 / r0 <= 1E-10) l = l2; else l = l1;
            //    cout << "    l=" << l << endl;
            if (l <= 1E-6 * r0) { /*cout << "NICHT genommen" << endl;*/   return false; }
            pout = Ps + l * K;
            //  cout << "GENOMMEN: pout=" << pout << endl;
            return true;
        }

        GOAT::maths::Vector<double> Ellipsoid::norm(const GOAT::maths::Vector<double>& p)
        {
            GOAT::maths::Vector<double> ps, n;
            ps = H * (p - P);

            n = GOAT::maths::Vector<double>(2.0 * ps[0] / (r[0] * r[0]), 2.0 * ps[1] / (r[1] * r[1]), 2.0 * ps[2] / (r[2] * r[2]));
            n = R * n / abs(n);
            return n;
        }

        /* Ellipsoid & Ellipsoid::operator = (const Ellipsoid& e)
         {
          if (this == &e) return *this;
          H=e.H;
          R=e.R;
          P=e.P;
          r=e.r;
          return *this;
         }*/

        bool Ellipsoid::isInside(const GOAT::maths::Vector<double>& p)
        {
            GOAT::maths::Vector<double> h;
            h = ediv(p, r);
            return abs2(h) < 1.0;
        }

        /*Ellipsoid* Ellipsoid::copy()
        {
          return new Ellipsoid (P,r,n,r0,alpha,e[0],e[1],e[2]);
        }

        void Ellipsoid::copy(Ellipsoid *E)
        {
          E->H=H;
          E->R=R;
          E->P=P;
          E->P2=P2;
          E->alpha=alpha;
          E->Ealpha=Ealpha;
          E->Ebeta=Ebeta;
          E->Egamma=Egamma;
          for (int i=0; i<3; i++)
          E->e[i]=e[i];
          E->n=n;
          E->r=r;
          E->r_2=r_2;
          E->r0=r0;
          E->type=type;
          E->initQuad();
         }

         */
        Ellipsoid& Ellipsoid::operator = (Ellipsoid& f)
        {
            if (this == &f) return *this;
            H = f.H;
            P = f.P;
            R = f.R;
            alpha = f.alpha;
            for (int i = 0; i < 3; i++)
                e[i] = f.e[i];
            n = f.n;
            por = f.por;
            pul = f.pul;
            r = f.r;
            r_2 = f.r_2;
            r0 = f.r0;
            P2 = f.P2;
            type = f.type;
            Ealpha = f.Ealpha;
            Ebeta = f.Ebeta;
            Egamma = f.Egamma;
            return *this;
        }

        Ellipsoid& Ellipsoid::operator = (Ellipsoid f)
        {
            if (this == &f) return *this;
            H = f.H;
            P = f.P;
            R = f.R;
            alpha = f.alpha;
            for (int i = 0; i < 3; i++)
                e[i] = f.e[i];
            n = f.n;
            por = f.por;
            pul = f.pul;
            r = f.r;
            r0 = f.r0;
            r_2 = f.r_2;
            P2 = f.P2;
            type = f.type;
            Ealpha = f.Ealpha;
            Ebeta = f.Ebeta;
            Egamma = f.Egamma;
            return *this;
        }

        void Ellipsoid::setr0(double r0)
        {
            GOAT::maths::Vector<double> h[8];
            /* r=r/this->r0*r0;
             r_2=GOAT::maths::Vector<double> (1/(r[0]*r[0]),1/(r[1]*r[1]),1/(r[2]*r[2]));
             P=P/this->r0*r0;*/
             // por=GOAT::maths::Vector<double> (r[0],r[1],r[2]);
             // pul=-por;
             // por=R*por+P;
             // pul=R*pul+P;
            initQuad();
            this->r0 = r0;
        }

        void Ellipsoid::initQuad()
        {
            GOAT::maths::Vector<double> h[8];

            for (int k = 0; k <= 1; k++)
                for (int l = 0; l <= 1; l++)
                    for (int m = 0; m <= 1; m++)
                    {
                        k == 0 ? h[k * 4 + l * 2 + m][0] = -r[0] : h[k * 4 + l * 2 + m][0] = r[0];
                        l == 0 ? h[k * 4 + l * 2 + m][1] = -r[1] : h[k * 4 + l * 2 + m][1] = r[1];
                        m == 0 ? h[k * 4 + l * 2 + m][2] = -r[2] : h[k * 4 + l * 2 + m][2] = r[2];
                    }

            pul = GOAT::maths::dzero;
            por = GOAT::maths::dzero;
            for (int i = 0; i < 8; i++)
            {
                h[i] = R * h[i];
                if (h[i][0] < pul[0]) pul[0] = h[i][0];
                if (h[i][1] < pul[1]) pul[1] = h[i][1];
                if (h[i][2] < pul[2]) pul[2] = h[i][2];
                if (h[i][0] > por[0]) por[0] = h[i][0];
                if (h[i][1] > por[1]) por[1] = h[i][1];
                if (h[i][2] > por[2]) por[2] = h[i][2];
            }
            pul = pul + P;
            por = por + P;
        }


        /*void Ellipsoid::setMatrix (double alpha, double beta, double gamma)
        {
         // Der Einschluss wird gedreht: alpha,beta,gamma sind die neuen Eulerwinkel !
         double ca,cb,cg;
         double sa,sb,sg;
         GOAT::maths::Vector<double> h[8];

         ca=cos(alpha); cb=cos(beta); cg=cos(gamma);
         sa=sin(alpha); sb=sin(beta); sg=sin(gamma);

         H(0,0)= cg*cb;   H(0,1)= sa*sb*cg+ca*sg; H(0,2)=-ca*sb*cg+sa*sg;
         H(1,0)=-cb*sg;   H(1,1)=-sa*sb*sg+ca*cg; H(1,2)= ca*sb*sg+cg*sa;
         H(2,0)=    sb;   H(2,1)=-sa*cb;          H(2,2)= ca*cb;

         R(0,0)=cb*cg;             R(0,1)=-cb*sg;          R(0,2)= sb;
         R(1,0)=ca*sg+sa*sb*cg;    R(1,1)=ca*cg-sa*sb*sg;  R(1,2)=-sa*cb;
         R(2,0)=sa*sg-ca*sb*cg;    R(2,1)=sa*cg+ca*sb*sg;  R(2,2)=ca*cb;


         for (int k=0; k<=1; k++)
           for (int l=0; l<=1; l++)
             for (int m=0; m<=1; m++)
             {
                k==0 ? h[k*4+l*2+m][0]=-r[0] : h[k*4+l*2+m][0]=r[0];
                l==0 ? h[k*4+l*2+m][1]=-r[1] : h[k*4+l*2+m][1]=r[1];
                m==0 ? h[k*4+l*2+m][2]=-r[2] : h[k*4+l*2+m][2]=r[2];
             }

         pul=zero;
         por=zero;

         for (int i=0; i<8; i++)
         {
          h[i]=R*h[i];

          if (h[i][0]<pul[0]) pul[0]=h[i][0];
          if (h[i][1]<pul[1]) pul[1]=h[i][1];
          if (h[i][2]<pul[2]) pul[2]=h[i][2];
          if (h[i][0]>por[0]) por[0]=h[i][0];
          if (h[i][1]>por[1]) por[1]=h[i][1];
          if (h[i][2]>por[2]) por[2]=h[i][2];
         }
         por=por+P;
         pul=pul+P;
         e[0]=H*ex;
         e[1]=H*ey;
         e[2]=H*ez;
         Ealpha=alpha;
         Ebeta=beta;
         Egamma=gamma;
        } */

        void Ellipsoid::binWrite(std::ofstream& os)
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
            r.binWrite(os);
            r_2.binWrite(os);
            P2.binWrite(os);
        }

        void Ellipsoid::binRead(std::ifstream& is)
        {
            type = OBJECTSHAPE_ELLIPSOID;
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
            r.binRead(is);
            r_2.binRead(is);
            P2.binRead(is);
        }

        std::ostream& operator << (std::ostream& os, Ellipsoid E)
        {
            os << "%r=" << E.r << std::endl;
            os << "%P=" << E.P << std::endl;
            os << "%Ealpha=" << E.Ealpha / M_PI * 180.0
                << "°   Ebeta=" << E.Ebeta / M_PI * 180.0
                << "°   Egamma=" << E.Egamma / M_PI * 180.0 << "°" << std::endl;
            return os;
        }

        double Ellipsoid::volume()
        {
            return 4.0 / 3.0 * M_PI * r[0] * r[1] * r[2];
        }

        void Ellipsoid::setr(GOAT::maths::Vector<double>& r)
        {
            this->r = r;
            r_2 = GOAT::maths::Vector<double>(1.0 / (r[0] * r[0]), 1.0 / (r[1] * r[1]), 1.0 / (r[2] * r[2]));
        }
        void Ellipsoid::setr(double a, double b, double c)
        {
            r = GOAT::maths::Vector<double>(a, b, c);
            r_2 = GOAT::maths::Vector<double>(1.0 / (r[0] * r[0]), 1.0 / (r[1] * r[1]), 1.0 / (r[2] * r[2]));
        }

        void Ellipsoid::seta(double a, bool VConst)
        {
            GOAT::maths::Vector<double> h;
            a = a * r0;
            if (VConst)
            {
                double t = sqrt(r[0] / a);
                h = GOAT::maths::Vector<double>(a, r[1] * t, r[2] * t);
            }
            else h = GOAT::maths::Vector<double>(a, r[1], r[2]);
            setr(h);
        }

        void Ellipsoid::setb(double b, bool VConst)
        {
            GOAT::maths::Vector<double> h;
            b = b * r0;
            if (VConst)
            {
                double t = sqrt(r[1] / b);
                h = GOAT::maths::Vector<double>(r[0] * t, b, r[2] * t);
            }
            else h = GOAT::maths::Vector<double>(r[0], b, r[2]);
            setr(h);
        }

        void Ellipsoid::setc(double c, bool VConst)
        {
            GOAT::maths::Vector<double> h;
            c = c * r0;
            if (VConst)
            {
                double t = sqrt(r[2] / c);
                h = GOAT::maths::Vector<double>(r[0] * t, r[1] * t, c);
            }
            else h = GOAT::maths::Vector<double>(r[0], r[1], c);
            setr(h);
        }

        GOAT::maths::Matrix<double> Ellipsoid::computeInertia()
        {
            GOAT::maths::Matrix<double> I;
            I(0, 0) = 1.0 / 5.0 * (r[1] * r[1] + r[2] * r[2]);
            I(1, 1) = 1.0 / 5.0 * (r[0] * r[0] + r[2] * r[2]);
            I(2, 2) = 1.0 / 5.0 * (r[0] * r[0] + r[1] * r[1]);
            return I / 1E-12; // 1E-12, da r hier in µm angegeben wird, ich will aber I in m²
        }
    }
}