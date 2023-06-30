/***************************************************************************
                          form.cpp  -  description
                             -------------------
    begin                : Wed Oct 24 2001
    copyright            : (C) 2001 by Thomas Weigel
    email                : weigel@lat.ruhr-uni-bochum.de
 ***************************************************************************/

#include "objectshape.h"
#include "misc.h"
#include "matrix.h" 
#ifdef WITH_OCTREE
#include "octree.h"
#endif
namespace GOAT
{
    namespace raytracing
    {
        ObjectShape::ObjectShape() {
            r0 = 1.0;
            H = maths::unity();
            R = maths::unity();
            sf = 1.0;
            Ealpha = 0;
            Ebeta = 0;
            Egamma = 0;
            Active = true;
        }

        ObjectShape::ObjectShape(const maths::Vector<double>& P,
            std::complex<double>  n,
            maths::Matrix<std::complex<double> >  alpha,
            const maths::Vector<double>& Ex,
            const maths::Vector<double>& Ey,
            const maths::Vector<double>& Ez,
            const int type)
        {
            Ealpha = 0.0;
            Ebeta = 0.0;
            Egamma = 0.0;

            Active = true;
            sf = 1.0;
            r0 = 1.0;
            this->type = type;
            this->P = P;
            this->alpha = alpha;
            this->n = n;
            ninel = n;
            e[0] = Ex;
            e[1] = Ey;
            e[2] = Ez;
            R(0, 0) = Ex[0]; R(0, 1) = Ex[1]; R(0, 2) = Ex[2];
            R(1, 0) = Ey[0]; R(1, 1) = Ey[1]; R(1, 2) = Ey[2];
            R(2, 0) = Ez[0]; R(2, 1) = Ez[1]; R(2, 2) = Ez[2];

            H(0, 0) = Ex[0]; H(0, 1) = Ey[0]; H(0, 2) = Ez[0];
            H(1, 0) = Ex[1]; H(1, 1) = Ey[1]; H(1, 2) = Ez[1];
            H(2, 0) = Ex[2]; H(2, 1) = Ey[2]; H(2, 2) = Ez[2];
        }

        ObjectShape::ObjectShape(const ObjectShape& F)
        {
            Active = true;
            sf = 1.0;
            type = F.type;
            P = F.P;
            alpha = F.alpha;
            Ealpha = F.Ealpha;
            Ebeta = F.Ebeta;
            Egamma = F.Egamma;
            H = F.H;
            R = F.R;
            por = F.por;
            pul = F.pul;
            r0 = F.r0;
            type = F.type;
            n = F.n;
            ninel = F.ninel;
            for (int i = 0; i < 3; i++)
                e[i] = F.e[i];

            // setMatrix (Ealpha,Ebeta,Egamma);
        }


        /*ObjectShape & ObjectShape::operator = (ObjectShape &f)
        {
         if (this == &f) return *this;
         H=f.H;
         P=f.P;
         R=f.R;
         alpha=f.alpha;
         e=f.e;
         n=f.n;
         por=f.por;
         pul=f.pul;
         //r=f.r;
         r0=f.r0;
         type=f.type;
         return *this;
        } */

        void ObjectShape::setMatrix(maths::Matrix<double> H)
        {
            this->H = H;
            R = invert(H);
            // R=invert(H);
            e[0] = H * maths::ex;
            e[1] = H * maths::ey;
            e[2] = H * maths::ez;
            initQuad();
        }

        void ObjectShape::rotate(maths::Vector<double> A, double phi)
        {
            R = rotMatrix(A, phi);
            H = invert(R);
            e[0] = H * maths::ex;
            e[1] = H * maths::ey;
            e[2] = H * maths::ez;
            initQuad();
        }

        void ObjectShape::setMatrix(double alpha, double beta, double gamma)
        {
            maths::Matrix <double> Dx, Dy, Dz;
            Dx = rotMatrix(maths::ex, alpha);
            Dy = rotMatrix(maths::ey, beta);
            Dz = rotMatrix(maths::ez, gamma);

            R = Dx * Dy * Dz;
            H = invert(R);

            e[0] = H * maths::ex;
            e[1] = H * maths::ey;
            e[2] = H * maths::ez;

            /*
             double ca,cb,cg;
             double sa,sb,sg;
             ca=cos(alpha); cb=cos(beta); cg=cos(gamma);
             sa=sin(alpha); sb=sin(beta); sg=sin(gamma);

             R(0,0)= cg*cb;   R(0,1)= sa*sb*cg+ca*sg; R(0,2)=-ca*sb*cg+sa*sg;
             R(1,0)=-cb*sg;   R(1,1)=-sa*sb*sg+ca*cg; R(1,2)= ca*sb*sg+cg*sa;
             R(2,0)=    sb;   R(2,1)=-sa*cb;          R(2,2)= ca*cb;

             H(0,0)=cb*cg;             H(0,1)=-cb*sg;          H(0,2)= sb;
             H(1,0)=ca*sg+sa*sb*cg;    H(1,1)=ca*cg-sa*sb*sg;  H(1,2)=-sa*cb;
             H(2,0)=sa*sg-ca*sb*cg;    H(2,1)=sa*cg+ca*sb*sg;  H(2,2)=ca*cb;

             // R=invert(H);
             e[0]=H*ex;
             e[1]=H*ey;
             e[2]=H*ez;*/
            Ealpha = alpha;
            Ebeta = beta;
            Egamma = gamma;
            initQuad();
        }

        /*
        void ObjectShape::initQuad()
        {
         initInc(this);
        }*/


        /*void ObjectShape::setr0(double r0)
        {
         setR0(this,r0);
        }*/

        maths::Matrix<double> computeInertia(ObjectShape* F)
        {
            switch (F->type)
            {
            case OBJECTSHAPE_ELLIPSOID: return ((Ellipsoid*)F)->computeInertia();
            case OBJECTSHAPE_SURFACE: return ((surface*)F)->computeInertia();
            }
        }

        void ObjectShape::setCenter2CoM()
        {
            maths::Vector<double> CoM = calcCoM();
            setCenter(CoM);
        }

        void ObjectShape::setCenter(maths::Vector<double> P)
        {
            switch (type)
            {
            case OBJECTSHAPE_SURFACE: ((surface*)this)->setCenter(P); break;
            default: this->P = P;
            }
        }

        bool ObjectShape::isOutsideWorld()
        {
            bool result = (pul[0] < -r0) || (por[0] > r0) ||
                          (pul[1] < -r0) || (por[1] > r0) ||
                          (pul[2] < -r0) || (por[2] > r0);
            return result;
        }

        bool intersectionTest(ObjectShape& A, ObjectShape& B)
        {
         bool result = (A.pul[0] <= B.por[0]) && (A.por[0] >= B.pul[0]) &&
                       (A.pul[1] <= B.por[1]) && (A.por[1] >= B.pul[1]) &&
                       (A.pul[2] <= B.por[2]) && (A.por[2] >= B.pul[2]);
         return result;             
        }       
    }
}