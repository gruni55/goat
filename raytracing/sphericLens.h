#include "lens.h"
#include "objectshape.h"
#pragma once
namespace GOAT
{
    namespace raytracing
    {
        /**
        * @brief This class provides a lens with spherical surfaces 
        * The lens consists of two side surfaces and one lateral surface. The properties of the lens is described in the
        * structure lp
        * Each side surface can either be convex, concave 
        * or flat. The center of the lens is situate in the center between the to side surfaces. The distance between this 
        * surfaces is the offset parameter in lp. 
        * Left side (convex):
        * \image html lens_convex_left.svg
        * 
        * Left side (concave):
        * \image html lens_concave_left.svg
        * 
        *  Left side (convex):
        * \image html lens_convex_right.svg
        * 
        * Right side (concave):
        * \image html lens_concave_right.svg
        */
        class sphericLens :
            public ObjectShape
        {
        public:
            sphericLens(const maths::Vector<double>& P,
                std::complex<double>  n,
                lensParms lp,
                GOAT::maths::Matrix<std::complex<double> >  alpha=maths::CUNITY,
                const maths::Vector<double>& Ex = maths::ex,
                const maths::Vector<double>& Ey = maths::ey,
                const maths::Vector<double>& Ez = maths::ez,
                const int type = OBJECTSHAPE_SPHERIC_LENS
            );

            void binWrite(std::ofstream& os);
            void binRead(std::ifstream& is);
            void scale(double sf);
            void initQuad();
            void setr0(double r0);
            /**
            * @brief returns surface normal at P
            * The arbitrary calculation of the normal at P is complex and time consuming, therefore P is used here as
            * a dummy variable. It is assumed, that before calling norm, the next() routine was used. The routine returns
            * the value of the normal at the last hit point calculated by next()
            */
            maths::Vector<double> norm(const maths::Vector<double>& P) { return currentnorm; } 
            /**
            * @brief returns the next intersection point of a ray with this object
            * The ray is desribed by the starting point of the ray p and the directional vector k (which is assumed to 
            * be normalized). Here, the current normal to the surface at the intersection point is calculated (which can be 
            * read by the norm() function.
            * The function returns true if there is an intersection point
            */
            bool next(const maths::Vector<double>& p, const maths::Vector<double>& k,
                maths::Vector<double>& pout);
            bool isInside(const maths::Vector<double>& p) { return false; } ///< not yet implemented
            double volume(); ///< calculates the volume of the lens

            // ------- Set methods -----------
            void setPos(maths::Vector<double> r); ///< set the position of the lens
            void setPos(double x, double y, double z); ///< set the position of the lens
            void setParms(lensParms lp) {this->lp=lp; init(); } ///< set the lens parameter

            maths::Vector<double> calcCoM() { return maths::dzero; }; ///< not yet implemented
            lensParms getParms() { return lp; }  ///< returns the parameters, which describe the lens

            
        private:
            void init();
            lensParms lp;
            maths::Vector<double> currentnorm;
        };
    }
}
