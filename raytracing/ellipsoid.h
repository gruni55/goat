#ifndef ELLIPSOID_H
#define ELLIPSOID_H

#include <complex>
#include "vector.h"
#include "matrix.h"
#include "objectshape.h"
#include <iostream>
#include <fstream>
/**  
  @brief This class represents an ellipsoid, defined by its half axis a, b and c according to the formula:
   \f$\frac{x^2}{a^2}+\frac{y^2}{b^2}+\frac{z^2}{c^2}=1\f$
*/
namespace GOAT
{
    namespace raytracing
    {
        class Ellipsoid : public ObjectShape
        {
        public:
            // Ellipsoid &operator = (const Ellipsoid& e);
            ~Ellipsoid() {};
            Ellipsoid();
            Ellipsoid(const ObjectShape&);
            /**
             * @brief Copy Constructor
             *
             */
            Ellipsoid(const Ellipsoid& E);
            /**
             * @brief Main constructor for class Ellipsoid.
             *
             * \param P Position of the ellipsoid (center of the ellipsoid)
             * \param r The components of this vector represents the length of the half axis along the x-, y- and z-axis
             * \param n (Complex) refractive index
             * \param r0 Radius of the calculation sphere (optional, needed by next()). If the Scene class is used, the radius of the Scene class will be used
             * for the calculation
             * \param alpha Polarisability (for inelastic calculations only)
             * Directions of the local coordinate system
             * \param Ex
             * \param Ey
             * \param Ez
             */
            Ellipsoid(
                const maths::Vector<double>& P,
                const maths::Vector<double>& r,
                std::complex<double>  n,
                double r0 = 1.0,
                const maths::Matrix<std::complex<double> > alpha = maths::CUNITY,
                const maths::Vector<double>& Ex = maths::ex,
                const maths::Vector<double>& Ey = maths::ey,
                const maths::Vector<double>& Ez = maths::ez
            );
            void scale(double sf); ///< Set the scaling factor (the half axis will be multiplied by this factor)
            /**
             * @brief Calculates the next intersection point with a ray, which is represented by a reference point and the direction vector.
             *
             * \param Ps Reference point of the ray
             * \param K Direction of the ray
             * \param pout Intersection point (zero if not found)
             * \return true, if an intersection point was found, otherwise false
             */
            bool next(const maths::Vector<double>& Ps, const maths::Vector<double>& K,
                maths::Vector<double>& pout);
            Ellipsoid& operator =  (Ellipsoid& f);
            Ellipsoid& operator = (Ellipsoid f);
            /**
             * @brief calculates the surface normal of the ellipsoid at a certain point (p must be on the surface)
             * @param p Point at which the normal will be calculated
             * @return surface normal
             */
            maths::Vector<double> norm(const maths::Vector<double>& p);
            /**
             * @brief Checks if a certain position is inside the object.
             *
             * \param p Position to check
             * \return true, if the position is inside otherwise false.
             */
            bool isInside(const maths::Vector<double>& p);
            double volume();  ///< Calculates the volume of the ellipsoid
            void initQuad(); ///< Sets the circumscribing cuboid (for use in inelastic calculations)
            friend std::ostream& operator << (std::ostream& os, Ellipsoid E);
            void binWrite(std::ofstream& os); ///< Writes information about the ellipsoid in a binary file.
            void binRead(std::ifstream& os); ///< Reads informatopn about an ellipsoid from a binary file. 
            double a() { return r[0]; } ///< Returns the length of the first half axis 
            double b() { return r[1]; }///< Returns the length of the second half axis 
            double c() { return r[2]; }///< Returns the length of the third half axis 
            maths::Vector<double> getr() { return r; }  ///< Returns a vector, where its components are the lengths of the three half axis
            void setPos(maths::Vector<double> r) { P = r; } ///< Set the position of the ellipsoid (center) to r
            void setPos(double x, double y, double z) { setPos(maths::Vector<double>(x, y, z)); } ///< Sets the position of the ellipsoid(center) to a vector, represented by its components x,y and z.
            void setr(maths::Vector<double>& r); ///< Sets the lengths of the three half axis, represented by the components of r
            void setr(double a, double b, double c); ///< Sets the lengths of the three half axis, represented by a, b and c
            /**
             * Changes the length of the first half axis. If VConst is set true, the lengths of the other half axis are changed by the factor \f$\sqrt{\frac{a_{old}}{a_{new}}}\f$, otherwise the other half axis remain the same
             */
            void seta(double a, bool VConst = false);
            /**
             * Changes the length of the second half axis. If VConst is set true, the lengths of the other half axis are changed by the factor \f$\sqrt{\frac{b_{old}}{b_{new}}}\f$ to keep the volume constant, otherwise the other half axis remain the same
             */
            void setb(double b, bool VConst = false);
            /**
             * Changes the length of the third half axis. If VConst is set true, the lengths of the other half axis are changed by the factor \f$\sqrt{\frac{c_{old}}{c_{new}}}\f$ to keep the volume constant, otherwise the other half axis remain the same
             */
            void setc(double c, bool VConst = false);
            void setr0(double r0);
            maths::Vector<double> calcCoM() { return P; } ///< Calculates the center of mass, which is in the case of an ellipsoid simply its center, P.
            maths::Matrix<double> computeInertia();   ///< Calculates the inertia matrix of the ellipsoid

            // protected :
            maths::Vector<double> r;   ///< GOAT::maths::Vector, where its components represent the lengths of the corresponding half axis
            maths::Vector<double> r_2;      ///< GOAT::maths::Vector, where its components represent the square of lengths of the corresponding half axis (for internal use only)
            maths::Vector<double> P2;      ///< Square of the position vector P: P2=(P[0]^2,P[1]^2,P[2]^2) (for internal use only)
        };
    }
}
#endif
