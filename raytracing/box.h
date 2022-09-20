#pragma once
#include "objectshape.h"
#include <iostream>


/**  
* @brief class which represents a box (cuboid). It is derived by class ObjectShape    
* This class is mainly used for the octree calculation. 
*/


namespace GOAT
{
	namespace raytracing
	{
		class Box :
			public ObjectShape
		{

		private:
			bool isOctree = false; ///< used for Octree calculation

		public:
			Box();
			Box(const ObjectShape& F);
			Box(const Box& B); ///< copy constructor
			/**
			 * Main constructor.
			 *
			 * \param P Position of the box (mid point)
			 * \param d Vector with the edge lengths in x-, y- and z-direction
			 * \param n (Complex) refractive index
			 * \param r0 Radius of the calculation sphere (needed by the next routine)
			 * \param alpha Polarisability matrix
			 * \param Ex Local coordinate system (first axis, usually x)
			 * \param Ey Local coordinate system (second axis, usually y)
			 * \param Ez Local coordinate system (third axis, usually z)
			 */
			Box(const maths::Vector<double>& P,
				const maths::Vector<double>& d,
				std::complex<double> n,
				double r0 = 1.0,
				const maths::Matrix<std::complex<double> > alpha = maths::CUNITY,
				const maths::Vector<double> Ex = maths::ex,
				const maths::Vector<double> Ey = maths::ey,
				const maths::Vector<double> Ez = maths::ez);
			~Box(); ///< destructor
			void binWrite(std::ofstream& os);  ///< writes object to a binary file
			void binRead(std::ifstream& is);   ///< reads object from a binary file
			void scale(double sf); ///< sets scaling factor
			double distance(const maths::Vector<double>& p, const maths::Vector<double>& k); ///< calculates the distance between a ray, described by the position vector p and its directional vector k 
			bool next(const maths::Vector<double>& p, const maths::Vector<double>& k, maths::Vector<double>& pout); ///< calculates the next crossing point between a ray represented by the position 
			maths::Vector<double> norm(const maths::Vector<double>& p); ///< normal vector at the position p
			bool isInside(const maths::Vector<double>& Ps); ///< returns true, if Ps is inside the box
			double volume() { return abs2(d); } ///< returns the volume of the box
			double getSize() { return abs(bounds[1] - bounds[0]); }  ///< returns a vector with the side lengths as its components
			void initQuad();  ///< calculates the circumferent cuboid
			void setr0(double r0); ///< sets the radius of the calculation sphere
			void calcDiag() ///< calculates the longest diagonal inside the box
			{
				diag[0] = bounds[1] - bounds[0];
				diag[1] = maths::Vector<double>(-diag[0][0], diag[0][1], diag[0][2]);
				diag[2] = maths::Vector<double>(-diag[0][0], -diag[0][1], diag[0][2]);
			}

			void setPos(maths::Vector<double> r) ///< sets the position P and the corresponding bounds
			{
				P = r;
				bounds[0] = P - d / 2.0;
				bounds[1] = P + d / 2.0;
			}
			void setPos(double x, double y, double z) { setPos(maths::Vector<double>(x, y, z)); } ///< sets the position P and the corresponding bounds
			void setD(maths::Vector<double> D) ///< sets the extensions in x-, y- and z-direction
			{
				d = D;
				bounds[0] = P - d / 2.0;
				bounds[1] = P + d / 2.0;
			}

			maths::Vector<double> bounds[2];  ///< positions of the two opposite corners  
			maths::Vector<double> d; ///< extensions of the box in x-, y- and z-direction
			maths::Vector<double> diag[3]; ///< diagonal 
			maths::Vector<double> calcCoM() { return maths::dzero; }
			void setOctree(bool isOctree) { this->isOctree = isOctree; } ///< box belongs to an octree calculation
		};

		std::ostream& operator<< (std::ostream& os, Box B); ///< output operator for the Box class
	}
}
