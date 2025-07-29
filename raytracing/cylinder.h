#pragma once
#include "objectshape.h"
namespace GOAT
{
	namespace raytracing
	{
		/**
		 * @brief This class provides a cylindrical object.
		 * Here, a cylindrical object is provided. The reference point is situated in the center of bottom face. The cylinder is described 
		 * by its height and radius. It is oriented in the z-direction
		 */
		class Cylinder :
			public ObjectShape
		{
		public:
			Cylinder();
			Cylinder(ObjectShape& os);
			Cylinder(Cylinder& c);
			Cylinder(
				const maths::Vector<double>& P,
				double r,
				double h,
				std::complex<double>  n,
				double r0 = 1.0,
				const maths::Matrix<std::complex<double> > alpha = maths::CUNITY,
				const maths::Vector<double>& Ex = maths::ex,
				const maths::Vector<double>& Ey = maths::ey,
				const maths::Vector<double>& Ez = maths::ez
			     );
			double radius() { return r; } 
			double height() { return h; }
			void setRadius(double r); ///<  
			void setHeight(double h); ///< set the height of the cylinder
			void setr0(double r0);
			void binWrite(std::ofstream& os);  
			void binRead(std::ifstream& is);
			maths::Vector<double> norm(const maths::Vector<double>& p);
			bool isInside(const maths::Vector<double>& p);
			double volume();  ///< Calculates the volume of the ellipsoid
			void initQuad(); ///< Sets the circumscribing cuboid (for use in inelastic calculations)
			maths::Vector<double> calcCoM() { return P; } ///< Calculates the center of mass, which is in the case of an ellipsoid simply its center, P.
			maths::Matrix<double> computeInertia();   ///< Calculates the inertia matrix of the ellipsoid
			void setPos(maths::Vector<double> r) { P = r; }; ///< sets reference point P 
			void setPos(double x, double y, double z) { P = maths::Vector<double>(x, y, z); };
			bool next(const maths::Vector<double>& Ps, const maths::Vector<double>& K,
				maths::Vector<double>& pout);
			Cylinder& operator =  (Cylinder& f);
			Cylinder& operator = (Cylinder f);
			void scale(double sf); ///< Set the scaling factor (the half axis will be multiplied by this factor)
		protected :
			double h = 1.0;
			double r = 1.0;

		};
	}
}