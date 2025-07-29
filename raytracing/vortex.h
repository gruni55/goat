#pragma once
#include "objectshape.h"
namespace GOAT
{
	namespace raytracing
	{
		/**
		 * @brief This class provides a vortex phase plate object.
		 * Here, a vortex phase plate object is provided. The reference point is situated in the center of bottom face. 
		 */
		class VortexPlate :
			public ObjectShape
		{
		public:
			VortexPlate();
			VortexPlate(ObjectShape& os);
			VortexPlate(VortexPlate& c);
			VortexPlate(
				const maths::Vector<double>& P,
				double r,
				double h,
				double dh,
				int m, 
				std::complex<double>  n,
				double r0 = 1.0,
				const maths::Matrix<std::complex<double> > alpha = maths::CUNITY,
				const maths::Vector<double>& Ex = maths::ex,
				const maths::Vector<double>& Ey = maths::ey,
				const maths::Vector<double>& Ez = maths::ez
			     );
			double radius() { return r; } 
			double height() { return h; }
			void setRadius(double r); 
			void setHeight(double h);
			void setr0(double r0); ///< set the radius of the calculation space
			void setm(int m); ///< set the topological charge
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
			VortexPlate& operator =  (VortexPlate& f);
			VortexPlate& operator = (VortexPlate f);
			void scale(double sf); ///< Set the scaling factor (the half axis will be multiplied by this factor)
			int order() {return m;} ///< returns the order (~topological charge)
			double vortexHeight() {return dh; } ///< returns the height of the vortex structure 
		protected :
			double h = 1.0; ///< height of the vortex (without vortex structure)
			double r = 1.0; ///< radius of the vortex
			int m = 1; ///< topological charge of the vortex
			double dh = 1.0; ///< height of the vortex structure

		};
	}
}