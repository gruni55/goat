/***************************************************************************
                          plane.h  -  description                              
                             -------------------                                         
    begin                : Sat Feb 19 2000                                           
    copyright            : (C) 2000 by Thomas Weigel                         
    email                : weigel@lat.ruhr-uni-bochum.de                                     
 ***************************************************************************/



#ifndef EBENE_H
#define EBENE_H


/**
  *@author Thomas Weigel
  */

#include<iostream>
#include<fstream>
#include "vector.h"
#include "matrix.h"
/**
 * @brief This class is used for the iray class. 
 * This class is intended for internal use only. It defines a plane, defined by a central position and two directional vectors. 
 */
namespace GOAT
{
	namespace raytracing
	{

		class Plane {
		public:
			Plane();
			/**
			 * @brief Constructor
			 *
			 * \param P Position (center) of the plane
			 * \param e1 first direction vector
			 * \param e2 second direction vector
			 */
			Plane(const maths::Vector<double>& P, const maths::Vector<double>& e1, const maths::Vector<double>& e2);
			/**
			 * @brief find intersection point with a sphere
			 * \param O center of the sphere
			 * \param r radius of the sphere
			 */
			void intersectSphere(maths::Vector<double> O, double r);
			void norm(); ///< normal of the plane
			void Normalenform();
			double distance(maths::Vector<double> R); ///< distance between the point R and the plane 
			void toString(char* S);
			void binWrite(std::ofstream& os); ///< writes plane into the binary file stream os
			void binRead(std::ifstream& is);  ///< reads plane from binary file stream is
			void rotate(double dr, double dtheta, double dphi);
			friend std::ostream& operator << (std::ostream& os, Plane E)
			{
				os << " P=" << E.P << std::endl;
				os << "n=" << E.n << std::endl;
				os << "e1=" << E.e1 << std::endl;
				os << "e2=" << E.e2 << std::endl;
				return os;
			}
			~Plane();
			maths::Vector<double> P, e1, e2, n;
		};

		const Plane Exz = Plane(maths::Vector<double>(0.0, 0.0, 0.0), maths::ex, maths::ez);
		const Plane Exy = Plane(maths::Vector<double>(0.0, 0.0, 0.0), maths::ex, maths::ey);
	}
}
#endif
