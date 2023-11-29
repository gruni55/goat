#pragma once
#include "vector.h"
#include "matrix.h"
namespace GOAT
{
	namespace raytracing
	{
		/**
		@brief This class describes a triangle, represented by its corner points. It is intented for internal purposes.
		 The triangle class is mainly used in class surface.
		 */
		class triangle
		{
		public:
			double u, v; ///< auxiliary variables for internal use
			maths::Vector<double> P[3];  ///< corner points of the triangle
			maths::Vector<double> n;     ///< normal of the triangle 
			maths::Vector<double> f[3];  ///< vectors which represents the sides of the triangle. For details refer to function calcSideVectors() 

			/**
			this function returns the surface area of the triangle
			*/
			double area()
			{
				maths::Vector<double> h1, h2;
				h1 = P[1] - P[0];
				h2 = P[2] - P[0];
				return abs(h1 % h2) / 2.0;
			}

			void calcSideVectors()  ///< This function calculates the three vectors, which represents represents the side vectors of the triangle: f0=P1-P0, f1=P2-P1, f2=P0-P1
			{
				f[0] = P[1] - P[0];
				f[1] = P[2] - P[1];
				f[2] = P[0] - P[2];
			}
			double distance(maths::Vector<double> p, maths::Vector<double> k);
			/**
			@brief Calculates the intersection point between the triangle and a straight line represented by a point r and the direction vector k
			@param r reference point on the straight line
			@param k direction vector of the straight line
			@param p intersection point (return parameter)
			@param eps accuracy (optional)
			@return 0: no intersection found  1: intersection point found
			*/
			int calcIntersectionPoint(maths::Vector<double> r, maths::Vector<double> k, maths::Vector<double>& p, double  eps = 1E-10);
			int calcIntersectionPoint(maths::Vector<double> r, maths::Vector<double> k, double& t, maths::Vector<double>& p, double  eps = 1E-10);
			void setnorm(void);   ///< calculates the surface normal n with help of the corner points  
			void setnorm(maths::Vector<double> n) { this->n = n; } ///< sets the surface normal n to the given value
			maths::Vector<double>  getnorm(void); ///< returns the surface normal n
			triangle();
			/**
			Constructor of the class triangle
			@param P1,P2,P3 : corner points of the triangle
			*/
			triangle(maths::Vector<double> P1, maths::Vector<double> P2, maths::Vector<double> P3);
			/**
			Constructor of the class triangle, the triangle is represented by a reference point and the three distance vectors
			between the reference point and the three corners of the triangle.
			@param P0: reference point
			@param ip1, ip2, ip3 : Vectors from the reference point P0 to the corner points
			*/
			triangle(maths::Vector<double> ip1, maths::Vector<double> ip2, maths::Vector<double> ip3, maths::Vector<double> P);
			triangle(const triangle& d); ///< Copy constructor
			/**
			* This operator returns the i-th corner point of the triangle
			* @param i
			* @return Vector of the i-th corner
			*/
			maths::Vector<double>& operator[](int i);
			/**
			* This operator returns the i-th corner point of the triangle
			* @param i
			* @return Vector of the i-th corner
			*/
			const maths::Vector<double>& operator[](int i) const; ///< Operator, gives back the Position of i-th corner
			triangle& operator=(const triangle& dr); ///< Assignment operator
			void binWrite(std::ofstream& os);
			void binRead(std::ifstream& is);
			~triangle(); ///< Destructor
			friend class surface;
		};

		std::ostream& operator << (std::ostream& os, const triangle& dr);
		/** @name Operators on triangle
		 *
		 */
		 ///@{
		triangle operator + (const triangle& dr, const maths::Vector<double>& v); ///< add vector v to the corners of the triangle dr
		triangle operator - (const triangle& dr, const maths::Vector<double>& v); ///< subtract vector v from the corners of the triangle dr
		triangle operator + (const maths::Vector<double>& v, const triangle& dr); ///< add vector v to the corners of the triangle dr
		triangle operator - (const maths::Vector<double>& v, const triangle& dr); ///< subtract vector v from the corners of the triangle dr
		triangle operator * (const maths::Matrix<double>& M, const triangle& dr); ///< multiply matrix  M to every corner of triangle dr
		triangle operator * (const triangle& dr, double a); ///< multiply the corners of the triangle dr with the scalar a
		triangle operator * (double a, const triangle& dr); ///< multiply the corners of the triangle dr with the scalar a
		triangle operator / (const triangle& dr, double a); ///< divide the corners of the triangle dr by the scalar a
		///@}
	}
}