#pragma once
#include "detector.h"
namespace GOAT
{
	namespace raytracing
	{
		/**
		* @brief This class makes a Kirchhoff calculation
		*/
		class Kirchhoff : public DetectorPlane
		{
		  public : 
			  /**
			  * @brief constructor
			  * 
			  */
			  Kirchhoff(double wvl, maths::Vector<double> P, maths::Vector<double> e1, maths::Vector<double> e2, int n1, int n2);
			  void calc(DetectorPlane* det);

		private:
			double k;
			double wvl;
			maths::Vector<std::complex<double> > point(DetectorPlane* det, maths::Vector<double> P);

		};
	}
}