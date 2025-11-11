#pragma once
#include "detector.h"
namespace GOAT
{
	namespace raytracing
	{
		/**
		* @brief This class makes a Kirchhoff calculation
		* This class is directly connect with a Detector. The calculation itself works as follows: 
		* At first, a normal raytracing step is performed to calculate the electric field at a detector. This detector is used as a 
		* source field for the next step, where the field at a given area is calculated with help of the Kirhhoff integral
		*/
		class Kirchhoff : public DetectorPlane
		{
		  public : 
			  /**
			  * @brief constructor
			  * \param wvl: Wavelength used for the caluclation of the electric field
			  * \param P: Position of the Kirchhoff plane
			  * \param e1: first direction vector of one edge of the detector
			  * \param e2: second direction vector of the edge of the detector perpendicular to the first one
			  * \param n1: number of cells in e1-direction
			  * \param n2: number of cells in e2-direction
			  */
			  Kirchhoff(double wvl, maths::Vector<double> P, maths::Vector<double> e1, maths::Vector<double> e2, int n1, int n2);

			  /**
			  *  @brief This method make the calculation
			  * With this method, the calculation of the Kirchhoff-integral will be performed for one detector. 
			  * \param det: a pointer to the detector, which acts as the source
			  */
			  void calc(DetectorPlane* det, bool clear=true);

			  /**
			  * @brief Do the Kirchhoff calculation with more than one detector as source
			  */
			  void calc(std::vector<DetectorPlane*> detList);

		private:
			double k;
			double wvl;
			maths::Vector<std::complex<double> > point(DetectorPlane* det, maths::Vector<double> P);

		};
	}
}