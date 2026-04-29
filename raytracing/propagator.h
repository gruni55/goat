#pragma once
#include <vector>
#include "detector.h"
namespace GOAT
{
	namespace raytracing
	{
		/**
		* @brief This class is the base class for all propagators. 
		* It inherits from DetectorPlane, so it can be used as a detector, but it also has a list of detectors as sources. 
		*  
		*/
		class Propagator : public DetectorPlane
		{
		public:
			/**
			* constructor 
			* \param wvl wavelength
			* \param P position of the propagator plane (center)
			* \param e1 first direction vector of the plane (length of the plane in this direction is determined by the length of e1)
			* \param e2 second direction vector of the plane (length of the plane in this direction is determined by the length of e2)
			* \param n1 number of points in the first direction
			* \param n2 number of points in the second direction
			*/
			Propagator(double wvl, maths::Vector<double> P, maths::Vector<double> e1, maths::Vector<double> e2, int n1, int n2);
			/**
			* constructor for square propagator, defined by the center position P, the normal n, the width d and the number of points in one direction N (so the array is N x N)
			* /param wvl wavelength
			* /param P position of the propagator plane (center)
			* /param n normal vector of the plane
			* /param d width of the plane in both directions
			* /param N number of points in one direction (so the array is N x N)
			*
			*/
			Propagator(double wvl, maths::Vector<double>P, maths::Vector<double> n, double d, int N);
			/**
			* Adds one detector to the list of sources. 
			*/
			void addDetector(DetectorPlane* det);
			/**
			* Adds a list of detectors to the list of sources.
			*/
			void addDetectorList(std::vector<DetectorPlane*> detList);
			/**
			* Deletes one detector from the list of sources.
			*/
			void delDetector(DetectorPlane* det);
			/**
			* Clears the list of sources, i.e. all detectors are removed from the list of sources.
			*/
			void clearSources();
			/**
			* Makes a calculation. clear determines if the content of the propagator should be cleared before the calculation. 
			* If clear is false, the result of the calculation will be added to the existing content.
			*/
			void calc(bool clear = true);
			/**
			* @brief Sets the number of threads for the calculation. 
			* If noThreads is larger than the maximum number of threads available, the maximum number of threads will be used.
			*/
			void setNumberOfThreads(int noThreads);
			/**
			* @brief Returns the number of threads used for the calculation.
			* This is the maximal number of threads that will be used for the calculation. It can be set with setNumberOfThreads().
			* 
			*/
			size_t numberOfThreads() { return noThreads; }
			/**
			* @brief Returns the number of sources, i.e. the number of detectors in the list of sources.
			*/
			size_t numberOfSources() { return sources.size(); }
			/**
			* @brief Returns the list of sources, i.e. the list of detectors in the list of sources.
			*/
			std::vector<DetectorPlane*> getSources() { return sources; }
			/**
			* @brief Returns the wavelength of the propagator.
			*/
			double getWavelength() { return wvl; }
			/**
			* @brief Sets the wavelength of the propagator.
			*/
			void setWavelength(double wvl) { this->wvl = wvl; }
			/**
			* @brief Do the Kirchhoff calculation with more than one detector as source
			*/
			void calc(std::vector<DetectorPlane*> detList);


		protected:
			/**
			 *  @brief This method make the calculation
			 * With this method, the calculation of the Kirchhoff-integral will be performed for one detector.
			 * \param det: a pointer to the detector, which acts as the source
			 */
			virtual void calcOne(DetectorPlane* det, bool clear = true) = 0;
			double k;
			double wvl;

			std::vector<DetectorPlane*> sources;
			int noThreads = 8;
		};
	}
}