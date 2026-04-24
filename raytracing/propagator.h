#pragma once
#include <vector>
#include "detector.h"
namespace GOAT
{
	namespace raytracing
	{
		class Propagator : public DetectorPlane
		{
		public:
			Propagator(double wvl, maths::Vector<double> P, maths::Vector<double> e1, maths::Vector<double> e2, int n1, int n2);
			Propagator(double wvl, maths::Vector<double>P, maths::Vector<double> n, double d, int N);
			void addDetector(DetectorPlane* det);
			void addDetectorList(std::vector<DetectorPlane*> detList);
			void delDetector(DetectorPlane* det);
			void clearSources();
			void calc(bool clear = true);
			void setNumberOfThreads(int noThreads);
			size_t numberOfThreads() { return noThreads; }
			size_t numberOfSources() { return sources.size(); }
			std::vector<DetectorPlane*> getSources() { return sources; }
			double getWavelength() { return wvl; }
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