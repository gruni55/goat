#include "propagator.h"
#include <omp.h>
namespace GOAT
{
	namespace raytracing
	{
		Propagator::Propagator(double wvl, maths::Vector<double> P, maths::Vector<double> e1, maths::Vector<double> e2, int n1, int n2)
			: DetectorPlane(P, e1, e2, n1, n2), wvl(wvl)
		{
			k = 2.0 * M_PI / wvl;			
			this->wvl = wvl;
		//	type = DETECTOR_PROPAGATOR;
		}

		Propagator::Propagator(double wvl, maths::Vector<double>P, maths::Vector<double> n, double d, int N)
			: DetectorPlane(P, n, d, N), wvl(wvl)
		{
			k = 2.0 * M_PI / wvl;
			this->wvl = wvl;
		//	type = DETECTOR_PROPAGATOR;
		}

		void Propagator::addDetector(DetectorPlane* det)
		{
			sources.push_back(det);
		}

		void Propagator::addDetectorList(std::vector<DetectorPlane*> detList)
		{
			sources.insert(sources.end(), detList.begin(), detList.end());
		}

		void Propagator::delDetector(DetectorPlane* det)
		{
			sources.erase(std::remove(sources.begin(), sources.end(), det), sources.end());
		}

		void Propagator::clearSources()
		{
			sources.clear();
		}

		void Propagator::setNumberOfThreads(int noThreads)
		{
			int maxThreads = omp_get_max_threads();
			if (noThreads > maxThreads) this->noThreads = maxThreads;
			else this->noThreads = noThreads;
		}

		void Propagator::calc(std::vector<DetectorPlane*> detList)
		{
			for (auto det : detList)
				calcOne(det, false);
		}

		void Propagator::calc(bool clear)
		{
			if (clear) clean();
			for (auto det : sources)
				calcOne(det, false);
		}
	}
}