#pragma once
#include "superarray.h"
#include "raytrace_usp_rt.h"
#include <vector>
#include "fft.h"

namespace GOAT
{
	namespace raytracing
	{
		class pulseCalculation_rt
		{
		public:
			pulseCalculation_rt(Scene S);
			void setDefaults();
			TrafoParms trafoparms;
			void calcTrafoParms();
			void setCenterWavelength(double wvl);
			void setBandwidth(double dWvl);
			void setRepetitionRate(double rep);
			void setSpectralRanges(int nI);
			void setSpatialResolution(double dx);
			void setRefractiveIndexFunctions(std::vector<std::function<std::complex<double>(double)>> nList);
			void setNumReflex(int numReflex);
			void setPulseWidth(double dt);
			void field(double t); ///< This function calculates the fields at time. Keep in mind, that it works only if the class has the list with the refractive index functions
			/*
				* Returns the number of cells per direction, nn along each axis of the calculation space.
				* The spatial resolution is then 2*r0/nn. This number is changed by calling #setSpatialResolution
				*/
			INDEX_TYPE getNumCellsPerDirection() { return nn; }
			Raytrace_usp_rt rt;
			Scene S;
			INDEX_TYPE nn = 0;
			double dWvl = 0.02;  ///< spectral width of the light (default 20nm)
			int numReflex = INEL_MAX_NREFLEX;
		};
	}
}