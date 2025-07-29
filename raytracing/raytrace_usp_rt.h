#pragma once
#include "raytrace.h"
#include "superarray.h"

namespace GOAT
{
	/** 
	* @brief Raytracer used for ultrashort pulse calculation with raytracing only
	*/
	namespace raytracing
	{
		class Raytrace_usp_rt :
			public Raytrace
		{
		public :
			Raytrace_usp_rt();
			Raytrace_usp_rt(const Scene& S, INDEX_TYPE n);
			~Raytrace_usp_rt();
			void clear(); ///< Clears the SuperArray grid (SA) 
			void init(); ///< Initialises the SuperArray grid (SA) according to the objects hold in Scene S
			void trace(double omega, std::complex<double> weight);  ///< Perform the raytracing process
			void storeData();
			void setRefractiveIndexFunctions(std::vector<std::function<std::complex<double>(double) > > nList); ///< sets the list of functions, which describe the wavelength dependend refractive
			void traceEnterObject(); ///< Function internally called when ray enters an object
			void traceLeaveObject(); ///< Function internally called when ray leaves an object
			
			std::vector<SuperArray <maths::Vector<std::complex<double> > > > SA;
			INDEX_TYPE n = 1;
			double omega, k0;
			std::complex<double> weight;
			std::vector<std::function<std::complex<double>(double) > > nList;
			int iR = 0; ///< Number of reflections to consider


		};
	}
}