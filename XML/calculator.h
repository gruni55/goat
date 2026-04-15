#pragma once
#include "raytrace.h"
#include "pulsecalculation_rt.h"
#include "xml.h"

/*****************************************************************//**
 * \file   calculator.h
 * \brief  This class processes calculations defined in a calculation job
 * This class processes calculations defined in a calculation job, which is part of the XML input. 
 * It uses the scene defined in the XML input and executes the calculation defined in the job.  
 *The class can be extended to support different types of calculations, e.g. pure raytracing calculations, Kirchhoff calculations, etc.
 * 
 * \author Thomas Weigel
 * \date   June 2024
 *********************************************************************/
namespace GOAT
{
	namespace XML
	{ 
		class Calculator
		{
		public:
			Calculator(raytracing::Scene& S, calculationJob& job);
			void exec();

			virtual ~Calculator() {}
		private:
			void pulseCalculation();
			void pureRaytraceCalculation();
			calculationJob job;
			raytracing::Scene S;
		};
	}
}