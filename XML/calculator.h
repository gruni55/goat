#pragma once
#include "raytrace.h"
#include "pulsecalculation_rt.h"
#include "xml.h"

/*****************************************************************//**
 * \file   calculator.h
 * \brief  This class processes calculations defined in a calculation job
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
			calculationJob job;
			raytracing::Scene S;
		};
	}
}