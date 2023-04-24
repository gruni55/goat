#include "pulsecalculation.h"
#include "refractive_index_functions.h"
/*****************************************************************//**
 * \file   axicon.cpp
 * \brief Example for the intensity distribution after an axicon with 10° angle, height 1mm, 
 * and 30mm diameter  
 * 
 * Light source: gaussian beam with waist 8mm
 * \author Thomas
 * \date   December 2021
 *********************************************************************/


int main(int argc, char** argv)
{
	// Light source definitions
        double r0=50000.0; 
	GOAT::maths::Vector<double> LSPos = -1000.0 * GOAT::maths::ex;      // Position of the light source
	int numRays = 1000000;                       // number of rays (per direction)  
	double wvl = 1.0;                        // wavelength
	double LSsize = 20000.0;                  // size of the light source in µm 

//	GOAT::raytracing::LightSrcGauss_mc LS(LSPos, numRays, wvl,waist,focusPos,LSsize); // initialize gaussian light source
        GOAT::raytracing::LightSrcPlane_mc LS (LSPos, numRays, wvl, LSsize);
        LS.setk(GOAT::maths::ex);  
	LS.setPol(GOAT::maths::Vector<std::complex<double> >(1.0, 0.0, 0.0));

	// Object (axicon) definitions
	GOAT::maths::Vector<double> objPos = GOAT::maths::dzero;           // Position of the axicon
	std::complex<double> nObj = 1.5;         // Refractive index 
	GOAT::raytracing::surface surf(objPos, nObj);
	surf.importBinSTL("axicon_10.stl");           


	// Scene definition 
	GOAT::raytracing::Scene S;
	S.setnS(1.0);                            // refractive index surrounding medium 
	S.setr0(r0);                           // radius of the calculation sphere (in µm) 
	S.addLightSource(&LS);                   // add light source to scene
	S.addObject(&surf);                      // add object (axicon) to scene

        // Set the refractive index functions
       std::vector<std::function<std::complex<double>(double) > > nList;
       nList.push_back(GOAT::raytracing::n_BK7);
       nList.push_back(GOAT::raytracing::n_Air);
       
       GOAT::raytracing::pulseCalculation pc(S);
       pc.setPulseWidth (100);
       pc.setRefractiveIndexFunctions(nList);
       pc.setSpatialResolution(1.0);
       pc.setReferenceTime(0.0);
       pc.field(3500);
   GOAT::raytracing::saveabsE(pc.trafo.SAres, "fieldabs.dat");

	
	return 0;
}
