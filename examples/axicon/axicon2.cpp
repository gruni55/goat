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
	int numRays = 10;                       // number of rays (per direction)  
	double wvl = 1.0;							 // wavelength
	double LSsize = 20000.0;                     // size of the light source in µm 
    
	GOAT::raytracing::LightSrcPlane_mc LS (LSPos, numRays, wvl, LSsize);
    LS.setk(GOAT::maths::ex);  
	LS.setPol(GOAT::maths::Vector<std::complex<double> >(1.0, 0.0, 0.0));

	// Object (axicon) definitions
	GOAT::maths::Vector<double> objPos = GOAT::maths::dzero;           // Position of the axicon
	std::complex<double> nObj = 1.5;         // Refractive index 
	GOAT::raytracing::surface surf(objPos, nObj);
	surf.importBinSTL("axicon_10.stl");           
	surf.setActive(false);

	GOAT::maths::Vector<double> boxPos(247.7, 0, 0);
	GOAT::maths::Vector<double> boxDim(100, 100, 10);
	GOAT::raytracing::Box box(boxPos, boxDim, 1.0);
	box.setActive(true);


	// Scene definition 
	GOAT::raytracing::Scene S;
	S.setnS(1.0);                            // refractive index surrounding medium 
	S.setr0(r0);                           // radius of the calculation sphere (in µm) 
	S.addLightSource(&LS);                   // add light source to scene
	S.addObject(&surf);                      // add object (axicon) to scene
	S.addObject(&box);

        // Set the refractive index functions
       std::vector<std::function<std::complex<double>(double) > > nList;
       nList.push_back(GOAT::raytracing::n_BK7);
       nList.push_back(GOAT::raytracing::n_Air);
	   nList.push_back(GOAT::raytracing::n_Air);
       
       GOAT::raytracing::pulseCalculation pc(S);
       pc.setPulseWidth (100);
       pc.setRefractiveIndexFunctions(nList);
       pc.setSpatialResolution(1.0);
       pc.setReferenceTime(0.0);
   
   

   double alpha = 10.0 / 180.0 * M_PI;
   double h = 20.0;
   double D = 20.0;
   double beta = asin(1.5 * sin(alpha));
   double d = h * sin(alpha);
   double l = h / tan(beta - alpha);
   double f = D + d + l;
   double t = (f + 250) * GOAT::raytracing::C_LIGHT_MU_FS;
   std::cout << "The electric field will be stored in " << argv[1] << std::endl;
   
   pc.field(t);
   GOAT::raytracing::saveabsE(pc.trafo.SAres, argv[1],1);
	return 0;
}
