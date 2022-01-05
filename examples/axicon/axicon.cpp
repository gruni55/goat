#include "raytrace.h"
/*****************************************************************//**
 * \file   axicon.cpp
 * \brief Example for the intensity distribution after an axicon with 10� angle, height 1mm, 
 * and 30mm diameter  
 * 
 * Light source: gaussian beam with waist 8mm
 * \author Thomas
 * \date   December 2021
 *********************************************************************/

int main(int argc, char** argv)
{
	// Light source definitions
	Vector<double> LSPos = -10000.0 * ez;      // Position of the light source
	Vector<double> focusPos(0, 0, 70000.0);  // Position of the focus
	int numRays = 801;                       // number of rays (per direction)  
	double wvl = 1.0;                        // wavelength
	double waist = 8000.0;                   // waist diameter in �m 
	double LSsize = 15000.0;                  // size of the light source in �m 

	LightSrcGauss LS(LSPos, numRays, wvl,waist,focusPos,LSsize); // initialize gaussian light source
	LS.setPol(Vector<std::complex<double> >(1.0, 0.0, 0.0));
	// Object (axicon) definitions
	Vector<double> objPos = dzero;           // Position of the axicon
	std::complex<double> nObj = 1.5;         // Refractive index 
	surface surf(objPos, nObj);
	surf.createsurface("axicon_2000.srf");           

	// Detector definition
	Vector<double> detPos(0, 0, 35000.0);    // Position of the detector
	Vector<double> detNorm(0, 0, -1);        // Surface normal of the detector
	double detSize = 20000.0;                // Size of the detector
	int detGridsize = 2001;                  // Number of Pixels (per direction) of the detector
	DetectorPlane dp( detPos, -detNorm, detSize, detGridsize);
	DetectorPlane dp2(-5000.0*ez, -detNorm, detSize, detGridsize);

	// Scene definition 
	Scene S;
	S.setnS(1.0);                            // refractive index surrounding medium 
	S.setr0(1E+8);                           // radius of the calculation sphere (in �m) 
	S.addLightSource(&LS);                   // add light source to scene
	S.addObject(&surf);                      // add object (axicon) to scene
	S.addDetector(&dp);                      // add detector to scene
	S.addDetector(&dp2);                      // add detector to scene
	
	 Raytrace_pure rp(S);                    // define raytracer without interfaces, just for the detectors
	 rp.setNumReflex(0);                     // consider no reflexions
	 // rp.trace();                             // start raytracing
	
	 dp.saveabs("axicon_intensity.dat");     // save detector contents as absolute value of the electric field
	 dp2.saveabs("axicon_intensity_before.dat");     // save detector contents as absolute value of the electric field

	 S.cleanAllDetectors();
	

	 // calculate the paths of the rays
	 LS.setAnzRays(21);                       // reduce the number of rays
	 LS.setD(15000);
	 Raytrace_Path rpt(S);                   // chose another raytracer  
	 rpt.setNumReflex(0);                    // once again, without reflexions
	 rpt.trace("axicon_rays.dat");           // make the raytracing  
	return 0;
}
