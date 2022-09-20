#include "raytrace.h"
/*****************************************************************//**
 * \file   layers.cpp
 * \brief  This program calculates the interference through a thin film layer as a function of the layer thickness.
 * 
 * \author Thomas
 * \date   January 2022
 *********************************************************************/
int main(int argc, char** argv)
{
	/* First, define the light source */
	GOAT::raytracing::LightSrcPlane LS(-200 * GOAT::maths::ex, 300, 1.0,50.0,I*GOAT::maths::ey); // plane wave with a diameter of 50µm and a wavelength 1µm
	LS.setk(GOAT::maths::ex);  // travels in positive x-direction
	
	/* Now, define the layer, described by a box (i.e. a cuboid) */
	double l = 5;
	GOAT::maths::Vector<double> d(l, 300, 300);
	GOAT::raytracing::Box Layer(GOAT::maths::dzero, d, 1.5); // define the layer with a thickness of l and an expansion in y- and z-direction  of 300µm
	double alpha = 30.0 / 180.0 * M_PI; // the layer is rotated by 30° around the y-axis
	Layer.setBeta(alpha);

	/* define the scene */
	GOAT::raytracing::Scene S;
	S.addLightSource(&LS);  // add the previously defined light source 
	S.setRaytype(LIGHTSRC_RAYTYPE_IRAY); // select the ray type
    S.addObject(&Layer);    // add the layer
	S.setr0(1000.0);        // define the diameter of the calculation sphere
	S.setnS(1.0);           // set the refractive index of the surrounding medium to 1
	
	int N = 1;
	/* Define  two detectors */
	GOAT::raytracing::DetectorPlane D(200.0 * GOAT::maths::ex, -GOAT::maths::ex,500, N);    //  one after the layer... 
	GOAT::raytracing::DetectorPlane D0(-180.0 * GOAT::maths::ex, -GOAT::maths::ex, 400, N); //  one before the layer 
	// Add the detectors to the scene
	S.addDetector(&D);     
	S.addDetector(&D0);

	GOAT::raytracing::Raytrace_pure RP(S); // Calculate with the detectors only
	RP.setNumReflex(2);  // Set the number of internal reflections to 2
	
	GOAT::maths::Vector<std::complex<double> > E1, E2;
	std::ofstream os("layers.dat"); 

	// loop over the different layer thicknesses
	for (double d = 1.0; d <= 3.0; d += 0.01)
	{
		((GOAT::raytracing::Box *)RP.S.Obj[0])->setD(GOAT::maths::Vector<double> (d, 300, 300)); // set the new thickness
		RP.S.cleanAllDetectors(); // Clean the detectors
		RP.trace();               // perform raytracing  

		// Sum up the intensities of all detectors
		E1 = GOAT::maths::czero;
		E2 = GOAT::maths::czero;
		for (int ix = 0; ix < N; ix++)
			for (int iy = 0; iy < N; iy++)
			{
				E1 = E1 + RP.S.Det[0]->D[ix][iy];
				E2 = E2 + RP.S.Det[1]->D[ix][iy];
			}
		os << d << "\t" << abs2(E1) << "\t" << abs2(E2) << std::endl;
		std::cout << d << "\t" << abs2(E1)  << "\t" << abs2(E2) << std::endl;
	}
	os.close();


	return 0;
}
