#include"raytrace.h"

/**
 * Example for elastic scattering by a spherical particle, illuminated by a plane wave.
 * (Experimental stage)
 */

int main(int argc, char** argv)
{
	Scene S; 
	S.setnS(1.0);			// set refractive index of surrounding medium to 1.0

	double Rsphere = 50.0;	// Radius of the sphere in µm

    LightSrcPlane LS (-300.0*ez, 12000, 1.0, 2.0*Rsphere,I*ey,LIGHTSRC_RAYTYPE_RAY);  // Setup of the incident plane wave
	LS.setk(ez);			// direction of the wave
	S.addLightSource(&LS);	// add light source to scene
	DetectorPlane D1 = DetectorPlane(0.0 * ez, -ez, 2.0 * Rsphere, 500);
	DetectorPlane D2 = DetectorPlane(10.0 * ez, -ez, 2.0 * Rsphere, 500);
	
	S.addDetector(&D1);
	S.addDetector(&D2);

	Ellipsoid obj(dzero, Vector<double>(Rsphere, Rsphere, Rsphere), 1.5); // define sphere with refractive index 1.5
	S.addObject(&obj);		// add object to scene
	S.setr0(1000.0);		// set radius of calculating sphere to 1000µm
	// S.setRaytype(LIGHTSRC_RAYTYPE_IRAY);
	
	/*Raytrace_Path rp(S);
	rp.setNumReflex(0);
	rp.trace("test.dat");
	*/
	 Raytrace_scattering R(S, 1800, 180);	// initialise raytracer for elastic scattering
	R.setNumReflex(10);		// set number of internal reflections to 5
	R.trace();				// start raytracing 
	R.save("scatt50.dat",SAVE_ABS2);	// save results to file (squared absolute value of the electric field |E|²
	D1.save("D1.dat");
	D2.save("D2.dat");
	//S.Det[0]->savePhase("D1.dat", SAVE_PHASE_Y);
	//S.Det[1]->savePhase("D2.dat", SAVE_PHASE_Y);
	
	return 0;
}
