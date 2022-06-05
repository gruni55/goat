#include"raytrace.h"

/**
 * Example for elastic scattering by a spherical particle, illuminated by a plane wave.
 * (Experimental stage)
 */

int main(int argc, char** argv)
{
	Scene S; 
	S.setnS(1.0);			// set refractive index of surrounding medium to 1.0

	double Rsphere = 10.0;	// Radius of the sphere in µm

    LightSrcPlane LS (-300.0*ez, 3000, 1.0, 2.0*Rsphere,I*ey,LIGHTSRC_RAYTYPE_RAY);  // Setup of the incident plane wave
	LS.setk(ez);			// direction of the wave
	S.addLightSource(&LS);	// add light source to scene
	

	Ellipsoid obj(dzero, Vector<double>(Rsphere, Rsphere, Rsphere), 1.5); // define sphere with refractive index 1.5
	S.addObject(&obj);		// add object to scene
	S.setr0(1000.0);		// set radius of calculating sphere to 1000µm
	
	
	Raytrace_scattering R(S,360,180);	// initialise raytracer for elastic scattering
	R.setNumReflex(5);		// set number of internal reflections to 5
	R.trace();				// start raytracing 
	R.save("scatt10.dat",SAVE_ABS2);	// save results to file (squared absolute value of the electric field |E|²
	return 0;
}
