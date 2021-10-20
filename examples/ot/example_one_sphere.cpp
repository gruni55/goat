#include "raytrace.h" // include file for the raytracing process

int main()
{
	double wvl = 1.064;
	double n_medium = 1.33;
	double n_particle = 1.59;
	double wvl_medium = wvl / n_medium;
	double radius = 5.0 * wvl_medium;
	
	
	// define object (ellipsoid with radius "radius" and refractive index n_particle positioned at the origin of the coordinate system
	Ellipsoid E(dzero, Vector<double>(radius, radius, radius), n_particle); 

	// define light source (gaussian beam starting position at -200ez, wavelength wvl, waist 1.0 (not interesting here) and fokus position at the origin)
	LightSrcGauss LS(-200.0 * ez, 500, wvl, 1.0, dzero);

	LS.P0 = 1.0; // Power of the light source is set to 1W

    Scene S; // define scene
	S.setr0(1000.0);        // set radius of the calculation sphere
	S.addObject(&E);        // add object to scene  
	S.addLightSource(&LS);  // add light source to scene
	S.setnS(n_medium);      // change the refractive index of the surrounding medium

	LS.setNA(1.02); // set light source's numerical aperture

	Raytrace_OT rot(S);     // define raytracer and initialize it with scene 

	std::ofstream os("example_one_sphere.dat"); // open file for writing
	for (double z = -10 * wvl_medium; z <= 10 * wvl_medium; z += 20 * wvl_medium / 299.0) // loop over position
	{
		E.P[2] = z;  // change position of the object (z-component)
		rot.trace(); // start raytracing 
		std::cout << z << "\t" << rot.F[0] << "\t" << rot.L[0] << std::endl;
		os << z << "\t" << rot.F[0] << "\t" << rot.L[0] << std::endl; // put result (force and torque) into file 
	}
	os.close(); // close file
}
