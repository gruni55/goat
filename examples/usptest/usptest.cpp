#include "raytrace_usp.h"
#include "fft.h"
#include <complex>
#include <functional>
#include <sstream>

std::complex<double> none(double wvl)
{
	return 1.0;
}

int main (int argc, char **argv)
{
	int nrays = 1000;
	/* ---  Let's first define the Scene ----
	    with one lightsource and one box, where the field is stored 
	*/
	
	double wvl = 1.0E-6; // peak wavelength
	
	// use light source with arbitrary ray distribution 
	GOAT::maths::Vector<double> P(-200E-6, 0, 0); // Position of the light source
	GOAT::raytracing::LightSrcPlane LS(P, nrays, 1E-6, 50E-6); // define Plane wave at P with nrays rays wavelength 1E-6 and size 50E-6
	GOAT::maths::Vector<double> k(1, 0, 0); // direction of the wave
	k = k / abs(k); 
	LS.setk(k);
	
	// Region of interest
	GOAT::raytracing::Box Box(GOAT::maths::dzero, GOAT::maths::Vector<double>(2E-6, 2E-6, 2E-6), 1.0);
	GOAT::raytracing::Scene S;
	
	S.setr0(1000.0E-6); 
	S.addLightSource(&LS);
	S.addObject(&Box);
	S.setnS(1.0);

	double dt = 1E-12;  // width of the pulse (in seconds)
	double sigma = dt / (2.0 * M_LN2);
	

	double dwvl= wvl * wvl * M_LN2 / (M_PI * M_PI * GOAT::raytracing::C_LIGHT * dt);// spectral FWHM
	std::cout << "dwvl=" << dwvl*1E+9 << "nm" << std::endl;
	std::vector<std::function<std::complex<double>(double) > > nList;
    
	std::cout << "Do raytracing...";
	GOAT::raytracing::Raytrace_usp rt(S,1000);
	rt.setNumReflex(0);
	rt.trace();
	std::cout << "done." << std::endl;
	
	GOAT::raytracing::TrafoParms parms;	
	parms.nI = 2;
	parms.nR = 1;
	parms.nS = 10;
	parms.lstart = wvl - dwvl;
	parms.lstop = wvl + dwvl;
    parms.nList.push_back(none);
	parms.nList.push_back(none);
	
	GOAT::raytracing::Trafo<GOAT::maths::Vector<double> > T(parms);
	std::string fname = "H:\\testsa.dat";
	// save(rt.SA[0], fname);

    
	dt = 1E-13;
	std::string numStr;
	std::stringstream sstr;
	double t;

	for (int i = 0; i < 20; i++)
	{
		t = i * dt;
		T.calc(rt.SA, t);
		sstr  << "h:\\data\\test" << i << ".dat";
		fname = sstr.str();		
		sstr.str("");
		std::cout << "  current filename:" << fname <<  "   t=" << t << std::endl;
		GOAT::raytracing::saveabsE(T.SAres, fname);
	}
 
  return 0;
}
