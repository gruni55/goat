#include "raytrace_usp.h"
#include "fft.h"
#include <complex>
#include <functional>

std::complex<double> none(double wvl)
{
	return 1.0;
}

int main (int argc, char **argv)
{
	int nrays = 100000;
	/* ---  Let's first define the Scene ----
	    with one lightsource and one box, where the field is stored 
	*/
	
	double wvl = 1.0E-6; // peak wavelength
	// use light source with arbitrary ray distribution 
	GOAT::raytracing::LightSrcPlane_mc LS(-GOAT::maths::ex * 200E-6, nrays, 1E-6, 100E-6);
	LS.setk(GOAT::maths::ex);
	GOAT::raytracing::Box Box(GOAT::maths::dzero, GOAT::maths::Vector<double>(100E-6, 100E-6, 100E-6), 1.0);
	GOAT::raytracing::Scene S;
	S.setr0(500.0E-6); 
	S.addLightSource(&LS);
	S.addObject(&Box);

	double dt = 1E-12;  // width of the pulse (in seconds)
	double sigma = dt / (2.0 * M_LN2);
	

	double dwvl= wvl * wvl * M_LN2 / (M_PI * M_PI * GOAT::raytracing::C_LIGHT * dt);// spectral FWHM
	std::cout << "dwvl=" << dwvl*1E+9 << "nm" << std::endl;
	std::vector<std::function<std::complex<double>(double) > > nList;
	
	std::cout << "Do raytracing...";
	GOAT::raytracing::Raytrace_usp rt(S,1000);
	rt.trace();
	std::cout << "done." << std::endl;

	GOAT::raytracing::TrafoParms parms;	
	parms.nI = 10;
	parms.lstart = wvl - dwvl;
	parms.lstop = wvl + dwvl;
    parms.nList.push_back(none);
	parms.nList.push_back(none);
	
	GOAT::raytracing::Trafo<GOAT::maths::Vector<double> > T(parms);
	std::string fname = "c:\\users\\weigetz9\\data\\testsa.dat";
	save(rt.SA[0], fname);


	
	T.calc(rt.SA, 1.0E-9);
 	std::cout << "Calculating done." << std::endl;
	fname="c:\\users\\weigetz9\\data\\test.dat";
	GOAT::raytracing::saveabsE(T.SAres, fname);
	
 
  return 0;
}
