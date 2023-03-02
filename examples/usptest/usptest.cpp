#include "raytrace_usp.h"
#include "fft.h"
#include <complex>
#include <functional>
#include <sstream>

std::complex<double> none(double wvl)
{
	return 1.0;
}

/**
* Wavelength dependent refractive index for glass (BK7) taken from refractiveindex.info
*/

std::complex<double> nGlass_BK7(double wvl)
{
	double wvl2 = wvl * wvl;
	std::complex<double> n = sqrt(1.03961212*wvl2/(wvl2- 0.00600069867)+ 0.231792344*wvl2/(wvl2- 0.0200179144)+ 1.01046945*wvl2/(wvl2- 103.560653)+1.0);
	return n;
}

std::complex<double> nGlass(double wvl)
{
	return 1.5075;
}

int main (int argc, char **argv)
{
	int nrays = 1000;
	double r0 = 1E+6;
	int nCells = 1000000;
	/* ---  Let's first define the Scene ----
	    with one lightsource and one box, where the field is stored 
	*/
	
	double wvl = 1.0; // peak wavelength
	
	// use light source with arbitrary ray distribution 
	GOAT::maths::Vector<double> P(-0.8E+6, 0, 0); // Position of the light source
	GOAT::raytracing::LightSrcPlane LS(P, nrays, 1.0, 50); // define Plane wave at P with nrays rays wavelength 1E-6 and size 50E-6
	GOAT::maths::Vector<double> k(1, 0, 0); // direction of the wave
	k = k / abs(k); 
	LS.setk(k);
	
	// Region of interest
	GOAT::raytracing::Box Box(GOAT::maths::dzero, GOAT::maths::Vector<double>(2.0, 2.0, 2.0), 1.0);
	GOAT::raytracing::Scene S;
	
	S.setr0(r0); 
	S.addLightSource(&LS);
	S.addObject(&Box);
	S.setnS(1.0);

	double dt = 1E-13;  // width of the pulse (in seconds)
	double sigma = dt / (2.0 * M_LN2);
	

	double dwvl= wvl * wvl * M_LN2 / (M_PI * M_PI * GOAT::raytracing::C_LIGHT_MU * dt);// spectral FWHM
	std::cout << "dwvl=" << dwvl*1E+3 << "nm" << std::endl;
	std::vector<std::function<std::complex<double>(double) > > nList;
    
	std::cout << "Do raytracing...";
	GOAT::raytracing::Raytrace_usp rt(S,nCells);
	rt.setNumReflex(0);
	rt.trace();
	std::cout << "done." << std::endl;
	
	GOAT::raytracing::TrafoParms parms;	
	parms.nI = 2;
	parms.nR = 1;
	parms.nS = 5000;
	parms.lstart = wvl - 20.0*dwvl;
	parms.lstop = wvl + 20.0*dwvl;
	/*
	parms.nList.push_back(nGlass);
	parms.nList.push_back(nGlass);
	*/
    parms.nList.push_back(nGlass_BK7);
	parms.nList.push_back(nGlass_BK7);
	
	parms.dt = dt;
	parms.wvl = wvl;
	/*
	parms.nList.push_back(none);
	parms.nList.push_back(none);
	*/
	
	GOAT::raytracing::Trafo<GOAT::maths::Vector<double> > T(parms);
	std::string fname = "H:\\testsa.dat";
	// save(rt.SA[0], fname);

    
	
	std::string numStr;
	std::stringstream sstr;
	double t;
	int n = 20000;
	std::ifstream is;
	double h;
	std::string str;
	std::ofstream os ("h:\\data\\timeb.dat");
 
	for (int i = 0; i < n; i++)
	{
		t = 4.464e-9- i * dt / 10.0;
		os.precision(7);
		os << t*1E+9 << std::endl;
		// t =  4.4615e-9-(double)n/2.0*dt/10.0+i * dt/10.0;
		T.calc(rt.SA, t);
		sstr  << "h:\\data\\testb" << i << ".dat";
		fname = sstr.str();		
		sstr.str("");
		std::cout << "current filename:" << fname << "\t t=" << t*1E+9;
		GOAT::raytracing::saveabsE(T.SAres, fname);
		is.open(fname);
		std::getline(is, str);
		is >> h;
		is.close();
		std::cout.precision(7);
		std::cout << "\t" << h << std::endl;
	}
	os.close();
  return 0;
}
