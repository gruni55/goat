#include "raytrace_usp.h"
#include "fft.h"
#include <complex>
#include <functional>
#include <sstream>
#include "pulsecalculation.h"

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

std::complex<double> nTest(double wvl)
{
	return 1.5065;
}


int main (int argc, char **argv)
{
	int nrays = 10000;  
	double r0 = 2E+6;      // Radius 2m
	int nCells = 2000000;  // resolution 1µm
	/* ---  Let's first define the Scene ----
	    with one lightsource and one box, where the field is stored 
	*/
	
	double wvl = 1.0; // peak wavelength
	
	// use light source with arbitrary ray distribution 
	GOAT::maths::Vector<double> P(-(0.5E+6+1.0), -1.0, -1.0);          // Position of the light source
	//GOAT::raytracing::LightSrcPlane LS(P, nrays, 1.0, 1);      // define Plane wave at P with nrays rays wavelength 1E-6 and size 50E-6
	GOAT::raytracing::LightSrcPlane_mc LS(P, nrays, 1.0, 1);
	GOAT::maths::Vector<double> k(1, 0, 0);					   // direction of the wave
	k = k / abs(k); 
	LS.setk(k);
	
	// Region of interest
	GOAT::raytracing::Box Box(GOAT::maths::dzero, GOAT::maths::Vector<double>(2.0, 2.0, 2.0), 1.0); // distance light source <-> "box" (edge) 1m
	GOAT::raytracing::Scene S;
	
	S.setr0(r0); 
	S.addLightSource(&LS);
	S.addObject(&Box);
	S.setnS(1.0);

	std::string fname;

    
	std::vector<std::function<std::complex<double>(double) > > nList; // List of the refractive index functions
	/* nList.push_back(none);
	nList.push_back(none);
	*/

	nList.push_back(nGlass_BK7);
	nList.push_back(nGlass_BK7);

	std::string numStr;
	std::stringstream sstr;
	double t;
	int n = 400;
	std::ifstream is;
	double h;
	std::string str;
	std::ofstream os ("h:\\data\\timec.dat");     // here, the time is stored
	std::ofstream oserg("h:\\data\\erg.dat");

	double dt = 100.0;									    	// width of the pulse (in femtoseconds)
	double sigma = dt / (2.0 * sqrt(2.0 * M_LN2));

	double dwvl = wvl * wvl * 8.0 * M_LN2 / (M_PI * GOAT::raytracing::C_LIGHT_MU_FS * dt);  // spectral FWHM
	std::cout << "dwvl=" << dwvl << "µm" << std::endl;
	GOAT::raytracing::pulseCalculation pc(S);     // The calculation class
	pc.setRefractiveIndexFunctions(nList);         
	pc.setPulseWidth(dt);
	pc.setSpatialResolution(2.0);
    double tref =1E6;
	pc.setReferenceTime(tref);
		
	int l = 0;
	//for (double t=5E6+40000; t<=5E6+60000;t+=10)
	// for (double t=5.06E+6; t<=5.1E+6; t+=10)
	for (double t = 2.535E+6; t <= 2.54E+6; t += 10)
	//for (double t = -5000; t <= 5000; t += 10)
	{		
		os.precision(12);
		os << std::scientific << t << std::endl;		
		pc.field(t); // calculate fields at time t 
		sstr  << "h:\\data\\testc" << l << ".dat";
		fname = sstr.str();		
		sstr.str("");
		os.precision(12);
		std::cout << std::scientific << t;
		GOAT::raytracing::saveabsE(pc.trafo.SAres, fname);
		
		is.open(fname);
		std::getline(is, str);
		is >> h;
		is.close();
		
		std::cout << "\t" << h << std::endl;
	    oserg << t << "\t" << h << std::endl;
		l++;
	}
	os.close();
	oserg.close();
  return 0;
}
