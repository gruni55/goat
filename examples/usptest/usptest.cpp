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

int main (int argc, char **argv)
{
	int nrays = 3;  
	double r0 = 2E+6;      // Radius 2m
	int nCells = 200000;  // resolution 1µm
	/* ---  Let's first define the Scene ----
	    with one lightsource and one box, where the field is stored 
	*/
	
	double wvl = 1.0; // peak wavelength
	
	// use light source with arbitrary ray distribution 
	GOAT::maths::Vector<double> P(-50-1.0, 0, 0);          // Position of the light source
	GOAT::raytracing::LightSrcPlane LS(P, nrays, 1.0, 1);      // define Plane wave at P with nrays rays wavelength 1E-6 and size 50E-6
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

	double dt = 1E-15;									    	// width of the pulse (in seconds)
	double sigma = dt / (2.0*sqrt(2.0 * M_LN2));                      

	double dwvl= wvl * wvl * 8.0 * M_LN2 / (M_PI *  GOAT::raytracing::C_LIGHT_MU * dt);  // spectral FWHM
	std::cout << "dwvl=" << dwvl << "µm" << std::endl;
	std::string fname;

    
	std::vector<std::function<std::complex<double>(double) > > nList; // List of the refractive index functions
	nList.push_back(none);
	nList.push_back(none);

	std::string numStr;
	std::stringstream sstr;
	double t;
	int n = 400;
	std::ifstream is;
	double h;
	std::string str;
	std::ofstream os ("h:\\data\\timec.dat");     // here, the time is stored
	GOAT::raytracing::pulseCalculation pc(S);     // The calculation class
	pc.setRefractiveIndexFunctions(nList);         
	pc.setPulseWidth(2E-15);
    double tref = 0.0;
	pc.setReferenceTime(tref);
		
	int l = 0;
	for (double t=-2E-14; t<=2E-14;t+=1E-16+1E-17)
	{		
		os.precision(12);
		os << std::scientific << (t-tref)*1E+9 << std::endl;		
		pc.field(t); // calculate fields at time t 
		sstr  << "h:\\data\\testc" << l << ".dat";
		fname = sstr.str();		
		sstr.str("");
		os.precision(12);
		std::cout << "current filename:" << fname << "\t t=" << std::scientific << (t-tref)*1E+9;
		GOAT::raytracing::saveabsE(pc.trafo.SAres, fname);
		
		is.open(fname);
		std::getline(is, str);
		is >> h;
		is.close();
		std::cout.precision(7);
		std::cout << "\t" << h << std::endl;
		l++;
	}
	os.close();
  return 0;
}
