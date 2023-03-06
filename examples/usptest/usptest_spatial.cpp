#include "raytrace_usp.h"
#include "raytrace_inel.h"
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
	int nrays = 100;
	double r0 = 1E+6;  // radius of the calculation sphere: 1m (diameter: 2m)
	int nCells = 2000000; //  number of cells in the virtual grid => resolution 0.1µm
	/* ---  Let's first define the Scene ----
	    with one lightsource and one box, where the field is stored 
	*/
	
	double wvl = 1.0; // peak wavelength 1µm
	
	GOAT::maths::Vector<double> P(-0.8E+6, 0, 0);		            // Position of the light source
	GOAT::raytracing::LightSrcPlane LS(P, nrays*nrays, 1.0, 20); // define Plane wave at P with nrays rays wavelength 1E-6 and size 10µm
	GOAT::maths::Vector<double> k(1, 0, 0);                         // direction of the wave	
	LS.setk(k);
	
	// Region of interest
	// GOAT::raytracing::Box Box(GOAT::maths::dzero, GOAT::maths::Vector<double>(50.0, 5.0, 5.0), 1.0);
	GOAT::raytracing::Ellipsoid Ell(GOAT::maths::dzero, GOAT::maths::Vector<double>(30, 20, 10), 1.5);
	GOAT::raytracing::Scene S;
	
	// Update the scene
	S.setr0(r0); 
	S.addLightSource(&LS);
	S.addObject(&Ell);
	S.setnS(1);

	double dt = 1E-13;                                     // width of the pulse (in seconds)
	

	double dwvl= wvl * wvl * M_LN2 / (M_PI * M_PI * GOAT::raytracing::C_LIGHT_MU * dt);// spectral FWHM
	std::cout << "dwvl=" << dwvl*1E+3 << "nm" << std::endl;
	std::vector<std::function<std::complex<double>(double) > > nList;
    
	std::cout << "Do raytracing...";
	// GOAT::raytracing::Raytrace_usp rt(S,nCells);
//	GOAT::raytracing::Raytrace_Path rt(S);
//	rt.trace("strahl.dat");
	GOAT::raytracing::Raytrace_Inel rt(S, nCells);
	GOAT::raytracing::RRTParms rrtparms;

	//rt.setNumReflex(0);
	rt.setExcitationFieldOnly();
	rt.trace(rrtparms);
	rt.exportExcitation("field.dat", GOAT::raytracing::INEL_EXPORT_EXCITATION_FIELD_ABS);
	saveabsE(rt.SGE[0], "feld.dat");

/*	GOAT::raytracing::TrafoParms parms;
	parms.nI = 2;
	parms.nR = 1;
	parms.nS = 500;
	parms.lstart = wvl - 20.0*dwvl;
	parms.lstop = wvl + 20.0*dwvl;
	*/
	/*
	parms.nList.push_back(nGlass);
	parms.nList.push_back(nGlass);
	
        parms.nList.push_back(nGlass_BK7);
	parms.nList.push_back(nGlass_BK7);
	
	parms.dt = dt;
	parms.wvl = wvl;
	
	parms.nList.push_back(none);
	parms.nList.push_back(none);
	*/
	/*
	GOAT::raytracing::Trafo<GOAT::maths::Vector<double> > T(parms);
	double t;
	int n = 20000;
        t=4.0015E-9;
        T.calc(rt.SA, t);
        saveabsE (T.SAres,"test.dat");
		*/
  return 0;
}
