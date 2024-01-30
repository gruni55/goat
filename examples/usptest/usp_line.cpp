#include "pulsecalculation.h"
#include "refractive_index_functions.h"
#include <fstream>
/**
* Example: Calculates the pulse after a glass box of 1m length
*  
*/


int main (int argc, char ** argv)
{
	double r0 = 1E+8;
	GOAT::raytracing::Scene S;
	S.setr0(r0);

	// ----------------- Light source ---------------
	GOAT::maths::Vector<double> LSPos(-(5E5 + 1E4), 0, 0);
	GOAT::raytracing::LightSrcPlane_mc LS(LSPos, 1, 1, 1E-3);
	LS.setk(GOAT::maths::ex);
	LS.setPol(GOAT::maths::vdc(GOAT::maths::ey));
	S.addLightSource(&LS);

	// ----------------- Glass ---------------------
	GOAT::maths::Vector<double> GlassBoxDim(1E5, 1, 1);
	GOAT::raytracing::Box GlassBox(GOAT::maths::dzero, GlassBoxDim, 1);
	GlassBox.setActive(false);
    S.addObject(&GlassBox);

	// ---------------- Detektor -------------------
	GOAT::maths::Vector<double> DetPos(5E5 + 10 + 2.5E4,0,0);
	// GOAT::maths::Vector<double> DetPos(0, 0, 0);
	GOAT::maths::Vector<double> DetDim(5E4, 0.5, 0.5);
	GOAT::raytracing::Box Det(DetPos, DetDim, 1);
	Det.setActive(true);
	S.addObject(&Det);

	// --------------- Wavelength dependent 
	std::vector< std::function< std::complex< double >(double) > > nList;

	nList.push_back(GOAT::raytracing::n_BK7);    // Glass box
	nList.push_back(GOAT::raytracing::n_Vacuum); // Detector
	nList.push_back(GOAT::raytracing::n_Vacuum); // Intermediate medium

	double pulseWidth = 100;
	double spatialRes = 0.25;
	
	GOAT::raytracing::pulseCalculation pc(S);
	pc.setSpatialResolution(spatialRes);
	pc.setRefractiveIndexFunctions(nList);
	pc.setSpectralRanges(500);
	pc.setNumWavelengthsPerRange(500);
	
	std::string FName = "C:\\users\\thomas\\data\\glasbox.dat";
	
	double time;
	std::ofstream os;
	os.open(FName);

	
	std::vector<double> d;
	std::vector<double> avgd;
	double fwhm;
	pc.setReferenceTime(1E6);
	double wvl = 0.5;
	// for (double wvl = 0.3; wvl <= 1.0; wvl += 0.001)
	{
		pc.setCenterWavelength(wvl);
		time=pc.findHitTime(1);
		pc.field(time+1E5);

	    GOAT::raytracing::saveabsE(pc.trafo.SAres, "C:\\users\\thomas\\data\\glasbox_field.dat", 1);

		for (int i = 0; i < pc.trafo.SAres.n[1][0]; i++)
			d.push_back(abs2(pc.trafo.SAres.G[1][i][1][1]));


		avgd = GOAT::maths::movingAvg(d, 50);
		std::size_t maxIndex;
		GOAT::maths::findmax(avgd, maxIndex);
		fwhm = GOAT::maths::FWHM(avgd, maxIndex) * spatialRes / GOAT::raytracing::C_LIGHT_MU_FS;
		
		os << wvl << "\t" << fwhm << "\t" << time << std::endl;
		std::cout << wvl << "\t" << fwhm << "\t" << time << std::endl;

		d.clear();
		avgd.clear();
	}

	os.close();



	return 0;
}