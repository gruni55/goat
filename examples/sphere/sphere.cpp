#include "refractive_index_functions.h"
#include "pulsecalculation.h"
/**
  This program calculates the pulse propagation within a spherical object. It gives the absolute values of electric field inside the sphere
*/

int main(int argc, char** argv)
{
	double r0 = 1000;  // radius of the calculation space

        /* Setup of the spherical object */
	double r = 100.0; 				        // radius of the sphere
	GOAT::maths::Vector<double> Pos(0, 0, 0);   		// Position of the sphere
	GOAT::maths::Vector<double> dim(50, r, r);		    // set its dimensions
	GOAT::raytracing::Ellipsoid E(Pos, dim, 1.5, r0);
	E.setActive(false);  // the electric field inside the sphere will be calculated

        GOAT::maths::Vector<double> BoxPos(-100, 0, 0);
        GOAT::maths::Vector<double> BoxDim (200,100,100);
        GOAT::raytracing::Box box(BoxPos,BoxDim, 1.0);     
        box.setActive(true);
		
		GOAT::maths::Vector<double> Box2Pos(250, 0, 0);
		GOAT::maths::Vector<double> Box2Dim(300, 300, 300);
		GOAT::raytracing::Box box2(Box2Pos, Box2Dim, 1.0);
		box2.setActive(true);
   
        /* Setup of the light source: Plane wave */
	GOAT::maths::Vector<double> LSPos(-510, 0, 0);       // Position of the light source
 	int numRays = 150000;                                      // number of rays (per direction => in total numRaysxnumRays ) 
        double LSdim = 150;                             // Edge length of the light source

	GOAT::raytracing::LightSrcPlane_mc LS(LSPos, numRays, 1.0, LSdim); // it's a plane wave
	LS.setk(GOAT::maths::ex);                               // Light source emitts in positive x-direction

        /* Put light source and object into the scene */
	GOAT::raytracing::Scene S;            
	S.addLightSource(&LS);
	S.addObject(&E);
        S.setnS(1.0);
//        S.addObject(&box);
	S.addObject(&box2);
	S.setr0(r0);
/*
	GOAT::raytracing::Raytrace_Path rp(S);
	rp.trace("test.dat"); 
*/	

        // Set the refractive index functions for the object and the surrounding medium 
	std::vector<std::function<std::complex<double>(double) > > nList;  // List to store the refractive index functions
	//nList.push_back(GOAT::raytracing::n_BK7);                          // function for the object (sphere)
    nList.push_back(GOAT::raytracing::n_BK7);                       // function for the box
	nList.push_back(GOAT::raytracing::n_Vacuum);                       // function for the surroundings
	nList.push_back(GOAT::raytracing::n_Vacuum);

        // Now, we prepare the calculation... 
	GOAT::raytracing::pulseCalculation pc(S);   
	pc.setPulseWidth(100);                                             // set pulse width to 100fs 
	pc.setRefractiveIndexFunctions(nList);                             // set the refractive index function list 
	pc.setSpatialResolution(1.0);                                      // set the spatial resolution of the grid to 100.0 µm 
	pc.setReferenceTime(0.0);                                          // set the reference time to 0.0fs
	                          
	double t = 700.0;                                                 // set the time (in fs) at which you want to calculate the fields
	pc.field(t); 	                                                   // make calculation
	GOAT::raytracing::saveabsE(pc.trafo.SAres, "fieldabs.dat", 0);     // store absolute values of the electric field inside the sphere (object number 0)  
	GOAT::raytracing::saveabsE(pc.trafo.SAres, "fieldabs1.dat", 1);

	return 0; 
}
