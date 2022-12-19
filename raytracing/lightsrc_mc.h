/**
* /file lightsrc_mc.h
* Here you can find the light sources with random ray distribution
*/
#pragma once
#include "lightsrc.h"
namespace GOAT
{
    namespace raytracing
    {
        /**
         * @brief This class provides a gaussian beam with arbitrary distributed rays 
        */
        class LightSrcGauss_mc : public LightSrcGauss
        {
            public: 
                LightSrcGauss_mc (const LightSrcGauss_mc & L);                
			/**
			 * @brief Main constructor for gaussian beam light source
			 * The direction of the beam is given by the starting point and the focus position.
			 * \param Pos Starting position of the beam (center)
			 * \param N number of rays (per direction)
			 * \param wvl wavelength
			 * \param w0 waist radius (virtual, only needed to calculate the electric field distribution correctly)
			 * \param focuspos Position of the focus
			 * \param D Size of starting area (side length)
			 * \param Pol Polarisation (on the starting area, at the axis of the beam)
			 * \param raytype (type of ray used for the calculation)
			 * \param r0 Radius of the calculation sphere
			 */
			LightSrcGauss_mc(maths::Vector<double> Pos, int N, double wvl, double w0, maths::Vector<double> focuspos, 
                          double D = 1.0, maths::Vector<std::complex<double> > Pol = maths::Vector<std::complex<double> >(0.0, 1.0, 0.0),
                         int raytype = LIGHTSRC_RAYTYPE_PRAY, double r0 = 1.0);
                int next(Ray_pow& S);
                int next (IRay& S);
                int next(tubedRay& ray);
                void reset();                
                // int next (tubedRay &S);
                GOAT::maths::Vector<double> genStartingPos ();
                double stddev;
                int rayCounter=0;
        };

        class LightSrcPlane_mc : public LightSrcPlane
        {
            public:
                LightSrcPlane_mc (const LightSrcPlane_mc & L);
                LightSrcPlane_mc (maths::Vector<double> Pos, int N, double wvl, double D = 100.0, 
                                  maths::Vector<std::complex<double> > Pol = maths::Vector<std::complex<double> >(0.0, 1.0, 0.0), 
                                  int raytype = LIGHTSRC_RAYTYPE_IRAY, double r0 = 100.0);
                int next(IRay& S);
                int next(tubedRay& S);
                int next(Ray_pow& S);
                GOAT::maths::Vector<double> genStartingPos ();
                int rayCounter=0;
                void reset();
        };

        
    }    
}