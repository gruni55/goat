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
        };

        /**
         * @brief Plane wave with random ray distribution.
         * This class provides a plane (square sized) wave with a given width.
         * The rays are arbitrarily but uniformly distributed inside the light source area. 
         */
        class LightSrcPlane_mc : public LightSrcPlane
        {
            public:
                LightSrcPlane_mc (const LightSrcPlane_mc & L);
                /**
                 * Constructor for arbitrary distributed plane wave.
                 * @param Pos Position of the source (center)
                 * @param N Total number of rays 
                 * @param wvl Wavelength
                 * @param Pol Direction of the polarization
                 * @param raytype Type of ray used for the raytracing process
                 * @param r0 Radius of the calculation space 
                 */
                LightSrcPlane_mc (maths::Vector<double> Pos, int N, double wvl, double D = 100.0, 
                                  maths::Vector<std::complex<double> > Pol = maths::Vector<std::complex<double> >(0.0, 1.0, 0.0), 
                                  int raytype = LIGHTSRC_RAYTYPE_IRAY, double r0 = 100.0);
                int next(IRay& S);
                int next(tubedRay& S);
                int next(Ray_pow& S);
                GOAT::maths::Vector<double> genStartingPos ();
                void reset();
        };


        /** @brief This class provides a light source, which emitts along a straight line (random ray distribution)
        */
        class LightSrcLine_mc : public LightSrc
        {
        public:
            LightSrcLine_mc();
            /**
            *  @brief constructor
            * @param Pos Position of the light source
            * @param N number of rays along straight line
            * @param size length of the light source
            * @param direction direction of light source (not the direction of emission !)
            */
            LightSrcLine_mc(maths::Vector<double> Pos, int N, double wvl, double size, maths::Vector<double> k = maths::ez, maths::Vector<double> direction = maths::ey);
            int next(RayBase* ray);
            int next(IRay& S);
            int next(Ray_pow& S);
            int next(tubedRay& S);
            void binWriteItem(std::ofstream& os) { /* to be implemented !!! */ }
            void binReadItem(std::ifstream& os) { /* to be implemented !!! */ }
            maths::Vector<double> direction;
            GOAT::maths::Vector<double> genStartingPos();
            double size;
        };

        /**
         * @brief  Ring shaped light source.
         * This class provides a ring shaped light source with arbitrary, uniform 
         * distributed rays within the light source area. 
         */

        class LightSrcRing_mc : public LightSrcPlane
        {
          public:
            LightSrcRing_mc(const LightSrcRing_mc& L);
            LightSrcRing_mc( maths::Vector<double> Pos, int N, double wvl,double rmin, double rmax,
                maths::Vector<std::complex<double> > Pol = maths::Vector<std::complex<double> >(0.0, 1.0, 0.0),
                int raytype = LIGHTSRC_RAYTYPE_IRAY, double r0 = 100.0);
            int next(IRay& S);
            int next(tubedRay& S);
            int next(Ray_pow& S);
            double getRmin () { return rmin; }
            double getRmax () { return rmax; }
            void setRmin(double rmin);
            void setRmax(double rmax); 
         private:
            double rmin = 0.0; ///< inner radius of the ring
            double rmax = 1.0; ///< outer radius of the ring
            GOAT::maths::Vector<double> genStartingPos();
            void reset(); ///< sets the internal ray counter to zero (intended for internal use only)
        };       


        /**
        * @brief The LightSrcRingGauss_mc provides a ring shaped light source with a gaussian distribution
        * This class describes a light source which cuts a ring out of a gaussian shape. The phase distribution is flat,
        * i.e. the phase remains constant along the whole light source surface
        */
       class LightSrcRingGauss_mc : public LightSrcPlane
        {
          public:
            LightSrcRingGauss_mc(const LightSrcRingGauss_mc& L);
            /**
             * @brief LightSrcRingGauss_mc
             * @param Pos Position of the light source
             * @param N number of rays
             * @param wvl wavelength
             * @param rmin inner radius
             * @param rmax outer radius
             * @param Pol polarisation
             * @param raytype type of ray
             * @param r0 radius of the calculation space
             */
            LightSrcRingGauss_mc( maths::Vector<double> Pos, int N, double wvl,double rmin, double rmax,
                maths::Vector<std::complex<double> > Pol = maths::Vector<std::complex<double> >(0.0, 1.0, 0.0),
                int raytype = LIGHTSRC_RAYTYPE_IRAY, double r0 = 100.0);
            int next(IRay& S);
            int next(tubedRay& S);
            int next(Ray_pow& S);
            void setFWHM (double r); ///< set width of the gaussian distribution to r
            double getRmax() { return rmax;} ///< returns outer radius
            double getRmin() { return rmin;} ///< returns inner radius
            void setRmin(double rmin); ///< set the inner radius to rmin
            void setRmax(double rmax); ///< set the outer radius to rmax

        private:
            double rmin = 0.0; ///< inner radius of the ring
            double rmax = 1.0; ///< outer radius of the ring          
            GOAT::maths::Vector<double> genStartingPos(); ///< generate a position within the ring (needed by the next() method)
            void reset();
            double sigma2; ///< sigma^2, used for internal purposes
        }; 
    }    
}
