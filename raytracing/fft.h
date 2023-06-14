#pragma once
#include <functional>
#include <vector>
#include "superarray.h"
#include "raytrace_usp.h"
#include "constants.h"

namespace GOAT
{
	namespace raytracing
	{
#ifndef I 
#define I std::complex<double> (0.0,1.0)
#endif
        constexpr int N_INTEGRAL = 500;

        /*! \brief Structure, which acts as a container for all informations needed to process the calculation.       
        * 
        * This structure is needed for the Fourier transformation which is used to calculate the temporal development 
        * of the electric field in space. Since dispersion is considered, for all objects and also for the surrounding medium
        * functions are necessary which calculate the complex-valued refractive index from the wavelength. These functions are
        * stored in the element nList. 
        */
        typedef struct TrafoParms
        {
            int nI=4;                     ///< defines the start of the integration range
            int nS=200;                ///< number of subdivision in the spectral range   
            int nR= INEL_MAX_NREFLEX;   ///< number of reflections considered in the raytracing part
            double omegaStart;              ///< lowest wavelength considered in the calculation 
            double omegaEnd;               ///< highest wavelength considered in the calculation 
            double omega0;                 ///< main frequency (corresponds to wvl)
            double wvl=1.0;                 ///< main wavelength   
            double dt=1E-15;                  ///< width of the pulse (in seconds) 
            std::vector<std::function<std::complex<double>(double) > > nList; ///< list of functions which describe the refractive index dependence on the wavelength (for each object one has to give one function) additionally one function for the surrounding medium
        };

        
        /*! \brief This class calculates the time dependence of a field which calculated before.
        * 
        * The calculation of the time dependence is made in two steps. In the first step, the fields are calculated 
        * for different wavelengths. In this class, the second step will be made. Here, the time dependence is calculated
        * with a kind of Fourier transformation. The class requires an array of SuperArray elements which contain the 
        * the information over all steps for every ray (length of the step and material). The process assumes, that the 
        * wavelength range can be subdivided into subranges in which the paths of the rays remain nearly the same. 
        * Since the phase is much more sensitive to changes of the refractive index, the phase change within these 
        * subranges is considered. The structure TrafoParms carries all information needed for the calculation.
        * The spectral resolution, needed to prevent "ghost" peaks is a function of the time difference dt to a certain 
        * reference time tref: dt=|tref-t|. By default tref is set to 0, but if you are interested of a range around 
        * a certain time, set tref to this time to reduce the amount of needed frequencies
        */
        class Trafo
        {
        public:
            Trafo();
            Trafo(TrafoParms);  ///< Constructor for initialization          
           
            /**
            * Perform a calculation at the time t
            * @param SA Array of SuperArray elements. Contains the fields for the different wavelengths  
            * defined by the parameters in the corresponding TrafoParms structure, given in the constructor
            */
            void calc(std::vector < std::vector<SuperArray <std::vector<gridEntry> > > >& SA, double t);
            SuperArray<maths::Vector<std::complex<double> > >SAres; ///< Container for the last result     
            void setRefractiveIndexFunctions(std::vector<std::function<std::complex<double>(double) > > nList);
            void setReferenceTime(double tref);
            void setTrafoParms(TrafoParms trafoparms); 

        private:
            void initResult(SuperArray<maths::Vector<std::complex<double> >>& SA);
            void initResult(double r0, int nx, int ny, int nz, ObjectShape** Obj, int numObjs);
            double pulseWeight(double omega);
        //    void createLTexpo();
            std::complex<double> calcPhase(std::vector<stepEntry> steps, double k0);
            GOAT::maths::Vector<std::complex<double> > calcOne(std::vector<stepEntry> steps, double t);
            GOAT::maths::Vector<std::complex<double> > integrate(double t, std::vector<gridEntry> ge, double omegastart, double omegastop);
            double twoSigma2; 
            double sigma2;
            double prefactor; // 1/(sigma * sqrt (2pi))
          //  double omegastart, omegastop;
            std::vector<double> freq;
            double tref = 0.0;
            TrafoParms tp;
            std::vector<std::function<std::complex<double>(double) > > nList; ///< List of the refractive index functions      
        };
	}
}
