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
        typedef struct 
        {
            int nI=4;                     ///< defines the number of spectral ranges
            int nS=200;                ///< number of subdivision per spectral range   
            int nR= INEL_MAX_NREFLEX;   ///< number of reflections considered in the raytracing part
            double omegaStart;              ///< lowest wavelength considered in the calculation 
            double omegaEnd;               ///< highest wavelength considered in the calculation 
            double omega0;                 ///< main frequency (corresponds to wvl)
            double wvl=1.0;                 ///< main wavelength (in µm)  
            double dt=100;                  ///< width of the pulse (in femto seconds) 
            std::vector<std::function<std::complex<double>(double) > > nList; ///< list of functions which describe the refractive index dependence on the wavelength (for each object one has to give one function) additionally one function for the surrounding medium
        } TrafoParms;

        
        /*! \brief This class calculates the time dependence of a field which calculated before.
        * 
        * The calculation of the time dependence is made in two steps. In the first step, the fields are calculated with help of 
        * raytracing. This class uses the information from the raytracing (path of the rays) to make the concrete calculation of the electric field 
        * for different wavelength equivalent to a Fourier transform. The class requires an array of SuperArray elements which contain the 
        * the information over all steps for every ray (length of the step and material).  The structure TrafoParms carries all information needed 
        * for the calculation. The spectral resolution, needed to prevent "ghost" peaks is a function of the time difference dt to a certain 
        * reference time tref: dt=|tref-t|. By default tref is set to 0, but if you are interested of a range around 
        * time, set tref to an approbriate time to reduce the amount of needed frequencies
        */
        class Trafo
        {
        public:
            Trafo();
            /**
            * @brief Constructor for initialization
            * Perform a calculation at the time t
            * @param SA Array of SuperArray elements. Contains the fields for the different wavelengths
            * defined by the parameters in the corresponding TrafoParms structure, given in the constructor
            */
            Trafo(TrafoParms);  
            void calc(std::vector < std::vector<SuperArray <std::vector<gridEntry> > > >& SA, double t); 
            /**
             * @brief calculation of the Fourier transformation.
             * The calculation of the Fourier transformation is made in a certain frequency intervall at the time t
             * @param SA this SuperArray structure contains the information gathered in the raytracing part. 
             * @param omegaStart Start frequency (in \f$fs^{-1}\f$) of the frequency intervall 
             * @param omegaStop End frequency (in \f$fs^{-1}\f$) of the frequency intervall 
             */
            void calc(std::vector<SuperArray <std::vector<gridEntry> > > & SA, double omegaStart, double omegaEnd, double t, bool do_clear=true);
            SuperArray<maths::Vector<std::complex<double> > >SAres; ///< Container for the last result     
            /**
             * @brief Set the refractive index functions.
             * Here, the refractive index functions were given with help of a std::vector. Each function must return a complex number 
             * (the refractive index) and has a double value (the wavelength) as argument
             */
            void setRefractiveIndexFunctions(std::vector<std::function<std::complex<double>(double) > > nList); 
            /**
             * @brief Sets the reference time (in fs).
             */
            void setReferenceTime(double tref);
            /**
             * @brief Sets the parameters used in the calculation.
             */
            void setTrafoParms(TrafoParms trafoparms);  ///< set the current transformation parameters (for details see: TrafoParms)
            void clear(); ///< Deletes all arrays (frees the memory)
            void empty(); ///< fills the entire result array with zero vectors 
            /**
             * @brief Resets and initialize the result array.
             * The result array will be cleared and renewed according to the settings of the SuperArray SA
             */
            void initResult(SuperArray<maths::Vector<std::complex<double> >>& SA); 
            /**
             * @brief Resets and initialize the result array.
             * The result array will be cleared and renewed according to the parameters
             * @param r0 Radius of the calculation space
             * @param nx Number of grid elements in x-direction
             * @param ny Number of grid elements in y-direction
             * @param nz Number of grid elements in z-direction
             * @param Obj List of all objects in the scene
             * @param numObjs Numnber of objects in List
             */
            void initResult(double r0, INDEX_TYPE nx, INDEX_TYPE ny, INDEX_TYPE nz, ObjectShape** Obj, int numObjs);
            std::vector<std::function<std::complex<double>(double) > > nList; ///< List of the refractive index functions      
        private:
            double pulseWeight(double omega);
            /**
             * @brief sets current values for the refractive indices
             * This method calculates the refractive indices of the different objects for a
             * given wavelength.
             * @param wvl current wavelength
             *
             */
            void setCurrNList(double wvl)
            {
                for (int i = 0; i < nList.size(); i++)
                    currNList[i] = nList[i](wvl);
            }
        //    void createLTexpo();
            std::complex<double> calcPhase(std::vector<stepEntry> steps, double k0);
            GOAT::maths::Vector<std::complex<double> > calcOne(std::vector<stepEntry> steps, double t);
            GOAT::maths::Vector<std::complex<double> > integrate(double t, std::vector<gridEntry> ge, double omegastart, double omegastop);
            double twoSigma2; 
            double sigma2; ///< for a gaussian pulse, the electric field, the temporal behavior is: \f$E(t)=E0 \cdot e^{-t^2/(2*\sigma^2)}\f$
            double prefactor; // 1/(sigma * sqrt (2pi))
            std::vector<double> freq;
            double tref = 0.0;
            TrafoParms tp;
            std::vector<std::complex<double > > currNList; ///< here, the refractive indices for the current wavelength are stored (for faster calculation)            
        };
	}
}
