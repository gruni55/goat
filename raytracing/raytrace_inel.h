/** \file raytrace_inel.h
 * Here you can find raytracing routines used for inelastic scattering
 * The file provides classes for inelastic scattering, but also for calculating the internal fields inside objects or in the whole 
 * calculation sphere 
 */
#pragma once
#include "raytrace.h"
#include "superarray.h"

namespace GOAT
{
    namespace raytracing
    {
        constexpr int INEL_CALCPHASE_EXCITATION=1;
        constexpr int INEL_CALCPHASE_RRT=2;
        constexpr int INEL_CALCPHASE_EXCITATION_ONLY=3;

        constexpr int INEL_MAX_NREFLEX=5;
        constexpr int INEL_SAVE_ABSE=1;

        constexpr int INEL_RADIATION_TYPE_FLOURESCENCE=1;
        constexpr int INEL_RADIATION_TYPE_RAMAN=2;

        constexpr int INEL_RADIATION_COHERENT=1;
        constexpr int INEL_RADIATION_INCOHERENT=2;

        constexpr int INEL_EXPORT_EXCITATION_FIELD_ABS=0;
        constexpr int INEL_EXPORT_EXCITATION_FIELD_VECTOR=1;


        /**   
        * This structure defines the detector for the inelastic scattering with a reference point P, the direction of detection n and the polarisation
        * directions (e1,e2) used for the calculation of the inelastic scattering. Default values: n=P/|P|,   
        */

        typedef struct
        {
            maths::Vector<double> P; ///< reference point of the detector
            maths::Vector<double> n; ///< direction of detection (backwards), by default set to n=P/|P|, but can be changed
            maths::Vector<double> e1, e2; ///< considered direction of polarisation, e1 is by default the projection of the z-axis onto a plane perpendicular to n with the reference point P. e2 is set to  
            double wvlinel=1.0; ///< Wavelength of the inelastic scattering
            int radiationType = INEL_RADIATION_TYPE_RAMAN; ///< Type of inelastic scattering: INEL_RADIATION_TYPE_FLOURESCENCE for flourescence, i.e. unpolarized or INEL_RADIATION_TYPE_RAMAN (default) for Raman scattering (polarized) 
            int coherency = INEL_RADIATION_INCOHERENT; ///< Scattering is coherent INEL_RADIATION_COHERENT or incoherent INEL_RADIATION_INCOHERENT
        } RRTParms;


        /**
         * @brief Calculates the parameters for the detector with help of the spherical coordinates theta and phi
         * This function generates a RRTParms structure for a scattering detection in the direction determined by the spherical coordinates theta and phi. 
         * 
         */
        RRTParms calcDetDirParms (double theta, double phi, double wvlinel); 

        /**   
        * @brief Class for calculating inelastic scattering (Raman) 
        * 
        * This class calculates the inelastic (Raman-) scattering of the active particles within the scene. The calculation is separated into two steps. 
        * At first, the exciting electric field is calculated and stored. In a second step the inelastic scattered field is calculated by following the rays
        * backwards from the detector towards the scattering objects. There it is weighted by the exciting field. Also the dipole characteristic of the 
        * scattering process is considered. The scattering behaviour is described by the polarizability tensor alpha. The theoretical background is described
        * in more detail in the articles below. Detailed information about the inelastic scattering calculation parameters are provided by the structure RRTParms. 
        *
        * 
        * 
        * -#  N. Velesco and G. Schweiger, �Geometrical optics calculation of inelastic scattering on large particles,� Appl. Opt. 38, 1046�1052 (1999)
        * -#  T. Weigel, J. Schulte, and G. Schweiger, �Inelastic scattering on particles with inclusions,� J. Opt. Soc. Am. A22, 1048�1052 (2005)
        * 
        */
        class Raytrace_Inel :
            public Raytrace
        {
        public:
            Raytrace_Inel();
            Raytrace_Inel(const Scene &S, int n);
            /**
            * @brief Calculates inelastic scattering
            * 
            * This function is the starting point for the inelastic scattering calculation. By the first call of this function, the excited field is calculated and then
            * the inelastic scattered field (RRT-field) is calculated. Here, the calculation in one direction is done. Since the calculation of the excited field is somewhat
            * time consuming, this field is not recalculated if the function is called a second time. If the scene has changed, just call the function resetCalculation() or 
            * change the currrent scene by calling sceneChanged(Scene &S) to start the whole calculation together from the beginning  (together with the excitation field). 
            * If only the excitation field is to be calculated in this function, then the function setExcitationFieldOnly() must be called beforehand. The exciting field is stored 
            * inside the object @see SGE (For details refer to #SuperArray). 
            * @param D Parameters for the inelastic scattering calculation
            */
            void trace(RRTParms D); 
            void traceEnterObject(); ///< Function internally called when ray enters an object
            void traceLeaveObject(); ///< Function internally called when ray leaves an object
            /**
             * @brief Change Scene to S
             * Changes the Scene to S and changes the calculation phase so that the calculation starts from the beginning when the trace() function is called, 
             * i.e. the excitation field is also recalculated.  
             */
            void sceneChanged(const Scene &S); ///< Change Scene to S (resets the calculation phase
            /**
             * @brief Export excitation field to (ASCII-) file
             * This function exports excitation field inside the objects. For each object, a separate file is created with name "<no of the object>_<fname>.dat". 
             * Each file contains the information about the field inside an array which contains the object. In the first line, the dimensions in x,y and z direction is given.
             * The kind of information, stored in the file is determined by the parameter savetype. In can have the following values: 
             *   INEL_EXPORT_EXCITATION_FIELD_ABS : the absolute value of the electric field is stored  (default value)
             *   INEL_EXPORT_EXCITATION_FIELD_VECTOR : the full maths::Vector is stored in the file. Each line stores the following values:
             *   real(Ex) imag(Ex) real(Ey) imag(Ey) real(Ez) imag(Ez)
             * 
             */
            void exportExcitation(std::string fname, int savetype = INEL_EXPORT_EXCITATION_FIELD_ABS); 
            double inel1, inel2; ///< result of the inelastic scattering for the two given polarisations
            void resetCalculation() { calcphase = INEL_CALCPHASE_EXCITATION; } ///< forces the calculation of the excited field when calling trace(RRTParms D)
            void setExcitationFieldOnly() { calcphase = INEL_CALCPHASE_EXCITATION_ONLY; } ///< Sets calculation phase so, that for the next call of trace(RRTParms D), only the excitation field is calculated.
            void unsetExcitationFieldOnly() { calcphase = INEL_CALCPHASE_RRT;  } ///< sets calculation phase in the way, that also the inelastic calculation will be done
            SuperArray* SGE;  ///< Here, the exciting field is stored 
            
        private:
            std::complex<double>  gewichte(maths::Vector<std::complex<double> > E, maths::Vector<std::complex<double> > p);
            void initExcitation();
            void initRRT();

            void traceExcitation();
            void saveExcitation();
            void saveRRT();
            void traceRRT();
            int calcphase;
            SuperArray* SGRRT1; 
            SuperArray* SGRRT2;
            
            bool* active;
            int n;
            int iR;
            LightSrc** LSExcit;
            int nLSExcit;
            RRTParms parms;
        };
    }
}


