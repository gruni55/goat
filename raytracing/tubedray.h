/***************************************************************************
                          tubed.h  -  description                              
                             -------------------                                         
    begin                : Fri Oct 15 1999                                           
    copyright            : (C) 1999 by Thomas Weigel                         
    email                : weigel@lat.ruhr-uni-bochum.de                                     
 ***************************************************************************/



#pragma once
#include "raybase.h"
#include <complex>
#include "fresnel.h"
#include "resutil.h"
#include "vector.h"
#include "matrix.h"
#include "plane.h"
#include "objectshape.h"

namespace GOAT {
    namespace raytracing {
        
        /*#ifndef KORR
        #define KORR 1E-8
        #endif */

#define KEIN_SCHNITTPUNKT     -1
#define GESCHNITTEN            0
#define IN_Plane               1

#ifndef EPS_WINKEL 
#define EPS_WINKEL 0 /// Bei Brechung werden nur Strahlen weiterverfolgt, die einen Winkel > EPS_WINKEL gegenüber der Oberfläche aufweisen
#endif 


        /**
          *@author Thomas Weigel
          */


          /**
             This structure is used for the tubedRay constructor to simplify the interchange of data
                  */
        typedef struct {
            std::complex<double>  n; /// current refractive index 
            maths::Vector<double> P[5];  /// current position of the edge ray (index 0-3) and of the center ray (index 4)
            maths::Vector<double> k[5];  /// current direction of the edge ray (index 0-3) and of the center ray (index 4)
            maths::Vector<std::complex<double> > E[5]; /// current electric field of the edge ray (index 0-3) and of the center ray (index 4)
            int numObj; /// number of objects within the scene
            ObjectShape** Obj; /// list of objects in the scene
            double r0; /// radius of the calculation space
            maths::Vector<std::complex<double> >*** Gitter;
        } tubedRayBuffer;



        /**
         * This structure is used for the gaussian beam to simplify the data transfer.
         */
        typedef struct
        {
            double w0;          /// (virtual) beam waist 
            maths::Vector<double> F;   /// position of the beam
            maths::Vector<double> k;   /// direction of the beam
            bool isGauss; /// is it really a gaussian beam
        } Gauss;


        void binWrite(Gauss gs, std::ofstream& os); /// write the Gauss structure to a binary file
        void binRead(Gauss& Gs, std::ifstream& is); /// read the Gauss structure from a binary file


        std::ostream& operator << (std::ostream& os, Gauss gs); /// output operator for the Gauss structure

#ifndef OHNE_GAUSS
#define OHNE_GAUSS
        const Gauss ohne_gauss = { 0,maths::Vector<double>(0,0,0),maths::Vector<double>(0,0,0),false };
#endif

        /**
        *  @brief This class represents a ray with a finite cross section.
        *  class tubedRay:
        *   This class represents a ray with a finite cross section. Therefore it is defined by four edge rays and one probe ray in the middle. Therefore,
        *   selffocusing of the ray can be considered. Each of these rays has its own electric field and directional vector. The probe ray is used to decide wether an object is
        *   hidden or not.
        */
        class tubedRay : public RayBase {
        public:
            tubedRay();
            tubedRay(const tubedRay& ray);
            tubedRay(tubedRayBuffer& B);
            /**
             * Constructor of tubedRay which is described by the
             *
             * \param g
             */
            tubedRay(Plane E, const maths::Vector<double>& p, double dy, double dz, const
                maths::Vector<std::complex<double> >& Pol, std::complex<double>  n0, double r0, double
                k0, const int Anzein = 0, ObjectShape** Einschluss = NULL, bool logRay = false);
            tubedRay(const maths::Vector<double>& p, double dy, double dz, const maths::Vector<std::complex<double> >& Pol, const maths::Vector<double>& K, std::complex<double>  n0, double r0, double k0, const int Anzein = 0, ObjectShape** Einschluss = NULL, bool logRay = false);
            void setGauss(Gauss g) { this->g = g; }
            bool next(); ///< do the next step
            maths::Vector<std::complex<double> > getE() { return E[4]; }
            int objectIndex() { return objIndex; }
            void reflectRay(RayBase*& tray, maths::Vector<double> n, std::complex<double> n1, std::complex<double> n2);
            tubedRay reflect(maths::Vector<double> n, std::complex<double>  n1, std::complex<double>  n2); /// Strahl wird reflektiert
            void refract(maths::Vector<double> n, std::complex<double>  n1, std::complex<double>  n2); /// Strahl wird gebrochen
            void tunnel(maths::Vector<std::complex<double> > Pol, std::complex<double>  n1, std::complex<double>  n2);
            void tunnel(maths::Vector<std::complex<double> > Pol, std::complex<double>  na, std::complex<double>  np, int l);
            double flaeche(); /// 
            void setRefract(std::complex<double>  n) { this->n = n; } ///< set current refractive index
            std::complex<double>  getRefract() { return n; } ///< returns the current refractive index
            void setk(int i, const maths::Vector<double>& K) { k[i] = K; }
            void setphi(int i, const std::complex<double>& p) { phi[i] = p; }
            void setP(int i, const maths::Vector<double>& p) { P[i] = p; }
            void setiR(int i) { iR = i; }
            void setGetunnelt(bool v) { getunnelt = v; }
            void setN0(std::complex<double> n) { n0 = n; }
            maths::Vector<double> getk(int i) { return k[i]; } ///< returns the directional vector of the i-th sub-ray
            maths::Vector<double> getk() { return k[4]; } ///< return the directional vector of the central sub-ray (same as getk(4) )
            maths::Vector<double> getP(int i) { return P[i]; } ///< returns the position of the i-th sub-ray
            maths::Vector<double> getP() { return P[4]; } ///< returns the position of the central sub-ray (same as getP(4) )
            std::complex<double>  getphi(int i) { return phi[i]; } 
            void getP(maths::Vector<double> p[5]) { p = P; } ///< returns an array whichs holds the current positions of the sub-rays
            maths::Vector<std::complex<double> > getAmp(int i) { return E[i]; } ///< returns the value of the electric field of the i-th sub-ray
            void setAmp(int i, const maths::Vector<std::complex<double> >& A) { E[i] = A; } ///< sets the value of the electric field of the i-th sub-ray
            ObjectShape* Einschluss(int i) { return Obj[i]; } 
            ~tubedRay();
            bool Einschluss() { return inObject; } ///< like isInObject (for compatibilty reasons, will be removed in future) 
            bool isInObject() { return inObject; } ///< Is the ray currently inside an object? 
            //  maths::Vector<double> checkEinschluss(int s, const maths::Vector<double>& Ps,int& Index);
            void checkObjectIntersection(int Index[5], maths::Vector<double>* Pmin);

            bool Getunnelt() { return getunnelt; }
            int currentObjectIndex() { return objIndex; }  ///< returns index of the currently hitted object
            double Intensity(int i) { return abs2(E[i]); } ///< returns the squared magnitude of the electric field
            friend std::ostream& operator << (std::ostream& os, tubedRay S);
            int reflections() { return iR; }
            /**
            * @brief This method checks if there is a phase jump due to self focusing
            * The method checks if neighboring edge rays cross each other between the start position and the end position. In this case a phase jump of -pi has to be considered.
            * The parameter epsilon is there for compatibility reasons and will not be considered (the parameter will be removed in future)
            * @param[in] P1 Array which holds the start position of the edge rays 
            * @param[in] P2 Array which holds the end position of the edge rays
            * @return value of the phase jump (0, \f$-\frac{\pi}{2}\f$ or \f$-\pi\f$) 
            */
            double pjump(maths::Vector<double> P1[5], maths::Vector<double> P2[5], const double epsilon = 1E-10);
            /**
            * @brief This method checks if there is a phase jump due to self focusing
            * The method checks if neighboring edge rays cross each other. In this case a phase jump of -pi has to be considered.
            * @param[in] P1 Array which holds the start position of the edge rays 
            * @param[in] P2 Array which holds the end position of the edge rays
            * @param[out] S Array of two vectors which hold the position of the foci for self focusing, if it exists otherwise it is set to a zero vector 
            * @return value of the phase jump (0, \f$-\frac{\pi}{2}\f$ or \f$-\pi\f$) 
            */
            double pjump(maths::Vector<double> P1[5], maths::Vector<double> P2[5], maths::Vector<double>* S);

            double pjump(void); ///< for internal use only
            double normVol(maths::Vector<double> P[5], maths::Vector<double> k, maths::Vector<double> n);
            bool schneidePlane(maths::Vector<double>* Erg, const Plane& E);
            int schneidePlane(const Plane& E, double d, maths::Vector<double>& S1, maths::Vector<double>& S2);
            maths::Vector<double> schneidePlane(const Plane& E, bool& found);
        public:
            /**
            * @brief This method reflects the beam
            * @param[in] n surface normal (angles of incidence will be calculated with help of the directional vectors k)             
            * @param[in] n1 refractive index first medium
            * @param[in] n2 refractive index second medium
            * @return refracted ray
            */
            tubedRay reflect(maths::Vector<double>* n, std::complex<double>  n1, std::complex<double>  n2); 
            
            /**
            * @brief This method refracts the beam
            * @param[in] n array which holds surface normals of the edge rays (angles of incidence will be calculated with help of the directional vectors k)
            * @param[in] n1 refractive index first medium
            * @param[in] n2 refractive index second medium
            */
            void refract(maths::Vector<double>* n, std::complex<double>  n1, std::complex<double>  n2); 

            /**
            * @brief calculates the Fresnel matrix for the reflected beam (used in method reflect)
            * @param[in] alpha angle of incidence (in radiants)
            * @param[in] n1 refractive index first medium
            * @param[in] n2 refractive index second medium
            * @return Fresnel matrix
            */
            maths::Matrix<std::complex<double> > Fresnel_reflect(double alpha, std::complex<double>  n1, std::complex<double>  n2);

            /**
            * @brief calculates the Fresnel matrix for the refracted beam (used in method reflect)
            * @param[in] alpha angle of incidence (in radiants)
            * @param[in] n1 refractive index first medium
            * @param[in] n2 refractive index second medium
            * @return Fresnel matrix
            */
            maths::Matrix<std::complex<double> > Fresnel_trans(double alpha, std::complex<double>  beta, std::complex<double>  n1, std::complex<double>  n2);


            void initElectricField(const maths::Vector<std::complex<double> >& Pol, const int numRays, double dx, Plane Eb); 
            void initElectricFieldGauss(int i, const Plane& Eb, const maths::Vector<std::complex<double> >& Pol);

            void initElectricField(const maths::Vector<std::complex<double> >& Pol, int AnzRays = 1);
            maths::Vector<double> crossPlane(const maths::Vector<double> Pe, const maths::Vector<double> n);
            double cross(const maths::Vector<double> P10, const maths::Vector<double> P11,
                const maths::Vector<double> P20, const maths::Vector<double> P21);
            double crossXAxis(const maths::Vector<double>& P1, const maths::Vector<double>& P2, const maths::Vector<double>& k);
            maths::Vector <double> nextCaustic(double& l);

            tubedRay& operator = (const tubedRay& S);
            // project(int i){};
            maths::Vector<double> OK;
            maths::Vector<double> P[5];
            maths::Vector<double> k[5];
            maths::Vector<std::complex<double> > E[5];

            std::complex<double>  phi[5];
            //EinschlussInfo *Ein;
            int iR;
            int pol;
            bool getunnelt;
            bool logRay;
            bool isValid;
            Gauss g;
            maths::Vector<double> ka;
            double KORR;
        };
    }
}
