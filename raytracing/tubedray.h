/***************************************************************************
                          tubed.h  -  description                              
                             -------------------                                         
    begin                : Fri Oct 15 1999                                           
    copyright            : (C) 1999 by Thomas Weigel                         
    email                : weigel@lat.ruhr-uni-bochum.de                                     
 ***************************************************************************/



#pragma once
#include <complex>
#include "fresnel.h"
#include "resutil.h"
#include "vector.h"
#include "matrix.h"
#include "plane.h"
#include "objectshape.h"


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
        Vector<double> P[5];  /// current position of the edge ray (index 0-3) and of the center ray (index 4)
        Vector<double> k[5];  /// current direction of the edge ray (index 0-3) and of the center ray (index 4)
        Vector<std::complex<double> > E[5]; /// current electric field of the edge ray (index 0-3) and of the center ray (index 4)
        int numObj; /// number of objects within the scene
        ObjectShape** Obj; /// list of objects in the scene
        double r0; /// radius of the calculation space
        Vector<std::complex<double> >*** Gitter;
    } tubedRayBuffer;



      /**
       * This structure is used for the gaussian beam to simplify the data transfer.
       */
    typedef struct
    {
        double w0;          /// (virtual) beam waist 
        Vector<double> F;   /// position of the beam
        Vector<double> k;   /// direction of the beam
        bool isGauss; /// is it really a gaussian beam
    } Gauss;


    void binWrite(Gauss gs, std::ofstream& os); /// write the Gauss structure to a binary file
    void binRead(Gauss& Gs, std::ifstream& is); /// read the Gauss structure from a binary file


    std::ostream& operator << (std::ostream& os, Gauss gs); /// output operator for the Gauss structure

#ifndef OHNE_GAUSS
#define OHNE_GAUSS
    const Gauss ohne_gauss = { 0,Vector<double>(0,0,0),Vector<double>(0,0,0),false };
#endif

#include "raybase.h"

    /**
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
        tubedRay(Plane E, const Vector<double>& p, double dy, double dz, const
            Vector<std::complex<double> >& Pol, std::complex<double>  n0, double r0, double
            k0, const int Anzein = 0, ObjectShape** Einschluss = NULL, bool logRay = false);
        tubedRay(const Vector<double>& p, double dy, double dz, const Vector<std::complex<double> >& Pol, const Vector<double>& K, std::complex<double>  n0, double r0, double k0, const int Anzein = 0, ObjectShape** Einschluss = NULL, bool logRay = false);
        void setGauss(Gauss g) { this->g = g; }
        bool next(); /// naechster Schritt
        Vector<std::complex<double> > getE() { return E[4]; }
        int objectIndex() { return objIndex; }
        void reflectRay(RayBase*& tray, Vector<double> n, std::complex<double> n1, std::complex<double> n2);
        void reflectRay(RayBase*& tray, Vector<double> n[5], std::complex<double> n1, std::complex<double> n2);
        tubedRay reflect(Vector<double> n, std::complex<double>  n1, std::complex<double>  n2); /// Strahl wird reflektiert
        void refract(Vector<double> n, std::complex<double>  n1, std::complex<double>  n2); /// Strahl wird gebrochen
        void tunnel(Vector<std::complex<double> > Pol, std::complex<double>  n1, std::complex<double>  n2);
        void tunnel(Vector<std::complex<double> > Pol, std::complex<double>  na, std::complex<double>  np, int l);
        double flaeche(); /// Flaeche zwischen den Randstrahlen
        void setRefract(std::complex<double>  n) { this->n = n; } /// Momentanen Brechungsindex setzen
        std::complex<double>  getRefract() { return n; } /// Gibt Brechungsindex zurueck
        void setk(int i, const Vector<double>& K) { k[i] = K; }
        void setphi(int i, const std::complex<double>& p) { phi[i] = p; }
        void setP(int i, const Vector<double>& p) { P[i] = p; }
        void setiR(int i) { iR = i; }
        void setGetunnelt(bool v) { getunnelt = v; }
        void setN0(std::complex<double> n) { n0 = n; }
        Vector<double> getk(int i) { return k[i]; }
        Vector<double> getk() { return k[4]; }
        Vector<double> getP(int i) { return P[i]; }
        Vector<double> getP() { return P[4]; }
        std::complex<double>  getphi(int i) { return phi[i]; }
        void getP(Vector<double> p[5]) { p = P; }
        Vector<std::complex<double> > getAmp(int i) { return E[i]; }
        void setAmp(int i, const Vector<std::complex<double> >& A) { E[i] = A; }
        ObjectShape* Einschluss(int i) { return Obj[i]; }
        ~tubedRay();
        bool Einschluss() { return inObject; }
        bool isInObject() { return inObject; }
        //  Vector<double> checkEinschluss(int s, const Vector<double>& Ps,int& Index);
        void checkObjectIntersection(int Index[5], Vector<double>* Pmin);

        bool Getunnelt() { return getunnelt; }
        int currentObjectIndex() { return objIndex; }
        double Intensity(int i) { return abs(E[i]); }
        friend std::ostream& operator << (std::ostream& os, tubedRay S);
        int reflections() { return iR; }
        double pjump(Vector<double> P1[5], Vector<double> P2[5], const double epsilon = 1E-10);
        double pjump(Vector<double> P1[5], Vector<double> P2[5], Vector<double>* S);
        double pjump(void);
        double normVol(Vector<double> P[5], Vector<double> k, Vector<double> n);
        bool schneidePlane(Vector<double>* Erg, const Plane& E);
        int schneidePlane(const Plane& E, double d, Vector<double>& S1, Vector<double>& S2);
        Vector<double> schneidePlane(const Plane& E, bool& found);
    public:
        /*!
        \param n Oberflaechennormale (in Richtung gegen den Strahl)
        \param n1 Brechungsindex Medium 1
        \param n2 Brechungsindex Medium 2
        Rueckgabewert: gebrochener Strahl
        */
        tubedRay reflect(Vector<double> n[5], std::complex<double>  n1, std::complex<double>  n2); /// Reflexion des Strahls
        /*!
        \param n Oberflaechennormale (in Richtung gegen den Strahl)
        \param n1 Brechungsindex Medium 1
        \param n2 Brechungsindex Medium 2
        */
        void refract(Vector<double> N[5], std::complex<double>  n1, std::complex<double>  n2); /// Brechung des Strahls


        Matrix<std::complex<double> > Fresnel_reflect(double alpha, std::complex<double>  n1, std::complex<double>  n2);
        Matrix<std::complex<double> > Fresnel_trans(double alpha, std::complex<double>  beta, std::complex<double>  n1, std::complex<double>  n2);
        void initElectricField(const Vector<std::complex<double> >& Pol, const int AnzRays, double dx, Plane Eb);
        void initElectricFieldGauss(int i, const Plane& Eb,
            const Vector<std::complex<double> >& Pol);

        void initElectricField(const Vector<std::complex<double> >& Pol, int AnzRays = 1);
        Vector<double> crossPlane(const Vector<double> Pe, const Vector<double> n);
        double cross(const Vector<double> P10, const Vector<double> P11,
            const Vector<double> P20, const Vector<double> P21);
        double crossXAxis(const Vector<double>& P1, const Vector<double>& P2, const Vector<double>& k);
        Vector <double> nextCaustic(double& l);

        tubedRay& operator = (const tubedRay& S);
        // project(int i){};
        Vector<double> OK;
        Vector<double> P[5];
        Vector<double> k[5];
        Vector<std::complex<double> > E[5];

        std::complex<double>  phi[5];
        //EinschlussInfo *Ein;
        int iR;
        int pol;
        bool getunnelt;
        bool logRay;
        bool isValid;
        Gauss g;
        Vector<double> ka;
        double KORR;
    };


