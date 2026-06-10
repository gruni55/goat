/***************************************************************************
                          objectshape.h  -  description
                             -------------------
    begin                : Wed Oct 24 2001
    copyright            : (C) 2001 by Thomas Weigel
    email                : weigel@lat.ruhr-uni-bochum.de
 ***************************************************************************/


 /**
  * @file This File is intended for internal use. e
  */
#pragma once

#ifndef INF
#define INF 1.0/0.0
#endif
#include "matrix.h"
#include "vector.h"

#include <fstream>
// #include <time.h>
namespace GOAT {
    namespace raytracing {
#define OBJECTSHAPE_NO_SHAPE    -1  //F/< No shape defined
#define OBJECTSHAPE_ELLIPSOID     10000  ///< Shape is an ellipsoid 
#define OBJECTSHAPE_SURFACE       10001 ///< Shape is triangulated surface
#define OBJECTSHAPE_CONE          10002 ///< Shape is a cone
#define OBJECTSHAPE_ASPHERIC_LENS 10003 ///< Shape is an aspheric lens
#define OBJECTSHAPE_SPHERIC_LENS  10004 ///< Shape is a spheric lens
#define OBJECTSHAPE_BOX           10005 ///< Shape is a box  
#define OBJECTSHAPE_CYLINDER      10006 ///< Shape is a cylinder
#define OBJECTSHAPE_VORTEX_PLATE  10007 ///< Shape is a vortex plate
#define OBJECTSHAPE_ROUGH_OBJECT     10008 ///< Shape is a rough object (surface roughness is considered)

#define FUNSURF          2  
#define SUPERELLIPSOID_D 17 
#define SUPERELLIPSOID   4
#define ZYLINDER     5
#define KREISKEGEL   6
#define KEGELSTUMPF  7
#define COMPOUND     8
#define SPIEGEL      12
#define SUPERELLIPSOID_N 10
#define ERYTHROCYTE      11
#define KEGELSTUMPF_HOHL 9 
#define HOHLFASER    13
#define NINCTYPES    14
#define LINSE	     15
#define ZYLINDER_HEXAGONAL 16
#define BOX         3
#define EPS 1E-10*r0

        enum class NFUNCTYPE { vacuum, air, bk7, silica, lasf5, pmma };

        /**
        * @brief Abstract base class for all volume objects
        * This abstract class provides a template for all volume objects. The refractive index is complex to be able to consider absorption.
        */
        class ObjectShape {
        public:
            ObjectShape();
            ObjectShape(const ObjectShape& F);


            /**
            @brief Constructor, as template for all derived classes.
            @param P position of the object (reference point)
            @param n refractive index (complex)
            @alpha polarizability matrix
            @param Ex, Ey, Ez direction of the object's coordinate system (default values: ex, ey and ez)
            @type  type of the object, defines the shape (default value: -1, no shape)
            */
            ObjectShape(const maths::Vector<double>& P,
                std::complex<double>  n,
                GOAT::maths::Matrix<std::complex<double> >  alpha,
                const maths::Vector<double>& Ex = maths::ex,
                const maths::Vector<double>& Ey = maths::ey,
                const maths::Vector<double>& Ez = maths::ez,
                const int type = -1
            );


            virtual void binWrite(std::ofstream& os) = 0;                    ///< binary writing to file   
            virtual void binRead(std::ifstream& os) = 0;                     ///< binary reading from file           
            virtual bool next(const maths::Vector<double>& p, const maths::Vector<double>& k,
                maths::Vector<double>& pout) = 0; ///< searches for the next (nearest) intersection of a ray with the object, p: current position of the ray, k: direction of the ray, pout: position of the crossing point. Returns true, if a crossing point was found.

            virtual maths::Vector<double> norm(const maths::Vector<double>& P) = 0;        ///< surface normal at the point P
            virtual bool isInside(const maths::Vector<double>& p) = 0;              ///< checks if point P is inside the object
            virtual double volume() = 0;                                      ///< returns the volume of the object



            // ----------- Setter ---------
            void scale(double sf);                                ///< sets scaling of the shape by the factor sf
            void setCenter(maths::Vector<double> P);                                ///< sets Center to P (check, if function is necessary)
            void setCenter2CoM(); ///< Calculates the center of mass (CoM) and sets the object's reference point to the CoM
           
            void setn(std::complex<double> n) { this->n = n; }                       ///< sets refractive index
            void setninel(std::complex<double> ninel) { this->ninel = ninel; }      ///< sets refractive index for inelastic (RRT) calculation
            void setActive(bool active) { Active = active; }                      ///< sets flag if the object is inelastic active, i.e. it will be considered for inelastic calculation   
            void setAlpha(double Alpha) { setMatrix(Alpha, Ebeta, Egamma); }   ///< sets rotation angle around x-axis
            void setBeta(double Beta) { setMatrix(Ealpha, Beta, Egamma); }     ///< sets rotation angle around y-axis 
            void setGamma(double Gamma) { setMatrix(Ealpha, Ebeta, Gamma); }   ///< sets rotation angle around z-axis
            void setPolMatrix(maths::Matrix<std::complex<double> >alpha) { this->alpha = alpha; }   ///< sets polarisability matrix
            void setMatrix(maths::Matrix<double> H);                                ///< sets the matrix for the transformation between object's coordinate system and outer coordinate system, H: transformation matrix
            void setMatrix(double alpha, double beta, double gamma);         ///< sets the matrix for the transformation between object's coordinate system and outer coordinate system, alpha, beta, gamma: angles (rotation around x-, y- and z-axis) to calculate transformation matrix
            void setVisible(bool visible) { this->visible = visible; } ///< set visiblity (used in GOATvis)
            void setID(const std::string& ID) { this->ID = ID; } ///< sets the ID string
			void setNfunc(std::function <std::complex<double>(double)> nfunc, NFUNCTYPE nfuncType = NFUNCTYPE::vacuum) { this->nfunc = nfunc; this->nfuncType = nfuncType; } ///< sets the function for the refractive index (used for inelastic (RRT) calculation)
			void setNfuncType(NFUNCTYPE nfuncType) { this->nfuncType = nfuncType; } ///< sets the type of the function for the refractive index (used for inelastic (RRT) calculation)
            virtual  void setr0(double r0) = 0;                                 ///< defines the radius of the calculation sphere
            virtual void setPos(maths::Vector<double> r) = 0; ///< sets reference point P 
            virtual void setPos(double x, double y, double z) = 0; ///< sets reference point P 

            // ----------- Getter ---------
            bool isOutsideWorld() const; ///< Test if bounding box is (partly) outside the calculation space
            int Type() const { return type; }                                         ///< returns the object's type
            std::complex<double> getninel() const { return ninel; }                      ///< returns refractive index
            std::complex<double> getn() const { return n; }                              ///< returns refractive index for inelastic (RRT) calculation
            bool getVisible() const { return visible; } ///< show the visiblity state (used in GOATvis)
			bool isVisible() const { return visible; } ///< show the visiblity state (used in GOATvis)
            double getr0() const { return r0; } ///< returns the radius of the calculation space stored in the object (only for internal use, it is set when the object is added to the scene)
			double getAlpha() const { return Ealpha; } ///< returns the rotation angle around x-axis
			double getBeta() const { return Ebeta; } ///< returns the rotation angle around y-axis
			double getGamma() const { return Egamma; } ///< returns the rotation angle around z-axis
            void getBBcorners(maths::Vector<double>& pul, maths::Vector<double>& por) const { pul = this->pul; por = this->por; } ///< returns the corners of the circumferent cuboid (bounding box)
		    maths::Vector<double> getBBoxMin() const { return pul; } ///< returns the corner of the circumferent cuboid with the lowest x,y and z coordinates
		    maths::Vector<double> getBBoxMax() const { return por; } ///< returns the corner of the circumferent cuboid with the highest x,y and z coordinates
            bool isActive() const { return Active; }                                 ///< returns true if the object should be considered for inelastic calculation
            std::string getID() const { return ID; } ///< returns the ID string
            maths::Vector<double>  getPos() const { return P; } ///< returns the position of the object
			maths::Matrix<double> getH() const { return H; } ///< returns the matrix for the transformation between object's coordinate system and outer coordinate system
			maths::Matrix<double> getR() const { return R; } ///< returns the matrix for the transformation back to the calculation system    
			double getScale() const { return sf; } ///< returns the scaling factor
			bool isRough() const { return rough; } ///< returns true if surface roughness is considered for the object
			NFUNCTYPE getFuncType() const { return nfuncType; } ///< returns the type of the function for the refractive index (used for inelastic (RRT) calculation)


            // ----------- Other functions ---------
            void rotate(maths::Vector<double> A, double phi); ///< sets the matrix for the transformation between object's coordinate system and outer coordinate system, A: rotation axis, phi: angle for rotation around A
            virtual void initQuad() = 0;                                      ///< calculates the circumferent cuboid (needed e.g. for the inelastic scattering calculations)
            virtual maths::Vector<double> calcCoM() = 0;     ///< calculates center of mass (needed by setCenter2CoM () )
            maths::Vector<double>& pos()
            {
                return P;
            }

            const maths::Vector<double>& pos() const
            {
                return P;
            }

            std::function <std::complex<double>(double)>& nFunc()
            {
                return nfunc;
            }

            const std::function <std::complex<double>(double)>& nFunc() const
            {
                return nfunc;
            }


        protected: 
            maths::Vector<double> P;                       ///< position of the object
            maths::Matrix<double> H, R;                     ///< matrices for the transformation in the local coordinate system (H) and back to the calculation system (R)
            std::complex<double>  n;                ///< refractive index of the object
            std::complex<double> ninel;             ///< refractive index of the object, used for inelastic (RRT) calculation
            maths::Matrix<std::complex<double> > alpha;    ///< polarisability matrix
            int type;                               ///< type of the object
            maths::Vector<double> pul, por;     ///< corners of the circumferent cuboid (lower left corner and upper right corner)
            maths::Vector<double> e[3];        ///< unity vectors, describing the directions of the local coordinate system
            double Ealpha, Ebeta, Egamma; ///< angles through which the object was rotated (around the x- (Ealpha), then the y- (Ebeta) and finally the z-axis (Egamma))
            double r0;                  ///< radius of the calculation sphere 
            double sf=1;         ///< scaling factor, it is used to scale the shape of the object     
            bool Active;   ///< should the object be considered for inelastic (RRT) calculations?
			bool rough=false; ///< is the surface of the object rough? (if true, surface roughness is considered in the calculations)
            double rho;        ///< mass density in \f$ kg/m^3 \f$
            /*
            * @brief Used in the visualization part (GOATvis) => if true object will be visualized
            * If this parameter is true, GOATvis will show the full representation of the object in the scene dialog otherwise only the bounding box is shown.
            * This parameter can be used e.g. for very heavy surface object to safe memory.
            */
            bool visible=true;
            std::function <std::complex<double>(double)>  nfunc;
            NFUNCTYPE nfuncType=NFUNCTYPE::vacuum;
			std::string ID = "Object"; ///< ID string, used for visualization (GOATvis) and for the user to identify the object
        };

        maths::Matrix<double> computeInertia(ObjectShape* F); ///< calculates inertia matrix
        bool intersectionTest(ObjectShape& A, ObjectShape& B); ///< Test if object A and object B may intersect each other (i.e. the bounding boxes around the objects intersect each other) 
    }
}

