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

#include "vector.h"
#include "matrix.h"
#include "resutil.h"
#include <fstream>
// #include <time.h>


#define OBJECTSHAPE_NO_SHAPE    -1  //F/< No shape defined
#define OBJECTSHAPE_ELLIPSOID    0  ///< Shape is ellipsoid 
#define OBJECTSHAPE_SURFACE      1 ///< Shape is triangulated surface
#define FUNSURF      2  
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
#define OBJECTSHAPE_BOX         3
#define EPS 1E-10*r0


/**  
* @brief Abstract base class for all volume objects 
* This abstract class provides a template for all volume objects. The refractive index is complex to be able to consider absorption. 
*/
class ObjectShape {
public:
	ObjectShape();
  ObjectShape (const ObjectShape &F);
 
  
  /**   
  @brief Constructor, as template for all derived classes. 
  @param P position of the object (reference point)
  @param n refractive index (complex)
  @alpha polarisabilty matrix
  optional parameters: 
  @param Ex, Ey, Ez direction of the object's coordinate system (default values: ex, ey and ez)
  @type  type of the object, defines the shape (default value: -1, no shape)
  */
  ObjectShape (const Vector<double> &P,                                 
        std::complex<double>  n,
        Matrix <std::complex<double> > alpha,
        const Vector<double> &Ex=ex,
        const Vector<double> &Ey=ey,
        const Vector<double> &Ez=ez,

        const int type=-1
       );


  virtual void binWrite (std::ofstream &os) = 0;                    ///< binary writing to file   
  virtual void binRead (std::ifstream &os) = 0;                     ///< binary reading from file
  virtual void scale(double sf) = 0;                                ///< sets scaling of the shape by the factor sf
  virtual bool next(const Vector<double> &p, const Vector<double> &k,
                     Vector<double> &pout) = 0; ///< searches for the next (nearest) intersection of a ray with the object, p: current position of the ray, k: direction of the ray, pout: position of the crossing point. Returns true, if a crossing point was found.
  
  virtual Vector<double> norm (const Vector<double> &P) = 0;        ///< surface normal at the point P
  virtual bool isInside (const Vector<double> &p) = 0;              ///< checks if point P is inside the object
  virtual double volume() = 0;                                      ///< returns the volume of the object
  
  void setCenter(Vector<double> P);                                ///< sets Center to P (check, if function is necessary)
  void setCenter2CoM(); ///< Calculates the center of mass (CoM) and sets the object's reference point to the CoM
 
  int Type() {return type;}                                         ///< returns the object's type
  virtual void initQuad ()= 0;                                      ///< calculates the circumferent cuboid (needed e.g. for the inelastic scattering calculations)
  virtual  void setr0(double r0)=0;                                 ///< defines the radius of the calculation sphere
  void setMatrix (Matrix<double> H);                                ///< sets the matrix for the transformation between object's coordinate system and outer coordinate system, H: transformation matrix
  void setMatrix (double alpha, double beta, double gamma);         ///< sets the matrix for the transformation between object's coordinate system and outer coordinate system, alpha, beta, gamma: angles (rotation around x-, y- and z-axis) to calculate transformation matrix
  void rotate (Vector<double> A, double phi); ///< sets the matrix for the transformation between object's coordinate system and outer coordinate system, A: rotation axis, phi: angle for rotation around A
  
  virtual void setP (Vector<double> r) = 0; ///< sets reference point P 
  virtual void setP (double x, double y, double z)= 0; ///< sets reference point P 
  virtual Vector<double> calcCoM() = 0;     ///< calculates center of mass (needed by setCenter2CoM () )
  void setn (std::complex<double> n) {this->n=n; }                       ///< sets refractive index
  void setninel (std::complex<double> ninel)  {this->ninel=ninel; }      ///< sets refractive index for inelastic (RRT) calculation
  std::complex<double> getninel () {return ninel; }                      ///< returns refractive index
  std::complex<double> getn () {return n; }                              ///< returns refractive index for inelastic (RRT) calculation
  void setPolMatrix (Matrix<std::complex<double> >alpha) {this->alpha=alpha;}   ///< sets polarisability matrix
  bool isActive () { return inelactive;}                                 ///< returns true if the object should be considered for inelastic calculation
  void setActive (bool active) {inelactive=active;}                      ///< sets flag if the object is inelastic active, i.e. it will be considered for inelastic calculation   
  void setAlpha(double Alpha) {setMatrix( Alpha,Ebeta,Egamma); }   ///< sets rotation angle around x-axis
  void setBeta(double Beta) {setMatrix( Ealpha,Beta,Egamma); }     ///< sets rotation angle around y-axis 
  void setGamma(double Gamma) {setMatrix( Ealpha,Ebeta,Gamma); }   ///< sets rotation angle around z-axis
  Vector<double> P;                       ///< position of the object
  Matrix<double> H,R;                     ///< matrices for the transformation in the local coordinate system (H) and back to the calculation system (R)
  std::complex<double>  n;                ///< refractive index of the object
  std::complex<double> ninel;             ///< refractive index of the object, used for inelastic (RRT) calculation
  Matrix<std::complex<double> > alpha;    ///< polarisability matrix
  int type;                               ///< type of the object
  Vector<double> pul,por;     ///< corners of the circumferent cuboid (lower left corner and upper right corner)
  Vector<double> e[3];        ///< unity vectors, describing the directions of the local coordinate system
  double Ealpha,Ebeta,Egamma; ///< angles through which the object was rotated (around the x- (Ealpha), then the y- (Ebeta) and finally the z-axis (Egamma))
  double r0;                  ///< radius of the calculation sphere 
  double sf;         ///< scaling factor, it is used to scale the shape of the object     
  bool inelactive;   ///< should the object be considered for inelastic (RRT) calculations?
  double rho;        ///< mass density in \f$ kg/m^3 \f$
};

Matrix<double> computeInertia(ObjectShape *F); ///< calculates inertia matrix



