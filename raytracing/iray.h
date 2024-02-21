/***************************************************************************
                          iray.h  -  description                              
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
#include "tubedray.h"
#include "misc.h"
#include "raybase.h"


 
namespace GOAT
{
	namespace raytracing
	{
		/**
		* @brief This class represents a single ray.
		*
		* This ray is intended to calculate with a single ray with two different polarisations. This class should not be used directly.
		* It is intended for use together with classes derived by LightSrc (like LightSrcPlane or LightSrcGauss)
		*/

		class IRay : public RayBase {
		public:
			IRay();
			/**
			 * @brief Contructor
			 *
			 * \param p Position of the ray
			 * \param Pol Polarisation (in the moment not used)
			 * \param K direction of the ray
			 * \param n0 refractive index
			 * \param r0 radius of the calculation sphere
			 * \param k0 wavenumber (=2*pi/wavelength in vacuum)
			 * \param numObj number of objects
			 * \param obj list of the objects
			 */
			IRay(const maths::Vector<double>& p,
				const maths::Vector<std::complex<double> >& Pol, const maths::Vector<double>& K,
				std::complex<double>  n0, double r0, double k0,
				const int numObj = 0, ObjectShape** obj = NULL);
			IRay(const IRay& r)
			{
				this->E1 = r.E1;
				this->E2 = r.E2;

				this->numObj = r.numObj;
				this->Obj = r.Obj;
				this->objIndex = r.objIndex;
				this->inObject = r.inObject;
				this->iR = r.iR;
				this->isValid = r.isValid;
				this->k = r.k;
				this->k0 = r.k0;
				this->KORR = r.KORR;
				this->n = r.n;
				this->OK = r.OK;
				this->P = r.P;
				this->r0 = r.r0;
				this->status = r.status;
				this->suppress_phase_progress = r.suppress_phase_progress;
			}


			bool next();

			maths::Vector<std::complex<double> > getE() { return E1; }
			/**
			 * @brief This function reflects ray on the surface.
			 * Here, the ray is reflected and the transmitted ray is created.
			 * \param n Surface normal
			 * \param n1 Refractive index (incident side)
			 * \param n2 Refractive index (transmitted side)
			 * \return Transmitted ray
			 */
			IRay reflect(maths::Vector<double> n, std::complex<double>  n1, std::complex<double>  n2);
			/**
			 * @brief This function refracts the ray on the surface.
			 *
			 *
			 * \param N Surface normal
			 * \param n1 Refractive index (incident side)
			 * \param n2 Refractive index (transmitted side)
			 */
			void refract(maths::Vector<double> N, std::complex<double>  n1, std::complex<double>  n2);
			void setRefract(std::complex<double>  n) { this->n = n; }  ///< Sets the current refractive index
			std::complex<double>  getRefract() { return n; } ///< Returns the current refractive index
			void setk(const maths::Vector<double>& K) { k = K; } ///< Sets the direction of the ray 
			void setP(const maths::Vector<double>& p) { P = p; } ///< Current position of the ray
			void setiR(int i) { iR = i; } ///< Sets the current reflexion counter (for internal use only)
			// void setGetunnelt (bool v) { getunnelt=v; } 
			maths::Vector<double> getk() { return k; } ///< Returns the direction of the ray
			maths::Vector<double> getP() { return P; } ///< Returns the current position of the ray
			int objectIndex() { return objIndex; } ///< Returns the index of the last hidden object (or -1 if no object was hidden)
			/**
			 * @brief Reflects ray on the surface
			 * /param tray Transmitted (=refracted) ray
			 * /param n Surface normal
			 * /param n1 Refractive index (incident side)
			 * /param n2 Refractive index (transmitted side)
			 */
			void reflectRay(RayBase*& tray, maths::Vector<double> n, std::complex<double> n1, std::complex<double> n2);
			ObjectShape* getObject(int i) { return Obj[i]; } ///< Returns the i-th object (for internal use only) 
			~IRay();
			bool isInObject() { return inObject; } ///< Returns true if ray is inside an object
			/**
			 * Checks wether there is an intersection between the ray and an object or not.
			 * In case of an intersection the function returns true and Index is set to the index of the object within the object list and Pmin gives the (nearest) intersection point
			 */
			bool checkObjectIntersection(int& Index, maths::Vector<double>& Pmin);
			int currentObjectIndex() { return objIndex; } ///< returns the index of the hidden object or -1
			friend std::ostream& operator << (std::ostream& os, IRay S);
			int reflections() { return iR; } ///< Returns the number of reflections the beam has already passed through

		public:
			maths::Vector<std::complex<double> > Pol1() { return E1 / abs(E1); } ///< direction of the first polarisation
			maths::Vector<std::complex<double> > Pol2() { return E2 / abs(E2); } ///< direction of the second polarisation

			maths::Matrix<std::complex<double> > Fresnel_reflect(double alpha, std::complex<double>  n1, std::complex<double>  n2); ///< returns Fresnel matrix for the reflection calculation (alpha: angle of incidence, n1,n2 are the refractive indices)
			maths::Matrix<std::complex<double> > Fresnel_trans(double alpha, std::complex<double>  beta, std::complex<double>  n1, std::complex<double>  n2); ///< returns Fresnel matrix for the transmission calculation (alpha: angle of incidence, n1,n2 are the refractive indices)
			void initElectricField(const Plane& Eb, const maths::Vector<std::complex<double> >& Pol, const int numOfRays = 1); ///< initialises the electric field with help of the Plane Eb and the polarisation Pol, numOfRays can be omitted. (Here, only one polarisation is given, the other is skipped)
			void initElectricField(const maths::Vector<std::complex<double> >& PolS, const maths::Vector<std::complex<double> >& PolP, const int AnzRays); ///< initialises the electric field vectors with help of polarisation vectors PolS and PolP, numOfRays can be omitted
			void initElectricField(const Plane& Eb, const maths::Vector<std::complex<double> >& Pol1, const maths::Vector<std::complex<double> >& Pol2, const int AnzRays);
			/**
			 * Initialises the ray for a gaussian beam according to
			 *
			 *
			 *
			 * \f$ \vec E=\vec{Pol} \cdot e^{-ik_0h} \cdot e^{-(r^2-h^2)/\sigma_2} \f$ where \f$ r=(\vec{F}-\vec{P})\cdot\vec{k},~h=r-|\vec{F}-\vec{P}| \f$ (\f$\vec{F}\f$: focal point)
			 *
			 *  E1 and E2 are set to the same value !
			 */
			void initElectricFieldGauss(double sigma2, maths::Vector<double> focuspos, maths::Vector<std::complex<double> > Pol);
			void initElectricFieldGauss(const Plane& Eb,
				const maths::Vector<std::complex<double> >& PolS,
				const maths::Vector<std::complex<double> >& PolP,
				Gauss g); ///< Initialises the ray for gaussian beam with help of the plane Eb and the polarisation vectors PolS and PolP. The parameters for the gaussian beam are stored in g

			void initElectricFieldGauss(maths::Vector<std::complex<double> >& Pol, Gauss g); ///< Initialises the ray for gaussian beam with help of the plane Eb and the polarisation vector Pol (second polarisation is not used, E2 is set to E1). The parameters for the gaussian beam are stored in g
			double cross(const maths::Vector<double> P10, const maths::Vector<double> P11, const maths::Vector<double> P20, const maths::Vector<double> P21);

			maths::Vector<double> crossPlane(const maths::Vector<double> Pe, const maths::Vector<double> n); ///< calculates the intersection point between the ray and a plane, which is defined by the vector P and the normal n
			/**
			 * Calculates the intersection point between the ray and a rectangle defined by the edge point P and the vectors e1 and e2. If there is no intersection point, a NaN valued vector is returned.
			 */
			maths::Vector<double> intersectRect(const maths::Vector<double> P, const maths::Vector<double> e1, const maths::Vector<double> e2);



			maths::Vector<double> P; ///< current position 
			maths::Vector<double> k; ///< current direction
			maths::Vector<std::complex<double> > E1, E2; ///< electric fields 

		protected:
			maths::Vector<double> OK;
			double KORR;
			bool isValid;
		};
	}
}
