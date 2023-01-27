#ifndef Lightsrc_H
#define Lightsrc_H

#include "raybase.h"
#include "vector.h"
#include "tubedray.h"
#include "iray.h"
#include "objectshape.h"
#include "ray_pow.h"

namespace GOAT
{
	namespace raytracing
	{
#define LIGHTSRC_RAYTYPE_RAY 1    ///< Ray class : tubedRay
#define LIGHTSRC_RAYTYPE_IRAY 2   ///< Ray class : IRay
#define LIGHTSRC_RAYTYPE_PRAY 3   ///< Ray class : Pow_Ray
constexpr int LIGHTSRC_SRCTYPE_PLANE=1;  ///< Light source is a plane wave
constexpr int LIGHTSRC_SRCTYPE_GAUSS=2;  ///< Light source is a gaussian wave
constexpr int LIGHTSRC_SRCTYPE_TOPHAT=3; ///< Light source is a top hat


#define LIGHTSRC_NOT_LAST_RAY 0  ///< Created ray is not the last ray 
#define LIGHTSRC_IS_LAST_RAY 1   ///< Created ray is the last ray 
#define LIGHTSRC_ERROR -1        ///< Error occurs within the ray creation
#define Z0 376.730313461         ///< Wave impedance of free space

#define LIGHTSRC_POL_X 0
#define LIGHTSRC_POL_Y 1
#define LIGHTSRC_POL_Z 2
#define LIGHTSRC_POL_USER_DEFINED 0 




		/**
		* @brief This abstract class is the basic class for all light sources used in this raytracing library. It provides all necessary interfaces.
		*
		* LightSrc: the virtual base class for describing light sources. Each light source is defined by a quadratic area as
		* starting area for the corresponding rays and a reference point Pos in the middle of this area. The direction has to be
		* defined separately in the derived class. It provides a standard programming interface for light sources. The virtual function next(...)
		* gives the next ray with the light source dependent parameters (position, direction, strength of the electric field etc.).
		* The return parameter has three possible values: LIGHTSRC_NOT_LAST_RAY: the created ray is not the last ray,
		* LIGHTSRC_IS_LAST_RAY: the last ray was created (i.e. the raytracing is at its end) and LIGHTSRC_ERROR: an error occurs.
		* The routine reset() resets the ray counter to its initial value, so the next call of the routine next() gives the first
		* ray. Since the rays need information about the objects to be able to calculate the crossing points with the object surface,
		* LightSrc provides routines to pass these informations to the rays. These routines are only needed for calculations without
		* the raytracing process, especially without the Scene class.
		*/
		class LightSrc

		{
		public:
			void reset();  ///< Reset everything (counting starts from the beginning)
			~LightSrc(void);
			LightSrc();
			LightSrc(const LightSrc&); ///< Copy constructor
			void clearObjects(); ///< clear object list
			void addObject(ObjectShape* obj);  ///< add single object to the object list
			void ObjectList(int Anz, ObjectShape** Obj);  ///< import object list
			/// 

			/**
			  * These functions return the next ray of the light source for further calculations
			  * @param ray next ray
			  * @retval LIGHTSRC_NOT_LAST_RAY, if ray is not the last ray,  LIGHTSRC_IS_LAST_RAY, if the ray is the last ray or LIGHTSRC_ERROR if an error occured
			  */
			  ///@{
			virtual int next(RayBase* ray) = 0; // gives back next ray
			virtual int next(Ray_pow& ray) = 0; // gives back next ray
			virtual int next(IRay& ray) = 0; // gives back next ray
			virtual int next(tubedRay& ray) = 0;  // gives back next ray
			///@}
			void binRead(std::ifstream& is); ///< writes content of LightSrc in a binary file, represented by is
			void binWrite(std::ofstream& os); ///< reads content of LightSrc from  a binary file, represented by os
			virtual void binWriteItem(std::ofstream& os) = 0; ///< writes content of LightSrc in a binary file, represented by is (has to be specified by the derived classes)focuspos
			virtual void binReadItem(std::ifstream& os) = 0; ///< reads content of LightSrc from a binary file, represented by os (has to be specified by the derived classes)
			int getNumObjs() { return numObjs; } ///< returns the number of Objects (only needed, when used seperately outside scene)
			ObjectShape* getObject(int i) { if ((i < 0) || (i > numObjs)) return NULL; return Obj[i]; }  ///< returns i-th item in the object list
			void setObject(ObjectShape* O, int i = -1); ///< exchanges i-th object in the object list
			int rayType() { return raytype; }  ///< returns the current ray type
			void setR0(double r0); ///< sets the radius of the calculation sphere
			double getDensity() { return density; } ///< returns the ray density, i.e. the number of rays per unit length (=D/N)
			void setD(double D) ///< sets the width of the light source
			{
				this->density = D / ((double)N);
				this->D = D;
				reset();
			} ///< sets the width of the light source (this resets also the ray counter)
			void setk(const maths::Vector<double>& k); ///< sets the main direction of the light source
			maths::Vector<double> getk() { return k; } ///< returns the main direction of the light source
			int getNumRays() { return N; } ///< returns the number of rays (per direction in space)
			void setNumRays(int N) ///< sets the number of rays (per direction in space)
			{
				this->N = N;
				density = D / ((double)N);
				reset();
			}
			void setPol(maths::Vector<std::complex<double> > pol) { Pol = pol; } ///< sets the polarisation 
			void setPos(maths::Vector<double> P); ///< sets the position of the light source. This is the center of the square area of the light source
			maths::Vector<double> getPos() { return Pos; }  ///< returns the position of the light source. This is the center of the square area of the light source
			void setN0(std::complex<double> n0) { this->n0 = n0; } ///< sets the complex valued refractive index of the intermediate medium
			ObjectShape** Obj;         ///< list of all objects
			// protected :
			maths::Vector<double> Pos; ///< position of the light source (center of the square area of the light source)
			int type;           ///< type of the light source
			double P0=1.0;        ///< power
			double density;     ///< ray density, i.e. distance between two neighboring rays
			maths::Vector<double> k;   ///< main direction of the light source   
			int N;  ///< number of rays (per direction)
			int i1; ///< first index of the ray inside the starting area (for internal use, -1 if the calculation has not yet been started)
			int	i2; ///< second index of the ray inside the starting area (for internal use, -1 if the calculation has not yet been started)
			maths::Vector<std::complex<double> > Pol; ///< polarisation (default: (0.0, 1.0, 0.0)
			maths::Vector<std::complex<double> > Pol2; ///< second polarisation (used by IRay)
			double r0;          ///< radius of the calculation sphere
			double wvl;         ///< wavelength
			int numObjs;        ///< number of objects
			std::complex<double> n0; ///< refractive index of the intermediate medium

			double D;           ///< width of the square light source area 
			maths::Vector<double> e1, e2;  ///< unit vectors that span the light source area 
			int raytype;        ///< Strahltyp : ray oder ISTRAHL (=RAY oder IRAY)
			int polType;        ///< Polarisationsrichtung (s.o.)  
			double Pall;        ///< overall Power  
			friend class LightSrcPlane;
			friend class LightSrcGauss;
			friend std::ostream& operator << (std::ostream& os, LightSrc* ls);
			bool suppress_phase_progress = false;
		};




		/**
		* @brief class derived from LightSrc. It represents a plane wave described by the electric field \f$ \vec{E}(\vec{P})=\vec{E}_0 \cdot e^{i\vec{k}\cdot\Delta\vec{P}} \f$
		*/
		class LightSrcPlane : public LightSrc

		{
		public:
			LightSrcPlane(void);
			/**
			 * Main constructor
			 * @param Pos Center of the start area
			 * @param N Number of rays (per direction, total number of rays:N*N)
			 * @param wvl Wavelength
			 * @param D Edge length of the quadratic start area (=size of the light source)
			 * @param Pol direction of the polarisation
			 */
			LightSrcPlane(maths::Vector<double> Pos, int N, double wvl, double D = 100.0, maths::Vector<std::complex<double> > Pol = maths::Vector<std::complex<double> >(0.0, 1.0, 0.0), int raytype = LIGHTSRC_RAYTYPE_IRAY, double r0 = 100.0); ///<constructor with the following parameters : Pos: position of the light source, k : main direction of the light source, N : number of rays per direction(total number of all rays is N* N)
			LightSrcPlane(const LightSrcPlane&); ///< copy Constructor
			~LightSrcPlane(void) {}; ///< Destructor, not needed
			int next(GOAT::raytracing::RayBase* ray);  ///< gives the next ray for the following calculations 
			int next(IRay& S);///< gives the next ray for the following calculations 
			int next(tubedRay& S);///< gives the next ray for the following calculations 
			int next(Ray_pow& S);///< gives the next ray for the following calculations 
			void binWriteItem(std::ofstream& os) { /* to be implemented !!! */ }
			void binReadItem(std::ifstream& os) { /* to be implemented !!! */ }


			// void turnSrc // to be done !!!
		};

		/**
		* @brief This class describes a focused gaussian beam.
		*
		 Class which describes a (focused) gaussian light source. The main direction is given by the source position and
		 the focal position. The electric field is calculated by
		 \f$\vec E(r,z)=\vec E_0 \frac{w_0}{w(z)}\cdot e^{\frac{r^2}{w^2(z)}}\cdot e^{-ik\frac{r^2}{2R(z)}}\cdot e^{i(\zeta(z)-kz)}\f$
		 The waist of the beam is only used for the correct electric field distribution within the starting area. Since we are working with geometrical optics
		 the rays follow straight lines inside the medium.
		*/
		class LightSrcGauss : public LightSrc
		{
		public:
			void binWriteItem(std::ofstream& os);
			void binReadItem(std::ifstream& is);
			int next(GOAT::raytracing::RayBase* ray);
			int next(Ray_pow& S);
			int next(IRay& S);
			int next(tubedRay& S);
			// void initElectricFieldGauss(GOAT::maths::Vector<double> &P, GOAT::maths::Vector<double> &k, GOAT::maths::Vector<std::complex<double> > &E); ///< Initialisierung eines Teilstrahls
			LightSrcGauss(void);
			LightSrcGauss(const LightSrcGauss&);
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
			LightSrcGauss(maths::Vector<double> Pos, int N, double wvl, double w0, maths::Vector<double> focuspos, double D = 1.0, maths::Vector<std::complex<double> > Pol = maths::Vector<std::complex<double> >(0.0, 1.0, 0.0), int raytype = LIGHTSRC_RAYTYPE_IRAY, double r0 = 1.0);
			void setW0(double w0) { this->w0 = w0; calcz0(); reset(); } ///< sets the beam's waist 
			maths::Vector<double> getFocuspos() { return focuspos; } ///< returns the position of the focus
			void setFocuspos(maths::Vector<double> fp) { focuspos = fp; f = abs(Pos - focuspos); reset(); } ///< sets the focus position to fp
			void setPos(maths::Vector<double> pos) { Pos = pos; f = abs(Pos - focuspos); reset(); } ///< sets the position of the light source to pos

			void setNA(double na); ///< sets the numerical aperture NA and recalculates the width D, the focal beam waist w0 and the Rayleigh-length z0
			void setWvl(double wvl); ///< sets the vacuum wavelength
			void setk(maths::Vector<double> k) { this->k = k; reset(); }	///< sets the main direction of light source 
			double calcz0() { z0 = M_PI * w0 * w0 / wvl; return z0; } ///< recalculates Rayleigh-length z0
			void reset()
			{
				k = focuspos - Pos;
				k = k / abs(k);
				e1 = k % maths::ez;
				if (abs(e1) < 1E-10) e1 = maths::ex;
				e2 = k % e1;
				e1 = e1 / abs(e1);
				e2 = e2 / abs(e2);

				i1 = 0;
				i2 = 0;
				calcNormfak();
			}
			double calcw(double z) ///< calculates the beam waist of the light beam at the distance z from the focal point, returns the value and sets the corrsponding local variable w (needed for next(), only for internal use)
			{
				w = w0 * sqrt(1.0 + z * z / (z0 * z0));
				return w;
			}

			void calcNormfak() ///< needed for next()
			{
				double l = abs(Pos - focuspos);
				calcz0();
				calcw(l);
				Normfak = Z0 * sqrt(2.0 / M_PI) / (n0 * l * M_PI) * w / w0 / w0 * exp(-2.0 * l * l / w / w);
				Normfak = Normfak / ((double)N * N) * l * l;
			}
			//protected : 
			double w0;        ///< Waist diameter (fictitious !), only used for the correct calculation of the electric field distribution within the starting area
			double f;         ///< distance between light source area and the focal point 
			double getNA() { return NA; } ///< returns numerical aperture

			std::complex<double> Normfak;
			maths::Vector<double> focuspos; ///< focal position	
			double z0;				 ///< Rayleigh length (for internal use only)
			maths::Vector<double> k;		 ///< direction of the gaussian beam
			double w;				 ///< radius of the beam (for internal use only)
			double NA;				 ///< numerical aperture (normalized by the intermediate refractive index)
		};

		/**
		 * @brief Writes a list of light sources into a binary file.
		 * @param os Output file stream to write the data
		 * @param nLS Number of light sources to write
		 * @param ls List of the light sources
		 */
		void binWriteLSList(std::ofstream& os, int nLS, LightSrc** ls);
		/**
		 * @brief Reads list of light sources from a binary file and allocates the memory for the light source list.
		 * @param is Input file stream to read the data
		 * @param nLS Number of light sources to read
		 * @param ls List of the light sources read from the input file
		 */
		void binReadLSList(std::ifstream& is, int nLS, LightSrc**& ls);
		/**
		 * @brief Copy the given light source list.
		 * @param d Destination light source list (memory will be allocated within this function)
		 * @param s Source light source list (this list will be copied)
		 * @param nLS Number of light sources to copy
		 */
		void copyLightSrcList(LightSrc**& d, LightSrc** s, int nLS);
	}
}
#endif
