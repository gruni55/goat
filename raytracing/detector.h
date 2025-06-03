#ifndef DETECTOR_H
#define DETECTOR_H

#include "vector.h"
#include "box.h"

#include <iostream>
#include <vector>

namespace GOAT
{
	namespace raytracing
	{
		#define DETECTOR_PLANE 20000
		#define DETECTOR_ANGLE 20001


		/**
		 * @brief The abstract Detector class provides an interface to a detector to store the information about the electric field into any kind of an array.
		 */
	class Detector
	{
	public:
		Detector(void);
		Detector(int n1, int n2); ///< Constructor which is called by the size of the array (n1 x n2)
		Detector(const Detector& Det); ///< copy constructor
		Detector& operator= (const Detector& D);
		maths::Vector<std::complex<double> >& operator () (int i1, int i2) { return D[i1][i2]; } ///< bracket operator which gives the content inside the array determined by the indices i1 and i2
		~Detector(void);

		void clean(); ///< Clean all data, i.e. all data is set to zero (Data array is not removed)
		void clear(); ///< clear 
		/**
		 * The virtual function checks if the ray hits the detector.
		 * @return true, if ray hits the detector otherwise false
		 * @param[in] P current position of the ray
		 * @param[in] k direction of the ray
		 * @param[out] i1 first index of the hit point within the detector array
		 * @param[out] i2 second index of the hit point within the detector array
		 * @param[out] d distance between ray and detector (along the ray)
		 */
		virtual bool cross(maths::Vector<double> P, maths::Vector<double> k, int& i1, int& i2, double& l) = 0;
		int N1(); ///< returns the dimension of the array in the first direction
		int N2(); ///< returns the dimension of the array in the second direction
		double D1(); ///< return the length in the first direction
		double D2(); ///< return the length in the second direction

		int Type() { return type; } ///< returns kind of detector
		/**
		* @load detector data from file
		* This function loads detector data from file. It expects a complex vector format. The first line has to start with "%n1 " followed by the number
		* of items in e1-direction. The second row must start with "%n2" followed by the number of items in e2-direction.  
		*/
		bool load(const char* fn); 
		void save(const char* fn); ///< stores the content (the whole vector)
		void saveabs(const char* fn);  ///< stores the content (absolute value of the electric field) of the detector array in the file determined by its filename fn 
		void savePhase(const char* fn, int coord);  ///< stores the content (phase of one component of the electric field, coord determines the coordinate 0,1,2 for x,y,z) of the detector array in the file fn 
		void savereal(const char* fn, int coord); ///< stores the content (real part of one component of the the electric field, coord determines the coordinate 0,1,2 for x,y,z) of the detector array in the file fn 
		void saveimag(const char* fn, int coord); ///< stores the content (imaginary part of one component of the the electric field, coord determines the coordinate 0,1,2 for x,y,z) of the detector array in the file fn 
			maths::Vector<std::complex<double> >** D; ///< Here, the data will be stored
			maths::Vector<double> position() { return P; } ///< returns the position of the detector 
			maths::Vector<double> norm() { return n; } ///< returns the normal vector of the detectors surface
		friend std::ostream& operator << (std::ostream& os, Detector& D);
		std::string fname;
		/**
		 * @brief Multiply with factor.
		 * This functions multiplies all elements of the detector with the factor fac
		 */
		void mult(double fac); 
		maths::Vector<double> gete1() { return e1; } ///< returns the direction of the first axis of the detector
		maths::Vector<double> gete2() { return e2; } ///< returns the direction of the second axis of the detector
	protected:
		maths::Vector<double> e1; ///< unit vector in the first direction 
		maths::Vector<double> e2; ///< unit vector in the second direction
		maths::Vector<double> P; ///< Position of the detector
		maths::Vector<double> n;  
		void init(int n1, int n2); ///< initialise array (for internal use only)
		double d1=0, d2=0;
		int n1=0, n2=0;
		int type=-1;
		friend class DetectorPlane;
	};


	/**
	 * @brief This class provides a plane detector, defined by its center point and its side length.
	 * The detector plane is spanned by the vectors e1 and e2. The absolute value of the vectors
	 * e1 and e2 are the cell widths in the corresponding direction.
	 */
	class DetectorPlane : public Detector
	{
	public:
		DetectorPlane(void);
		/**
		 * Constructor which defines a square detector defined by the center Position P, the surface normal n, the width d and the number of cells in one direction N (so the array is N x N).
		 */
		DetectorPlane(maths::Vector<double> P, maths::Vector<double> n, double d, int N);
		/**
		 * Constructor which defines a more general detector.
		 *
		 * \param P  Position vector
		 * \param e1 first direction vector
		 * \param e2 second direction vector
		 * \param n1 number of cells in e1 direction
		 * \param n2 number of cells in e2 direction
		 */
		DetectorPlane(maths::Vector<double> P, maths::Vector<double> e1, maths::Vector<double> e2, int n1, int n2);
		bool cross(maths::Vector<double> P, maths::Vector<double> k, int& i1, int& i2, double& l);	///< implementation of the intersection checking function for the plane detector		
	};

	/*class DetectorBox : public Detector
	{
	public:
		DetectorBox(maths::Vector<double> P, maths::Vector<double> d, maths::Vector<int> n);
		bool cross(maths::Vector<double> P, maths::Vector<double> k, int& i1, int& i2, double& l);
	     
		std::vector<std::vector<std::vector<std::complex<double> > > > data;
		Box box;
	};*/
#define SAVE_X 0
#define SAVE_Y 1
#define SAVE_Z 2
#define SAVE_PHASE_X 3
#define SAVE_PHASE_Y 4
#define SAVE_PHASE_Z 5
#define SAVE_ABS 6
}
}
#endif
