#ifndef DETECTOR_H
#define DETECTOR_H

#include "vector.h"
#define DETECTOR_ANGLE 1
#define DETECTOR_PLANE 2
#include <iostream>

/**
 * The abstract Detector class provides an interface to a detector to store the information about the electric field into any kind of an array. 
 */
class Detector 
{
public:
	Detector(void);
	Detector(int n1, int n2); ///< Constructor which is called by the size of the array (n1 x n2)
	Detector(const Detector &Det); ///< copy constructor
	Detector& operator= (const Detector& D); 
	Vector<std::complex<double> >& operator () (int i1,int i2) { return D[i1][i2]; } ///< bracket operator which gives the content inside the array determined by the indices i1 and i2
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
	virtual bool cross(Vector<double> P, Vector<double> k, int &i1, int &i2, double &l) = 0; 
	int N1(); ///< returns the dimension of the array in the first direction
	int N2(); ///< returns the dimension of the array in the second direction
	int Type() {return type; } ///< returns kind of detector
	void saveabs (char *fn);  ///< stores the content (absolute value of the electric field) of the detector array in the file determined by its filename fn 
	void savePhase (char *fn, int coord);  ///< stores the content (phase of one component of the electric field, coord determines the coordinate 0,1,2 for x,y,z) of the detector array in the file fn 
    void savereal (char *fn, int coord); ///< stores the content (real part of one component of the the electric field, coord determines the coordinate 0,1,2 for x,y,z) of the detector array in the file fn 
    void saveimag (char *fn, int coord); ///< stores the content (imaginary part of one component of the the electric field, coord determines the coordinate 0,1,2 for x,y,z) of the detector array in the file fn 
	Vector<std::complex<double> >** D; ///< Here, the data will be stored

friend std::ostream& operator << (std::ostream &os, Detector& D);  
	protected: 
		Vector<double> e1; 
		Vector<double> e2;
		Vector<double> P;
		Vector<double> n;
		void init(int n1, int n2); ///< initialise array (for internal use only)
		double d1, d2;
int n1,n2;
int type;
friend class DetectorPlane;
}; 


/**
 * This class provides a plane detector, defined by its center point and its side length. The detector plane is spanned by the vectors e1 and e2. The absolute value of the vectors
 * e1 and e2 are the cell widths in the corresponding direction. 
 */
class DetectorPlane : public Detector
{
public:
  DetectorPlane (void);
  /**
   * Constructor which defines a square detector defined by the center Position P, the surface normal n, the width d and the number of cells in one direction N (so the array is N x N).
   */
  DetectorPlane(Vector<double> P, Vector<double> n, double d, int N); 
  /**
   * Constructor which defines a more general detector.
   * 
   * \param P  Position vector
   * \param e1 first direction vector
   * \param e2 second direction vector
   * \param n1 number of cells in e1 direction
   * \param n2 number of cells in e2 direction
   */
  DetectorPlane (Vector<double> P, Vector<double> e1, Vector<double> e2, int n1, int n2);
  bool cross(Vector<double> P, Vector<double> k, int &i1, int &i2, double &l);	///< implementation of the intersection checking function for the plane detector
};

#define SAVE_X 0
#define SAVE_Y 1
#define SAVE_Z 2
#define SAVE_PHASE_X 3
#define SAVE_PHASE_Y 4
#define SAVE_PHASE_Z 5
#define SAVE_ABS 6

#endif
