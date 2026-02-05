 #pragma once
#include "detector.h"
#include "objectshape.h"
#include "superarray.h"
#include "vector.h"
namespace GOAT
{
	namespace raytracing
	{

		maths::Vector<std::complex<double> > point(DetectorPlane* det, maths::Vector<double> P, double wvl);

		/**
		* @brief This class makes a Kirchhoff calculation
		* This class is directly connect with a Detector. The calculation itself works as follows: 
		* At first, a normal raytracing step is performed to calculate the electric field at a detector. This detector is used as a 
		* source field for the next step, where the field at a given area is calculated with help of the Kirhhoff integral
		*/
		class Kirchhoff : public DetectorPlane
		{
		  public : 
			  /**
			  * @brief constructor
			  * \param wvl: Wavelength used for the caluclation of the electric field
			  * \param P: Position of the Kirchhoff plane
			  * \param e1: first direction vector of one edge of the detector
			  * \param e2: second direction vector of the edge of the detector perpendicular to the first one
			  * \param n1: number of cells in e1-direction
			  * \param n2: number of cells in e2-direction
			  */
			  Kirchhoff(double wvl, maths::Vector<double> P, maths::Vector<double> e1, maths::Vector<double> e2, int n1, int n2);
			  void addDetector(DetectorPlane* det);
			  void addDetectorList(std::vector<DetectorPlane*> detList);
			  void delDetector(DetectorPlane* det);
			  void clearSources();
			  void calc(bool clear = true);
			  void setNumberOfThreads(int noThreads);
			  int numberOfThreads() { return noThreads; }

		private:
			/**
			 *  @brief This method make the calculation
			 * With this method, the calculation of the Kirchhoff-integral will be performed for one detector.
			 * \param det: a pointer to the detector, which acts as the source
			 */
			void calc(DetectorPlane* det, bool clear = true);

			/**
			* @brief Do the Kirchhoff calculation with more than one detector as source
			*/
			void calc(std::vector<DetectorPlane*> detList);

			double k;
			double wvl;
			
			std::vector<DetectorPlane*> sources;
			int noThreads = 8;
		};

		/**
		*  @brief This class makes a 3D Kirchhoff calculation
		* This class makes a 3D Kirchhoff calculation especially for pulsed calculations. It needs a Box object to define the volume, where the field is calculated.
		*/
		class Kirchhoff3D
		{
		  public:
			  /**
			  * @brief constructor
			  * Constructs the Kirchhoff3D object with the box defining the calculation volume
			  * \param box: Box object defining the calculation volume
			  * \param numCellsPerDir: number of cells per direction (the world koordinate system will be divided in numCellsPerDir x numCellsPerDir x numCellsPerDir cells)
			  */
			  explicit Kirchhoff3D(Box* box, int numCellsPerDir); ///< constructor with the box defining the calculation volume
			  void addDetector(DetectorPlane* det); ///< adds one detector as source
			  void addDetectorList(std::vector<DetectorPlane*> detList); ///< adds a list of detectors as sources
			  void clean() { field3D.fill(maths::czero); }; ///< cleans the calculated field (sets all values to zero)
			  void calc(double wvl, int noThreads = 8); ///< performs the calculation for the given wavelength wvl with noThreads threads
			  /**
			  * @brief bracket operator to access the calculated field 
			  * \param ix: index in x-direction
			  *	\param iy: index in y-direction
			  * \param iz: index in z-direction
			  */
			  maths::Vector<std::complex<double>>& operator () (INDEX_TYPE ix, INDEX_TYPE iy, INDEX_TYPE iz) { return field3D(0, ix, iy, iz); }
			  const SuperArray<maths::Vector<std::complex<double>>>& field() const {
				  return field3D;
			  }
			  SuperArray<maths::Vector<std::complex<double>>> field3D; ///< 3D array storing the calculated field

		  private:
			  /**
			 *  @brief This method make the calculation
			 * With this method, the calculation of the Kirchhoff-integral will be performed for one detector.
			 * \param det: a pointer to the detector, which acts as the source
			 * \param clear: if true, the field3D array will be cleared before calculation
			 */
			  void calc(DetectorPlane* det, double wvl, int numThreads, bool clear = true);
			  Box *box; ///< list of boxes defining the calculation volume
			  std::vector<DetectorPlane*> sources; ///< list of detectors acting as sources
			 
		};
	}
}