#pragma once
#include "raytrace_inel.h"
#include "raytrace_field.h"
#include "gridentry.h"
#include <vector>


namespace GOAT
{
    namespace raytracing
    {

        /**
        * @brief This class provides a raytracer to calculate the field distribution in a box detector for a ultrashort pulsed light source
        * Unlike the class Raytrace_usp, the field in more than one  object or in parts of objects can stored. 
        * As in Raytrace_usp, in the SuperArray grid the information (step length, index of the object involved in the step) for each step from the light source to the corresponding point in the grid is stored. For each reflection order one grid will be created. The field will be stored within one or more boxes (called "box detector")
       */
        class Raytrace_Field_usp :
            public Raytrace        {
        public: 
            Raytrace_Field_usp();
            Raytrace_Field_usp(const Scene& S, INDEX_TYPE n);
            ~Raytrace_Field_usp();
            void traceOneRay(RayBase* ray, int& Reflexions, int& recur);
            void storeData(maths::Vector<double> PStart, maths::Vector<double> Pen, maths::Vector<std::complex<double> > EStart);
	    void clean(); ///< removes all content from SuperArray grid (SA)
            void clear(); ///< Clears the SuperArray grid (SA) 
            void init(); ///< Initialises the SuperArray grid 
            void trace(); ///< make the raytracing
            void storeStack(maths::Vector<double> PStart, maths::Vector<double> PStop); ///< 
            void traceEnterObject() {} ///< has to be added since this class is a child of Raytrace (makes nothing)
            void traceLeaveObject() {} ///< has to be added since this class is a child of Raytrace (makes nothing)

            void addBoxDetector(Box* box) { BoxDetector.push_back(box); } ///< add one box detector
            void addBoxDetectorList(std::vector<Box*> BoxDetector); ///< add a list of box detectors
            int indexCurrentDetector = -1; ///< index of the current box detector (for internal use only)
            std::vector<Box*> BoxDetector; ///< list of the box detectors
            std::vector<SuperArray <std::vector<gridEntry> > > SA; ///< Here, the information for the steps is stored (needed for Fourier transform in pulsecalculation_field)       
            int findBoxDetectorIntersection(maths::Vector<double> P, maths::Vector<double> k, maths::Vector<double>& pStart, maths::Vector<double>& pStop); ///< searches for the intersection point with the next box detector
            INDEX_TYPE n = 1; ///< Number of cells in one direction
            int iR = 0; ///< Number of reflections to consider
            gridEntry stack; ///< here, the information from the light source until the region of interest (=box) is reached
            maths::Vector<double> pDetStart, pDetStop; ///< 
            int oldObjIndex=-1;
        };
    }
}
