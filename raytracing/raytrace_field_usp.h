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
        * 
        */
        class Raytrace_Field_usp :
            public Raytrace        {
        public: 
            Raytrace_Field_usp();
            Raytrace_Field_usp(const Scene& S, INDEX_TYPE n);
            ~Raytrace_Field_usp();
            void traceOneRay(RayBase* ray, int& Reflexions, int& recur);
            void storeData(maths::Vector<double> PStart, maths::Vector<double> Pen, maths::Vector<std::complex<double> > EStart);
            void clear(); ///< Clears the SuperArray grid (SA) 
            void init();
            void trace();
            void storeStack(maths::Vector<double> PStart, maths::Vector<double> PStop);
            void traceEnterObject() {}
            void traceLeaveObject() {}
            void addBoxDetector(Box* box) { BoxDetector.push_back(box); }
            void addBoxDetectorList(std::vector<Box*> BoxDetector);
            int indexCurrentDetector = -1;
            std::vector<Box*> BoxDetector;
            std::vector<SuperArray <std::vector<gridEntry> > > SA;       
            int findBoxDetectorIntersection(maths::Vector<double> P, maths::Vector<double> k, maths::Vector<double>& pStart, maths::Vector<double>& pStop);
            INDEX_TYPE n = 1; ///< Number of cells in one direction
            int iR = 0; ///< Number of reflections to consider
            gridEntry stack; ///< here, the information from the light source until the region of interest (=box) is reached
            maths::Vector<double> pDetStart, pDetStop;
        };
    }
}
