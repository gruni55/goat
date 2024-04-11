#pragma once
#include "raytrace.h"
#include "superarray.h"
#include <vector>
namespace GOAT
{
    namespace raytracing
    {
        /**
        * @brief Class which stores the electric field inside a box 
        * This raytracer can store the electric field even if there are objects in at least one box (we will call this box a "box detector").
        * The different box detectors must not interfere ! (In the moment, this will not be checked!)
        */
        class Raytrace_Field :
            public Raytrace
        {
        public:
            Raytrace_Field();
            Raytrace_Field(Scene& S);
            void traceOneRay(RayBase* ray, int& Reflexions, int& recur);            
            void addBoxDetector(Box* box); ///< add a box as detector
            void init(); ///< do some initialisation (e.g. clear the superarray)           
            SuperArray<maths::Vector<std::complex<double> > > SE;            
            int iR = 0; ///< Number of reflections to consider
            void trace(); ///< Start the raytracing process
            void setResolution(double res); ///< set the resolution (which will be used for all detectors)

        private:         
            void traceEnterObject(); 
            void traceLeaveObject();
            void storeData(maths::Vector<double> PStart, maths::Vector<double> Pen, maths::Vector<std::complex<double> > EStart);
            int findBoxDetectorIntersection(maths::Vector<double> P, maths::Vector<double> k, maths::Vector<double>& pout); 
            std::vector<Box *> BoxDetector;  
            double resolution=0.1;
            INDEX_TYPE numCellsPerDirection;
            int indexCurrentDetector=-1;
        };
    }
}
