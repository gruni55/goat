#pragma once
#include "raytrace.h"
#include "raytrace_inel.h"
#include "vector.h"
#include <vector>

namespace GOAT
{
    namespace raytracing
    {
        typedef struct stepEntry ///< describes the data stored for one ray step
        {
            double l; ///< step size
        int matIndex; ///< material index, stores which object (=object number) is hidden. -1 if ray moves in the surroundings
        };

        typedef struct gridEntry ///< Type of one entry in the grid
        {
            std::vector<stepEntry> step; ///< holds info for one step
            maths::Vector<std::complex<double> > E; ///<electric field (without phase advance along the path) 
        };


        /**
        * @brief This class performs ray tracing in preparation for a later calculation of short pulses
        */
        class Raytrace_usp :
            public Raytrace
        {
        public:
            Raytrace_usp();
            Raytrace_usp(const Scene& S, int n);
            void init();
            void trace();
            void storeData();
            void traceEnterObject(); ///< Function internally called when ray enters an object
            void traceLeaveObject(); ///< Function internally called when ray leaves an object
            

            std::vector<SuperArray <gridEntry> > SA;
            int n = 0;
            int iR = 0;
        };
    }
}

