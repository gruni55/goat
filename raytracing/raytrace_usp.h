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


        /**
        * @brief Structure which holds the 
        */
        typedef struct gridEntry ///< Type of one entry in the grid
        {
            std::vector<stepEntry> step; ///< holds info for one step
            maths::Vector<std::complex<double> > E; ///<electric field (without phase advance along the path) 
        };


        /**
        * @brief This class performs ray tracing in preparation for a later calculation of short pulses
        * This class performs ray tracing and stores the length and the refractive index as an index of each step into a 
        * SuperArray class named SA. This information is necessary for the calculation of the short pulses.
        */
        class Raytrace_usp :
            public Raytrace
        {
        public:
            Raytrace_usp();
            /**
            * @brief Constructor 
            * @param S: Scene which holds all information about the scene (light sources, objects...)
            * @param n: Number of grid elements per coordinate across the entire width (i.e. 2*r0) used in the SuperArray SA
            */
            Raytrace_usp(const Scene& S, int n); 
            void init(); ///< Initialises the SuperArray grid (SA) according to the objects hold in Scene S
            void trace(); ///< Perform the raytracing process
            void storeData(); ///< Store the data in the grid
            void traceEnterObject(); ///< Function internally called when ray enters an object
            void traceLeaveObject(); ///< Function internally called when ray leaves an object
            

            std::vector<SuperArray <gridEntry> > SA; ///< Grid, that holds the information, needed to further calculate the electric field for short pulses
            int n = 0; ///< Number of cells in one direction
            int iR = 0; ///< Number of reflections to consider
        };
    }
}

