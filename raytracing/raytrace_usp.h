#pragma once
#include "raytrace.h"
#include "raytrace_inel.h"
#include "vector.h"
#include "gridentry.h"
#include <vector>

namespace GOAT
{
    namespace raytracing
    {
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
            Raytrace_usp(const Scene& S, INDEX_TYPE n);
            void clear(); ///< Clears the SuperArray grid (SA) 
            void init(); ///< Initialises the SuperArray grid (SA) according to the objects hold in Scene S
            void trace(); ///< Perform the raytracing process
            void storeData(); ///< Store the data in the grid
            void traceEnterObject(); ///< Function internally called when ray enters an object
            void traceLeaveObject(); ///< Function internally called when ray leaves an object
            void setn(INDEX_TYPE n); ///< changes the number of (virtual) cells in the calculation space
            /**
            * @brief Grid, that holds the information, needed to further calculate the electric field for short pulses
            * The vector holds the arrays for the different reflection orders (each reflection order is stored in a separate array)
            * In each element of the SuperArray, the steps to this point are stored (step length, index of the material)
            */
            std::vector<SuperArray <std::vector<gridEntry> > > SA; 
            INDEX_TYPE n = 1; ///< Number of cells in one direction
            int iR = 0; ///< Number of reflections to consider
            gridEntry stack; ///< here, the information from the light source until the region of interest (=box) is reached
        };
    }
}

