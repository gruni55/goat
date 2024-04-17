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
        class Raytrace_field_usp :
            public Raytrace_Field
        {
        public: 
            Raytrace_field_usp();
            Raytrace_field_usp(Scene &S);

        protected:
            void clear(); ///< Clears the SuperArray grid (SA) 
            void init(); ///< Initialises the SuperArray grid (SA) according to the objects hold in Scene S
            void traceEnterObject();
            void traceLeaveObject();
            void storeData(maths::Vector<double> PStart, maths::Vector<double> Pen, maths::Vector<std::complex<double> > EStart);
            void storeEntry(maths::Vector<double> PStart, maths::Vector<double> PStop);
            std::vector<SuperArray <std::vector<gridEntry> > > SA;
            INDEX_TYPE n = 1; ///< Number of cells in one direction
            int iR = 0; ///< Number of reflections to consider
            gridEntry stack; ///< here, the information from the light source until the region of interest (=box) is reached
        };
    }
}
