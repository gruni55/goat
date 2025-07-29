#pragma once
#include<vector>
#include"vector.h"
#include <iostream>
namespace GOAT
{
	namespace raytracing
	{
        /**
        * @brief Structure to store the information about one raytracing step
        * In this structure, the information about a single raytracing step is stored. It is
        * used within the ultrashort pulse calculation process
        */
        typedef struct  ///< describes the data stored for one ray step
        {
            double l; ///< step size
            int matIndex; ///< material index, stores which object (=object number) is hidden. -1 if ray moves in the surroundings
        } stepEntry; 

        /**
         * @brief operator << Output operator for the stepEntry structure
         * @param os stream to write on
         * @param se stepEntry object to write
         * @return
         */
        std::ostream & operator << (std::ostream &os, const stepEntry &se);
  

        /**
        * @brief Structure which holds all steps from the light source to the grid point
        * In this structure every step is stored from the light source to the grid point and additonally
        * the electric field at the grid point. For the electric field determination, only the Fresnel matrices,
        * which describes the interaction between the ray and the object's surface, will be considered 
        * This structure is needed within the ultrashort pulse calculation
        */
        typedef struct 
        {
            std::vector<stepEntry> step; ///< holds info for one step
            maths::Vector<std::complex<double> > E; ///<electric field (without phase advance along the path) 
        } gridEntry;             
                 
        std::ostream & operator << (std::ostream &os, const gridEntry &ge);
  }
}
