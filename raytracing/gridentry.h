#pragma once
#include<vector>
#include"vector.h"
#include <iostream>
namespace GOAT
{
	namespace raytracing
	{
        typedef struct  ///< describes the data stored for one ray step
        {
            double l; ///< step size
            int matIndex; ///< material index, stores which object (=object number) is hidden. -1 if ray moves in the surroundings
        } stepEntry; 

        std::ostream & operator << (std::ostream &os, const stepEntry &se);
  

        /**
        * @brief Structure which holds the
        */
        typedef struct ///< Type of one entry in the grid
        {
            std::vector<stepEntry> step; ///< holds info for one step
            maths::Vector<std::complex<double> > E; ///<electric field (without phase advance along the path) 
        } gridEntry;             
         
        std::ostream & operator << (std::ostream &os, const gridEntry &ge);
  }
}
