#include "vector.h"
#include "matrix.h"

#pragma once
#include <vector>
namespace GOAT
{
	namespace raytracing
	{
        enum Curvature {convex,concave,flat};
        /**
       * @brief Structure to describe one side of an spheric lens 
       * For this type of lens one side can either be spherical (convex or  concave) or flat 
       */
        typedef struct
        {
            /**
            * @brief center of the sphere, which describes the corresponding surface of the lens
            * Vector from the position vector of the lens to the center of the sphere
            */ 
            maths::Vector<double> P; 
            double R; ///< Radius of the sphere
            double shift = 0; ///< how much must the curve be shifted (with consideration of the offset)
            Curvature curvature; ///< Curvature (convex, concave or flat)
        } lensSide;


        /**
        * @brief Structure, which holds the full information about the spheric lens
        * The side surfaces of the lens are described by left and right. The thickness of the lens is described by offset
        * and the height/radius of the lens by radius (for details, see also the documentation of the class sphericLens)
        */
        typedef struct
        {
            lensSide left, right; ///< descriptions of the left (towards negative-z) and the right side of the lens
            double offset; ///< distance between the surfaces 
            double radius; ///
        } lensParms;

     
        void binWrite(lensSide ls, std::ofstream& os); ///< write a lensSide structure in a binary file
        void binWrite(lensParms lp, std::ofstream& os); ///< write a lensParms structure in a binary file
        void binRead(lensSide ls, std::ifstream& is); ///< read a lensSide structure from a binary file       
        void binRead(lensParms lp, std::ifstream& is);///< read a lensParms structure from a binary file       
	}
}