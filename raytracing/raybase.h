#pragma once
#include "vector.h"
#include "objectshape.h"
/**
 * Abstract base class for all rays used for the raytracing process. This is the parent class, from which all ray classes are derived. 
 */
namespace GOAT
{
    namespace raytracing
    {
#define RAYBASE_STATUS_NONE       0
#define RAYBASE_STATUS_FIRST_STEP 1
        class RayBase {
        public:
            virtual bool next() = 0;
            virtual maths::Vector<double> getk() = 0; ///< gives back direction Vector
            virtual maths::Vector<std::complex<double> > getE() = 0; ///< gives back electric field strength
            virtual maths::Vector<double> getP() = 0; ///< gives back position
            virtual bool isInObject() = 0; ///< true if ray is in an object otherwise false
            virtual int objectIndex() = 0; ///< index of the current object or -1
            virtual void reflectRay(RayBase*& tray, maths::Vector<double> n, std::complex<double> n1, std::complex<double> n2) = 0; ///< reflects ray
            RayBase* tray=0; ///< transmitted ray (used by the raytracer, for internal use only)
            bool inObject=false; ///< is in an object
            ObjectShape** Obj=0; ///< list of all objects
            int numObj=0; ///< number of objects
            int objIndex=0; ///< index of the current object
            std::complex<double>  n, n0; ///< current refractive index and refractive index of the host material
            double k0=2.0*M_PI; ///< wave number
            double r0=1.0, rc=1.0; ///< radius of the calculation sphere
            int iR=0; ///< number of reflections already done
            bool suppress_phase_progress=false; ///< suppress phase change in next(), needed for short pulse consideration
            int status = RAYBASE_STATUS_NONE;
        };
    }
}