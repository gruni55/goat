#pragma once
#include "vector.h"
/**
 * Abstract base class for all rays used for the raytracing process. This is the parent class, from which all ray classes are derived. 
 */
class RayBase {
public:
    virtual bool next() = 0;
    virtual Vector<double> getk() = 0; ///< gives back direction Vector
    virtual Vector<std::complex<double> > getE() = 0; ///< gives back electric field strength
    virtual Vector<double> getP() = 0; ///< gives back position
    virtual bool isInObject() = 0; ///< true if ray is in an object otherwise false
    virtual int objectIndex() = 0; ///< index of the current object or -1
    virtual void reflectRay(RayBase *&tray, Vector<double> n, std::complex<double> n1, std::complex<double> n2) = 0; ///< reflects ray
    RayBase* tray; ///< transmitted ray (used by the raytracer, for internal use only)
    bool inObject; ///< is in an object
    ObjectShape** Obj; ///< list of all objects
    int numObj; ///< number of objects
    int objIndex; ///< index of the current object
    std::complex<double>  n, n0; ///< current refractive index and refractive index of the host material
    double k0; ///< wave number
    double r0, rc; ///< radius of the calculation sphere
    int iR; ///< number of reflections already done
};
