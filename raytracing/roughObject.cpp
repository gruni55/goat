#include "roughObject.h"
#include "cone.h"
#include "sphericLens.h"    
#include "cylinder.h"
#include "vortex.h"
namespace GOAT
{
	namespace raytracing
	{
       
        ObjectShape* makeRough(ObjectShape* obj, double sigma)
        {
            if (!obj)
                return nullptr;

            if (obj->isRough())
                return obj;

            switch (obj->Type())
            {
            case OBJECTSHAPE_SURFACE:
            {
                auto robj = new roughObject<surface>(static_cast<surface*>(obj), sigma);
                delete obj;
				return robj;
            }
            case OBJECTSHAPE_ELLIPSOID:
            {
                auto robj=new roughObject<Ellipsoid>(
                    static_cast<Ellipsoid*>(obj), sigma);
                delete obj;
                return robj;
            }

            case OBJECTSHAPE_BOX:
            {
                auto robj = new roughObject<Box>(
                    static_cast<Box*>(obj), sigma);
                delete obj;
                return robj;
            }
			case OBJECTSHAPE_CONE:
			{
                auto robj = new roughObject<Cone>(
                    static_cast<Cone*>(obj), sigma);
                delete obj;
                return robj;
			}
			case OBJECTSHAPE_CYLINDER:
			{
                auto robj = new roughObject<Cylinder>(
                    static_cast<Cylinder*>(obj), sigma);
                delete obj;
                return robj;
			}
			case OBJECTSHAPE_SPHERIC_LENS:
            {
                auto robj = new roughObject<sphericLens>(
                    static_cast<sphericLens*>(obj), sigma);
                delete obj;
                return robj;
			}
			case OBJECTSHAPE_VORTEX_PLATE:
			{
                auto robj = new roughObject<VortexPlate>(
                    static_cast<VortexPlate*>(obj), sigma);
                delete obj;
                return robj;
			}
            default:
                return nullptr;
            }
        }
	}
}