#include "roughObject.h"
namespace GOAT
{
	namespace raytracing
	{
       
        ObjectShape* makeRough(ObjectShape* obj, double sigma)
        {
            switch (obj->Type())
            {
            case OBJECTSHAPE_SURFACE:
                return new roughObject<surface>(
                    static_cast<surface*>(obj), sigma);

            case OBJECTSHAPE_ELLIPSOID:
                return new roughObject<Ellipsoid>(
                    static_cast<Ellipsoid*>(obj), sigma);

            case OBJECTSHAPE_BOX:
                return new roughObject<Box>(
                    static_cast<Box*>(obj), sigma);

            default:
                return nullptr;
            }
        }
	}
}