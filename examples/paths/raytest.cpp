#include <iostream>
#include "ellipsoid.h"
#include "lightsrc.h"
#include "raytrace.h"
#include "vector.h"
#include <fstream>
#include "surface.h"
#include "octree.h"
#include "detector.h"


int main()
{
    Scene S; 
   
    // define object
    const int numberObjects = 1;
    Ellipsoid obj(Vector<double> (0,0,-10),Vector<double>(30,10,10),1.5);
    obj.setBeta(50.0 / 180.0 * M_PI);
    S.setr0(100);
    S.addObject(&obj);

    Ellipsoid obj2(Vector<double>(-60, 0, -20), Vector<double>(10, 10, 20), 1.6);
    obj2.setBeta(-30.0 / 180.0 * M_PI);
    S.addObject(&obj2);

    LightSrcPlane ls = LightSrcPlane(-50.0 * ez, 30, 1.0, 30.0 );
    S.addLightSource(&ls);

    LightSrcGauss ls2(Vector<double>(-60.0, 0.0, 20.0), 20, 1.0, 0.1, Vector<double>(-60.0, 0.0, 0.0));
    ls2.setNA(0.9);
    S.addLightSource(&ls2);
    

    Raytrace_Path rp;
    rp.setShowOutgoingRays(false);
    rp.setScene(S);
    rp.trace("test.dat");
}
