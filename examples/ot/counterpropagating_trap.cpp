/*****************************************************************//**
 * \file   ot_test.cpp
 * \brief  Example for optical tweezers calculation
 *
 * This program calculates the optical forces on a spherical particle with radius 5µm in the field of a focused
 * gaussian beam with an NA of 1.2 in water (refractive index: 1.33)
 *
 * \author Thomas
 * \date   July 2021
 *********************************************************************/
#include <fstream>
#include <Raytrace.h>

int main()
{
    std::ofstream os;
    Scene S;

    /* Define object, sphere as ellipsoid with radius 5µm and refractive index 1.59 */
    double wvl = 1.064;
    double wvlm = wvl / 1.33;
    Ellipsoid obj(dzero, Vector<double>(10,5,5), 1.5);

    /* Defining Scene with calculation radius 1000µm and a refractive index of 1.33 */
    S.setnS(1.33);
    S.setr0(1000.0);
    S.addObject(&obj);  // add the object to the scene

   /* Define light source. Here, a gaussian light source is used with focus point at (0,0,0) positioned at (0,0,-200)
       with 200 rays per direction and a wavelength of 1.064µm and a (fictious) waist diameter of 0.5µm. The waist diameter
       will later be overwritten by the definition of the NA. The x-position of the particle is varied from -10µm to 10µm */

    double d = 5;
    LightSrcGauss ls1 = LightSrcGauss(-200.0 * ez, 200, 1.064, 0.5, dzero+d*ez);
    LightSrcGauss ls2 = LightSrcGauss(200.0 * ez, 200, 1.064, 0.5, dzero-d*ez);

    /* set the power to 1W */
    ls1.P0 = 1.0;
    ls2.P0 = 1.0;

    /* add light source to scene */
    S.addLightSource(&ls1);
    S.addLightSource(&ls2);

    /* set the NA to 1.2 */
    ls1.setNA(1.2);
    ls2.setNA(1.2);

    /* initialize raytracer for optical tweezers */
    Raytrace_OT ot(S);

    os.open("cp5.dat");

    /* variation loop */
    for (double z = -10; z <= 10; z += 0.05)
    {
        obj.P[2] = z; // change object's position
        ot.trace();            // do raytracing 

        /* The result is stored in the array F for the force and in L for the angular momenta
           of the optical tweezers raytracer ot. There is an entry for each object. */

        os << z << "\t" << ot.F[0] << "\t" << ot.L[0] << std::endl;
        std::cout << z << "\t" << ot.F[0] << "\t" << ot.L[0] << std::endl;
    }
    os.close();
}

