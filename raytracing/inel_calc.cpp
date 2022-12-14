#include "inel_calc.h"
#include "plane.h"
#ifdef _MSC_VER
  #ifndef copysign
    #define copysign _copysign
  #endif 
#endif 
namespace GOAT
{
  namespace raytracing
  {
    maths::Vector<double> startpunkt(int ix, int iy, grid& git, 
                              double dr, double theta, double phi)
    {
    maths::Vector<double> start, exrot, eyrot, r0;


    exrot = drehvektor(maths::ex, theta, phi);
    eyrot = drehphivektor(maths::ey, phi);

    r0 = ursprungrot(dr, theta, phi, exrot, eyrot, git);
    start = r0 + ix*exrot + iy*eyrot;

    start[0] = start[0] + git.xmax/2.0;
    start[1] = start[1] + git.ymax/2.0;
    start[2] = start[2] + git.zmax/2.0;

    return start;
    }

    maths::Vector<double> startpunkt(int ix, int iy, grid& git, double theta, double phi)
    {
    maths::Vector<double> start, exrot, eyrot, r0;

    exrot = drehvektor(maths::ex, theta, phi);
    eyrot = drehphivektor(maths::ey, phi);
    r0 = ursprungrot(0.0, theta, phi, exrot, eyrot, git);
    start = r0 + ix*exrot + iy*eyrot;

    start[0] = start[0] + git.xmax/2.0;
    start[1] = start[1] + git.ymax/2.0;
    start[2] = start[2] + git.zmax/2.0;

    return start;
    }

    maths::Vector<double> drehvektor(maths::Vector<double> vein, double theta, double phi)
    {
    maths::Vector <double> xyzergebnis,kugergebnis;

    kugergebnis=cart2sphere(vein);

    kugergebnis[1]=kugergebnis[1]+theta;
    kugergebnis[2]=kugergebnis[2]+phi;

    xyzergebnis = sphere2cart(kugergebnis);

    return xyzergebnis;
    }

    maths::Vector<double> drehphivektor(maths::Vector<double> vein, double phi)
    {
    maths::Vector <double> xyzergebnis,kugergebnis;

    kugergebnis=cart2sphere(vein);

    kugergebnis[2]=kugergebnis[2]+phi;
    xyzergebnis = sphere2cart(kugergebnis);

    return xyzergebnis;
    }

    maths::Vector<double> ursprungrot(double dr, double theta, double phi, 
                              maths::Vector<double>& exrot, maths::Vector<double>& eyrot, 
                              grid& git)
    {
    maths::Vector<double> rm, rmrot, dm;
    double nvx, nvy;
    rm[0] = 0.0; rm[1] = 0.0; rm[2] = git.zmax/2.0+dr;

    rmrot = drehvektor(rm,theta,phi);
    nvx = git.xmax/(2.0*git.dx);
    nvy = git.ymax/(2.0*git.dy);
    dm =  (2*nvx-1)*git.dx/2.0* exrot + 
          (2*nvy-1)*git.dy/2.0* eyrot;
    return (rmrot - dm);
    }

    #ifdef WITH_SUPERGITTER
    SuperArray verfolgung(maths::Vector<double> p0, maths::Vector<double> k0, SuperArray &git)
    {
    SuperArray partikel;

    partikel = git;

    // Berechnung des Eintrittspunktes i. d. Partikel
    return partikel;
    }

    maths::Vector<double> pnext(maths::Vector<double> p0, maths::Vector<double> k0, SuperArray &git)
    {

    double lambdax, lambday, lambdaz, lambda;
    double signx, signy, signz;

    signx = copysign(1,k0[0]); signy = copysign(1,k0[1]); signz = copysign(1,k0[2]);

    lambdax = ((floor(p0[0]/git.d[0]) + signx + (signx<0)*(fmod(p0[0],git.d[0])!=0))*git.d[0]-p0[0])/k0[0];
    lambday = ((floor(p0[1]/git.d[1]) + signy + (signy<0)*(fmod(p0[1],git.d[1])!=0))*git.d[1]-p0[1])/k0[1];
    lambdaz = ((floor(p0[2]/git.d[2]) + signz + (signz<0)*(fmod(p0[2],git.d[2])!=0))*git.d[2]-p0[2])/k0[2];

    if (lambdax<lambday)
      lambda = lambdax;
    else 
      lambda = lambday;

    if (lambda>lambdaz)
      lambda = lambdaz; 

    return  p0 + lambda*k0;

    }

    maths::Vector<double> pnext(maths::Vector<double> p0, maths::Vector<double> k0, SuperArray &git,double eps) 
    {
    double lambdax, lambday, lambdaz, lambda;
    double signx, signy, signz;
    double sx, sy, sz;

    /***********************************************************************************************
     * Algorithmus zur Bestimmung der Schnittpunkte des Strahls mit den gridebenen               *
     *                                                                                             *
     * Prinzip:                                                                                    *
     * Ausgehend vom Punkt p_0 werden die Schnittpunkte mit Ebenen x=const.,                       *
     * y=const. und z=const. gebildet. Hierzu werden zunaechst die Stuetzpunkte s_x, s_y und s_z,  *
     * welche den Abstand dieser Ebenen vom Ursprung des grids angeben, berechnet. Der           * 
     * Schnittpunkt mit dem kleinsten Abstand ist dann der gesuchte Punkt.                         *
     *                                                                                             *
     * Die Stuetzpunkte s_i bestimmen sich folgendermassen:                                        *
     *                                                                                             *  
     *                                                                                             *
     * Die i-Komponente des Punktes p_0 ist von der Form                                           *
     *                                                                                             *
     *  p_{0,i} = m*di + a,  m ganz, a >=0                                                         *
     *                                                                                             *
     * wobei di die Schrittweite des grids in i-Richtung ist.                                    *
     *                                                                                             *
     * Ist k_i > 0 und a = 0, so ist der Stuetzpunkt s_i gegeben durch (m+1)*di=(m+sign(k_i)+0)*di *
     * Ist k_i > 0 und a > 0, so ist der Stuetzpunkt s_i gegeben durch (m+1)*di=(m+sign(k_i)+0)*di *
     * Ist k_i < 0 und a = 0, so ist der Stuetzpunkt s_i gegeben durch (m-1)*di=(m+sign(k_i)+0)*di *
     * Ist k_i < 0 und a > 0, so ist der Stuetzpunkt s_i gegeben durch   m*di  =(m+sign(k_i)+1)*di *
     *                                                                                             *
     * Der Schnittpunkt mit dem kleinsten Abstand ist dann der gesuchte Punkt.                     *
     ***********************************************************************************************/

    signx = copysign(1.0,k0[0]); signy = copysign(1.0,k0[1]); signz = copysign(1.0,k0[2]);

    sx = (floor((p0[0]+signx*2*eps)/git.d[0]) + signx + (signx<0)*(fmod(p0[0]+signx*2*eps,git.d[0])!=0))*git.d[0];
    sy = (floor((p0[1]+signy*2*eps)/git.d[1]) + signy + (signy<0)*(fmod(p0[1]+signy*2*eps,git.d[1])!=0))*git.d[1];
    sz = (floor((p0[2]+signz*2*eps)/git.d[2]) + signz + (signz<0)*(fmod(p0[2]+signz*2*eps,git.d[2])!=0))*git.d[2];

    lambdax = (sx - p0[0])/k0[0];
    lambday = (sy - p0[1])/k0[1];
    lambdaz = (sz - p0[2])/k0[2];

    if (lambdax<=lambday)
      lambda = lambdax;
    else 
      lambda = lambday;

    if (lambda>=lambdaz)  lambda = lambdaz;


    return  p0 + lambda*k0;

    }

    #endif
    maths::Vector<double> pnext(maths::Vector<double> p0, maths::Vector<double> k0,
    maths::Vector<double> d, double eps)
    {
    double lambdax, lambday, lambdaz, lambda;
    double signx, signy, signz;
    double sx, sy, sz;

    signx = copysign(1.0,k0[0]); signy = copysign(1.0,k0[1]); signz = copysign(1.0,k0[2]);

    sx = (floor((p0[0]+signx*2*eps)/d[0]) + signx + (signx<0)*(fmod(p0[0]+signx*2*eps,d[0])!=0))*d[0];
    sy = (floor((p0[1]+signy*2*eps)/d[1]) + signy + (signy<0)*(fmod(p0[1]+signy*2*eps,d[1])!=0))*d[1];
    sz = (floor((p0[2]+signz*2*eps)/d[2]) + signz + (signz<0)*(fmod(p0[2]+signz*2*eps,d[2])!=0))*d[2];

    lambdax = (sx - p0[0])/k0[0];
    lambday = (sy - p0[1])/k0[1];
    lambdaz = (sz - p0[2])/k0[2];

    if (lambdax<=lambday)
      lambda = lambdax;
    else
      lambda = lambday;

    if (lambda>=lambdaz)
      lambda = lambdaz;


    return  p0 + lambda*k0;

    }


    maths::Vector<double> pnext(Plane E, maths::Vector<double> p0s, maths::Vector<double> k0s,
    maths::Vector<double> d, double eps)
    {
    double lambdax, lambday, lambdaz, lambda;
    double signx, signy, signz;

    double sx, sy, sz;
    maths::Vector<double> p0,k0;

    p0=maths::Vector<double> (p0s*E.e1,p0s*E.e2,p0s*E.n);
    k0=maths::Vector<double> (k0s*E.e1,k0s*E.e2,k0s*E.n);


    signx = copysign(1.0,k0[0]); signy = copysign(1.0,k0[1]); signz = copysign(1.0,k0[2]);

    sx = (floor((p0[0]+signx*2*eps)/d[0]) + signx + (signx<0)*(fmod(p0[0]+signx*2*eps,d[0])!=0))*d[0];
    sy = (floor((p0[1]+signy*2*eps)/d[1]) + signy + (signy<0)*(fmod(p0[1]+signy*2*eps,d[1])!=0))*d[1];
    sz = (floor((p0[2]+signz*2*eps)/d[2]) + signz + (signz<0)*(fmod(p0[2]+signz*2*eps,d[2])!=0))*d[2];

    lambdax = (sx - p0[0])/k0[0];
    lambday = (sy - p0[1])/k0[1];
    lambdaz = (sz - p0[2])/k0[2];

    if (lambdax<=lambday)
      lambda = lambdax;
    else
      lambda = lambday;

    if (lambda>=lambdaz)
      lambda = lambdaz;


    return  p0s + lambda*k0s;

    }


  }
}
