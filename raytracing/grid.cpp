#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include "grid.h"
#include <complex>
#include "goodies.h"


namespace GOAT
{
  namespace raytracing
  {
    grid::grid()
    {
      gridarray = 0;
      dx = 0.0;  dy = 0.0;  dz = 0.0; 
      dx2 = 0.0; dy2 = 0.0; dz2 = 0.0;
      nxmax = 0; nymax = 0; nzmax = 0;
    }

    grid::grid(int nx, int ny, int nz, double x, double y,
                  double z)
    {
    xmax = x; ymax = y; zmax = z;
    nxmax = nx; nymax = ny; nzmax = nz;

    dx = xmax/nx; dx2 = dx/2.0;
    dy = ymax/ny; dy2 = dy/2.0;
    dz = zmax/nz; dz2 = dz/2.0;

    gridarray = newarray();
    }


    grid::grid(int nx, int ny, int nz, double r)
    {
    xmax = 2.0*r; ymax = xmax; zmax = ymax;
    nxmax = nx; nymax = ny; nzmax = nz;
    rP=r;

    dx = xmax/nx; dx2 = dx/2.0;
    dy = ymax/ny; dy2 = dy/2.0;
    dz = zmax/nz; dz2 = dz/2.0;

    
    gridarray = newarray();
    }

    grid::~grid()
    {
    delarray();
    }

    // Wird fuer die Konstruktoren der grid-Klasse ben�tigt,
    // allokiert gridarray 

    maths::Vector<std::complex<double> >& grid::operator()(int i, int j, int k)
    {
    if (gridarray[i][j][k]!=0)
    {
      return *gridarray[i][j][k];
    }
    else
      { 
        DUMMY=maths::Vector<std::complex<double> > (INF,INF,INF);
      return DUMMY;
      }
    }

    maths::Vector<std::complex<double> >& grid::operator()(maths::Vector<int>& P)
    {
    if (gridarray[P[0]][P[1]][P[2]]!=0)
    {
      return *gridarray[P[0]][P[1]][P[2]];
    }
    else
      { 
        DUMMY=maths::Vector<std::complex<double> > (INF,INF,INF);
      return DUMMY;
      }
    }
    grid& grid::operator=(const grid& g)
    {
    if (this != &g)
    {
      xmax = g.xmax; ymax =g.ymax; zmax = g.zmax;
      nxmax = g.nxmax; nymax = g.nymax; nzmax = g.nzmax;
      dx = g.dx; dy  = g.dy; dz = g.dz;
      dx2 = g.dx2; dy2  = g.dy2; dz2 = g.dz2;
      delete gridarray;
      gridarray=newarray();
      for (int ix=0; ix<nxmax; ix++)
      {
      for (int iy=0; iy<nymax; iy++)
      {
        for (int iz=0; iz<nzmax; iz++)
        {
        if (g.gridarray[ix][iy][iz]!=0)
        {
          gridarray[ix][iy][iz] = new maths::Vector<std::complex<double> >;
        *gridarray[ix][iy][iz] = *g.gridarray[ix][iy][iz];
        }
        else
          gridarray[ix][iy][iz] = 0;
        }
      }
      }
    }
    return *this;
    }

    // Wird fuer die Konstruktoren der grid-Klasse ben�tigt,
    // allokiert gridarray

    maths::Vector<std::complex<double> > ****grid::newarray()
    {
    maths::Vector<std::complex<double> > ****hilfarray;

    hilfarray = new maths::Vector<std::complex<double> >***[nxmax];

    for (int ix = 0; ix < nxmax; ix++)
    {
      hilfarray[ix] = new maths::Vector<std::complex<double> >**[nymax];
      for (int iy = 0; iy < nymax; iy++)
      {
      hilfarray[ix][iy] = new maths::Vector<std::complex<double> >*[nzmax];
      for (int iz = 0; iz < nzmax; iz++)
      {
        hilfarray[ix][iy][iz] = 0;
      }
      }
    }
    return hilfarray;
    }

    void grid::delarray()
    {
    for( int ix=0; ix < nxmax; ix++)
    {
      for( int iy=0; iy < nymax; iy++)
      {
      delete [] gridarray[ix][iy];
      }
      delete [] gridarray[ix];
    }
    delete [] gridarray;
    }


    int grid::isInside(int ix, int iy, int iz,
                          double x0, double y0, double z0, double r0)
    // Basisfunktion zur Bestimmung der Lage eines 
    // gridpunktes ix, iy, iz innerhal einer        
    // Kugel vom Radius r_0 mit Mittelpunkt bei
    // x_0, y_0, z_0
    {
    if (maths::sqr((2*ix+1)*dx2-x0)+maths::sqr((2*iy+1)*dy2-y0)+maths::sqr((2*iz+1)*dz2-z0)<=maths::sqr(r0))
    {
      return 1;
    }
    else
      return 0;
    } 

    int grid::isInside(int ix, int iy, int iz, maths::Vector<double> P, double r0)
    // Basisfunktion zur Bestimmung der Lage eines 
    // gridpunktes ix, iy, iz innerhal einer        
    // Kugel vom Radius r_0 mit Mittelpunkt bei P
    {
    if (maths::sqr((2*ix+1)*dx2-P[0])+maths::sqr((2*iy+1)*dy2-P[1])+maths::sqr((2*iz+1)*dz2-P[2])<=maths::sqr(r0))
    {
      return 1;      
    }
    else
      return 0;
    }

    void grid::init_grid(objectInfo *ein, int anzein)
    // Initialisierung des grids mit anzein Einschl�ssen
    {
    for(int nein=0;nein<anzein;nein++)
    {
    //  double x0 = ein[nein].P[0];
    //  double y0 = ein[nein].P[1];
    //  double z0 = ein[nein].P[2];
    //  double r0 = ein[nein].a;q
    
      for(int ix=0;ix<nxmax;ix++)
      {
        for(int iy=0;iy<nymax;iy++)
        {
          for(int iz=0;iz<nzmax;iz++)
          {
          if (gridarray[ix][iy][iz]==0)
          {
            if (inObject(ix,iy,iz,ein[nein]))
            {
              gridarray[ix][iy][iz] =  new maths::Vector<std::complex<double> >;
              *gridarray[ix][iy][iz] = maths::Vector<std::complex<double> >(1,1,1);
    //          cout << "Einschluss bei ix=" << ix << ", iy= " << iy << ", iz= " << iz << "\n";
    //          cout << "test:" <<  *gridarray[ix][iy][iz] << "\n";
            }
            else
              gridarray[ix][iy][iz] = 0;
          }
        }
        }
      }
    }
    }


    int grid::in_kugel(int ix, int iy, int iz)
    // Test, ob ein vorgegebener Punkt ix, iy, iz innerhalb des 
    // Partikels liegt
    {
    double x0 = xmax/2;
    double y0 = ymax/2;
    double z0 = zmax/2;
    double r0 = x0;
    int hilf = isInside(ix, iy, iz, x0, y0, z0, r0);
    return hilf;
    }

    int grid::inObject(int ix, int iy, int iz, objectInfo ein)
    {
    double r0 = ein.a*rP;
    int hilf = isInside(ix, iy, iz, ein.P*rP, r0);
    return hilf;
    }

    void grid::show_grid(char *fname)
    {
    std::ofstream os;
    os.open(fname);
    
    for (int ix=0; ix<nxmax; ix++)
    {
      for (int iy=0; iy<nymax; iy++)
      {
      for (int iz=0; iz<nzmax; iz++)
      {
      if (gridarray[ix][iy][iz]!=0)
        {

          os << (2*ix+1)*dx2 << " " << (2*iy+1)*dy2 << " " << (2*iz+1)*dz2 << "\n";
        }
      }
      }
    }
    os.close();
    }

    void grid::show_grid()
    {
    for (int ix=0; ix<nxmax; ix++)
    {
      for (int iy=0; iy<nymax; iy++)
      {
      for (int iz=0; iz<nzmax; iz++)
      {
    //   cout << "Gitt:" << gridarray[ix][iy][iz] << "\n";
      if (gridarray[ix][iy][iz]!=0)
        {
        std::cout << (2*ix+1)*dx2 << " " << (2*iy+1)*dy2 << " " << (2*iz+1)*dz2 << "\n";
        }
    //   else
    //     cout << "kein Einschluss" << "\n";
      }
      }
    }
    }

    maths::Vector<int> grid::gridpunkt(double x0, double y0, double z0)
    {
    maths::Vector<int> punkt;

    punkt[0] = (int)(x0/dx);
    punkt[1] = (int)(y0/dy);
    punkt[2] = (int)(z0/dz);

    return punkt;
    }

    maths::Vector<int> grid::gridpunkt(maths::Vector<double> &P)
    {
    maths::Vector<int> punkt;

    punkt[0] = (int)(P[0]/dx);
    punkt[1] = (int)(P[1]/dy);
    punkt[2] = (int)(P[2]/dz);

    return punkt;
    }
  }
}
