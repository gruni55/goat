#ifndef __grid_H__
#define __grid_H__

#include <complex>
#include "vector.h"
#include "resutil.h"

namespace GOAT
{
    namespace raytracing
    {

        class grid
        {

        maths::Vector <std::complex<double> > ****gridarray;
        maths::Vector <std::complex<double> > ****newarray();
        void delarray();

        // Testet allgemein, ob eine gridzelle ix, iy, iz innerhalb einer 
        // Kugel mit Radius r0 und Mittelpunkt x0, y0, z0 liegt
        // wird von in_kugel und in_einschluss benoetigt

        int ist_innen(int ix, int iy, int iz,
                    double x0, double y0, double z0, double r0);
        int ist_innen(int ix, int iy, int iz, maths::Vector<double> P, double r0);

        double dx2, dy2, dz2;

        public:

        double rP;
        int nxmax,nymax,nzmax;
        double xmax, ymax,  zmax, dx, dy, dz;
        maths::Vector<std::complex<double> > DUMMY;

        grid();
        grid(int nx, int ny, int nz, double xmax, double ymax, double zmax);
        grid(int nx, int ny, int nz, double r);

        ~grid();

        maths::Vector<std::complex<double> >& operator() (int  i, int j, int k);
        maths::Vector<std::complex<double> >& operator() (maths::Vector<int>& P);

                        grid& operator=  (const grid& g);

        // Initialisiert ein grid mit anzein Einschluessen, die in ein
        // stehen
        void init_grid(objectInfo *ein, int anzein);

        // Erzeugt leeres grid
        void init_grid();

        // Ausgabe des grids auf Standardausgabe
        void show_grid();

        // Ausgabe des grids in Datei fname
        void show_grid(char* fname);

        // Gibt die Koordinaten der gridzelle aus, in der der Punkt mit
        // den  Koordinaten x0, y0 und z0 liegt
        maths::Vector<int> gridpunkt(double x0, double y0, double z0);
        maths::Vector<int> gridpunkt(maths::Vector<double>& P);
        // Testet, ob die gridzelle ix, iy, iz innerhalb des Partikels liegt
        int  in_kugel(int ix, int iy, int iz);

        // Testet, ob die gridzelle ix, iy, iz innerhalb des Einschlusses ein
        // liegt
        int  in_einschluss(int ix, int iy, int iz, objectInfo ein);

        void set_parms(int nx, int ny, int nz, int xmax, int ymax, int zmax);
        };
    }
}
#endif /*  __grid_H__ */
