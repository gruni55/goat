
#ifndef DBL_EPSILON
#define DBL_EPSILON 2.22045E-16
#endif
//#define c 299792458          // Vakuum-Lichtgeschwindigkeit in m/s
#include <complex>
#include <math.h>
#include <stdlib.h>
#include "vector.h"
#include "resutil.h"
#include <iostream>
namespace GOAT
{
    namespace raytracing
    {

        std::ostream& operator << (std::ostream& os, GlobalParms& parms)
        {
            os << parms.nx << "  " << parms.ny << "   " << parms.alpha << "  " << parms.db << std::endl;
            os << parms.dx << "  " << parms.dy << "   " << parms.dxy << "  " << parms.n0 << std::endl;
            os << parms.bmax << "  " << parms.r0 << "   " << parms.r0end << "  " << parms.l0 << std::endl;
            os << parms.k0 << "  " << parms.pol << "   " << parms.ResRad << "   " << parms.ResAzi << std::endl;
            os << parms.numObj << "  " << parms.AngleTol << "  " << parms.evan << std::endl;
            os << parms.phase << "  " << parms.logscale << "  " << parms.tunneln << "  " << parms.ColMax << "  " << parms.ColMin
                << "  " << parms.EinX << std::endl;
            os << parms.AnzReflex << "  " << parms.AnzRays << std::endl;
            return os;
        }

        std::istream& operator >> (std::istream& is, GlobalParms& parms)
        {
            is >> parms.nx >> parms.ny >> parms.alpha >> parms.db;
            is >> parms.dx >> parms.dy >> parms.dxy >> parms.n0;
            is >> parms.bmax >> parms.r0 >> parms.r0end >> parms.l0;
            is >> parms.k0 >> parms.pol >> parms.ResRad >> parms.ResAzi;
            //is >> parms.numObj   >> parms.AngleTol   >> parms.evan >> endl;
            is >> parms.phase >> parms.logscale >> parms.tunneln >> parms.ColMax >> parms.ColMin;
            is >> parms.EinX >> parms.AnzReflex >> parms.AnzRays;
            return is;
        }

        double absc(const std::complex<double>& z)
        {
            return sqrt(real(z) * real(z) + imag(z) * imag(z));
        }
        /*
        complex<double> sqrt (const complex<double>& z)
        {
         double phi;
         double absz=absc(z);
         complex<double> Erg;
         if (absz==0) return complex<double>(0.0,0.0);
         phi=asin(imag(z)/absz);

         Erg=sqrt(absz) * complex<double>(cos(phi/2),sin(phi/2));

         return Erg;
        }

        complex<double> Ln (const complex<double>& z)
        {
         double x,y;
         x=real(z); y=imag(z);
         return complex<double>(0.5*log(x*x+y*y),asin(y/abs(z)));
        }

        complex<double> log (const complex<double>& z)
        {
         double phi;
         double absz=abs(z);
         phi=asin(imag(z)/absz);

         return complex<double>(log(absz),phi);
        }

        complex<double> exp (const complex<double>& z)
        {
         double x=exp(real(z));
         return complex<double>(x*cos(imag(z)),x*sin(imag(z)));
        }

        complex<double> sin(const complex<double>& z)
        {
         return (exp(I*z)-exp(-I*z))/(2.0*I);
        }

        complex<double> cos(const complex<double>& z)
        {
         return (exp(I*z)+exp(-I*z))/2.0;
        }
        */

        /*complex<double>  asin (const complex<double> & z)
        {
         complex<double>  Erg;
         Erg=-I*log(I*z+sqrt(1.0-z*z));
         return Erg;
        }*/

        std::complex<double>  acos(const std::complex<double>& z)
        {
            std::complex<double>  Erg;
            Erg = -I * log(z + sqrt(z * z - 1.0));
            return Erg;
        }

        std::complex<double> atanh(const std::complex<double>& z)
        {
            std::complex<double> h;
            h = 0.5 * (log(1.0 + z) - log(1.0 - z));
            return h;
        }

        void output(int nx, int ny, maths::Vector<std::complex<double> >** G)
        {
            for (int i = 0; i <= nx; i++)
                for (int j = 0; j <= ny; j++)
                    std::cout << "[" << i << "," << j << "]=" << G[i][j] << std::endl;
        }


        void init_Strahl(GlobalParms Parms, StrahlArray* Strahl)
        {
            double dPy = 2.0 * Parms.r0 / ANZ_STRAHLEN;

            int iS = 0;
            for (double Py = -Parms.bmax; Py <= Parms.bmax; Py += dPy)
            {
                Strahl[iS].S.pol = Parms.pol;
                if (Parms.pol == SENKRECHT) Strahl[iS].S.EAmp = maths::Vector<double>(0, 0, 1);
                else Strahl[iS].S.EAmp = maths::Vector<double>(0, 1, 0);
                Strahl[iS].S.E = 1.0;
                Strahl[iS].S.k = maths::Vector<double>(1, 0, 0);
                //dreh(ray[iS].S.k,Parms.alpha);
                //dreh(ray[iS].S.EAmp,Parms.alpha);
                Strahl[iS].S.phi = 0.0;
                Strahl[iS].P = maths::Vector<double>(-Parms.bmax, Py, 0);
                iS++;
            }
        }


        void sub_Gitter(GlobalParms parms, maths::Vector<std::complex<double> >** Erg,
            maths::Vector<std::complex<double> >** Gitter1,
            maths::Vector<std::complex<double> >** Gitter2)
        {
            // Hier werden die Gitter Erg,Gitter1 und Gitter2 addiert
            // und in Erg abgespeichert: Erg=Erg+Gitter1+Gitter2 

            for (int j = 0; j < parms.nx; j++)
                for (int k = 0; k < parms.ny; k++)
                {
                    Erg[j][k] = Gitter1[j][k] - Gitter2[j][k];
                }
        }


        void add_Gitter(GlobalParms parms, maths::Vector<std::complex<double> >** Erg,
            maths::Vector<std::complex<double> >** Gitter1,
            maths::Vector<std::complex<double> >** Gitter2)
        {
            // Hier werden die Gitter Erg,Gitter1 und Gitter2 addiert
            // und in Erg abgespeichert: Erg=Erg+Gitter1+Gitter2 

            for (int j = 0; j <= parms.nx; j++)
                for (int k = 0; k <= parms.ny; k++)
                {
                    Erg[j][k] = Erg[j][k] + Gitter1[j][k] + Gitter2[j][k];
                }
        }

        void clear(GlobalParms parms, maths::Vector<std::complex<double> >** Gitter)
        {
            for (int j = 0; j < parms.ny; j++)
                for (int k = 0; k < parms.nx; k++)
                {
                    Gitter[j][k] = maths::Vector<std::complex<double> >(0.0, 0.0, 0.0);
                }
        }

        void Delete(int n, maths::Vector<std::complex<double> >** Gitter)
        {
            for (int i = n - 1; i >= 0; i--)
                delete Gitter[i];
            delete Gitter;
        }

        void copy(StrahlInfo& dest, StrahlInfo src)
        {
            dest.E = src.E;
            dest.k = src.k;
            dest.EAmp = src.EAmp;
            dest.phi = src.phi;
        }

        maths::Vector <std::complex<double> >** newGitter(int n, int m)
        {
            maths::Vector <std::complex<double> >** Erg;
            Erg = (maths::Vector<std::complex<double> > **) calloc(n, sizeof(maths::Vector<std::complex<double> >*));
            for (int i = 0; i <= n; i++)
                Erg[i] = (maths::Vector<std::complex<double> >*) calloc(m, sizeof(maths::Vector<std::complex<double> >));
            for (int i = 0; i < n + 1; i++)
                for (int j = 0; j < m + 1; j++)
                    Erg[i][j] = maths::Vector<std::complex<double> >(0.0, 0.0, 0.0);
            return Erg;
        }

        void checkObjectIntersection(maths::Vector<double>& anf, const maths::Vector<double>& end,
            StrahlInfo& S, int numObj, objectInfo* Obj,
            maths::Vector<double>& Ps, int& Index)
            /*
              Sucht nach dem naechsten Einschluss auf dem Weg des Strahls zwischen
              Punkt anf und Punkt end.

              Eingabeparameter:
              anf, end: Anfangs- und Endpunkte des Suchintervalls
              S       : enthaelt Infos ueber den Verlauf des Strahls
              numObj  : Anzahl der Einschluesse
              Ein     : enthaelt Infos ueber die Einschluesse

              Rueckgabewerte:
              Ps      : Schnittpunkt mit der Einschlusskante
              Index   : Nummer des Einschlusses (-1, falls nicht gefunden)
            */

        {
            double B, C, D, l1, l2;
            double lmin = abs(end - anf);
            Index = -1;
            for (int i = 0; i < numObj; i++)
            {
                B = S.k * (anf - Obj[i].P);
                C = abs(anf - Obj[i].P);
                C = C * C - Obj[i].a * Obj[i].a;
                D = B * B - C;
                if (D >= 0.0)
                {
                    l1 = sqrt(D);
                    l2 = -B - l1;
                    l1 = -B + l1;
                    if (((l2 < l1) && (l2 > 0.0)) || (l1 < 0.0)) { l1 = l2; }
                    if ((lmin > l1) && (l1 > 0.0)) { lmin = l1; Index = i; }
                }
            }
            Ps = anf + lmin * S.k;
        }

        void checkObjectIntersection(double r0, maths::Vector<double>& anf, const maths::Vector<double>& end,
            StrahlInfo& S, int numObj, objectInfo* Obj,
            maths::Vector<double>& Ps, int& Index)
            /*
              Sucht nach dem naechsten Einschluss auf dem Weg des Strahls zwischen
              Punkt anf und Punkt end.

              Eingabeparameter:
              anf, end: Anfangs- und Endpunkte des Suchintervalls
              S       : enthaelt Infos ueber den Verlauf des Strahls
              numObj  : Anzahl der Einschluesse
              Ein     : enthaelt Infos ueber die Einschluesse
              (hier sind alle Werte fuer den Einschluss relativ zur Radius r0
               des Partikels angegeben !!!!!!!)

              Rueckgabewerte:
              Ps      : Schnittpunkt mit der Einschlusskante
              Index   : Nummer des Einschlusses (-1, falls nicht gefunden)
            */

        {
            double B, C, D, l1, l2;
            double lmin = abs(end - anf);
            Index = -1;
            for (int i = 0; i < numObj; i++)
            {
                B = S.k * (anf - Obj[i].P * r0);
                C = abs(anf - Obj[i].P * r0);
                C = C * C - Obj[i].a * Obj[i].a * r0 * r0;
                D = B * B - C;
                if (D >= 0.0)
                {
                    l1 = sqrt(D);
                    l2 = -B - l1;
                    l1 = -B + l1;
                    if (((l2 < l1) && (l2 > 0.0)) || (l1 < 0.0)) { l1 = l2; }
                    if ((lmin > l1) && (l1 > 0.0)) { lmin = l1; Index = i; }
                }
            }
            Ps = anf + lmin * S.k;
        }
        void checkObjectIntersection(double r0, maths::Vector<double>& anf, const maths::Vector<double>& end,
            const maths::Vector<double> k, int numObj, objectInfo* Obj,
            maths::Vector<double>& Ps, int& Index)
            /*
              Sucht nach dem naechsten Einschluss auf dem Weg des Strahls zwischen
              Punkt anf und Punkt end.

              Eingabeparameter:
              anf, end: Anfangs- und Endpunkte des Suchintervalls
              S       : enthaelt Infos ueber den Verlauf des Strahls
              numObj  : Anzahl der Einschluesse
              Ein     : enthaelt Infos ueber die Einschluesse
              (hier sind alle Werte fuer den Einschluss relativ zur Radius r0
               des Partikels angegeben !!!!!!!)

              Rueckgabewerte:
              Ps      : Schnittpunkt mit der Einschlusskante
              Index   : Nummer des Einschlusses (-1, falls nicht gefunden)
            */

        {
            double B, C, D, l1, l2;
            double lmin = abs(end - anf);
            Index = -1;
            for (int i = 0; i < numObj; i++)
            {
                B = k * (anf - Obj[i].P * r0);
                C = abs(anf - Obj[i].P * r0);
                C = C * C - Obj[i].a * Obj[i].a * r0 * r0;
                D = B * B - C;
                if (D >= 0.0)
                {
                    l1 = sqrt(D);
                    l2 = -B - l1;
                    l1 = -B + l1;
                    if (((l2 < l1) && (l2 > 0.0)) || (l1 < 0.0)) { l1 = l2; }
                    if ((lmin > l1) && (l1 > 0.0)) { lmin = l1; Index = i; }
                }
            }
            Ps = anf + lmin * k;
        }
        double minmax(double a, double b)
        {
            return (a > b) ? a : b;
        }

        maths::Vector<double> sph2kart(maths::Vector<double> v)
        {
            /*
               Trafo: sphaerisches -> kartesisches Koord.-System
            */

            maths::Vector <double> Erg;
            Erg[0] = sin(v[1]) * cos(v[2]);
            Erg[1] = sin(v[1]) * sin(v[2]);
            Erg[2] = cos(v[1]);
            Erg *= v[0];
            return Erg;
        }

        maths::Vector<double> kart2sph(maths::Vector<double> v)
        {
            /*
             Trafo: kartesisches  -> sphaerisches Koord.-System
            */

            maths::Vector<double> Erg;
            Erg[0] = abs(v);
            Erg[1] = asin(v[2] / Erg[0]);
            Erg[2] = atan(v[1] / v[0]);
            return Erg;
        }

        Point get_grid_point(const GlobalParms& parms, maths::Vector<double> P)
        {
            Point Erg;
            Erg.x = (int)((P[0] + parms.r0) / parms.dx + 0.5);
            Erg.y = (int)((P[1] + parms.r0) / parms.dy + 0.5);
            return Erg;
        }

        maths::Vector <double> set_grid_point(const GlobalParms& parms, Point P)
        {
            maths::Vector <double> Erg;
            Erg[0] = P.x * parms.dx - parms.r0;
            Erg[1] = P.y * parms.dy - parms.r0;
            return Erg;
        }

        double grad(const GlobalParms& parms, maths::Vector<std::complex<double> >** Gitter, const
            maths::Vector<double>& P)
        {
            Point G;
            double Erg;
            maths::Vector<std::complex<double> > r, rx, ry;
            G = get_grid_point(parms, P);
            //cout << "Gitter-> " << Gitter[G.x][G.y] << endl;
            r = Gitter[G.x][G.y];
            rx = Gitter[G.x + 1][G.y];
            ry = Gitter[G.x][G.y + 1];
            Erg = abs(r[0] - rx[0] + r[0] - ry[1] + r[1] - rx[0] + r[1] - ry[1] + r[2] - rx[0] + r[2] - ry[1]);
            return Erg;
        }

        /*Vector<double> next (const GlobalParms& parms,
                                       const Vector<double>& P0,
                                       const  Vector<double>& k)
        {

            sucht den nächsten Gitterpunkt, wenn P0 nicht (sicher) auf einer
            Gitterstrebe sitzt

            Wird von <next> verwendet.
            Parameter:
            parms : Infos ueber die Kugel (außen)
            P0    : Punkt an dem man sich gerade befindet
            k     : Einheitsvektor in Ausbreitungsrichtung


         double x,y;
         double x1,y1,x2,y2;
         double m,c;
         double abs1,abs2;
         Vector<double> Erg;

         x=P0[0]+parms.r0;
         y=P0[1]+parms.r0;

         x=x/(2.0*parms.r0)*parms.nx;
         y=y/(2.0*parms.r0)*parms.ny;


         cout << "x=" << x << "    y=" << y << "k=   " << k << endl;
         m=k[1]*parms.ny / (k[0] * parms.nx);
         c=y-m*x;

         if ((fabs(x-floor(x+0.5)<1E-10))||(fabs(y-floor(y+0.5))<1E-10))
         {
          if (k[0]<0.0) x-=0.1; else x+=0.1;
          y=m*x+c;
         }

         if (k[0]>0.0)
         {
            if (k[1]>0.0) { x1=ceil(x); y1=m*x1+c; y2=ceil(y); x2=(y2-c)/m; } // kx,ky>0
            else { x1=ceil(x); y1=m*x1+c;  y2=floor(y); x2=(y2-c)/m; } // kx>0,ky<0
         }
         else
         {
           if (k[1]>0.0) { x1=floor(x); y1=m*x1+c; y2=ceil(y); x2=(y2-c)/m; } // kx<0,ky>0
           else { x1=floor(x); y1=m*x1+c; y2=floor(y); x2=(y2-c)/m; } // kx<0,ky<0
         }

         abs1=sqrt(x1*x1+y1*y1);
         abs2=sqrt(x2*x2+y2*y2);
         if (abs2>abs1)
         {
          Erg= Vector<double> (x1,y1,0.0);
          if (abs(Erg-P0)<1E-10) Erg=Vector<double> (x2,y2,0.0);
         }
         else
          {
          Erg= Vector<double> (x2,y2,0.0);
          if (abs(Erg-P0)<1E-10) Erg=Vector<double> (x1,y1,0.0);
         }
         Erg[0]=Erg[0]/parms.nx*2.0*parms.r0-parms.r0;
         Erg[1]=Erg[1]/parms.ny*2.0*parms.r0-parms.r0;
         cout << "Erg=" << Erg << "    r0=" << parms.r0 << endl;
         return Erg;
        }
        */
        maths::Vector<double> next(const GlobalParms& parms, const maths::Vector<double>& P0,
            const maths::Vector<double>& k)
        {

            /*    berechnet den naechsten Schnittpunkt des Strahls mit einer Gitterlinie
                P0 : derzeitiger Ort
                k  : Ausbreitungsvektor

            */
            double m, c;
            double hx, hy;
            double x0, y0, x1, y1, x2, y2;
            maths::Vector<double> Erg;

            m = k[1] / k[0] * parms.dx / parms.dy;

            x0 = (P0[0] + parms.r0) / parms.dx;
            y0 = (P0[1] + parms.r0) / parms.dy;

            c = y0 - m * x0;

            hx = x0;
            hy = y0;

            if (fabs(hx - floor(hx + 0.5)) <= parms.dx / 1000.0)
            {
                // liegt auf einer Spalte
                if (fabs(hy - floor(hy + 0.5)) <= parms.dy / 1000.0)
                {
                    // liegt aber auch auf einer Zeile -> Kreuzpunkt getroffen 
                    if (k[0] >= 0.0)
                    {
                        if (k[1] >= 0.0) { x1 = x0 + 1.0; y1 = m * x1 + c; y2 = y0 + 1.0; x2 = (y2 - c) / m; }
                        else { x1 = x0 + 1.0; y1 = m * x1 + c; y2 = y0 - 1.0; x2 = (y2 - c) / m; }
                    }
                    else
                    {
                        if (k[1] >= 0.0) { x1 = x0 - 1.0; y1 = m * x1 + c; y2 = y0 + 1.0; x2 = (y2 - c) / m; }
                        else { x1 = x0 - 1.0; y1 = m * x1 + c; y2 = y0 - 1.0; x2 = (y2 - c) / m; }
                    }
                }
                else // Spalte
                {
                    if (k[0] >= 0.0)
                    {
                        if (k[1] >= 0.0) { y1 = ceil(y0);  x1 = (y1 - c) / m; x2 = x0 + 1.0; y2 = m * x2 + c; }
                        else { y1 = floor(y0); x1 = (y1 - c) / m; x2 = x0 + 1.0; y2 = m * x2 + c; }
                    }
                    else
                    {
                        if (k[1] >= 0.0) { y1 = ceil(y0);  x1 = (y1 - c) / m; x2 = x0 - 1.0; y2 = m * x2 + c; }
                        else { y1 = floor(y0); x1 = (y1 - c) / m; x2 = x0 - 1.0; y2 = m * x2 + c; }
                    }
                }
            }
            else
            {
                if (fabs(hy - floor(hy + 0.5)) <= parms.dy / 1E+7)
                {
                    // liegt auf einer Zeile
                    if (k[0] >= 0.0)
                    {
                        if (k[1] >= 0.0) { x1 = ceil(x0); y1 = m * x1 + c; y2 = y0 + 1.0; x2 = (y2 - c) / m; }
                        else { x1 = ceil(x0); y1 = m * x1 + c; y2 = y0 - 1.0; x2 = (y2 - c) / m; }
                    }
                    else
                    {
                        if (k[1] >= 0.0) { x1 = floor(x0); y1 = m * x1 + c; y2 = y0 + 1.0; x2 = (y2 - c) / m; }
                        else { x1 = floor(x0); y1 = m * x1 + c; y2 = y0 - 1.0; x2 = (y2 - c) / m; }
                    }
                }

                else
                {
                    // liegt irgendwo in der Zelle

                    if (k[0] >= 0.0)
                    {
                        if (k[1] >= 0.0) { x1 = ceil(x0);    y1 = m * x1 + c; y2 = ceil(y0); x2 = (y2 - c) / m; }
                        else { y1 = floor(y0); x1 = (y1 - c) / m; x2 = ceil(x0);   y2 = x2 * m + c; }
                    }
                    else
                    {
                        if (k[1] >= 0.0) { y1 = ceil(y0);  x1 = (y1 - c) / m; x2 = floor(x0); y2 = m * x2 + c; }
                        else { y1 = floor(y0); x1 = (y1 - c) / m; x2 = floor(x0); y2 = m * x2 + c; }
                    }
                }
            }
            if ((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0) < (x2 - x0) * (x2 - x0) + (y2 - y0) * (y2 - y0))
                Erg = maths::Vector<double>(x1 * parms.dx - parms.r0, y1 * parms.dy - parms.r0, 0.0);
            else
                Erg = maths::Vector<double>(x2 * parms.dx - parms.r0, y2 * parms.dy - parms.r0, 0.0);
            return Erg;
        }

        void minmax(double x, double dx, int& min, int& max)
        {
            int i = int(x / dx);
            if (fabs((double)(i + 1.0) * dx - x) <= 1.0E-10) i++;
            if (fabs((double)i * dx - x) <= 1.0E-10)
            {
                min = i - 1;
                max = i + 1;
            }
            else
            {
                min = i;
                max = i + 1;
            }
        }
        maths::Vector<double> nextP(maths::Vector<double> P, maths::Vector<double> k, maths::Vector<double> OK, double rK, bool& found)
            /*
              Berechnet den nächsten Schnittpunkt eines Strahls der Richtung k mit einer Kugel, die sich am Ort OK
              befindet und einen Radius rK hat. P ist der Ort an dem man sich gerade befindet
            */
        {
            maths::Vector<double> Ps, Erg;
            double l1, l2, l;
            double P2, Pk, k2;
            Ps = P - OK;
            k2 = k * k;
            Pk = Ps * k;
            P2 = Ps * Ps;
            l1 = (Pk * Pk - P2 * k2 + rK * rK * k2);
            if (l1 < 0.0) { found = false; return maths::Vector<double>(INF, INF, INF); }
            l2 = (-Pk + sqrt(l1)) / k2;
            l1 = (-Pk - sqrt(l1)) / k2;
            if (l1 <= 1E-10) l = l2; else l = l1;
            if (l <= 0.0) { found = false; return maths::Vector<double>(INF, INF, INF); }
            Erg = P + l * k;
            found = true;
            return Erg;
        }

        std::ostream& operator << (std::ostream& os, objectInfo E)
        {
            os << "P=" << E.P << std::endl;
            os << "n=" << E.n << "     a=" << E.a << std::endl;
            return os;
        }

        bool operator == (objectInfo a, objectInfo b)
        {
            return (a.P == b.P) && (a.n == b.n) && (a.a == b.a) && (a.alpha == b.alpha);
        }

        std::ostream& savebinGlobalParms(std::ostream& os, GlobalParms parms)
        {
            os.write((char*)&parms.nx, (char)sizeof(parms.nx));
            os.write((char*)&parms.ny, (char)sizeof(parms.ny));
            os.write((char*)&parms.alpha, (char)sizeof(parms.alpha));
            os.write((char*)&parms.db, (char)sizeof(parms.db));
            os.write((char*)&parms.dx, (char)sizeof(parms.dx));
            os.write((char*)&parms.dy, (char)sizeof(parms.dy));
            os.write((char*)&parms.dxy, (char)sizeof(parms.dxy));
            os.write((char*)&parms.AnzReflex, (char)sizeof(parms.AnzReflex));
            os.write((char*)&parms.AnzRays, (char)sizeof(parms.AnzRays));
            os.write((char*)&parms.bmax, (char)sizeof(parms.bmax));
            os.write((char*)&parms.r0, (char)sizeof(parms.r0));
            os.write((char*)&parms.r0end, (char)sizeof(parms.r0end));
            os.write((char*)&parms.l0, (char)sizeof(parms.l0));
            os.write((char*)&parms.k0, (char)sizeof(parms.k0));
            os.write((char*)&parms.n0, (char)sizeof(parms.n0));
            os.write((char*)&parms.pol, (char)sizeof(parms.pol));
            os.write((char*)&parms.ResRad, (char)sizeof(parms.ResRad));
            os.write((char*)&parms.ResAzi, (char)sizeof(parms.ResAzi));
            os.write((char*)&parms.numObj, (char)sizeof(parms.numObj));
            os.write((char*)&parms.AngleTol, (char)sizeof(parms.AngleTol));
            os.write((char*)&parms.evan, (char)sizeof(parms.evan));
            os.write((char*)&parms.PolAngle, (char)sizeof(parms.PolAngle));
            os.write((char*)&parms.phase, (char)sizeof(parms.phase));
            os.write((char*)&parms.logscale, (char)sizeof(parms.logscale));
            os.write((char*)&parms.tunneln, (char)sizeof(parms.tunneln));
            os.write((char*)&parms.ColMax, (char)sizeof(parms.ColMax));
            os.write((char*)&parms.ColMin, (char)sizeof(parms.ColMin));
            os.write((char*)&parms.EinX, (char)sizeof(parms.EinX));
            return os;
        }

        std::istream& loadbinGlobalParms(std::istream& is, GlobalParms& parms)
        {
            is.read((char*)&parms.nx, (char)sizeof(parms.nx));
            std::cout << "parms.nx=" << parms.nx << std::endl;
            is.read((char*)&parms.ny, (char)sizeof(parms.ny));
            is.read((char*)&parms.alpha, (char)sizeof(parms.alpha));
            is.read((char*)&parms.db, (char)sizeof(parms.db));
            is.read((char*)&parms.dx, (char)sizeof(parms.dx));
            is.read((char*)&parms.dy, (char)sizeof(parms.dy));
            is.read((char*)&parms.dxy, (char)sizeof(parms.dxy));
            is.read((char*)&parms.AnzReflex, (char)sizeof(parms.AnzReflex));
            is.read((char*)&parms.AnzRays, (char)sizeof(parms.AnzRays));
            is.read((char*)&parms.bmax, (char)sizeof(parms.bmax));
            is.read((char*)&parms.r0, (char)sizeof(parms.r0));
            is.read((char*)&parms.r0end, (char)sizeof(parms.r0end));
            is.read((char*)&parms.l0, (char)sizeof(parms.l0));
            is.read((char*)&parms.k0, (char)sizeof(parms.k0));
            is.read((char*)&parms.n0, (char)sizeof(parms.n0));
            is.read((char*)&parms.pol, (char)sizeof(parms.pol));
            is.read((char*)&parms.ResRad, (char)sizeof(parms.ResRad));
            is.read((char*)&parms.ResAzi, (char)sizeof(parms.ResAzi));
            is.read((char*)&parms.numObj, (char)sizeof(parms.numObj));
            is.read((char*)&parms.AngleTol, (char)sizeof(parms.AngleTol));
            is.read((char*)&parms.evan, (char)sizeof(parms.evan));
            is.read((char*)&parms.PolAngle, (char)sizeof(parms.PolAngle));
            is.read((char*)&parms.phase, (char)sizeof(parms.phase));
            is.read((char*)&parms.logscale, (char)sizeof(parms.logscale));
            is.read((char*)&parms.tunneln, (char)sizeof(parms.tunneln));
            is.read((char*)&parms.ColMax, (char)sizeof(parms.ColMax));
            is.read((char*)&parms.ColMin, (char)sizeof(parms.ColMin));
            is.read((char*)&parms.EinX, (char)sizeof(parms.EinX));
            return is;
        }

        RRTParmsInfo readRRTParms(bool old, std::ifstream* is)
        {
            RRTParmsInfo erg;
            is->read((char*)&erg.Ebene, (char)sizeof(erg.Ebene));
            is->read((char*)&erg.Pol, (char)sizeof(erg.Pol));
            is->read((char*)&erg.angmin, (char)sizeof(erg.angmin));
            is->read((char*)&erg.angmax, (char)sizeof(erg.angmax));
            is->read((char*)&erg.nang, (char)sizeof(erg.nang));
            is->read((char*)&erg.wave, (char)sizeof(erg.wave));
            is->read((char*)&erg.isKoherent, (char)sizeof(erg.isKoherent));
            return erg;
        }

        RRTParmsInfo readRRTParms(bool old, std::ifstream& is)
        {
            RRTParmsInfo erg;
            is.read((char*)&erg.Ebene, (char)sizeof(erg.Ebene));
            is.read((char*)&erg.Pol, (char)sizeof(erg.Pol));
            is.read((char*)&erg.angmin, (char)sizeof(erg.angmin));
            is.read((char*)&erg.angmax, (char)sizeof(erg.angmax));
            is.read((char*)&erg.nang, (char)sizeof(erg.nang));
            is.read((char*)&erg.wave, (char)sizeof(erg.wave));
            is.read((char*)&erg.isKoherent, (char)sizeof(erg.isKoherent));
            return erg;
        }

        void writeRRTParms(std::ofstream& os, RRTParmsInfo erg)
        {
            os.write((char*)&erg.Ebene, (char)sizeof(erg.Ebene));
            os.write((char*)&erg.Pol, (char)sizeof(erg.Pol));
            os.write((char*)&erg.angmin, (char)sizeof(erg.angmin));
            os.write((char*)&erg.angmax, (char)sizeof(erg.angmax));
            os.write((char*)&erg.nang, (char)sizeof(erg.nang));
            os.write((char*)&erg.wave, (char)sizeof(erg.wave));
            os.write((char*)&erg.isKoherent, (char)sizeof(erg.isKoherent));
        }

        void readRRTParms(std::ifstream& is, RRTParmsInfo& erg)
        {
            is.read((char*)&erg.Ebene, (char)sizeof(erg.Ebene));
            is.read((char*)&erg.Pol, (char)sizeof(erg.Pol));
            is.read((char*)&erg.angmin, (char)sizeof(erg.angmin));
            is.read((char*)&erg.angmax, (char)sizeof(erg.angmax));
            is.read((char*)&erg.nang, (char)sizeof(erg.nang));
            is.read((char*)&erg.wave, (char)sizeof(erg.wave));
            is.read((char*)&erg.isKoherent, (char)sizeof(erg.isKoherent));
        }

        GlobalParms readGlobalParms(bool old, std::ifstream& is)
        {
            char dummy;
            GlobalParms erg;
            is.read((char*)&erg.nx, (char)sizeof(erg.nx));
            is.read((char*)&erg.ny, (char)sizeof(erg.ny));
            is.read((char*)&erg.alpha, (char)sizeof(erg.alpha));
            is.read((char*)&erg.db, (char)sizeof(erg.db));
            is.read((char*)&erg.dx, (char)sizeof(erg.dx));
            is.read((char*)&erg.dy, (char)sizeof(erg.dy));
            is.read((char*)&erg.dxy, (char)sizeof(erg.dxy));
            is.read((char*)&erg.AnzReflex, (char)sizeof(erg.AnzReflex));
            is.read((char*)&erg.AnzRays, (char)sizeof(erg.AnzRays));
            is.read((char*)&erg.bmax, (char)sizeof(erg.bmax));
            is.read((char*)&erg.r0, (char)sizeof(erg.r0));
            is.read((char*)&erg.r0end, (char)sizeof(erg.r0end));
            is.read((char*)&erg.l0, (char)sizeof(erg.l0));
            is.read((char*)&erg.k0, (char)sizeof(erg.k0));
            is.read((char*)&erg.n0, (char)sizeof(erg.n0));
            is.read((char*)&erg.pol, (char)sizeof(erg.pol));
            is.read((char*)&erg.ResRad, (char)sizeof(erg.ResRad));
            is.read((char*)&erg.ResAzi, (char)sizeof(erg.ResAzi));
            is.read((char*)&erg.numObj, (char)sizeof(erg.numObj));
            is.read((char*)&erg.AngleTol, (char)sizeof(erg.AngleTol));
            is.read((char*)&erg.evan, (char)sizeof(erg.evan));
            is.read((char*)&erg.PolAngle, (char)sizeof(erg.PolAngle));
            is.read((char*)&erg.phase, (char)sizeof(erg.phase));
            is.read((char*)&erg.logscale, 1);
            if (old) { is.read((char*)&dummy, 1); is.read((char*)&dummy, 1); }
            is.read((char*)&erg.tunneln, 1);
            if (old) { is.read((char*)&dummy, 1); is.read((char*)&dummy, 1); }
            is.read((char*)&erg.ColMax, (char)sizeof(erg.ColMax));
            is.read((char*)&erg.ColMin, (char)sizeof(erg.ColMin));
            is.read((char*)&erg.EinX, (char)sizeof(erg.EinX));
            return erg;
        }

        GlobalParms readGlobalParms(bool old, std::ifstream* is)
        {
            char dummy;
            GlobalParms erg;
            is->read((char*)&erg.nx, (char)sizeof(erg.nx));
            is->read((char*)&erg.ny, (char)sizeof(erg.ny));
            is->read((char*)&erg.alpha, (char)sizeof(erg.alpha));
            is->read((char*)&erg.db, (char)sizeof(erg.db));
            is->read((char*)&erg.dx, (char)sizeof(erg.dx));
            is->read((char*)&erg.dy, (char)sizeof(erg.dy));
            is->read((char*)&erg.dxy, (char)sizeof(erg.dxy));
            is->read((char*)&erg.AnzReflex, (char)sizeof(erg.AnzReflex));
            is->read((char*)&erg.AnzRays, (char)sizeof(erg.AnzRays));
            is->read((char*)&erg.bmax, (char)sizeof(erg.bmax));
            is->read((char*)&erg.r0, (char)sizeof(erg.r0));
            is->read((char*)&erg.r0end, (char)sizeof(erg.r0end));
            is->read((char*)&erg.l0, (char)sizeof(erg.l0));
            is->read((char*)&erg.k0, (char)sizeof(erg.k0));
            is->read((char*)&erg.n0, (char)sizeof(erg.n0));
            is->read((char*)&erg.pol, (char)sizeof(erg.pol));
            is->read((char*)&erg.ResRad, (char)sizeof(erg.ResRad));
            is->read((char*)&erg.ResAzi, (char)sizeof(erg.ResAzi));
            is->read((char*)&erg.numObj, (char)sizeof(erg.numObj));
            is->read((char*)&erg.AngleTol, (char)sizeof(erg.AngleTol));
            is->read((char*)&erg.evan, (char)sizeof(erg.evan));
            is->read((char*)&erg.PolAngle, (char)sizeof(erg.PolAngle));
            is->read((char*)&erg.phase, (char)sizeof(erg.phase));
            is->read((char*)&erg.logscale, 1);
            if (old) is->read((char*)&dummy, 1);
            is->read((char*)&erg.tunneln, 1);
            if (old) is->read((char*)&dummy, 1);
            is->read((char*)&erg.ColMax, (char)sizeof(erg.ColMax));
            is->read((char*)&erg.ColMin, (char)sizeof(erg.ColMin));
            is->read((char*)&erg.EinX, (char)sizeof(erg.EinX));
            return erg;
        }

        GlobalParms::GlobalParms()
        {
            r0 = 1.0;
            nx = 100;
            ny = 100;
            alpha = 0.0;
            dx = r0 / ((double)nx);
            dy = r0 / ((double)ny);
            dxy = sqrt(dx * dx + dy * dy);
            db = 0.0;
            AnzReflex = 3;
            AnzRays = 300;
            bmax = 0;
            r0end = r0;
            l0 = 1.0;
            k0 = 2.0 * M_PI / l0;
            n0 = 1.0;
            pol = SENKRECHT;
            ResRad = 1;
            ResAzi = 1;
            numObj = 0;
            AngleTol = 1E-5;
            evan = 0;
            PolAngle = 0;
            phase = 0;
            logscale = false;
            tunneln = false;
            ColMax = 0;
            ColMin = 0;
            EinX = 0;
            nOrientAvgAlpha = 2;
            nOrientAvgBeta = 2;
            nOrientAvgGamma = 2;
        }

        RRTParmsInfo::RRTParmsInfo()
        {
            Ebene = VAR_XZ;
            Pol = SENKRECHT;
            Strahlungsart = RAMAN;
            wave = 1.0;
            angmin = 0;
            angmax = M_PI;
            nang = 91;
            isKoherent = false;
        }
    }
}
