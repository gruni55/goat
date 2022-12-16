
#pragma once

#include <stdlib.h>
#include <fstream>
#include "vector.h"
#include "resutil.h"
#include "objectshape.h"

namespace GOAT
{
   namespace raytracing 
   {
      /** 
       * This class provides a virtual 3D-grid of complex vectors which can be used to store e.g. the electric field in a volume. 
       *  It circumscribes a sphere with radius r0 and virtual means here, that only  the parts which are needed were allocated, 
       *  but for the programmer it seems like a whole array or grid. 
       *  By adding an object to the grid, the memory for the grid around the object is allocated, so one can store the electric field 
       * inside the object. The cells of the grid can be addressed with help of the bracket()-operator with different arguments 
       **/ 
      class SuperArray
      {
         public :
         SuperArray ();
         /**
          * @brief Constructor for SuperArray
          * @param r0 Radius of the circumscribed sphere
          * @p nx, @p ny, @p nz number of cells in the x-, y-, z- direction
          * @param type describes wether the whole space of the grid is allocated (IN_HOST) or in the objects only (IN_OBJECT) 
         */
         SuperArray(double r0,int nx,int ny,int nz, const int typ=IN_HOST); 
         SuperArray(double r0,int nx, int ny, int nz,ObjectShape **Obj, int anzEin, const bool isAbsolute=false, const int typ=0);
         ~SuperArray () { clear();};

         /**
          * 
         */
         void copy (const SuperArray &S);
         /**
          * 
          *  
          */ 
         void add (const SuperArray& S);
         void sub (const SuperArray& S);
         friend double sumabs (const SuperArray &S);
         friend double sumabs2 (const SuperArray &S);
         friend double abs2sum (const SuperArray &S);
         bool addInc(ObjectShape **Obj,int anzEin, const bool isAbsolute=false); // anzEin Einschluesse hinzufuegen
                                           // isAbsolute=true => Orts-/Groessenangaben in absoluten Koordinaten 
         bool addInc (ObjectShape *E,const bool isAbsolute=false);  // Einschluss hinzufuegen (isAbsolute s.o.) 
         maths::Vector<int> gitterpunkt (maths::Vector<double> P);
         bool inObject (maths::Vector<double> P, int i);
         bool inObject (maths::Vector<int> Pi, int i);
         bool inObject (int ix, int iy, int iz, int i);
         maths::Vector<std::complex<double> >& operator () (int ix, int iy, int iz); ///< gives the content of the cell[ix][iy][iz]
         /**
          * @brief returns the content of the i-th object, with the grid-coordinates ix,iy,iz
          * @param i number of the object
          * @param ix 
         */
         maths::Vector<std::complex<double> >& operator () (int i, int ix, int iy, int iz, bool isObjectCoordinate=true);
         maths::Vector<std::complex<double> >& operator () (maths::Vector<int> Pi);
         maths::Vector<std::complex<double> >& operator () (int i,maths::Vector<int> Pi);
         maths::Vector<std::complex<double> >& operator () (maths::Vector<double> P);
         maths::Vector<std::complex<double> >& operator () (int i, maths::Vector<double> P);
         void saveExPhase (char* FName,int i=0);
         void saveEyPhase (char * FName,int i=0);
         void saveEzPhase (char* FName,int i=0);
         void saveExPol (char* FName,int i=0);
         void saveEyPol (char* FName,int i=0);
         /**
          * @brief store the absolute value of the field inside the i-th object
          * @param FName name of the file where the data should be stored
          * @param i number of the object 
         */
         void saveEzPol (char* FName,int i=0);
         /**
          * @brief store the absolute value of the field inside the i-th object
          * @param FName name of the file where the data should be stored
          * @param i number of the object 
         */
         void saveabsE (const char* FName, int i=0); 
         /**
          * @brief store the full electric field inside the i-th object
          * 
          * The full electric field is stored in the following way:
          * real(Ex(x0,y0,z0)) imag(Ex(x0,y0,z0)) real(Ey(x0,y0,z0)) imag(Ey(x0,y0,z0)) real(Ez(x0,y0,z0)) imag(Ez(x0,y0,z0)) 
          * real(Ex(x0,y0,z1)) imag(Ex(x0,y0,z1)) ...
          *    :                      :  ...
          * real(Ex(x0,y0,zN)) imag(Ex(x0,y0,zN)) ...
          * real(Ex(x0,y1,z0)) imag(Ex(x0,y1,z0)) ...
          *    :                      :
          * real(Ex(x0,y1,zN)) imag(Ex(x0,y1,zN)) ...
          *    :                      :
          * 
          * @param FName name of the file where the data should be stored
          * @param i number of the object 
         */
         void saveFullE(const char* FName, int i=0);
         void makeReal ();
         void fill(const maths::Vector<std::complex<double> > &x);
         SuperArray& operator = (const SuperArray &S);
         void clear();
         void allockugel();
         maths::Vector<int> kugelindex(maths::Vector<int> Pi);
         maths::Vector<std::complex<double> > kugelwert(maths::Vector<int> Pi);
         maths::Vector<std::complex<double> > kugelwert(int ix, int iy, int iz);
         friend std::ostream&   operator << (std::ostream &os, const SuperArray &S);
 
 /* 
    Hier kommen die eigentlichen Daten: 
    Ein : Array mit den Einschluessen (Orte und Ausdehnungen der Einschluesse werden absolut angegeben)
    anzEin : Anzahl Einschluesse
    type: Verteilung der aktiven Molek�le
          0: nur im Einschluss
          1: nur im Host
          2: �berall
    ywerte,
    zwerte : Arrays zur Bestimmung der Dimensionen von K
    r0  : Radius des Host-Partikels=halbe Breite des Uebergitters
    G   : Array mit den Untergittern, wird ben�tigt falls type 0 ist; ObjectShape : G[Einschlussindex][ix][iy][iz]
    K   : Kugelgitter, wird ben�tigt falls type 1 oder 2 ist
    Pul : Array mit den Basispunkten der Untergitter 
    n   : Array mit den Ausdehnungen der    "   n[Einschlussindex][nx][ny][nz] (nx,ny,nz sind als Vektor zusammengefasst)
    nges : Ausdehnung des Uebergitters als Vektor (nx,ny,nz) zusammengefasst 
    d   : Gitterkonstanten als Vektor (dx,dy,dz) zusammengefasst
 */ 
  


         ObjectShape **Obj;
         int anzEin,type;
         int* ywerte;
         int **zwerte;
         maths::Vector<std::complex<double> >  ****G;
         maths::Vector<std::complex<double> >  ***K;
         maths::Vector<std::complex<double> > dummy;
         maths::Vector<int> *Pul,*n,nges;
         maths::Vector<double> d;    
         double r0;
         bool isequal;
         bool iscleared;
         maths::Vector<std::complex<double> > pc;
         //  int Fehler;
         maths::Matrix<double> H, R;
};
   }
}