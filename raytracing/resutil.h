#include <complex>
#pragma once
#ifndef complex_I
#define complex_I
//std::complex<double> I=std::complex<double>(0.0,1.0);    // imaginaere Einheit
#define I std::complex<double>(0.0,1.0)
#endif

#ifndef INF
#define INF 1.0/0.0
#endif
#ifndef RESUTIL_H
#define RESUTIL_H 

#define SENKRECHT  1
#define PARALLEL   0
#define MULTI      1
#define ANZ_STRAHLEN 500     // Anzahl der einfallenden Strahlen

#define IS_NOTHING 0
#define IS_ROW 1
#define IS_COL 2

#include "vector.h"
#include "matrix.h"
#define SQRT3  1.73205080756887729352744634150587 // Wurzel aus 3
#define RAMAN 0
#define FLUORESZENZ 1

#define IN_INC  1 
#define IN_HOST 2
#define IN_INC_AND_HOST 3

#ifndef VAR_XZ 
#define VAR_XZ   1
#endif
#ifndef VAR_YZ 
#define VAR_YZ   2
#endif
#ifndef VAR_XY
#define VAR_XY   3
#endif


typedef struct 
{
 int x,y;
} Point;

typedef struct 
{
 char Name[255];
 double n;
} OptProp;


typedef struct 
{
 /*
    E    : Richtungsvektor des Felds
    EAmp : Amplitude
    k    : Ausbreitungsvektor (normiert)
    phi  : Phase
    pol  : Polarisationsrichtung (senkrecht oder parallel)
    alpha : Winkel von einem Reflexionspunkt zur naechsten  (vom
            Partikelmittelpkt. gerechnet
    m    : Anzahl der bisher durchgefuehrten Strahlumlaeufe
     */


 bool getunnelt;
 int pol;
 Vector<double> k,EAmp;
 std::complex<double> E;
 std::complex<double> phi;
 int m; 
 double b; 
} StrahlInfo;

typedef struct 
{
 
 Vector <double> P;
 StrahlInfo      S; 
} StrahlArray;

typedef struct 
{
 /*
    P  : Ort des Einschlusses 
    n  : Brechungsindex  
    a  : Radius (in Einheiten des Partikelradius´)
    alpha : Polarisierbarkeit;
 */
 Vector <double> P;
 std::complex<double> n;
 double a;
 Matrix<double> alpha;  
} EinschlussInfo;

class RRTParmsInfo
{
public :
	RRTParmsInfo();
	int Ebene;
	int Pol;
	int Strahlungsart;
	double angmin, angmax;
	int nang;
	double wave;
	bool isKoherent;
};


class GlobalParms
{
	/*
	nx,ny     : Anzahl der Gitterpunkte in x/y-Richtung
	alpha     : Einfallswinkel (relativ zur x-Achse)
	AnzReflex : Anzahl der max. durch zufuehrenden Reflexionen an der Oberflaeche
	pro Strahl
	AnzRays   : Anzahl Strahlen
	dx, dy    : Breite und Hoehe einer Gitterzelle
	dxy       : Diagonale einer Gitterzelle
	r0        : Groesse des Partikels
	l0        : (Vakuum-)Wellenlaenge
	n0        : Brechungsindex des Partikels
	pol       : Polarisation des einfallenden Strahls
	ResRad    : radiale Modenzahl
	ResAzi    : azimutale Modenzahl (phi-Richtung)
	numObj    : Anzahl Einschluesse
	ColMin, ColMax : Minimal-/Maximalwert bei Farbeinteilung (in % des Maximums)
	*/

 public :
	 GlobalParms();

	 /// Radius der "Weltkugel"
	 double r0;

	 /// Anzahl der Gitterpunkte in x/y-Richtung
	 int nx, ny;
	 /// Einfallswinkel (relativ zur x-Achse)
	 double alpha;
	 /// Breite einer Gitterzelle  
	 double dx;
	 /// Hoehe einer Gitterzelle  
	 double dy;
	 /// Diagonale einer Gitterzelle
	 double dxy;
	 double db;
	 /// Anzahl der max. durch zufuehrenden Reflexionen an der Oberflaeche pro Strahl
	 int AnzReflex;
	 /// Anzahl Strahlen
	 int AnzRays;
	 double bmax, r0end, l0, k0;
	 /// Brechungsindex des Umgebungsmediums
	 std::complex<double> n0;
	 /// Polarisationsrichtung: kann die Werte SENKRECHT und PARALLEL annehmen
	 int pol;
	 int ResRad, ResAzi;
	 int numObj;
	 double AngleTol, evan;
	 double PolAngle;
	 int phase;
	 bool logscale;
	 bool tunneln;
	 double ColMax, ColMin;
	 int EinX;
	 int nOrientAvgAlpha, nOrientAvgBeta, nOrientAvgGamma;
};



std::ostream& operator << (std::ostream& os, GlobalParms& parms);
std::istream& operator >> (std::istream& is, GlobalParms& parms);

double abs2(double x);

void output (int nx, int ny, Vector<std::complex<double> > **G);
void init_Strahl (GlobalParms Parms, StrahlArray *Strahl);
void sub_Gitter (GlobalParms parms, Vector<std::complex<double> > **Erg, 
                Vector<std::complex<double> > **Gitter1,
                Vector<std::complex<double> > **Gitter2);
void add_Gitter (GlobalParms parms, Vector<std::complex<double> > **Erg, 
                Vector<std::complex<double> > **Gitter1,
                Vector<std::complex<double> > **Gitter2);
void clear (GlobalParms parms, Vector<std::complex<double> > **Gitter);
void Delete (int n, Vector<std::complex<double> > **Gitter);
void copy (StrahlInfo &dest, StrahlInfo src);
double minmax (double a, double b);
Vector<double> kart2sph (Vector<double> v);
Vector<double> sph2kart (Vector<double> v);

Point get_grid_point (const GlobalParms &parms,Vector<double> P);
Vector <double> set_grid_point (const GlobalParms& parms, Point P);
double grad (const GlobalParms& parms, Vector<std::complex<double> > **Gitter, const
      Vector<double> &P);
//std::complex<double>  asin(const std::complex<double> &);
std::complex<double>  acos(const std::complex<double> &);
Vector<double> next (const GlobalParms& parms, 
                               const Vector<double>& P0,
                               const  Vector<double>& k);

void minmax (double x, double dx, int &min, int &max);
void checkObjectIntersection (Vector<double>& anf, const Vector<double>& end,
                      StrahlInfo& S, int numObj, EinschlussInfo *Obj,
		      Vector<double>& Ps, int &Index);
void checkObjectIntersection (double r0, Vector<double>& anf, const Vector<double>& end,
                      StrahlInfo& S, int numObj, EinschlussInfo *Obj,
		      Vector<double>& Ps, int &Index);
void checkObjectIntersection (double r0, Vector<double>& anf, const Vector<double>& end,
                      const Vector<double> k, int numObj, EinschlussInfo *Obj,
		      Vector<double>& Ps, int &Index);

Vector<double> nextP (Vector<double> P, Vector<double> k, Vector<double> OK,double rK,bool &found); 
void toString (char *S, EinschlussInfo *E, int i);
std::ostream& operator << (std::ostream& os, EinschlussInfo E);
      	
bool operator == (EinschlussInfo a, EinschlussInfo b);
std::ostream &savebinGlobalParms (std::ostream &os, GlobalParms parms);
std::istream & loadbinGlobalParms (std::istream &os, GlobalParms &parms);
GlobalParms readGlobalParms (bool old , std::ifstream &is);
GlobalParms readGlobalParms (bool old , std::ifstream *is);
void writeRRTParms (std::ofstream &os, RRTParmsInfo erg);
 RRTParmsInfo readRRTParms (bool old, std::ifstream *is);
RRTParmsInfo readRRTParms (bool old, std::ifstream &is);
void readRRTParms (std::ifstream &is,RRTParmsInfo &erg);
#endif
