#include <complex>
#include <limits>
#include "vector.h"
#include "matrix.h"

#pragma once
namespace GOAT
{
	namespace raytracing
	{
#ifndef INF
		constexpr double INF = std::numeric_limits<double>::infinity();
#endif
#ifndef complex_I
#define complex_I
		//std::complex<double> I=std::complex<double>(0.0,1.0);    // imaginaere Einheit
#define I std::complex<double>(0.0,1.0)
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

#define SQRT3  1.73205080756887729352744634150587 // Wurzel aus 3
#define RAMAN 0
#define FLUORESZENZ 1

#define IN_OBJECT  1 
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
			int x, y;
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
			maths::Vector<double> k, EAmp;
			std::complex<double> E;
			std::complex<double> phi;
			int m;
			double b;
		} StrahlInfo;

		typedef struct
		{

			maths::Vector <double> P;
			StrahlInfo      S;
		} StrahlArray;

		typedef struct
		{
			/*
			   P  : Ort des Einschlusses
			   n  : Brechungsindex
			   a  : Radius (in Einheiten des Partikelradius�)
			   alpha : Polarisierbarkeit;
			*/
			maths::Vector <double> P;
			std::complex<double> n;
			double a;
			maths::Matrix<double> alpha;
		} objectInfo;


		/**
		* @brief Class used to set the parameters for inelastic scattering (may replaced later)
		*/
		class RRTParmsInfo
		{
		public:
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

		public:
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

		void output(int nx, int ny, maths::Vector<std::complex<double> >** G);
		void init_Strahl(GlobalParms Parms, StrahlArray* Strahl);
		void sub_Gitter(GlobalParms parms, maths::Vector<std::complex<double> >** Erg,
			maths::Vector<std::complex<double> >** Gitter1,
			maths::Vector<std::complex<double> >** Gitter2);
		void add_Gitter(GlobalParms parms, maths::Vector<std::complex<double> >** Erg,
			maths::Vector<std::complex<double> >** Gitter1,
			maths::Vector<std::complex<double> >** Gitter2);
		void clear(GlobalParms parms, maths::Vector<std::complex<double> >** Gitter);
		void Delete(int n, maths::Vector<std::complex<double> >** Gitter);
		void copy(StrahlInfo& dest, StrahlInfo src);
		double minmax(double a, double b);
		maths::Vector<double> kart2sph(maths::Vector<double> v);
		maths::Vector<double> sph2kart(maths::Vector<double> v);

		Point get_grid_point(const GlobalParms& parms, maths::Vector<double> P);
		maths::Vector <double> set_grid_point(const GlobalParms& parms, Point P);
		double grad(const GlobalParms& parms, maths::Vector<std::complex<double> >** Gitter, const maths::Vector<double>& P);
		//std::complex<double>  asin(const std::complex<double> &);
		std::complex<double>  acos(const std::complex<double>&);
		maths::Vector<double> next(const GlobalParms& parms,
			const maths::Vector<double>& P0,
			const maths::Vector<double>& k);

		void minmax(double x, double dx, int& min, int& max);
		void checkObjectIntersection(maths::Vector<double>& anf, const maths::Vector<double>& end,
			StrahlInfo& S, int numObj, objectInfo* Obj,
			maths::Vector<double>& Ps, int& Index);
		void checkObjectIntersection(double r0, maths::Vector<double>& anf, const maths::Vector<double>& end,
			StrahlInfo& S, int numObj, objectInfo* Obj,
			maths::Vector<double>& Ps, int& Index);
		void checkObjectIntersection(double r0, maths::Vector<double>& anf, const maths::Vector<double>& end,
			const maths::Vector<double> k, int numObj, objectInfo* Obj,
			maths::Vector<double>& Ps, int& Index);

		maths::Vector<double> nextP(maths::Vector<double> P, maths::Vector<double> k, maths::Vector<double> OK, double rK, bool& found);
		void toString(char* S, objectInfo* E, int i);
		std::ostream& operator << (std::ostream& os, objectInfo E);

		bool operator == (objectInfo a, objectInfo b);
		std::ostream& savebinGlobalParms(std::ostream& os, GlobalParms parms);
		std::istream& loadbinGlobalParms(std::istream& os, GlobalParms& parms);
		GlobalParms readGlobalParms(bool old, std::ifstream& is);
		GlobalParms readGlobalParms(bool old, std::ifstream* is);
		void writeRRTParms(std::ofstream& os, RRTParmsInfo erg);
		RRTParmsInfo readRRTParms(bool old, std::ifstream* is);
		RRTParmsInfo readRRTParms(bool old, std::ifstream& is);
		void readRRTParms(std::ifstream& is, RRTParmsInfo& erg);
	}
}
#endif
