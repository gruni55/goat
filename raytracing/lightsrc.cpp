#include "lightsrc.h"
#include "ray_pow.h"
#include "misc.h"

namespace  GOAT
{
	namespace raytracing
	{

		double erf(double x)
		{
			// constants    
			double a1 = 0.254829592;
			double a2 = -0.284496736;
			double a3 = 1.421413741;
			double a4 = -1.453152027;
			double a5 = 1.061405429;
			double p = 0.3275911;

			// Save the sign of x    
			int sign = 1;
			if (x < 0)
				sign = -1;
			x = fabs(x);    // A&S formula 7.1.26    
			double t = 1.0 / (1.0 + p * x);
			double y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * exp(-x * x);
			return sign * y;
		}

		LightSrc::LightSrc()
		{
			type = LIGHTSRC_ERROR;
			Obj = 0;
			numObjs = 0;
			P0 = 1.0;
			Pall = 0.0;
			i1 = 0;
			i2 = 0;
		}
		LightSrc::LightSrc(const LightSrc& L)
		{
			Pos = L.Pos;
			type = L.type;
			density = L.density;
			e1 = L.e1;
			e2 = L.e2;
			k = L.k;
			i1 = L.i1;
			i2 = L.i2;
			Pol = L.Pol;
			r0 = L.r0;
			wvl = L.wvl;
			P0 = L.P0;
			numObjs = L.numObjs;
			Obj = L.Obj;
			// copyFormList(Ein,L.Ein,numObjs);
			D = L.D;
			raytype = L.raytype;
			n0 = 1.0;
			N = L.N;
			Pall = 0;
		}

		void LightSrc::setk(const maths::Vector<double>& k)
		{
			this->k = k / abs(k);
			e1 = k % maths::ez;
			if (abs(e1) < 1E-10) e1 = k % maths::ex;
			e2 = k % e1;
			e1 = e1 / abs(e1);
			e2 = e2 / abs(e2);
		}

		void LightSrc::binRead(std::ifstream& is)
		{
			double nre, nim;
			Pos.binRead(is);
			is.read((char*)&type, (char)sizeof(type));
			is.read((char*)&density, (char)sizeof(density));
			k.binRead(is);
			is.read((char*)&N, (char)sizeof(N));
			Pol.binRead(is);
			is.read((char*)&r0, (char)sizeof(r0));
			is.read((char*)&wvl, (char)sizeof(wvl));
			is.read((char*)&numObjs, (char)sizeof(numObjs));
			is.read((char*)&nre, (char)sizeof(nre));
			is.read((char*)&nim, (char)sizeof(nim));
			n0 = std::complex<double>(nre, nim);
			is.read((char*)&D, (char)sizeof(D));
			e1.binRead(is);
			e2.binRead(is);
			is.read((char*)&raytype, (char)sizeof(raytype));
			switch (type)
			{
			case LIGHTSRC_SRCTYPE_PLANE: ((LightSrcPlane*)this)->binReadItem(is); break;
			case LIGHTSRC_SRCTYPE_GAUSS: ((LightSrcGauss*)this)->binReadItem(is); break;
			}
			binReadIncList(is, Obj, numObjs);
		}


		void LightSrc::setR0(double r0)
		{
			this->r0 = r0;
			if (numObjs > 0)
				for (int i = 0; i < numObjs; i++)
					Obj[i]->setr0(r0);
		}

		void binWriteLSList(std::ofstream& os, int nLS, LightSrc** ls)
		{
			for (int i = 0; i < nLS; i++)
			{
				os.write((char*)&ls[i]->type, sizeof(ls[i]->type));
				ls[i]->binWrite(os);
				switch (ls[i]->type)
				{
				case LIGHTSRC_SRCTYPE_PLANE: ((LightSrcPlane*)ls[i])->binWriteItem(os); break;
				case LIGHTSRC_SRCTYPE_GAUSS: ((LightSrcGauss*)ls[i])->binWriteItem(os); break;
				}
			}
		}

		void binReadLSList(std::ifstream& is, int nLS, LightSrc**& ls)
		{
			int type;
			ls = (LightSrc**)malloc(sizeof(LightSrc*) * nLS);
			for (int i = 0; i < nLS; i++)
			{

				is.read((char*)&type, sizeof(type));
				switch (type)
				{
				case LIGHTSRC_SRCTYPE_PLANE:
					ls[i] = new LightSrcPlane();
					ls[i]->binRead(is);
					((LightSrcPlane*)ls[i])->binReadItem(is);
					break;
				case LIGHTSRC_SRCTYPE_GAUSS:
					ls[i] = new LightSrcGauss();
					((LightSrcGauss*)ls[i])->binRead(is);
					((LightSrcGauss*)ls[i])->binReadItem(is);
					break;
				}
			}
		}


		void LightSrc::binWrite(std::ofstream& os)
		{
			double nre, nim;
			Pos.binWrite(os);
			os.write((char*)&type, (char)sizeof(type));
			os.write((char*)&density, (char)sizeof(density));
			k.binWrite(os);
			os.write((char*)&N, (char)sizeof(N));
			Pol.binWrite(os);
			os.write((char*)&r0, (char)sizeof(r0));
			os.write((char*)&wvl, (char)sizeof(wvl));
			os.write((char*)&numObjs, (char)sizeof(numObjs));
			nre = real(n0);
			nim = imag(n0);
			os.write((char*)&nre, (char)sizeof(nre));
			os.write((char*)&nim, (char)sizeof(nim));
			os.write((char*)&D, (char)sizeof(D));
			e1.binWrite(os);
			e2.binWrite(os);
			os.write((char*)&raytype, (char)sizeof(raytype));
			switch (type)
			{
			case LIGHTSRC_SRCTYPE_PLANE: ((LightSrcPlane*)this)->binWriteItem(os); break;
			case LIGHTSRC_SRCTYPE_GAUSS: ((LightSrcGauss*)this)->binWriteItem(os); break;
			}
			binWriteIncList(os, Obj, numObjs);
		}


		void LightSrc::addObject(ObjectShape* obj)
		{
			if (numObjs > 0) Obj = (ObjectShape**)realloc(Obj, sizeof(ObjectShape*) * (numObjs + 1));
			else Obj = (ObjectShape**)malloc(sizeof(ObjectShape*));
			Obj[numObjs] = obj;
			/*switch (obj->type)
			{
			case FUNSURF: Ein[numObjs] = new Funsurf(*((Funsurf*)obj)); break;
			case OBJECTSHAPE_ELLIPSOID: Ein[numObjs] = new Ellipsoid(*((Ellipsoid*)obj)); break;
			case ERYTHROCYTE: Ein[numObjs] = new Erythrocyte(*((Erythrocyte*)obj)); break;
			case SUPERELLIPSOID: Ein[numObjs] = new Superellipsoid(*((Superellipsoid*)obj)); break;
			case SUPERELLIPSOID_N: Ein[numObjs] = new Superellipsoid_n(*((Superellipsoid_n*)obj)); break;
			case SUPERELLIPSOID_D: Ein[numObjs] = new Superellipsoid_D(*((Superellipsoid_D*)obj)); break;
			case OBJECTSHAPE_SURFACE: Ein[numObjs] = new surface(*((surface*)obj)); break;
			case ZYLINDER: Ein[numObjs] = new Zylinder(*((Zylinder*)obj)); break;
			case ZYLINDER_HEXAGONAL: Ein[numObjs] = new ZylinderHexagonal(*((ZylinderHexagonal*)obj)); break;
			case KREISKEGEL: Ein[numObjs] = new Kegel(*((Kegel*)obj)); break;
			case KEGELSTUMPF: Ein[numObjs] = new Kegelstumpf(*((Kegelstumpf*)obj)); break;
			case KEGELSTUMPF_HOHL: Ein[numObjs] = new KegelstumpfHohl(*((KegelstumpfHohl*)obj)); break;
			case COMPOUND: Ein[numObjs] = new Compound(*((Compound*)obj)); break;
			case LINSE: Ein[numObjs] = new Linse(*((Linse*)obj)); break;
			case OBJECTSHAPE_BOX: Ein[numObjs] = new Box(*((Box*)obj)); break;
			}*/
			Obj[numObjs]->r0 = r0;
			Obj[numObjs]->initQuad();
			numObjs++;
		}

		void LightSrc::ObjectList(int Anz, ObjectShape** Obj)
		{
			for (int i = 0; i < Anz; i++)
				addObject(Obj[i]);
			/*	if (this->numObjs != 0) delete[] this->Ein;
				numObjs=Anz;
				copyFormList(this->Ein,Ein,Anz);*/
		}

		void LightSrc::setPos(maths::Vector<double> P)
		{

			Pos = P;
		}

		void LightSrc::setObject(ObjectShape* O, int i)
			/**
			   O : Pointer of the object
			   i : Index of the object, which will be changed. The object will be inserted at the end of the object list, if i<0 or larger than the number of objects.

			**/
		{
			if ((i < 0) || (i > numObjs - 1))
			{
				Obj = (ObjectShape**)realloc(Obj, sizeof(ObjectShape*) * (numObjs + 1));
				Obj[numObjs] = O;
				// copyInc(Ein[numObjs],O);
				numObjs++;
			}
			else
			{
				// delete Ein[i];
				Obj[i] = O;
				// copyInc(Ein[i],O);
			}
		}

		LightSrcPlane::LightSrcPlane(void)
		{
			type = LIGHTSRC_SRCTYPE_PLANE;
			k = maths::ez;
			density = 1;
			raytype = LIGHTSRC_RAYTYPE_IRAY;
			numObjs = 0;
			Obj = 0;
			e1 = maths::ex;
			e2 = maths::ey;

			reset();
		}



		LightSrcPlane::LightSrcPlane(maths::Vector<double> Pos, int N, double wvl, double D, maths::Vector<std::complex<double> > Pol, int raytype, double r0) : LightSrc()
		{
			e1 = maths::ex;
			e2 = maths::ey;

			this->Pos = Pos;
			this->density = D / ((double)N);
			this->type = LIGHTSRC_SRCTYPE_PLANE;
			this->k = maths::ez;
			this->raytype = raytype;
			this->Pol = Pol;
			this->r0 = r0;
			this->wvl = wvl;
			this->D = D;
			this->N = N;
			this->n0 = 1.0;
			numObjs = 0;
			Obj = 0;
			reset();
		}

		LightSrcPlane::LightSrcPlane(const LightSrcPlane& LS)
		{
			type = LIGHTSRC_SRCTYPE_PLANE;
			e1 = LS.e1;
			e2 = LS.e2;

			this->Pos = LS.Pos;
			this->density = LS.density;
			this->type = LIGHTSRC_SRCTYPE_PLANE;
			this->k = LS.k;
			this->raytype = LS.raytype;
			this->Pol = LS.Pol;
			this->r0 = LS.r0;
			this->wvl = LS.wvl;
			this->D = LS.D;
			this->N = LS.N;
			this->n0 = LS.n0;
			// copyFormList(this->Ein, LS.Ein, LS.numObjs);
			this->Obj = LS.Obj;
			this->polType = LS.polType;
			this->numObjs = LS.numObjs;
			Pall = 0;
		}

		int LightSrcPlane::next(RayBase* ray)
		{
			switch (raytype)
			{
			case LIGHTSRC_RAYTYPE_IRAY: return next(*(IRay*)ray); break;
			case LIGHTSRC_RAYTYPE_PRAY: return next(*(Ray_pow*)ray); break;
			case LIGHTSRC_RAYTYPE_RAY:
			default: return next(*(tubedRay*)ray);
			}
		}

		int LightSrcPlane::next(IRay& S)
		{

			Plane E;
			maths::Vector<double> P = Pos + (i1 * density - D / 2.0) * e1 + (i2 * density - D / 2.0) * e2;
			E.e1 = e1;
			E.e2 = e2;
			E.n = k;
			S = IRay(P, Pol, k, 1.0, r0, 2.0 * M_PI / wvl, numObjs, Obj);
			S.E1 = Pol / (N * N);
			S.E2 = Pol2 / (N * N);
			// S.init_Efeld(E,Pol);
			i1++;

			if (i1 * density > D) { i1 = 0; i2++; }
			if (i2 * density >= D) { return LIGHTSRC_IS_LAST_RAY; }
			return LIGHTSRC_NOT_LAST_RAY;
		}

		int LightSrcPlane::next(Ray_pow& S)
		{

			Plane E;
			double Pow;

			maths::Vector<double> P = Pos + (i1 * density - D / 2.0) * e1 + (i2 * density - D / 2.0) * e2;
			// maths::Vector<double> P=Pos+(i1*density-D/2.0)*e1; // NUR ZU TESTZWECKEN !!!!!!!!!!!!!!
			E.e1 = e1;
			E.e2 = e2;
			E.n = k;
			Pow = 1.0 / ((double)(N * N) * D * D);
			S = Ray_pow(Pow, P, Pol, k, 1.0, r0, 2.0 * M_PI / wvl, numObjs, Obj);

			S.initElectricField(E, Pol);
			S.P = P;
			S.E1 = Pol;
			S.E2 = sqrt(Pow) * Pol / (double)(N * N);
			S.k = k;
			i1++;
			Pall += abs2(S.E2);
			// if (i1*density>D) return LIGHTSRC_IS_LAST_RAY; // NUR ZU TESTZWECKEN !!!!!!!!!!!!!!
			if (i1 * density > D) { i1 = 0; i2++; }
			if (i2 * density >= D) { return LIGHTSRC_IS_LAST_RAY; }
			return LIGHTSRC_NOT_LAST_RAY;
		}


		int LightSrcPlane::next(tubedRay& S)
		{
			double Pow = 1.0 / (N * N * D * D);
			maths::Vector<double> P = Pos + (i1 * density - D / 2.0) * e1 + (i2 * density - D / 2.0) * e2;
			S = tubedRay(P, density, density, sqrt(Pow) * Pol, k, 1.0, r0, 2.0 * M_PI / wvl, numObjs, Obj);
			S.setN0(n0);
			// S.init_Efeld(Pol,1);
			i1++;
			if (i1 * density > D) { i1 = 0; i2++; }
			if (i2 * density >= D) return LIGHTSRC_IS_LAST_RAY;
			return LIGHTSRC_NOT_LAST_RAY;
		}

		void LightSrc::reset()
		{
			i1 = 0;
			i2 = 0;
			Pall = 0;
		}


		LightSrc::~LightSrc(void)
		{
			if (numObjs > 0) delete[] Obj;
			numObjs = 0;
			reset();
		}

		void LightSrc::clearObjects()
		{
			if (numObjs > 0)
			{
				/*	 for (int i=0; i<numObjs; i++)
						 delete Ein[i];
					 delete[] Ein;*/
				numObjs = 0;
			}
		}

		LightSrcGauss::LightSrcGauss(void)
		{
			type = LIGHTSRC_SRCTYPE_GAUSS;
			Pall = 0.0;
			k = maths::ez;
			density = 1;
			raytype = LIGHTSRC_RAYTYPE_IRAY;
			numObjs = 0;
			Obj = 0;
			e1 = maths::ex;
			e2 = maths::ey;
			n0 = 1.0;
			Pol = maths::Vector<std::complex<double> >(0, 1.0, 0);
			polType = LIGHTSRC_POL_Y;
			reset();
		}


		LightSrcGauss::LightSrcGauss(maths::Vector<double> Pos, int N, double wvl, double w0, maths::Vector<double> focuspos, double D, maths::Vector<std::complex<double> > Pol, int raytype, double r0) : LightSrc()
			/**
			  Konstruktor für Gauss-Strahlungsquelle
			  Strahlen laufen aus einem Rechteck am Ort Pos der Breite D in z-Richtung auf den Fokuspunkt zu
			  Parameter:
			  Pos : Position der Quelle (Mitte)
			  N   : Anzahl Strahlen je Raumrichtung
			  wvl : Wellenlänge
			  w0  : Taillendurchmesser
			  focuspos : Position des Fokus
			  D   : Breite der Lichtquelle
			  r0  : Radius der "Weltkugel"
			**/
		{
			/*
			e1=ex;
			e2=ey;
			*/
			// Erst mal die Vektoren e1 und e2 setzen

			Pall = 0.0;
			this->Pos = Pos;
			this->density = D / ((double)N);
			this->type = LIGHTSRC_SRCTYPE_GAUSS;
			k = focuspos - Pos;
			k = k / abs(k);
			e1 = k % maths::ez;
			if (abs(e1) < 1E-10) e1 = maths::ex;
			e2 = k % e1;
			e1 = e1 / abs(e1);
			e2 = e2 / abs(e2);

			this->raytype = raytype;
			this->Pol = Pol;
			this->r0 = r0;
			this->wvl = wvl;
			this->D = D;
			this->focuspos = focuspos;
			Pol = maths::Vector<std::complex<double> >(0, 1.0, 0);
			polType = LIGHTSRC_POL_Y;

			this->Pol = Pol;
			this->w0 = w0;
			this->N = N;
			this->P0 = P0;
			f = abs(Pos - focuspos);
			z0 = M_PI * w0 * w0 / wvl;
			numObjs = 0;
			Obj = 0;
			n0 = 1.0;
			double d = abs(Pos - focuspos);
			calcz0();
			calcw(d);
			double theta = atan(wvl / (M_PI * w0));
			NA = real(n0) * sin(theta);
			reset();
		}


		void LightSrcGauss::setWvl(double wvl)
		{
			this->wvl = wvl;
			double theta = asin(real(NA / this->n0));
			w0 = wvl / (M_PI * tan(theta));
		}

		int LightSrcGauss::next(IRay& S)
		{
			maths::Vector<double> fp, P = Pos + (i1 * density - D / 2.0) * e1 + (i2 * density - D / 2.0) * e2;
			// maths::Vector<double> P=Pos+(i1*density-D/2.0)*e1; // NUR ZU TESTZWECKEN !!!!!!!!!! 
			double x1, x2, x3;
			x1 = P * e1;
			x2 = P * e2;
			double r2 = x1 * x1 + x2 * x2;

			double s2 = w0 * w0 / log(2.0);
			double g;
			double L = 0.1;
			double R;

			maths::Vector<double> h, hk;
			maths::Matrix<double> DM;
			maths::Vector<std::complex<double> > E;
			std::complex<double> E0;
			double absh, gamma;
			double zeta;
			double absfp;

			fp = focuspos - P;  // Hier beginnt der Strahl
			hk = fp / abs(fp);  // Normierter Richtungsvektor vom Startpunkt zum Fokus


			h = k % hk;
			absh = abs(h);

			x3 = fp * k;
			R = x3 * (1.0 + z0 * z0 / (x3 * x3));
			if (z0 == 0) zeta = M_PI / 2.0;
			else zeta = atan(x3 / z0);


			if (absh == 0)
			{
				hk = k;
				gamma = 0;
			}
			else
			{
				h /= absh;
				gamma = std::acos(k * hk / (abs(k) * abs(hk)));
			}

			// Pol=maths::Vector<std::complex<double> >(0,1,0);  // y-Polarisation

			DM = rotMatrix(h, gamma);
			S = IRay(P, Pol, hk, 1.0, r0, 2.0 * M_PI / wvl, numObjs, Obj);
			S.k = k;
			S.n = n0;
			absfp = fp * k;
			g = gaussw(-fp * k, 2.0 * M_PI / std::real(S.k0), w0); // w(z)   
			r2 = S.P[0] * S.P[0] + S.P[1] * S.P[1];

			//E0=sqrt(1.0/M_PI*g/(absfp*w0*w0)*sqrt(2.0/M_PI)*exp(-2.0*absfp*absfp/(g*g))
			//		*1.0/(1.0-erf(sqrt(2.0)*absfp/g)))*w0/g*exp(-r2/g/g);

			E0 = sqrt(P0) * sqrt(2.0 / M_PI / g / g) * exp(-r2 / g / g); // ohne Korrektur
			// E0=1;
			maths::Vector<double> F = focuspos;

			S.k = focuspos - S.P;   // Richtungsvektor auf den Fokus gerichtet
			S.k /= abs(S.k);
			h = k % S.k;
			absh = abs(h);
			if (absh != 0)
			{
				h /= absh;
				gamma = std::acos(S.k * k);
				DM = rotMatrix(h, gamma);
				S.E1 = E0 * DM * Pol;
				S.E2 = E0 * DM * Pol2;
			}
			else
			{
				S.E1 = E0 * Pol;
				S.E2 = E0 * Pol2;
			}
			S.n = n0;
			i1++;
			if (i1 * density >= D) { i1 = 0; i2++; }
			if (i2 * density >= D) return LIGHTSRC_IS_LAST_RAY;
			return LIGHTSRC_NOT_LAST_RAY;

		}

		int LightSrcGauss::next(Ray_pow& S)
		{
			maths::Vector<double> fp, P = Pos + (i1 * density - D / 2.0) * e1 + (i2 * density - D / 2.0) * e2;
			//  maths::Vector<double> fp,P=Pos+(i1*density-D/2.0)*e1; // nur zu TESTZWECKEN !!!
				// P : Startort, density=Anzahl Strahlen/Längeneinheit, D: Breite der Lichtquelle
				//     i1=momentaner Index in e1-Richtung (normalerweise x-Richtung), i2=momentaner Index in e2-Richtung (normalerweise y-Richtung)

			double x1, x2, x3; // Hilfsgrößen
			x1 = P * e1;
			x2 = P * e2;
			double r2 = x1 * x1 + x2 * x2;  // Quadrat des Abstands von der Laserstrahlachse

			double s2 = w0 * w0 / log(2.0);  // w0 : Strahltaille 
			double g;
			double L = 0.1;
			double R;

			maths::Vector<double> h, hk;
			maths::Matrix<double> DM;
			maths::Vector<std::complex<double> > E;
			std::complex<double> E0;
			double absh, gamma;
			double zeta;
			double absfp;
			calcNormfak();
			fp = focuspos - P;  // Hier beginnt der Strahl
			hk = fp / abs(fp);  // Normierter Richtungsvektor vom Startpunkt zum Fokus

			h = k % hk;
			absh = abs(h);

			x3 = fp * k;
			R = x3 * (1.0 + z0 * z0 / (x3 * x3));
			if (z0 == 0) zeta = M_PI / 2.0;
			else zeta = atan(x3 / z0);


			if (absh == 0)
			{
				hk = k;
				gamma = 0;
			}
			else
			{
				h /= absh;
				gamma = std::acos(k * hk / (abs(k) * abs(hk)));
			}

			// Pol=maths::Vector<std::complex<double> >(0,1,0);  // y-Polarisation

			DM = rotMatrix(h, gamma);
			S = Ray_pow(1, P, Pol, hk, n0, r0, 2.0 * M_PI / wvl, numObjs, Obj);
			S.k = k;
			S.n = n0;
			absfp = fp * k;
			g = gaussw(-fp * k, 2.0 * M_PI / std::real(S.k0), w0); // w(z)  


			r2 = S.P[0] * S.P[0] + S.P[1] * S.P[1];

			S.Pow = P0 * k * hk * exp(-2.0 * r2 / g / g) / ((double)N * N) * D * D / g;//*sqrt(8.0/M_PI);// *real(Normfak);

			// cout << "%w=" << g << "   Pow=" << S.Pow << "   P0=" << P0 << endl;
#ifdef MIT_NORMIERUNG
			S.Pow = k * hk * exp(-2.0 * r2 / g / g) / ((double)N * N) * D * D / g;
#endif
			S.E1 = Pol;
			S.E2 = sqrt(S.Pow) * Pol;

			maths::Vector<double> F = focuspos;
			S.k = focuspos - S.P;   // Richtungsvektor auf den Fokus gerichtet
			S.k /= abs(S.k);
			h = k % S.k;
			absh = abs(h);
			if (absh != 0)
			{
				h /= absh;
				gamma = std::acos(S.k * k);
				DM = rotMatrix(h, gamma);
				S.E1 = DM * S.E1;
				S.E2 = DM * S.E2;
			}
			// else S.E1=E0*Pol;
			S.n = n0;
			i1++;
			if (type == LIGHTSRC_SRCTYPE_TOPHAT)
			{
				S.E2 = S.E2 / abs(S.E2);
				S.Pow = 1.0;
			}
			Pall += abs2(S.E2);
			//  if (i1*density>D) { return LIGHTSRC_IS_LAST_RAY; } // nur zu TESTZWECKEN !!!!
			if (i1 * density > D) { i1 = 0; i2++; }
			if (i2 * density > D) return LIGHTSRC_IS_LAST_RAY;
			return LIGHTSRC_NOT_LAST_RAY;
		}

		int LightSrcGauss::next(RayBase* ray)
		{
			switch (raytype)
			{
			case LIGHTSRC_RAYTYPE_IRAY: return next(*(IRay*)ray); break;
			case LIGHTSRC_RAYTYPE_PRAY: return next(*(Ray_pow*)ray); break;
			case LIGHTSRC_RAYTYPE_RAY:
			default: return next(*(tubedRay*)ray);
			}
		}

		int LightSrcGauss::next(tubedRay& S)
		{
			maths::Vector<double> Ph = (i1 * density - D / 2.0) * e1 + (i2 * density - D / 2.0) * e2;
			maths::Vector<double> P = Pos + Ph;
			maths::Vector<double> k = focuspos - P;  // Richtung des Strahles
			k = k / abs(k);
			double r = abs(Ph);
			double s = w0 / (2.0 * sqrt(log(2.0)));
			double E0 = exp(-r * r / s / s);
			S = tubedRay(P, density, density, Pol, k, 1.0, r0, 2.0 * M_PI / wvl, numObjs, Obj);

			S.setN0(n0);
			S.n = n0;
			for (int i = 0; i < 5; i++) S.E[i] = Pol * E0;
			i1++;
			if (i1 * density > D) { i1 = 0; i2++; }
			if (i2 * density >= D) return LIGHTSRC_IS_LAST_RAY;
			return LIGHTSRC_NOT_LAST_RAY;
		}

		/*int LightSrcGauss::next(ray &S)
		{
			maths::Vector<double> fp,P=Pos+(i1*density-D/2.0)*e1+(i2*density-D/2.0)*e2;
			 // maths::Vector<double> P=Pos+(i1*density-D/2.0)*e1; // NUR ZU TESTZWECKEN !!!!!!!!!!
		 double x1,x2,x3;
		 x1=P*e1;
		 x2=P*e2;
		 double r2=x1*x1+x2*x2;

		 double s2=w0*w0/log(2.0);
		 double g;
		 double L=0.1;
		 double R;

		 maths::Vector<double> h,hk;
		 maths::Matrix<double> DM;
		 maths::Vector<std::complex<double> > E;
		 std::complex<double> E0;
		 double absh,gamma;
		  double zeta;
		  double Eph;
		  double l;
		  double absfp;

		  fp=focuspos-P;  // Hier beginnt der Strahl
		  hk=fp/abs(fp);  // Normierter Richtungsvektor vom Startpunkt zum Fokus


		  h=k%hk;
		  absh=abs(h);

		  x3=fp*k;
		  R=x3*(1.0+z0*z0/(x3*x3));
		  if (z0==0) zeta=M_PI/2.0;
		  else zeta=atan(x3/z0);


		  if (absh==0)
		  {
			  hk=k;
			  gamma=0;
		  }
		  else
		  {
		  h/=absh;
		  gamma=acos(k*hk/(abs(k)*abs(hk)));
		  }

		  Pol=maths::Vector<std::complex<double> >(0,1,0);  // y-Polarisation

		  DM=rotMatrix(h,gamma);
		  S=ray(P,density,density,Pol,hk,1.0,r0,2.0*M_PI/wvl,numObjs,Ein);
		  S.setN0(n0);
		  S.n=n0;
		  absfp=fp*k;
		  g=w0/(2.0*sqrt(log(2.0)));

		   g=gaussw(-fp*k,2.0*M_PI/S.k0,w0); // w(z)
		r2=S.P[4][0]*S.P[4][0]+S.P[4][1]*S.P[4][1];

		E0=sqrt(1.0/M_PI*g/(absfp*w0*w0)*sqrt(2.0/M_PI)*exp(-2.0*absfp*absfp/(g*g))
				*1.0/(1.0-erf(sqrt(2.0)*absfp/g)))*w0/g*exp(-r2/g/g);

		//  E0=sqrt(2.0/M_PI/g/g)*exp(-r2/g/g); // ohne Korrektur
		// E0=1;
		maths::Vector<double> F=focuspos;
		double B,C2;

		for (int i=0; i<5; i++)
		{
			 S.k[i]=focuspos-S.P[i];   // Richtungsvektor auf den Fokus gerichtet
			 S.k[i]/=abs(S.k[i]);
			 h=k%S.k[i];
			 absh=abs(h);
			 if (absh!=0)
			 {
				 h/=absh;
				 gamma=acos(S.k[i]*k);
				 DM=rotMatrix(h,gamma);
				 S.E[i]=E0*DM*Pol;
			 }
			 else S.E[i]=E0*Pol;
			 S.n0=n0;
		}
		  i1++;
		//  cout << "E=" << S.E[4] << endl;
		 // if (i1*density>D) return LIGHTSRC_IS_LAST_RAY; // NUR ZU TESTZWECKEN !!!!!!!!!!!!!!!!!!
		  if (i1*density>D) {i1=0; i2++; }
		  if (i2*density>=D) return LIGHTSRC_IS_LAST_RAY;
		  return LIGHTSRC_NOT_LAST_RAY;

		}*/

		LightSrcGauss::LightSrcGauss(const LightSrcGauss& L) : LightSrc(L)
		{
			focuspos = L.focuspos;
			n0 = L.n0;
			wvl = L.wvl;
			w0 = L.w0;
			P0 = L.P0;
			k = L.k;
			numObjs = L.numObjs;
		}

		/*
		void LightSrcGauss::initElectricFieldGauss(maths::Vector<double> &P, maths::Vector<double> &n, maths::Vector<std::complex<double> > &E)
		{

		  Pol=E;
		  k=n;
		  focuspos=P;
		  calcz0();
		   double d=abs(Pos-focuspos);
		  calcw(d);
		}  */


		/*LightSrcGauss::LightSrcGauss(const LightSrc &L) : LightSrc(L)
		{

		};*/

		std::ostream& operator << (std::ostream& os, LightSrc* ls)
		{
			os << "% Anzahl Objekte:" << ls->numObjs << std::endl;
			for (int i = 0; i < ls->numObjs; i++)
			{
				os << "%------------------- Objekt " << i << " -----------------" << std::endl;
				os << "% P=" << ls->Obj[i]->P << std::endl;
				if (ls->Obj[i]->type == OBJECTSHAPE_ELLIPSOID)
				{
					Ellipsoid* e = (Ellipsoid*)ls->Obj[i];
					os << "% Halbachsen= " << e->r << std::endl;
				}
				os << "% Brechungsindex n=" << real(ls->Obj[i]->n) << "+ i*" << imag(ls->Obj[i]->n) << std::endl;
				os << "% -----------------  ENDE ----------------" << std::endl;

			}
			os << "% Pos=" << ls->getPos() << std::endl;
			os << "% D=" << ls->D << "   k=" << ls->getk() << "  N=" << ls->N << std::endl;
			os << "% density=" << ls->density << "  n0=" << ls->n0 << "   Pow=" << ls->P0 << std::endl;
			switch (ls->type)
			{
			case LIGHTSRC_SRCTYPE_GAUSS: os << "Fokuspos=" << ((LightSrcGauss*)ls)->getFocuspos() << std::endl; break;
			};
			return os;
		}


		void LightSrcGauss::binWriteItem(std::ofstream& os)
		{
			os.write((char*)&w0, (char)sizeof(w0));
			focuspos.binWrite(os);
			os.write((char*)&z0, (char)sizeof(z0));
			os.write((char*)&P0, (char)sizeof(P0));
			os.write((char*)&f, (char)sizeof(f));
		}

		void LightSrcGauss::binReadItem(std::ifstream& is)
		{
			is.read((char*)&w0, (char)sizeof(w0));
			focuspos.binRead(is);
			is.read((char*)&z0, (char)sizeof(z0));
			is.read((char*)&P0, (char)sizeof(P0));
			is.read((char*)&f, (char)sizeof(f));
			k = focuspos - Pos;
			k = k / abs(k);
		}

		void LightSrcGauss::setNA(double NA)
		{
			this->NA = NA;
			double theta = asin(NA / real(n0));
			double l = abs(Pos - focuspos);
			w0 = wvl / (M_PI * tan(theta));
			calcz0();
			w=calcw(l);
			// D = 2.0 * l * wvl / (M_PI * w0);
			D = 6 * w;
			density = D / ((double)N);
			
		}

		void copyLightSrcList(LightSrc**& d, LightSrc** s, int nLS)
		{
			if ((nLS > 0) && (s != 0))
			{
				d = (LightSrc**)malloc(sizeof(LightSrc*) * nLS);
				for (int i = 0; i < nLS; i++)
				{
					switch (s[i]->type)
					{
					case LIGHTSRC_SRCTYPE_PLANE: d[i] = new LightSrcPlane(*((LightSrcPlane*)s[i])); break;
					case LIGHTSRC_SRCTYPE_GAUSS: d[i] = new LightSrcGauss(*((LightSrcGauss*)s[i])); break;
					case LIGHTSRC_SRCTYPE_TOPHAT: d[i] = new LightSrcGauss(*((LightSrcGauss*)s[i]));
						d[i]->type = LIGHTSRC_SRCTYPE_TOPHAT;
						break;
					}
				}
			}
		}
	}
}