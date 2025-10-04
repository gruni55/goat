#include "lightsrc.h"
#include "lightsrc_mc.h"
#include "ray_pow.h"
#include "misc.h"
#include <complex>

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
			numObjs = 0;
			P0 = 1.0;
			Pall = 0.0;
			i1 = 0;
			i2 = 0;
			numRaysRT = 20;
		}
		LightSrc::LightSrc(const LightSrc& L)
		{
			Pos = L.Pos;
			type = L.type;
			density = L.density;
			e1 = L.e1;
			e2 = L.e2;
			setk(L.k);
			i1 = L.i1;
			i2 = L.i2;
			Pol = L.Pol;
			initPol = Pol;
			r0 = L.r0;
			wvl = L.wvl;
			k0 = L.k0;
			P0 = L.P0;
			numObjs = L.numObjs;
			numRaysRT = L.numRaysRT;
			Obj = L.Obj;
			// copyFormList(Ein,L.Ein,numObjs);
			D = L.D;
			D1 = L.D1;
			D2 = L.D2;

			raytype = L.raytype;
			n0 = 1.0;
			N = L.N;
			Pall = 0;
			suppress_phase_progress = L.suppress_phase_progress;

		}

		void LightSrc::setk(const maths::Vector<double>& k)
		{
			this->k = k / abs(k);
			e1 = k % maths::ez;
			if (abs(e1) < 1E-10) e1 = k % maths::ex;
			e2 = k % e1;
			e1 = e1 / abs(e1);
			e2 = e2 / abs(e2);
			rotVec = maths::cart2sphere(k);
			adjustDirection();
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
			is.read((char*)&D1, (char)sizeof(D1));
			is.read((char*)&D2, (char)sizeof(D2));
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
			os.write((char*)&D1, (char)sizeof(D1));
			os.write((char*)&D2, (char)sizeof(D2));
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
			Obj.push_back(obj);
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

		void LightSrc::ObjectList(int Anz, std::vector<ObjectShape*> Obj)
		{
			for (int i = 0; i < Anz; i++)
				addObject(Obj[i]);
			/*	if (this->numObjs != 0) delete[] this->Ein;
				numObjs=Anz;
				copyFormList(this->Ein,Ein,Anz);*/
		}

		void LightSrc::setPol(maths::Vector<std::complex<double>> pol)
		{
			initPol = pol;
			adjustDirection();
		}

		void LightSrc::setPos(maths::Vector<double> P)
		{

			Pos = P;
		}

		void LightSrc::adjustDirection()
		{
			// Adjust polarisation vector
			maths::Vector<double> polReal = maths::real(initPol);
			maths::Vector<double> polImag = maths::imag(initPol);


			maths::Vector<double> sCoord;
			maths::Vector<double> polRealnew, polImagnew;

			// note: (r,theta,phi) => sCoord
			//       (dr,dtheta,dphi) => rotVec 

			sCoord = maths::cart2sphere(polReal);
			polRealnew[0] = cos(rotVec[2]) * cos(rotVec[1]) * polReal[0] - sin(rotVec[2]) * cos(rotVec[1]) * polReal[1] + cos(sCoord[2] + rotVec[2]) * sin(rotVec[1]) * polReal[2];
			polRealnew[1] = sin(rotVec[2]) * cos(rotVec[1]) * sCoord[0] + cos(rotVec[2]) * cos(rotVec[1]) * polReal[1] + sin(sCoord[2] + rotVec[2]) * sin(rotVec[1]) * polReal[2];
			polRealnew[2] = cos(rotVec[1]) * polReal[2] - sCoord[0] * sin(sCoord[1]) * sin(rotVec[1]);

			sCoord = maths::cart2sphere(polImag);
			polImagnew[0] = cos(rotVec[2]) * cos(rotVec[1]) * polImag[0] - sin(rotVec[2]) * cos(rotVec[1]) * polImag[1] + cos(sCoord[2] + rotVec[2]) * sin(rotVec[1]) * polImag[2];
			polImagnew[1] = sin(rotVec[2]) * cos(rotVec[1]) * sCoord[0] + cos(rotVec[2]) * cos(rotVec[1]) * polImag[1] + sin(sCoord[2] + rotVec[2]) * sin(rotVec[1]) * polImag[2];
			polImagnew[2] = cos(rotVec[1]) * polImag[2] - sCoord[0] * sin(sCoord[1]) * sin(rotVec[1]);

			Pol = maths::Vector<std::complex<double> >(polRealnew[0] + I * polImagnew[0], polRealnew[1] + I * polImagnew[1], polRealnew[2] + I * polImagnew[2]);
		}

		void LightSrc::setObject(ObjectShape* O, int i)
			/**
			   O : Pointer of the object
			   i : Index of the object, which will be changed. The object will be inserted at the end of the object list, if i<0 or larger than the number of objects.

			**/
		{
			if ((i < 0) || (i > numObjs - 1))
			{
				Obj.push_back(O);
				numObjs++;
			}
			else
			{
				Obj[i] = O;
			}
		}

		LightSrcPlane::LightSrcPlane(void)
		{
			type = LIGHTSRC_SRCTYPE_PLANE;
			k = maths::ez;
			density = 1;
			raytype = LIGHTSRC_RAYTYPE_IRAY;
			numObjs = 0;
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
			setk(maths::ez);
			this->raytype = raytype;
			/*
			this->Pol = Pol;
			initPol = Pol;*/
			this->r0 = r0;
			this->wvl = wvl;
			this->D = D;
			this->D1 = D;
			this->D2 = D;
			this->N = N;
			this->n0 = 1.0;
			numObjs = 0;
			setPol(Pol);
			reset();
		}

		LightSrcPlane::LightSrcPlane(const LightSrcPlane& LS) : LightSrc(LS)
		{
			type = LIGHTSRC_SRCTYPE_PLANE;
			e1 = LS.e1;
			e2 = LS.e2;

			this->Pos = LS.Pos;
			this->density = LS.density;
			this->type = LIGHTSRC_SRCTYPE_PLANE;
			this->raytype = LS.raytype;
			this->initPol = LS.initPol;
			this->Pol = LS.Pol;
			this->r0 = LS.r0;
			this->wvl = LS.wvl;
			this->D = LS.D;
			this->D1 = LS.D1;
			this->D2 = LS.D2;
			this->N = LS.N;
			this->n0 = LS.n0;
			// copyFormList(this->Ein, LS.Ein, LS.numObjs);
			this->Obj = LS.Obj;
			this->polType = LS.polType;
			this->numObjs = LS.numObjs;
			suppress_phase_progress = LS.suppress_phase_progress;
			Pall = 0;
			setk(LS.k);
		}

		int LightSrcPlane::next(RayBase* ray)
		{
			ray->suppress_phase_progress = suppress_phase_progress;
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
			maths::Vector<double> P = Pos + (i1 * density - D1 / 2.0) * e1 + (i2 * density - D2 / 2.0) * e2;
			E.e1 = e1;
			E.e2 = e2;
			E.n = k;
			S = IRay(P, Pol * sqrt(P0), k, 1.0, r0, 2.0 * M_PI / wvl, numObjs, Obj);
			S.suppress_phase_progress = suppress_phase_progress;
			S.E1 = Pol;
			S.E2 = Pol2;
			// S.init_Efeld(E,Pol);
			i1++;

			if (i1 * density > D1) { i1 = 0; i2++; }
			if (i2 * density >= D2) { return LIGHTSRC_IS_LAST_RAY; }
			return LIGHTSRC_NOT_LAST_RAY;
		}

		int LightSrcPlane::next(Ray_pow& S)
		{

			Plane E;
			double Pow;

			maths::Vector<double> P = Pos + (i1 * density - D1 / 2.0) * e1 + (i2 * density - D2 / 2.0) * e2;
			// maths::Vector<double> P=Pos+(i1*density-D/2.0)*e1; // NUR ZU TESTZWECKEN !!!!!!!!!!!!!!
			E.e1 = e1;
			E.e2 = e2;
			E.n = k;
			Pow = P0 / ((double)(N * N) * D1 * D2);
			S = Ray_pow(Pow, P, Pol, k, 1.0, r0, 2.0 * M_PI / wvl, numObjs, Obj);
			S.suppress_phase_progress = suppress_phase_progress;
			S.initElectricField(E, Pol);
			S.P = P;
			S.E1 = Pol;
			S.E2 = sqrt(Pow) * Pol / (double)(N * N);
			S.k = k;
			i1++;
			Pall += abs2(S.E2);
			// if (i1*density>D) return LIGHTSRC_IS_LAST_RAY; // NUR ZU TESTZWECKEN !!!!!!!!!!!!!!
			if (i1 * density > D1) { i1 = 0; i2++; }
			if (i2 * density >= D2) { return LIGHTSRC_IS_LAST_RAY; }
			return LIGHTSRC_NOT_LAST_RAY;
		}


		int LightSrcPlane::next(tubedRay& S)
		{
			double Pow = P0 / (N * N * D2 * D1);
			maths::Vector<double> P = Pos + (i1 * density - D1 / 2.0) * e1 + (i2 * density - D2 / 2.0) * e2;
			S = tubedRay(P, density, density, sqrt(Pow) * Pol, k, 1.0, r0, 2.0 * M_PI / wvl, numObjs, Obj);
			S.suppress_phase_progress = suppress_phase_progress;
			S.setN0(n0);
			// S.init_Efeld(Pol,1);
			i1++;
			if (i1 * density > D1) { i1 = 0; i2++; }
			if (i2 * density >= D2) return LIGHTSRC_IS_LAST_RAY;
			return LIGHTSRC_NOT_LAST_RAY;
		}


		LightSrcRing::LightSrcRing() : LightSrc()
		{
			type = LIGHTSRC_SRCTYPE_RING;
			setk(maths::ez);
			density = 2.0 * rmax / ((double)N);
			raytype = LIGHTSRC_RAYTYPE_IRAY;
			numObjs = 0;
			e1 = maths::ex;
			e2 = maths::ey;
			D1 = 2.0 * rmax;
			D2 = 2.0 * rmax;
			D = 2.0 * rmax;
			reset();
		}

		LightSrcRing::LightSrcRing(maths::Vector<double> Pos, int N, double wvl, double rmin, double rmax, maths::Vector<std::complex<double> > Pol, int raytype, double r0) : LightSrc()
		{
			setk(maths::ez);
			e1 = maths::ex;
			e2 = maths::ey;

			this->Pos = Pos;
			this->density = 2.0 * rmax / ((double)N);
			this->type = LIGHTSRC_SRCTYPE_RING;
			D = 2.0 * rmax;
			D1 = 2.0 * rmax;
			D2 = 2.0 * rmax;

			this->raytype = raytype;
			this->Pol = Pol;
			this->initPol = Pol;
			this->r0 = r0;
			this->wvl = wvl;
			this->rmin = rmin;
			this->rmax = rmax;
			this->N = N;
			this->n0 = 1.0;
			numObjs = 0;
			reset();
		}

		void LightSrcRing::setRmin(double rmin)
		{
			if (rmin < rmax) this->rmin = rmin;
		}

		void LightSrcRing::setRmax(double rmax)
		{
			this->rmax = rmax;
			D = rmax / (double)N;
			D1 = 2.0 * rmax;
			D2 = D1;
			density = 2.0 * rmax / ((double)N);
		}

		int LightSrcRing::next(RayBase* ray)
		{
			ray->suppress_phase_progress = suppress_phase_progress;
			switch (raytype)
			{
			case LIGHTSRC_RAYTYPE_IRAY: return next(*(IRay*)ray); break;
			case LIGHTSRC_RAYTYPE_PRAY: return next(*(Ray_pow*)ray); break;
			case LIGHTSRC_RAYTYPE_RAY:
			default: return next(*(tubedRay*)ray);
			}
		}

		int LightSrcRing::next(IRay& S)
		{

			Plane E;
			maths::Vector<double> P;
			double absP;
			bool found = false;
			do
			{
				P = (i1 * density - D1 / 2.0) * e1 + (i2 * density - D2 / 2.0) * e2;
				absP = abs(P);
				if ((absP < rmin) || (absP > rmax))
				{
					i1++;
					if (i1 * density > D1) { i1 = 0; i2++; }
					if (i2 * density >= D2) { return LIGHTSRC_IS_LAST_RAY; }
				}
				else found = true;
			} while (!found);
			P = Pos + P;

			E.e1 = e1;
			E.e2 = e2;
			E.n = k;
			S = IRay(P, Pol * sqrt(P0), k, 1.0, r0, 2.0 * M_PI / wvl, numObjs, Obj);
			S.suppress_phase_progress = suppress_phase_progress;
			S.E1 = Pol / (N * N);
			S.E2 = Pol2 / (N * N);
			// S.init_Efeld(E,Pol);
			i1++;

			if (i1 * density > D1) {
				i1 = 0; i2++; std::cout << "% i2=" << i2 << std::endl;
			}
			if (i2 * density >= D2) { return LIGHTSRC_IS_LAST_RAY; }

			return LIGHTSRC_NOT_LAST_RAY;
		}

		int LightSrcRing::next(Ray_pow& S)
		{

			Plane E;
			double Pow;
			maths::Vector<double> P;
			double absP;
			bool found = false;
			do
			{
				P = (i1 * density - D1 / 2.0) * e1 + (i2 * density - D2 / 2.0) * e2;
				absP = abs(P);
				if ((absP < rmin) || (absP > rmax))
				{
					i1++;

					if (i1 * density > D1) { i1 = 0; i2++; }
					if (i2 * density >= D2) { return LIGHTSRC_IS_LAST_RAY; }
				}
				P = Pos + P;
			} while (!found);
			E.e1 = e1;
			E.e2 = e2;
			E.n = k;
			Pow = P0 / ((double)(N * N) * D1 * D2);
			S = Ray_pow(Pow, P, Pol, k, 1.0, r0, 2.0 * M_PI / wvl, numObjs, Obj);
			S.suppress_phase_progress = suppress_phase_progress;
			S.initElectricField(E, Pol);
			S.P = P;
			S.E1 = Pol;
			S.E2 = sqrt(Pow) * Pol / (double)(N * N);
			S.k = k;
			i1++;
			Pall += abs2(S.E2);
			// if (i1*density>D) return LIGHTSRC_IS_LAST_RAY; // NUR ZU TESTZWECKEN !!!!!!!!!!!!!!
			if (i1 * density > D1) { i1 = 0; i2++; }
			if (i2 * density >= D2) { return LIGHTSRC_IS_LAST_RAY; }
			return LIGHTSRC_NOT_LAST_RAY;
		}


		int LightSrcRing::next(tubedRay& S)
		{
			double Pow = P0 / (N * N * D2 * D1);

			Plane E;
			maths::Vector<double> P;
			double absP;
			bool found = false;
			do
			{
				P = (i1 * density - D1 / 2.0) * e1 + (i2 * density - D2 / 2.0) * e2;
				absP = abs(P);
				if ((absP < rmin) || (absP > rmax))
				{
					i1++;

					if (i1 * density > D1) { i1 = 0; i2++; }
					if (i2 * density >= D2) { return LIGHTSRC_IS_LAST_RAY; }
				}
				else found = true;
				P = Pos + P;
			} while (!found);

			S = tubedRay(P, density, density, sqrt(Pow) * Pol, k, 1.0, r0, 2.0 * M_PI / wvl, numObjs, Obj);
			S.suppress_phase_progress = suppress_phase_progress;
			S.setN0(n0);
			// S.init_Efeld(Pol,1);
			i1++;
			if (i1 * density > D1) { i1 = 0; i2++; }
			if (i2 * density >= D2) return LIGHTSRC_IS_LAST_RAY;
			return LIGHTSRC_NOT_LAST_RAY;
		}





		void LightSrc::reset()
		{
			i1 = 0;
			i2 = 0;
			Pall = 0;
			rayCounter = 0;
			Isum1 = 0;
			Isum2 = 0;
			switch (type)
			{
			case LIGHTSRC_SRCTYPE_PLANE_MC: ((LightSrcPlane_mc*)this)->reset(); break;
			case LIGHTSRC_SRCTYPE_GAUSS_MC: ((LightSrcGauss_mc*)this)->reset(); break;
			}
		}


		LightSrc::~LightSrc(void)
		{
			if (numObjs > 0) {
				Obj.clear();
				Obj.shrink_to_fit();
			}
			numObjs = 0;
			reset();
		}

		void LightSrc::clearObjects()
		{
			if (numObjs > 0)
			{
				Obj.clear();
				numObjs = 0;
			}
		}
	

		void LightSrc::removeObject(ObjectShape* obj)
		{
			for (std::vector<raytracing::ObjectShape*>::iterator it = Obj.begin(); it != Obj.end(); ++it)
				if (*it == obj)
				{
					delete* it;
					Obj.erase(it);
					break;
				}
			numObjs = Obj.size();
		}

		LightSrcGauss::LightSrcGauss(void)
		{
			type = LIGHTSRC_SRCTYPE_GAUSS;			
			Pall = 0.0;
			setk(maths::ez);
			density = 1;
			raytype = LIGHTSRC_RAYTYPE_IRAY;
			numObjs = 0;
			e1 = maths::ex;
			e2 = maths::ey;
			n0 = 1.0;
			Pol = maths::Vector<std::complex<double> >(0, 1.0, 0);
			polType = LIGHTSRC_POL_Y;
			reset();
		}


		LightSrcGauss::LightSrcGauss(maths::Vector<double> Pos, int N, double wvl, double w0, maths::Vector<double> focuspos, double D, maths::Vector<std::complex<double> > Pol, int raytype, double r0) : LightSrc()
			/*
			  Konstruktor f�r Gauss-Strahlungsquelle
			  Strahlen laufen aus einem Rechteck am Ort Pos der Breite D in z-Richtung auf den Fokuspunkt zu
			  Parameter:
			  Pos : Position der Quelle (Mitte)
			  N   : Anzahl Strahlen je Raumrichtung
			  wvl : Wellenl�nge
			  w0  : Taillendurchmesser
			  focuspos : Position des Fokus
			  D   : Breite der Lichtquelle
			  r0  : Radius der "Weltkugel"
			*/
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
			this->raytype = raytype;
			
			this->r0 = r0;
			this->wvl = wvl;
			this->D = D;
			this->D1 = D;
			this->D2 = D;
			this->focuspos = focuspos;		
			this->Pol = Pol;
			this->w0 = w0;
			this->N = N;
			this->P0 = P0;
			this->suppress_phase_progress = suppress_phase_progress;
			setPol(Pol);
			reset();
		}
        void LightSrcGauss::setk(maths::Vector<double> k) 
		{ 
			this->k = k; reset(); 
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

			S.suppress_phase_progress = suppress_phase_progress;

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
				// P : Startort, density=Anzahl Strahlen/L�ngeneinheit, D: Breite der Lichtquelle
				//     i1=momentaner Index in e1-Richtung (normalerweise x-Richtung), i2=momentaner Index in e2-Richtung (normalerweise y-Richtung)

			double x1, x2, x3; // Hilfsgr��en
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
			S.suppress_phase_progress = suppress_phase_progress;
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
			ray->suppress_phase_progress = suppress_phase_progress;
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
			double r = abs(Ph);
			double s = w0 / (2.0 * sqrt(log(2.0)));
			double E0 = exp(-r * r / s / s);
			S = tubedRay(P, density, density, Pol, k, 1.0, r0, 2.0 * M_PI / wvl, numObjs, Obj);
			S.suppress_phase_progress = suppress_phase_progress;

			S.setN0(n0);
			S.n = n0;
			double r2 = abs2(Pos - P);
			std::complex<double> phase = exp(I * (-k0 * n0 * f + zeta - k0 * r2 / (2.0 * R)));

			for (int i = 0; i < 5; i++) S.E[i] = Pol * E0 * phase;
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
			reset();
		}

		std::complex<double> LightSrcGauss::calcStartPhase(maths::Vector<double> P)
		{
			std::complex<double> phase;
			maths::Vector<double> hv = P - Pos;
			double r2 = hv * hv;
			double z = k * (focuspos - P);    // we don't use the absolute value of Pos-focuspos, to get the right sign of z 
			R = z * (1 + z0 * z0 / (z * z));  // calculate the curvature of the phase front at P
			double gouy = atan(z / z0); // Gouy phase
			phase = exp(-I * (k0 * n0 * r2 / (2.0 * R) + (k0 * n0 * z - gouy)));
			return phase;
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
		LightSrcLine::LightSrcLine() : LightSrc()
		{
			type = LIGHTSRC_SRCTYPE_LINE;
		}

		LightSrcLine::LightSrcLine(maths::Vector<double> Pos, int N, double wvl, double size, maths::Vector<double> k, maths::Vector<double> direction) : LightSrc ()
		{
			type = LIGHTSRC_SRCTYPE_LINE;
			this->Pos = Pos;
			this->N = N;
			this->wvl = wvl;
			this->size = size;
			this->density = size / ((double)N);
			Pol = maths::Vector<std::complex<double> >(0, 1, 0);
			this->D = size;
			this->D1 = size;
			this->direction = direction/abs(direction);
			setk(k);
		}

		int LightSrcLine::next(RayBase* ray)
		{
			ray->suppress_phase_progress = suppress_phase_progress;
			switch (raytype)
			{
			case LIGHTSRC_RAYTYPE_IRAY: return next(*(IRay*)ray); break;
			case LIGHTSRC_RAYTYPE_PRAY: return next(*(Ray_pow*)ray); break;
			case LIGHTSRC_RAYTYPE_RAY:
			default: return next(*(tubedRay*)ray);
			}

			return 0;
		}

		int LightSrcLine::next(IRay& S)
		{
			Plane E;
			maths::Vector<double> P = Pos + (i1 * density - D1 / 2.0) * direction;
			E.e1 =direction;
			E.n = k;
			S = IRay(P, Pol * sqrt(P0), k, 1.0, r0, 2.0 * M_PI / wvl, numObjs, Obj);
			S.suppress_phase_progress = suppress_phase_progress;
			S.E1 = Pol / (N * N);
			S.E2 = Pol2 / (N * N);
			// S.init_Efeld(E,Pol);
			i1++;

			if (i1 * density > D1) { return LIGHTSRC_IS_LAST_RAY; }			
			return LIGHTSRC_NOT_LAST_RAY;
		}

		int LightSrcLine::next(Ray_pow& S)
		{

			Plane E;
			double Pow;

			maths::Vector<double> P = Pos + (i1 * density - D1 / 2.0) * direction;
			// maths::Vector<double> P=Pos+(i1*density-D/2.0)*e1; // NUR ZU TESTZWECKEN !!!!!!!!!!!!!!
			E.e1 = direction;
			
			E.n = k;
			Pow = P0 / ((double)(N * N) * D1 * D2);
			S = Ray_pow(Pow, P, Pol, k, 1.0, r0, 2.0 * M_PI / wvl, numObjs, Obj);
			S.suppress_phase_progress = suppress_phase_progress;
			S.initElectricField(E, Pol);
			S.P = P;
			S.E1 = Pol;
			S.E2 = sqrt(Pow) * Pol / (double)(N * N);
			S.k = k;
			i1++;
			Pall += abs2(S.E2);
			// if (i1*density>D) return LIGHTSRC_IS_LAST_RAY; // NUR ZU TESTZWECKEN !!!!!!!!!!!!!!
			if (i1 * density > D1) { return LIGHTSRC_IS_LAST_RAY; }			
			return LIGHTSRC_NOT_LAST_RAY;
		}


		int LightSrcLine::next(tubedRay& S)
		{
			double Pow = P0 / (N * N * D2 * D1);
			maths::Vector<double> P = Pos + (i1 * density - D1 / 2.0) * direction;
			S = tubedRay(P, density, density, sqrt(Pow) * Pol, k, 1.0, r0, 2.0 * M_PI / wvl, numObjs, Obj);
			S.suppress_phase_progress = suppress_phase_progress;
			S.setN0(n0);
			// S.init_Efeld(Pol,1);
			i1++;
			if (i1 * density > D1) { return LIGHTSRC_IS_LAST_RAY; }
			return LIGHTSRC_NOT_LAST_RAY;
		}


}
}
