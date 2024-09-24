#include "raytrace_inel.h"
#include "inel_calc.h"

namespace GOAT
{
	namespace raytracing
	{
		RRTParms calcDetDirParms(double theta, double phi, double wvlinel)
		{
			RRTParms D;
			
			D.n = -maths::Vector<double>(cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta));

			if (fabs(D.n * maths::ez)-1<1E-10)
			{
				D.e1 = maths::ex;
				D.e2 = maths::ey;
			}
			else
			{
				D.e1 = maths::ez * (1.0 - maths::ez * D.n);
				D.e2 = D.e1 % D.n;
			}
			D.wvlinel = wvlinel;
			return D;
		}

		Raytrace_Inel::Raytrace_Inel()
		{
			active = 0;
			SGE = 0;
			SGRRT1 = 0;
			SGRRT2 = 0;
			nLSExcit = 0;
			LSExcit = 0;
			n = 0;
			iR = 0;
			calcphase = INEL_CALCPHASE_EXCITATION;
		}

		Raytrace_Inel::Raytrace_Inel(const Scene& S, int n = 500) : Raytrace((Scene)S)
		{
			SGRRT1 = 0;
			SGRRT2 = 0;
			nLSExcit = 0;
			LSExcit = 0;
			this->n = n;
			SGE = 0;
			iR = 0;
			calcphase = INEL_CALCPHASE_EXCITATION;
			active = new bool[S.nObj];
			for (int i = 0; i < S.nObj; i++) active[i] = true;
			sceneChanged(S);
		}

		std::complex<double>  Raytrace_Inel::gewichte(maths::Vector<std::complex<double> > E, maths::Vector<std::complex<double> > p)
		{
			/*
			Berechnet aus dem Dipolmoment p und den Eigenschaften des Strahls (Efeld,Richtung)
			die Gewichtung fuer das RRT
			obwohl p vom Typ "maths::Vector<complex<double> >" ist, wird angenommen, dass p ein reeller Vektor
			ist.
			S.k wird wie gewoehnlich als normiert betrachtet
			*/

			std::complex<double>  Erg;
			if (parms.radiationType == INEL_RADIATION_TYPE_FLOURESCENCE) return abs(p);
			if ((abs(E) == 0.0) || (abs(p) == 0.0)) return 0.0;
			Erg = p * E / abs(E);
			return Erg;
		}

		void Raytrace_Inel::trace(RRTParms D)
		{
			// Das anregende Feld muss nur einmal berechnet werden
			if ((calcphase == INEL_CALCPHASE_EXCITATION) || (calcphase == INEL_CALCPHASE_EXCITATION_ONLY))
			{
				std::cout << "starting calculation excitation field" << std::endl;
				std::cout << "start initialization...";
				initExcitation();
				std::cout << "done." << std::endl;
				S.setRaytype(LIGHTSRC_RAYTYPE_RAY);
				std::cout << "start raytracing...";
				Raytrace::trace();
				std::cout << "done." << std::endl;
				S.resetLS();
				std::cout << "adding fields...";
				for (int i = 1; i < INEL_MAX_NREFLEX; i++)
					SGE[0].add(SGE[i]);
				std::cout << "done." << std::endl;
				if (calcphase == INEL_CALCPHASE_EXCITATION) calcphase = INEL_CALCPHASE_RRT;
			}
			
			// Jetzt folgt die Berechnung des R�ckverfolgten (inelastischen) Felds
			if (calcphase != INEL_CALCPHASE_EXCITATION_ONLY)
			{
				std::cout << "starting inelastic calculation..." << std::endl;
				LightSrcPlane LSRRT(-D.n * S.r0, n*2.0, D.wvlinel, 2.0 * S.r0);
				LSRRT.setk(D.n);
				S.addLightSourceRRT(&LSRRT, maths::Vector<std::complex<double> >(D.e1[0], D.e1[1], D.e1[2]), maths::Vector<std::complex<double> >(D.e2[0], D.e2[1], D.e2[2]));
				initRRT();
				traceRRT();
				
				// Die unterschiedlichen Reflexionsordnungen m�ssen zusammengef�hrt werden
				for (int i = 1; i < INEL_MAX_NREFLEX; i++)
				{
					SGRRT1[0].add(SGRRT1[i]);
					SGRRT2[0].add(SGRRT2[i]);
				}
				
				// Jetzt muss das Ergebnis noch mit dem anregenden Feld und der Dipolcharakteristik gewichtet werden
				std::complex<double> g;
				for (int i = 0; i < SGRRT1[0].numObjs; i++)
					for (int ix = 0; ix < SGRRT1[0].n[i][0]; ix++)
						for (int iy = 0; iy < SGRRT1[0].n[i][1]; iy++)
							for (int iz = 0; iz < SGRRT1[0].n[i][2]; iz++)
							{
								g = gewichte(SGRRT1[0].G[i][ix][iy][iz], S.Obj[i]->alpha * SGE[0].G[i][ix][iy][iz]);
								SGRRT1[0].G[i][ix][iy][iz] *= g;
								g = gewichte(SGRRT2[0].G[i][ix][iy][iz], S.Obj[i]->alpha * SGE[0].G[i][ix][iy][iz]);
								SGRRT2[0].G[i][ix][iy][iz] *= g;
							}
				double  anzrays2 = (double)S.LSRRT->getNumRays() * (double)S.LSRRT->getNumRays();
				
				if (parms.coherency==INEL_RADIATION_COHERENT)
				{
					inel1 = abs2sum(SGRRT1[0]) / anzrays2 * SGRRT1[0].nges[0] * SGRRT1[0].nges[1] * SGRRT1[0].nges[2];
					inel2 = abs2sum(SGRRT2[0]) / anzrays2 * SGRRT2[0].nges[0] * SGRRT2[0].nges[1] * SGRRT2[0].nges[2];
				}
				else
				{
					inel1 = sumabs2(SGRRT1[0]) / anzrays2 * SGRRT1[0].nges[0] * SGRRT1[0].nges[1] * SGRRT1[0].nges[2];
					inel2 = sumabs2(SGRRT2[0]) / anzrays2 * SGRRT2[0].nges[0] * SGRRT2[0].nges[1] * SGRRT2[0].nges[2];
				}
			}
			else calcphase = INEL_CALCPHASE_RRT;
		}

		void Raytrace_Inel::traceEnterObject()
		{
		}

		void Raytrace_Inel::traceLeaveObject()
		{
		if ((calcphase == INEL_CALCPHASE_EXCITATION) || (calcphase == INEL_CALCPHASE_EXCITATION_ONLY)) saveExcitation();
			else saveRRT(); 
		}

		void Raytrace_Inel::sceneChanged(const Scene &S)
		{ 
			// Wenn die Anzahl der Objekte ge�ndert wird, werden die Listen mit den inelastischen Brechungsindizes 
			// und den Polarisierbarkeiten auf die Defaultwerte zur�ckgesetzt 
			
			calcphase = INEL_CALCPHASE_EXCITATION; 
			Raytrace::setScene(S);	
		}
		
		void Raytrace_Inel::exportExcitation(std::string fname, int savetype)
		{
			std::string full_fn;
			for (int i = 0; i < S.nObj; i++)
			{
				full_fn =  fname + "_" + std::to_string(i) + ".dat";
				switch (savetype)
				{
					case INEL_EXPORT_EXCITATION_FIELD_ABS: saveabsE(SGE[0],full_fn, i); break;
					case INEL_EXPORT_EXCITATION_FIELD_VECTOR: saveFullE(SGE[0],full_fn, i); break;
				}
			}
		}

		void Raytrace_Inel::initExcitation()
		{
			if (S.nObj > 0)
			{
				SGE = new SuperArray<maths::Vector<std::complex<double> > >[INEL_MAX_NREFLEX];
				for (int i = 0; i < INEL_MAX_NREFLEX; i++)
				{
					SGE[i] = SuperArray<maths::Vector<std::complex<double> > >(S.r0, n, n, n, IN_OBJECT);
					for (int j = 0; j < S.nObj; j++)
					{
						SGE[i].addInc(S.Obj[j]);
					}
				}
			}
			
		}

		void Raytrace_Inel::initRRT()
		{ 
			// Erst mal den Speicher (Supergitter) allozieren 
			if (S.nObj > 0)
			{
				SGRRT1 = new SuperArray<maths::Vector<std::complex<double> > >[INEL_MAX_NREFLEX]; // 2 wegen der beiden Reflexionsordnungen
				SGRRT2 = new SuperArray<maths::Vector<std::complex<double> > >[INEL_MAX_NREFLEX];

				for (int i = 0; i < INEL_MAX_NREFLEX; i++)
				{
					SGRRT1[i] = SuperArray<maths::Vector<std::complex<double> > >(S.r0, n, n, n, IN_OBJECT);
					SGRRT2[i] = SuperArray<maths::Vector<std::complex<double> > >(S.r0, n, n, n, IN_OBJECT);
					for (int j = 0; j < S.nObj; j++)
					{
						SGRRT1[i].addInc(S.Obj[j]);
						SGRRT2[i].addInc(S.Obj[j]);
					}
				}
			}
		}

		void Raytrace_Inel::traceExcitation()
		{
			S.setRaytype(LIGHTSRC_RAYTYPE_RAY);
			Raytrace::trace();
			for (int i = 1; i < INEL_MAX_NREFLEX; i++)
				SGE[0].add(SGE[i]);
		}

		void Raytrace_Inel::saveExcitation()
		{
			double l,s = 0.0;
			double L = abs(PStop - PStart);
			std::complex<double> phase;
			maths::Vector<double> P=PStart;
			maths::Vector<double> Pnew;
			maths::Vector<std::complex<double> > EG, Eh, Ef;
			maths::Vector<INDEX_TYPE> cell;
			double k0 = S.LS[currentLS]->getWavenumber();
			if ((S.Obj[currentObj]->Active) && (L < 2.0 * S.r0))
			{
				while (s < L)
				{
					SGE[iR].Error=NO_ERRORS;
					Pnew = pnext(P, kin, SGE[iR],currentIndex, 1E-5);
					l = abs(Pnew - P);
					s += l;
					phase = exp(I * (s - l / 2.0) * k0 * S.Obj[currentObj]->n);
					cell = SGE[iR].gitterpunkt((Pnew + P) / 2.0);
					/*if (SGE[iR].Error==SUPERGITTER)
					{
						std::cout << "cell: " << cell[0] << "\t" << cell[1] << "\t" << cell[2] << "\tP:" <<  (Pnew + P) / 2.0 << std::endl;
					}
					else*/
					if (SGE[iR].Error==NO_ERRORS)
					{
					EG = SGE[iR](currentObj, cell);
					Ef = EStart * phase;
					Eh = EG + Ef *sqrt(l);
					Eh = norm(Eh);
					SGE[iR](currentObj, cell) = sqrt(abs2(EG) + abs2(Ef) * l) * Eh;			
					}
					P = Pnew;
				}
			}
		}

		void Raytrace_Inel::traceRRT()
		{
			useRRTParms = true;
			Raytrace::trace();
			useRRTParms = false;
		}

		void Raytrace_Inel::saveRRT()
		{
			double l, s = 0.0;
			double L = abs(PStop - PStart);
			std::complex<double> phase;
			maths::Vector<double> P = PStart;
			maths::Vector<double> Pnew;
			maths::Vector<std::complex<double> > EG, Eh, Ef;
			maths::Vector<INDEX_TYPE> cell;
			double k0 = 2.0 * M_PI / parms.wvlinel;
			if ((S.Obj[currentObj]->Active) && (L < 2.0 * S.r0))
			{
				while (s < L)
				{
					Pnew = pnext(P, kin, SGRRT1[iR], currentIndex, 1E-5);
					l = abs(Pnew - P);
					s += l;
					phase = exp(I * (s - l / 2.0) * k0 * S.Obj[currentObj]->n);
					cell = SGRRT1[iR].gitterpunkt((Pnew + P) / 2.0);
					
					EG = SGRRT1[iR](currentObj, cell);
					Ef = EStart * phase;
					Eh = EG + Ef * sqrt(l);
					Eh = norm(Eh);
					SGRRT1[iR](currentObj, cell) = sqrt(abs2(EG) + abs2(Ef) * l) * Eh;

					EG = SGRRT2[iR](currentObj, cell);
					Ef = EStart2 * phase;
					Eh = EG + Ef * sqrt(l);
					Eh = norm(Eh);
					SGRRT2[iR](currentObj, cell) = sqrt(abs2(EG) + abs2(Ef) * l) * Eh;

					P = Pnew;
				}
			}
		}
	}
}
