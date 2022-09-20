#include "raytrace.h"
#include <float.h>

namespace GOAT
{
	namespace raytracing
	{
#define MAX_RECURSIONS 20
		Raytrace::Raytrace()
		{
			init();
		}

		Raytrace::Raytrace(const Scene& S)
		{
			init();
			this->S = S;
			numReflex = RAYTRACE_MAX_REFLEXIONS;
		}

		void Raytrace::setScene(const Scene& S)
		{
			this->S = S;
		}

		void Raytrace::setNumReflex(int numReflex)
		{
			this->numReflex = numReflex;
		}


		void Raytrace::trace()
		{
			/*RayBase* tray;
			RayBase* ray;
			*/
			int statusLS;
			int Reflexions = 0;
			int recursions = 0;
			lost = 0;
			switch (S.raytype)
			{
			case LIGHTSRC_RAYTYPE_IRAY: ray = new IRay; tray = new IRay;  break;
			case LIGHTSRC_RAYTYPE_PRAY: ray = new Ray_pow; tray = new Ray_pow; break;
			case LIGHTSRC_RAYTYPE_RAY:
			default: ray = new tubedRay;  tray = new tubedRay;
			}

			if (!useRRTParms)
				for (int i = 0; i < S.nLS; i++) // Schleife über die Lichtquellen
				{
					S.resetLS();
					do
					{
						currentLS = i;
						Abbruch = false;
						Reflexions = 0;
						statusLS = S.LS[i]->next(ray);
						recursions = 0;
						traceOneRay(ray, Reflexions, recursions); // Verfolgung eines Teilstrahls
					} while (statusLS != LIGHTSRC_IS_LAST_RAY);
					//	delete ray;
						//delete tray;
				}
			else
			{
				ray = new IRay; tray = new IRay;

				do
				{
					Abbruch = false;
					Reflexions = 0;
					statusLS = S.LSRRT->next(ray);
					recursions = 0;
					traceOneRay(ray, Reflexions, recursions); // Verfolgung eines Teilstrahls
				} while (statusLS != LIGHTSRC_IS_LAST_RAY);
				//	delete ray;
				//	delete tray;
			}
		}


		void Raytrace::traceOneRay(RayBase* ray, int& Reflexions, int& recur)
		{
			RayBase* tray = 0;
			double stepSize;
			bool Abbruch = false;
			int objIndex;
			Abbruch = recur > MAX_RECURSIONS;
			recur++;

			while ((Reflexions <= numReflex) && (!Abbruch))
			{

				int oldObjIndex = ray->objectIndex();

				EStart = ray->getE();
				PStart = ray->getP();


				if ((S.raytype == LIGHTSRC_RAYTYPE_IRAY) || useRRTParms) EStart2 = ((IRay*)ray)->E2;
				if (S.raytype == LIGHTSRC_RAYTYPE_PRAY) PowIn = ((Ray_pow*)ray)->Pow;
				Abbruch = !ray->next();
				objIndex = ray->objectIndex();
				EStop = ray->getE();
				PStop = ray->getP();


				if ((S.raytype == LIGHTSRC_RAYTYPE_IRAY) || useRRTParms) EStop2 = ((IRay*)ray)->E2;
				kin = ray->getk();

				if (S.nDet > 0)
				{
					int i1, i2;
					double l;
					stepSize = abs(PStop - PStart);
					std::complex<double> n;
					if (ray->isInObject() && (objIndex > -1)) n = S.Obj[objIndex]->n;
					else n = S.nS;
					for (int i = 0; i < S.nDet; i++)
					{
						if (S.Det[i]->cross(PStart, kin, i1, i2, l))
							if ((l <= stepSize) && (l > 0)) S.Det[i]->D[i1][i2] += EStart * exp(I * ray->k0 * n * l);
					}
				}

				if (abs(PStart - PStop) / S.r0 < 1E-10)
				{
					Abbruch = true;
					lost++;
				}

				if (!Abbruch)
				{
					if (ray->isInObject()) // Strahl befindet sich im Objekt
					{
						if (useRRTParms) ray->reflectRay(tray, -S.Obj[objIndex]->norm(PStop), S.Obj[objIndex]->ninel, S.nS);
						else
						{
							copyRay(tray, ray);
							ray->reflectRay(tray, -S.Obj[objIndex]->norm(PStop), S.Obj[objIndex]->n, S.nS);
						}
						kref = ray->getk();
						ktrans = tray->getk();
						if (S.raytype == LIGHTSRC_RAYTYPE_PRAY)
						{
							PowRef = ((Ray_pow*)ray)->Pow;
							PowTrans = ((Ray_pow*)tray)->Pow;
						}
						traceLeaveObject();
						int tReflexions = 0;
						currentObj = ray->objIndex;
						traceOneRay(tray, tReflexions, recur);

						Reflexions++;
						//delete tray;
					}
					else
						if (objIndex > -1)
						{
							maths::Vector<double> n = S.Obj[objIndex]->norm(PStop);
							if (useRRTParms)
							{
								copyRay(tray, ray);
								ray->reflectRay(tray, n, S.nS, S.Obj[objIndex]->n);
							}
							else
							{
								copyRay(tray, ray);
								ray->reflectRay(tray, n, S.nS, S.Obj[objIndex]->n);
							}
							kref = ray->getk();
							ktrans = tray->getk();
							currentObj = objIndex;
							if (S.raytype == LIGHTSRC_RAYTYPE_PRAY)
							{
								PowRef = ((Ray_pow*)ray)->Pow;
								PowTrans = ((Ray_pow*)tray)->Pow;
							}
							traceEnterObject();
							int tReflexions = 0;
							traceOneRay(tray, Reflexions, recur);
							//delete tray;
							Reflexions++;
							Abbruch = true;
						}
						else
						{
							traceStopObject();
							Abbruch = true;
						}
				}
				if (tray != 0) { delete tray; tray = 0; }
			}
		}

		void Raytrace::copyRay(RayBase*& dest, RayBase* src)
		{
			switch (S.raytype)
			{
			case LIGHTSRC_RAYTYPE_IRAY: dest = new IRay(*(IRay*)src);  break;
			case LIGHTSRC_RAYTYPE_PRAY: dest = new Ray_pow(*(Ray_pow*)src); break;
			case LIGHTSRC_RAYTYPE_RAY:
			default: dest = new tubedRay(*(tubedRay*)src);
			}

		}


		void Raytrace::init()
		{
			useRRTParms = false;
			ray = 0;
			tray = 0;
			Abbruch = false;
			PowIn = 0;
			PowRef = 0;
			currentLS = -1;
			currentObj = -1;
			lost = 0;
		}

		Raytrace_OT::Raytrace_OT()
		{
			F = 0;
			L = 0;
		}

		Raytrace_OT::Raytrace_OT(Scene S)
		{
			F = 0;
			L = 0;
			S.setRaytype(LIGHTSRC_RAYTYPE_PRAY);
			this->S = S;
			//	for (int i = 0; i < nLS; i++) LS[i]->raytype = LIGHTSRC_RAYTYPE_PRAY;
		}

		void Raytrace_OT::setScene(const Scene& S)
		{
			Raytrace::setScene(S);
			this->S = S;
			this->S.setRaytype(LIGHTSRC_RAYTYPE_PRAY);
		}

		void Raytrace_OT::trace()
		{
			if (F != 0) delete F;
			if (L != 0) delete L;
			F = new maths::Vector<double>[S.nObj];
			L = new maths::Vector<double>[S.nObj];

			f = new maths::Vector<double> *[S.nLS];
			l = new maths::Vector<double> *[S.nLS];
			for (int i = 0; i < S.nLS; i++)
			{
				f[i] = new maths::Vector<double>[S.nObj];
				l[i] = new maths::Vector<double>[S.nObj];
			}

			Raytrace::trace();

			double I = 0;
			for (int i = 0; i < S.nLS; i++)
				I += S.LS[i]->Pall;

			for (int j = 0; j < S.nObj; j++)
			{
				F[j] = maths::dzero;
				L[j] = maths::dzero;
				for (int i = 0; i < S.nLS; i++)
				{
					F[j] += f[i][j] / I * S.LS[i]->P0 * real(S.nS) / c_light;
					L[j] += l[i][j] * 1E-6 / I * S.LS[i]->P0 * real(S.nS) / c_light;
				}
			}
			for (int i = 0; i < S.nLS; i++)
			{
				delete[] f[i];
				delete[] l[i];
			}
			delete[] f;
			delete[] l;
		}

		void Raytrace_OT::traceLeaveObject()
		{
			// std::cout << PStart << "   " << PStop << std::endl;
			maths::Vector<double> fe, fr, ft, fg;
			maths::Vector<double> r;

			fe = (kin * PowIn) * real(S.nS / S.Obj[currentObj]->n);
			fr = (-kref * PowRef);
			ft = (-ktrans * PowTrans);
			fg = ft + fe + fr;
			r = PStop - S.Obj[currentObj]->P;
			f[currentLS][currentObj] += fg;
			l[currentLS][currentObj] += r % fg;
		}

		void Raytrace_OT::traceEnterObject()
		{
			// std::cout << PStart << "   " << PStop << std::endl;
			maths::Vector<double> fe, fr, ft, fg;
			maths::Vector<double> r;

			fe = (kin * PowIn);
			fr = (-kref * PowRef);
			ft = (-ktrans * PowTrans) * real(S.nS / S.Obj[currentObj]->n);
			fg = fe + fr + ft;
			r = PStop - S.Obj[currentObj]->P;
			f[currentLS][currentObj] += fg;
			l[currentLS][currentObj] += r % fg;
		}

		Scene::Scene()
		{
			nLS = 0;
			LS = 0;
			nObj = 0;
			Obj = 0;
			r0 = DBL_MAX;
			nS = 1.0;
			raytype = LIGHTSRC_RAYTYPE_RAY;
			LS = 0;
			LSRRT = 0;
			nDet = 0;
		}


		void Scene::addObject(ObjectShape* obj)
		{
			if (nObj == 0)
				Obj = (ObjectShape**)malloc(sizeof(ObjectShape*));
			else
				Obj = (ObjectShape**)realloc(Obj, sizeof(ObjectShape*) * (nObj + 1));
			Obj[nObj] = obj;
			obj->r0 = r0;
			obj->initQuad();
			nObj++;
		}

		void Scene::addObjectList(int nobj, ObjectShape** obj)
		{
			for (int i = 0; i < nobj; i++)
				addObject(obj[i]);
		}

		void Scene::removeAllObjects()
		{
			if (nObj > 0)
			{
				free(Obj);
				for (int i = 0; i < nLS; i++)  // remove objects from all light sources
					LS[i]->clearObjects();
				nObj = 0;
			}
		}

		void Scene::removeObject(int index)
		{
			if ((index < nObj) && (index >= 0))
			{
				for (int i = index; i < nObj - 1; i++)
					Obj[i] = Obj[i + 1];
				nObj--;
				if (nObj < 0) nObj = 0;
			}
		}


		void Scene::addLightSource(LightSrc* ls, int raytype)
		{
			if (nLS == 0)
				LS = (LightSrc**)malloc(sizeof(LightSrc));
			else
				LS = (LightSrc**)realloc(LS, sizeof(LightSrc) * (nLS + 1));
			/*
				switch (ls->type)
				{
				case LIGHTSRC_SRCTYPE_PLANE: LS[nLS] = new LightSrcPlane(*(LightSrcPlane*)ls); break;
				case LIGHTSRC_SRCTYPE_GAUSS: LS[nLS] = new LightSrcGauss(*(LightSrcGauss*)ls); break;
				}
				*/
			LS[nLS] = ls;
			LS[nLS]->clearObjects();
			if (nObj > 0) LS[nLS]->ObjectList(nObj, Obj);
			LS[nLS]->raytype = raytype;
			LS[nLS]->setR0(r0);
			LS[nLS]->setN0(nS);
			nLS++;
		}

		void Scene::removeLightSrc(int index)
		{
			if ((index < nLS) && (index >= 0))
			{
				for (int i = index; i < nLS - 1; i++)
					LS[i] = LS[i + 1];
				nLS--;
				if (nLS < 0) nLS = 0;
			}
		}

		void Scene::removeAllLightSources()
		{
			if (nLS > 0)
			{
				free(LS);
				nLS = 0;
			}
		}

		void Scene::addLightSourceRRT(LightSrc* ls, maths::Vector<std::complex<double> >Pol1, maths::Vector<std::complex<double> >Pol2)
		{
			if (LSRRT != 0) delete LSRRT;
			/*switch (ls->type)
			{
			case LIGHTSRC_SRCTYPE_PLANE: LSRRT = new LightSrcPlane(*(LightSrcPlane*)ls); break;
			case LIGHTSRC_SRCTYPE_GAUSS: LSRRT = new LightSrcGauss(*(LightSrcGauss*)ls); break;
			}
			*/
			LSRRT = ls;
			LSRRT->clearObjects();
			if (nObj > 0) LSRRT->ObjectList(nObj, Obj);
			LSRRT->raytype = LIGHTSRC_RAYTYPE_IRAY;
			LSRRT->setR0(r0);
			LSRRT->setN0(nS);
			LSRRT->Pol = Pol1;
			LSRRT->Pol2 = Pol2;
		}

		void Scene::addLightSourceList(int nls, LightSrc** ls)
		{
			for (int i = 0; i < nls; i++)
			{
				addLightSource(ls[i]);
				ls[i]->raytype = raytype;
			}
		}

		void Scene::addDetector(Detector* D)
		{
			if (nDet == 0)
				Det = (Detector**)malloc(sizeof(Detector));
			else
				Det = (Detector**)realloc(Det, sizeof(Detector) * (nDet + 1));
			Det[nDet] = D;
			nDet++;
		}

		void Scene::addDetectorList(int nDet, Detector** D)
		{
			for (int i = 0; i < nDet; i++) addDetector(D[i]);
		}

		void Scene::cleanAllDetectors()
		{
			if (nDet > 0)
				for (int i = 0; i < nDet; i++) Det[i]->clean();
		}

		void Scene::removeAllDetectors()
		{
			if (nDet > 0)
			{
				/*for (int i = 0; i < nDet; i++)
					delete Det[i];*/
				free(Det);
				nDet = 0;
			}
		}

		void Scene::removeDetector(int index)
		{
			if ((index < nDet) && (index >= 0))
			{
				for (int i = index; i < nDet - 1; i++)
					Det[i] = Det[i + 1];
				nDet--;
				if (nDet < 0) nDet = 0;
			}
		}


		void Scene::setr0(double r0)
		{
			this->r0 = r0;
			if (nLS > 0)
				for (int i = 0; i < nLS; i++)
					LS[i]->setR0(r0);

			if (nObj > 0)
				for (int i = 0; i < nObj; i++)
					Obj[i]->setr0(r0);
		}

		void Scene::setnS(std::complex<double> nS)
		{
			if (nLS > 0)
				for (int i = 0; i < nLS; i++)
				{
					LS[i]->setN0(nS);
				}
			this->nS = nS;
			this->nSRRT = nS;
		}

		void Scene::setnSRRT(std::complex<double> nS)
		{
			this->nSRRT = nS;
		}

		Scene::Scene(const Scene& S)
		{
			LSRRT = S.LSRRT;
			LS = S.LS;
			nLS = S.nLS;
			Obj = S.Obj;
			nObj = S.nObj;
			r0 = S.r0;
			nS = S.nS;
			raytype = S.raytype;
			Det = S.Det;
			nDet = S.nDet;
		}

		void Scene::setRaytype(int raytype)
		{
			if (nLS > 0)
				for (int i = 0; i < nLS; i++)
					LS[i]->raytype = raytype;
			this->raytype = raytype;
		}

		void Scene::resetLS()
		{
			if (nLS > 0)
				for (int i = 0; i < nLS; i++)
				{
					LS[i]->clearObjects();
					LS[i]->ObjectList(nObj, Obj);
					LS[i]->reset();
				}
		}

		int Scene::testLS()
		{
			return -1;
		}

		Raytrace_Path::Raytrace_Path() : Raytrace()
		{
		}

		Raytrace_Path::Raytrace_Path(Scene S) : Raytrace(S)
		{
		}

		void Raytrace_Path::trace(std::string FName)
		{
			//S.setRaytype(LIGHTSRC_RAYTYPE_IRAY);
			os.open(FName);
			storeInFile = true;
			Raytrace::trace();
			os.close();
		}

		void Raytrace_Path::trace()
		{
			storeInFile = false;
			if (numRays > 0)
			{
				if (P1 != 0) free(P1);
				if (P2 != 0) free(P2);
			}
			numRays = 0;
			Raytrace::trace();
			numRays--;
		}

		void Raytrace_Path::setShowOutgoingRays(bool show)
		{
			showOutgoingRays = show;
		}

		bool Raytrace_Path::getShowOutgoingRays()
		{
			return showOutgoingRays;
		}

		void Raytrace_Path::traceLeaveObject()
		{
			storeStep();
		}

		void Raytrace_Path::traceEnterObject()
		{
			storeStep();
		}

		void Raytrace_Path::storeStep()
		{
			if (storeInFile)
				os << PStart << "\t" << PStop << std::endl;
			else
			{
				if (numRays == 0)
				{
					P1 = (maths::Vector<double> *) malloc(sizeof(maths::Vector<double>));
					P2 = (maths::Vector<double> *) malloc(sizeof(maths::Vector<double>));
				}
				else
				{
					P1 = (maths::Vector<double> *) realloc(P1, sizeof(maths::Vector<double>) * (numRays + 1));
					P2 = (maths::Vector<double> *) realloc(P2, sizeof(maths::Vector<double>) * (numRays + 1));
				}
				P1[numRays] = PStart;
				P2[numRays] = PStop;
				numRays++;
			}
		}

		void Raytrace_Path::traceStopObject()
		{
			if (showOutgoingRays) storeStep();
		}

		Raytrace_pure::Raytrace_pure() : Raytrace()
		{

		}

		Raytrace_pure::Raytrace_pure(const Scene& S) : Raytrace(S)
		{}
	}
}