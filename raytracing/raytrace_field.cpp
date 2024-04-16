#include "raytrace_field.h"
#include "inel_calc.h"
namespace GOAT {
	namespace raytracing {
		Raytrace_Field::Raytrace_Field() : Raytrace()
		{
		}

		Raytrace_Field::Raytrace_Field(Scene& S) : Raytrace(S)
		{
		}

		void Raytrace_Field::traceOneRay(RayBase* ray, int& Reflexions, int& recur)
		{
			double l;
			std::complex <double> n;
			RayBase* tray = 0;
			double stepSize;
			bool Abbruch = false;
			int objIndex;
			Abbruch = recur > MAX_RECURSIONS;
			recur++;
			maths::Vector<double> pDet;
			while ((Reflexions < numReflex) && (!Abbruch))
			{

				// save first the infos at the beginning of the next step
				int oldObjIndex = ray->objectIndex();
				EStart = ray->getE();
				PStart = ray->getP();

				if ((S.raytype == LIGHTSRC_RAYTYPE_IRAY) || useRRTParms) EStart2 = ((IRay*)ray)->E2;
	//			if (S.raytype == LIGHTSRC_RAYTYPE_PRAY) PowIn = ((Ray_pow*)ray)->Pow;

				Abbruch = !ray->next();		// Is there no further crossing with an object ?		
				Abbruch = Abbruch || (abs(EStart) < 10.0 * std::numeric_limits<double>::min()); // Stop, if absolute value of the electric field is smaller than 10*smallest number

				// save now the infos about the state after the step
				objIndex = ray->objectIndex();
				EStop = ray->getE();
				PStop = ray->getP();
				
				if ((S.raytype == LIGHTSRC_RAYTYPE_IRAY) || useRRTParms) EStop2 = ((IRay*)ray)->E2;
				kin = ray->getk();

				// search a hit with a detector within the last step
		/*	if (S.nDet > 0) 
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
						{
							if ((l <= stepSize) && (l > 0)) S.Det[i]->D[i1][i2] += EStart * exp(I * ray->k0 * n * l);
						}
					}
				}
			   */

				if (abs(PStart - PStop) / S.r0 < 10.0 * std::numeric_limits<double>::min()) // if the step is less than 1E-10 times the world radius the program assumes, that the ray hasn't moved => stop calculation
				{
					Abbruch = true;
					lost++;
				}


				if (!Abbruch)
				{
					// search for next box detector 
					if (indexCurrentDetector < 0) // ray is not in a box detector
					{
						indexCurrentDetector = findBoxDetectorIntersection(PStart, kin, pDet);
						if (indexCurrentDetector >= 0) // intersection point with box detector found
						{
							l = abs(PStart - pDet);
							if ( l < abs(PStart - PStop)) // is the intersection point with the box detector before the intersection with the next object?
							{
								std::complex <double> n;
								if (currentObj > -1) n = S.Obj[currentObj]->n;
								else n = S.nS;
								storeData(pDet, PStop, EStart * exp(I * ray->k0 * n * l)); // Store the electric field								
							}
							else // intersection point with box detector is outside PStart and PStop
							{								
								indexCurrentDetector = -1;
							}
						}

					}
					else // ray is already in a box detector
					{
						BoxDetector[indexCurrentDetector]->next(PStart, kin, pDet);
						l = abs(PStart - pDet);
						if ( l < abs(PStart - PStop)) // the end of the box detector is inside the object
						{
							if (objIndex > -1) n = S.Obj[objIndex]->n;
							else n = S.nS;
							storeData(PStart, pDet, EStart);
							indexCurrentDetector = -1;
						}
						else
						{
							l = abs(PStop - PStart);
							if (currentObj > -1) n = S.Obj[currentObj]->n;
							else n = S.nS;
							storeData(PStart, PStop, EStart);
						}
					}
					if (ray->isInObject()) // Is the ray inside an object ?
					{
						if (useRRTParms) ray->reflectRay(tray, -S.Obj[objIndex]->norm(PStop), S.Obj[objIndex]->ninel, S.nS);
						else
						{
							ray->status = RAYBASE_STATUS_NONE;
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
						int tReflexions = -1;
						currentObj = ray->objIndex;
						traceOneRay(tray, tReflexions, recur);
						Reflexions++;
						//delete tray;
					}
					else
						if (objIndex > -1) // an object was hit
						{
							maths::Vector<double> n = S.Obj[objIndex]->norm(PStop);
							//                                  std::cout << "n=" << S.Obj[objIndex]->n << std::endl;
														// std::cout << PStop << "\t" << n << std::endl;
										//				std::cout << "PStart=" << PStart << "\tPStop=" << PStop << "\tn=" << n << std::endl;
							if (useRRTParms)
							{
								copyRay(tray, ray);
								ray->reflectRay(tray, n, S.nS, S.Obj[objIndex]->n);
							}
							else
							{
								copyRay(tray, ray);
								ray->reflectRay(tray, n, S.nS, S.Obj[objIndex]->n);
								//					std::cout << "*k=" << ((tubedRay*)tray)->k[4]  <<  std::endl;
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
							ray->status = RAYBASE_STATUS_NONE;
							tray->status = RAYBASE_STATUS_NONE;
							int tReflexions = -1;
							traceOneRay(tray, Reflexions, recur);
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

		void Raytrace_Field::traceEnterObject()
		{
		}

		void Raytrace_Field::traceLeaveObject()
		{
		}

		void Raytrace_Field::addBoxDetector(Box* box)
		{
			BoxDetector.push_back(box);			
		}

		void Raytrace_Field::storeData(maths::Vector<double> PStart, maths::Vector<double> PStop, maths::Vector<std::complex<double>> EStart)
		{
			double s = 0.0;
			double l;
			double L = abs(PStop - PStart);
			maths::Vector < std::complex<double> >E = EStart;
			maths::Vector<double> k = (PStop - PStart);
			k /= abs(k);
			maths::Vector<double> P = PStart;
			maths::Vector<double> Pnew;
			maths::Vector<INDEX_TYPE> cell;
			std::complex<double> n;
			if (currentObj > -1) n = S.Obj[currentObj]->n;
			else n = S.nS;
			bool cancel = false;
			while ((s < L) && (!cancel))
			{
				Pnew = pnext(P, k, SE, 1E-100);  // search next grid cell		
				l = abs(Pnew - P);					  // length of the last step  					
				cancel = (l < 1E-15); // cancel, if the step is less than 1E-15µm
				s += l;               // path inside the detector
				cell = SE.gitterpunkt((Pnew + P) / 2.0); // get cell index (global)			    
				SE(indexCurrentDetector, cell);
				if (SE.Error == NO_ERRORS)
				SE(indexCurrentDetector, cell) += E * exp(I * ray->k0 * n * s);
				else
					SE.Error = NO_ERRORS;
				P = Pnew;
			}
		}

		int Raytrace_Field::findBoxDetectorIntersection(maths::Vector<double> P, maths::Vector<double> k, maths::Vector<double>& pout)
		{
			double l = S.r0;
			double lh;
			maths::Vector<double> ph;
			int detectorFound=-1;
			pout = P;
			for (int i = 0; i < BoxDetector.size(); i++)
			{
				if (BoxDetector[i]->next(P, k, ph))
				{
					lh = abs(P - ph);
					if (lh < l)
					{
						detectorFound = i;
						pout = ph;
					}
				}
			}
			return detectorFound;
		}
		void Raytrace_Field::init()
		{
			setResolution(resolution);
			SE.clear();	
			SE = SuperArray<maths::Vector<std::complex<double> > >(S.r0, numCellsPerDirection, numCellsPerDirection, numCellsPerDirection);
			for (int i = 0; i < BoxDetector.size(); i++)
			{
				BoxDetector[i]->setActive(true);
				SE.addInc(BoxDetector[i]);
			}
		}
		void Raytrace_Field::trace()
		{
			SE.r0=S.r0;
			init();
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
				for (int i = 0; i < S.nLS; i++) // Schleife �ber die Lichtquellen
				{
					S.resetLS();
					do
					{
						//						std::cout << "%---------------------------------" << std::endl;
						currentLS = i;
						Abbruch = false;
						Reflexions = -1;
						statusLS = S.LS[i]->next(ray);
						recursions = 0;
						ray->status = RAYBASE_STATUS_FIRST_STEP;
						ray->suppress_phase_progress = S.suppress_phase_progress;
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
					ray->suppress_phase_progress = S.suppress_phase_progress;
					recursions = 0;
					traceOneRay(ray, Reflexions, recursions); // Verfolgung eines Teilstrahls
				} while (statusLS != LIGHTSRC_IS_LAST_RAY);
				//	delete ray;
				//	delete tray;
			}
		}

		void Raytrace_Field::setResolution(double res)
		{
			numCellsPerDirection = S.r0 / res;
			resolution = S.r0 / (INDEX_TYPE)numCellsPerDirection;
		}
	}
}
