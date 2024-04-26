 #include "inel_calc.h"
#include "raytrace_field_usp.h"


namespace GOAT
{
	namespace raytracing
	{		
		void Raytrace_Field_usp::clear()
		{
			if (SA.size() > 0)
			{
				for (int i = 0; i < INEL_MAX_NREFLEX; i++)
					SA[i].clear();
				SA.clear();
			}
		}

                void Raytrace_Field_usp::clean()
		{
                 if (SA.size() > 0)
                     for (int i = 0; i < INEL_MAX_NREFLEX; i++)
                      for (int iBox=0; BoxDetector.size(); iBox)
			for (INDEX_TYPE ix=0; ix<SA[i].n[iBox][0]; ix++) 
			  for (INDEX_TYPE iy=0; iy<SA[i].n[iBox][1]; iy++) 
			    for (INDEX_TYPE iz=0; iz<SA[i].n[iBox][2]; iz++) 
                             SA[i].G[iBox][ix][iy][iz].clear();
		} 

		void Raytrace_Field_usp::addBoxDetectorList(std::vector<Box*> BoxDetector) 
		{
			this->BoxDetector = BoxDetector;
			for (int i = 0; i < INEL_MAX_NREFLEX; i++)
				for (int j = 0; j < BoxDetector.size(); j++)
				SA[i].addInc(BoxDetector[j]);
		}

		void Raytrace_Field_usp::init()
		{
			clear();
			stack.step.clear();
			stack.E = GOAT::maths::Vector<std::complex<double> >(0, 0, 0);
			S.resetLS();
			currentIndex = GOAT::maths::Vector<INDEX_TYPE>(-1, -1, -1);
		//	if (S.nObj > 0)
			{
				SA = std::vector<SuperArray <std::vector<gridEntry>  > >(INEL_MAX_NREFLEX);
				for (int i = 0; i < INEL_MAX_NREFLEX; i++)
				{
					SA[i] = SuperArray<std::vector<gridEntry> >(S.r0, n, n, n, IN_OBJECT);				
				}
			}
		}

		void Raytrace_Field_usp::storeStack(maths::Vector<double> PStart, maths::Vector<double> PStop)
		{
			if (ray->status == RAYBASE_STATUS_FIRST_STEP)
				stack.step.clear();
			stepEntry se;
			se.l = abs(PStop - PStart);
			se.matIndex = S.nObj;  // We have to use the refractive index of the surrounding medium			
			stack.step.push_back(se);
		}


		void Raytrace_Field_usp::storeData(maths::Vector<double> PStart, maths::Vector<double> PStop, maths::Vector<std::complex<double> >EStart)
		{
			double s = 0.0;
			double l;
			double L = abs(PStop - PStart);
			maths::Vector < std::complex<double> >E = EStart;
			maths::Vector<double> P = PStart;
			maths::Vector<double> Pnew;
			maths::Vector<INDEX_TYPE> cell;
			stepEntry ge;
			gridEntry gridStack;
			bool cancel = false;
			// std::cout << PStart << "\t" << PStop << std::endl;
			currentIndex = GOAT::maths::Vector<INDEX_TYPE>(-1, -1, -1);

			if ((L < 2.0 * S.r0) )
			{
				// each cell entry consists of the stack, i.e. all steps until the detector was hidden and the 
				// length of the step from the surface to the cell (last crossing point)
				gridStack.step.insert(gridStack.step.end(), stack.step.begin(), stack.step.end());
				gridStack.step.push_back(ge);
				while ((s < L) && (!cancel))
				{
					Pnew = pnext(P, kin, SA[iR], currentIndex, 1E-100);  // search next grid cell					
					l = abs(Pnew - P);					  // length of the last step  					
					cancel = (l < 1E-15); // cancel, if the step is less than 1E-15µm
					if (cancel) std::cout << "% ABBRUCH !!!!  " << P << "," << l << std::endl;

					s += l;               // path inside the detector
					cell = SA[iR].gitterpunkt((Pnew + P) / 2.0); // get cell index (global)

					// prepare cell entry
					ge.l = s;


					// set the right material index 
					if (currentObj < 0) ge.matIndex = S.nObj;
					else ge.matIndex = currentObj;

					// put everything in the Array
					SA[iR](indexCurrentDetector, cell);
					if (SA[iR].Error == NO_ERRORS)
					{
						gridStack.step.back() = ge;
						SA[iR](indexCurrentDetector, cell).push_back(gridStack);
						SA[iR](indexCurrentDetector, cell).back().E = E;
						// gridStack.step.push_back(ge);						
					}
					else
					{
						SA[iR].Error = NO_ERRORS;
						//std::cout << "ERROR" << std::endl;
					}
					P = Pnew;
				}
			}
		}


		int Raytrace_Field_usp::findBoxDetectorIntersection(maths::Vector<double> P, maths::Vector<double> k, maths::Vector<double>& pStart, maths::Vector<double>& pStop)
		{
			double l = S.r0;
			double lh;
			maths::Vector<double> ph;
			int detectorFound = -1;
			pStart = P;
			if (indexCurrentDetector > -1) 
			{
                              detectorFound=indexCurrentDetector;
                              BoxDetector[indexCurrentDetector]->next(P,k,pStop);
			}
                        else
                        {
				for (int i = 0; i < BoxDetector.size(); i++)
				{
					if (BoxDetector[i]->next(P, k, ph))
					{
						lh = abs(P - ph);
						if (lh < l)
						{
							detectorFound = i;
							pStart = ph;
						}	
					}	
				}
				if (detectorFound > -1) BoxDetector[detectorFound]->next(pStart, k, pStop);

			}
			return detectorFound;
		}

		Raytrace_Field_usp::Raytrace_Field_usp()
		{
		}

		Raytrace_Field_usp::Raytrace_Field_usp(const Scene& S, INDEX_TYPE n) : Raytrace(S)
		{
			this->n = n;		
			init();
		}

		Raytrace_Field_usp::~Raytrace_Field_usp()
		{
			clear();
			BoxDetector.clear();
		}

		void Raytrace_Field_usp::trace()
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
						indexCurrentDetector = -1;
					        stack.step.clear();
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

		void Raytrace_Field_usp::traceOneRay(RayBase* ray, int& Reflexions, int& recur)
		{
			double l;
			std::complex <double> n;
			RayBase* tray = 0;
			double stepSize;
			bool Abbruch = false;
			int objIndex;
			Abbruch = recur > MAX_RECURSIONS;
			recur++;
			gridEntry hstack;
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
						indexCurrentDetector = findBoxDetectorIntersection(PStart, kin, pDetStart, pDetStop);
						if (indexCurrentDetector >= 0) // intersection point with box detector found
						{
							l = abs(PStart - pDetStart);
							if (l < abs(PStart - PStop)) // is the first intersection point with the box detector before the intersection with the next object?
							{
								std::complex <double> n;
								double l1 = abs(pDetStart - pDetStop);
								double l2 = abs(pDetStart - PStop);
								if (l1 > l2) pDetStop = PStop; // second intersection point with the box detector is after the next surface

								storeStack(PStart, pDetStart); // store the path from the starting point to the first intersection with the box detector								
								storeData(pDetStart, pDetStop, EStart); // Store the electric field								
								storeStack(pDetStart, PStop);
							}
							else // intersection point with box detector is outside PStart and PStop
							{
								storeStack(PStart, PStop);
								indexCurrentDetector = -1;
							}
						}

					}
					else // ray is already in a box detector
					{
						BoxDetector[indexCurrentDetector]->next(PStart, kin, pDetStop);						
						double l1 = abs(PStart - pDetStop);
						double l2 = abs(PStart - PStop);
						if (l1 < l2) // the end of the box detector is full in the object or in the surrounding medium
/**
* ATTENTION !!!: In the moment, we can't consider more than one box detector within one step !!!!!!
*/
						{											
							storeData(PStart, pDetStop, EStart);
//							storeStack(pDetStop, PStop); 
							indexCurrentDetector = -1;
						}
						else // step is full in the box detector
						{														
/*							storeStack(P)
							storeData(PStart, PStop, EStart);*/
							storeData(PStart, PStop, EStart);
							storeStack(PStart, PStop);
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
						hstack = stack;
						traceOneRay(tray, tReflexions, recur);
						stack = hstack;
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
							hstack = stack;
							traceOneRay(tray, Reflexions, recur);
							stack = hstack;
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
	}
}

