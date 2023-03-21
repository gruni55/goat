#include "raytrace_usp.h"
#include "inel_calc.h"
namespace GOAT
{
	namespace raytracing
	{
		Raytrace_usp::Raytrace_usp()
		{
		}

		Raytrace_usp::Raytrace_usp(const Scene& S, int n) : Raytrace(S)
		{
			this->n = n;
		}

		void Raytrace_usp::init()
		{
			currentIndex = GOAT::maths::Vector<int>(-1, -1, -1);
			if (S.nObj > 0)
			{			
				SA = std::vector<SuperArray <std::vector<gridEntry>  > >(INEL_MAX_NREFLEX);
				for (int i = 0; i < INEL_MAX_NREFLEX; i++)
				{
					SA[i] = SuperArray<std::vector<gridEntry> > (S.r0, n, n, n, IN_OBJECT);
					for (int j = 0; j < S.nObj; j++)
						SA[i].addInc(S.Obj[j]);
				}
			}
		}

		void Raytrace_usp::trace()
		{
			init();
			S.setRaytype(LIGHTSRC_RAYTYPE_RAY);
			S.suppress_phase_progress = true;
			Raytrace::trace();
		/*	for (int i = 1; i < INEL_MAX_NREFLEX; i++)
				SA[0].add(SA[i]);*/
		}
		
		void Raytrace_usp::storeData()
		{
			double s = 0.0;
			double l;
			double L = abs(PStop - PStart);
			maths::Vector < std::complex<double> >E=EStart;
			maths::Vector<double> P = PStart;
			maths::Vector<double> Pnew;
			maths::Vector<int> cell;
			stepEntry ge;
			gridEntry gridStack;		
			bool cancel = false;
			currentIndex = GOAT::maths::Vector<int>(-1, -1, -1);
			// std::cout << PStart << "/" << PStop << std::endl;
			if (L < 2.0 * S.r0)
			{	
				// each cell entry consists of the stack, i.e. all steps until the detector was hidden and the 
				// length of the step from the surface to the cell (last crossing point)
				gridStack.step.insert(gridStack.step.end(), stack.step.begin(), stack.step.end());
						
				while ( (s < L) && (!cancel) )
				{
					Pnew = pnext(P, kin, SA[iR],1E-100);  // search next grid cell					
					l = abs(Pnew - P); // length of the last step  					
					cancel = (l < 1E-15); // cancel, if the step is less than 1E-15µm
					if (cancel) std::cout << "ABBRUCH !!!!  " << P << "," << l << std::endl;
										
					s += l;               // path inside the detector
					cell = SA[iR].gitterpunkt((Pnew + P) / 2.0); // get cell index (global)

					// prepare cell entry
					ge.l = abs(P-Pnew); 

					// set the right material index 
					if (currentObj < 0) ge.matIndex = S.nObj;
					else ge.matIndex = currentObj;
					
					// put everything in the Array
					SA[iR](currentObj, cell);
				    if (SA[iR].Error == NO_ERRORS)
					{						
						SA[iR](currentObj, cell).push_back(gridStack);							
						SA[iR](currentObj, cell).back().E = E;
 						gridStack.step.push_back(ge);
					}		
					else
					{
						SA[iR].Error = NO_ERRORS;
						std::cout << "ERROR" << std::endl;
					}
					P = Pnew;					
				}
			}
		}
		void Raytrace_usp::traceEnterObject()
		{
			if (ray->status == RAYBASE_STATUS_FIRST_STEP)
				stack.step.clear();
			stepEntry se;
			se.l = abs(PStop - PStart);
			se.matIndex = S.nObj;  // We have to use the refractive index of the surrounding medium			
			stack.step.push_back(se);
		}

		void Raytrace_usp::traceLeaveObject()
		{
			stepEntry se;
			se.l = abs(PStop - PStart);
			se.matIndex = currentObj;			
			storeData();
			stack.step.push_back(se);			
		}
			
	}

}
