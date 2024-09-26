#include "raytrace_usp_rt.h"
#include "raytrace_inel.h"
#include "inel_calc.h"
namespace GOAT
{
	namespace raytracing
	{
		Raytrace_usp_rt::Raytrace_usp_rt()
		{}
		

		Raytrace_usp_rt::~Raytrace_usp_rt()
		{
			clear();
			nList.clear();
		}

		Raytrace_usp_rt::Raytrace_usp_rt(const Scene& S, INDEX_TYPE n) : Raytrace(S)
		{
			this->n = n;
			init();
		}

		void Raytrace_usp_rt::clear()
		{
			if (SA.size() > 0)
			{
				for (int i = 0; i < INEL_MAX_NREFLEX; i++)
					SA[i].clear();
				SA.clear();
			}
		}

		void Raytrace_usp_rt::init()
		{
			clear();
			
			S.resetLS();
//			currentIndex = GOAT::maths::Vector<INDEX_TYPE>(-1, -1, -1);
			if (S.nObj > 0)
			{
				SA = std::vector<SuperArray <maths::Vector<std::complex<double> > > >(INEL_MAX_NREFLEX);
				for (int i = 0; i < INEL_MAX_NREFLEX; i++)
				{
					SA[i] = SuperArray <maths::Vector<std::complex<double> > > (S.r0, n, n, n, IN_OBJECT);
					for (int j = 0; j < S.nObj; j++)
						SA[i].addInc(S.Obj[j]);
				}
			}
		}

		void Raytrace_usp_rt::trace(double omega, std::complex<double> weight)
		{
			this->weight = weight;
			double wvl = 2.0 * M_PI * C_LIGHT_MU_FS / omega;			
			k0 = 2.0 * M_PI / wvl;
			S.setRaytype(LIGHTSRC_RAYTYPE_IRAY);			
			for (int i = 0; i < S.nObj; i++)
				S.Obj[i]->setn(nList[i](wvl));

			for (int i = 0; i < S.nLS; i++)
				S.LS[i]->setWavelength(wvl);

			Raytrace::trace();			
		}

		void Raytrace_usp_rt::storeData()
		{
			double s = 0.0;
			double l;
			double L = abs(PStop - PStart);
			maths::Vector < std::complex<double> >E = EStart;
			maths::Vector<double> P = PStart;
			maths::Vector<double> Pnew;
			maths::Vector<INDEX_TYPE> cell;
			bool cancel = false;
			
//			currentIndex = GOAT::maths::Vector<INDEX_TYPE>(-1, -1, -1);

			if ((L < 2.0 * S.r0) && S.Obj[currentObj]->isActive())
			{				
				while ((s < L) && (!cancel))
				{
					Pnew = pnext(P, kin, SA[iR], 1E-100);  // search next grid cell					
					l = abs(Pnew - P);					  // length of the last step  					
					cancel = (l < 1E-15); // cancel, if the step is less than 1E-15µm
					if (cancel) std::cout << "% Abort !!!!  " << P << "," << l << std::endl;
					s += l;               // path inside the detector
					cell = SA[iR].gitterpunkt((Pnew + P) / 2.0); // get cell index (global)

					// put everything in the Array
					SA[iR](currentObj, cell);
					if (SA[iR].Error == NO_ERRORS)
					{
						SA[iR](currentObj, cell) += EStart * exp(I * k0 * s) * weight;
					}
					else
					{
						SA[iR].Error = NO_ERRORS;
					}
					P = Pnew;
				}
			}
		}

		void Raytrace_usp_rt::setRefractiveIndexFunctions(std::vector<std::function<std::complex<double>(double)>> nList)
		{
			this->nList = nList;
		}

		void Raytrace_usp_rt::traceEnterObject()
		{						
		}

		void Raytrace_usp_rt::traceLeaveObject()
		{
			storeData();			
		}


	}
}
