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
			if (S.nObj > 0)
			{			
				SA = std::vector<SuperArray <gridEntry>>(INEL_MAX_NREFLEX);
				for (int i = 0; i < INEL_MAX_NREFLEX; i++)
				{
					SA[i] = SuperArray<gridEntry>(S.r0, n, n, n, IN_OBJECT);
					for (int j = 0; j < S.nObj; j++)
						SA[i].addInc(S.Obj[j]);
				}
			}
		}

		void Raytrace_usp::trace()
		{
			S.setRaytype(LIGHTSRC_RAYTYPE_RAY);
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

			if (L < 2.0 * S.r0)
			{
				while (s < L)
				{
					Pnew = pnext(P, kin, SA[iR], 1E-5);
					l = abs(Pnew - P);
					s += l;
					cell = SA[iR].gitterpunkt((Pnew + P) / 2.0);
					ge.l = l;
					ge.matIndex = currentObj;
					SA[iR](currentObj, cell).step.push_back(ge);
					SA[iR](currentObj, cell).E = E;
					P = Pnew;
				}
			}
		}
		void Raytrace_usp::traceEnterObject()
		{
		}

		void Raytrace_usp::traceLeaveObject()
		{
			storeData();
		}
			
	}

}
