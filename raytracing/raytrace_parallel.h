#pragma once
#include "raytrace.h"
#include <queue>
#include <mutex>
#include <barrier>

namespace GOAT
{
	namespace raytracing
	{
		template <class T> class RaytraceParallel
		{
		public: 
			RaytraceParallel(Scene &S);
			void trace();
			Scene S;
			private:
				RayBase* getNextRay(int iLS);
				
				bool useRRTParms = false;
		};


		template <class T> RaytraceParallel<T>::RaytraceParallel(Scene& S) : S(S)
		{
			this->S = S;
		}

		template	<class T> RayBase* RaytraceParallel<T>::getNextRay(int iLS)
		{
			int statusLS;
			RayBase* ray = nullptr;
			if (!useRRTParms)
			{
				switch (S.raytype)
				{
				case LIGHTSRC_RAYTYPE_IRAY: ray = new IRay; break;
				case LIGHTSRC_RAYTYPE_PRAY: ray = new Ray_pow; break;
				case LIGHTSRC_RAYTYPE_RAY:
				default: ray = new tubedRay; 
				}

				statusLS = S.LS[iLS]->next(ray);
				if (statusLS == LIGHTSRC_IS_LAST_RAY)
					return nullptr;
				else
					return ray;

			}
			else
			{
				statusLS = S.LSRRT->next(ray);
				ray->status = RAYBASE_STATUS_FIRST_STEP;
				ray->suppress_phase_progress = S.suppress_phase_progress;
				return ray;
			}
			return nullptr;
		}

		template <class T> void RaytraceParallel<T>::trace()
		{
			S.resetLS();
			std::vector<std::thread> threads;
			
			// Create worker threads
			std::barrier startBarrier(S.getNumberOfThreads() + 1);

			for (int i = 0; i < S.getNumberOfThreads(); i++)
			{
				threads.emplace_back([this, &startBarrier]()
					{
						Scene Slocal = S;
						copyLightSrcList(Slocal.LS, S.LS, S.getNumberOfLightSources());

						for (int l = 0; l < Slocal.getNumberOfLightSources(); l++)
						{
							Slocal.LS[l]->setNumRays(S.LS[l]->getNumRays() / S.getNumberOfThreads());
						}
						T rt(Slocal);

						startBarrier.arrive_and_wait();  // warten

						rt.trace();                      // alle starten danach
					});
			}

			startBarrier.arrive_and_wait();          // Hauptthread gibt Start frei
			


			
			for (auto& t : threads) // Wait for all threads to finish
				t.join();
		}
	}
}