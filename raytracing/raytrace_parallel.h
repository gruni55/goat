#pragma once
#include "raytrace.h"
#include <queue>
#include <mutex>
#include <barrier>

namespace GOAT
{
	namespace raytracing
	{
		/** @brief This class implements a parallel raytracing algorithm using multiple threads.
		* This class is a template class that takes a raytracing class as a template parameter. 
		* All rays are generated in parallel by multiple threads, each thread processes a subset of the rays. 
		* The results are then combined to produce the final output. The class uses a barrier to synchronize the threads and ensure that all threads have completed 
		* their work before proceeding to the next step. Here an example, how to use this class:
		* (Suppose we already have a Scene object S)
		* @code
		*   S.setNumberOfThreads(4); // Set the number of threads to use for parallel raytracing
		*	GOAT::raytracing::RaytraceParallel<GOAT::raytracing::Raytrace_pure> rtParallel(S); // Create a RaytraceParallel object with the Scene S with the template parameter Raytrace_pure
		*   rtParallel.trace(); // Start the parallel raytracing process
		* @endcode
		* The result is provided in the detectors of the Scene object S. 
		* The class can be used with any raytracing class that implements the required interface, such as Raytrace_pure, Raytrace_Inel, Raytrace_Path, etc.
		* Up to now, the class can only be used together with MC-light sources. 
		*/

		template <class T> class RaytraceParallel
		{
		public: 
			RaytraceParallel(Scene &S);
			void trace();
			void requestStop() { stopFlag->store(true); } ///< Request to stop the raytracing process
			void setStopFlag(std::shared_ptr<std::atomic<bool>> flag) { stopFlag = flag; } ///< Sets the stop flag for the raytracing process
			Scene S;
			
			private:
				RayBase* getNextRay(int iLS);
				
				bool useRRTParms = false;
				std::shared_ptr<std::atomic<bool>> stopFlag =
				std::make_shared<std::atomic<bool>>(false); ///< flag to stop calculation, e.g. if the user wants to stop the calculation
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
			std::barrier<> startBarrier{ S.getNumberOfThreads() + 1 };

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
						rt.setStopFlag(stopFlag);
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