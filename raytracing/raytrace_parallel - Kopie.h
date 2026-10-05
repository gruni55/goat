#pragma once
#include "raytrace.h"
#include <queue>
#include <mutex>

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
			std::queue <RayBase*> rayQueue;
			std::vector<std::thread> threads;
			std::mutex queueMutex;
			std::condition_variable cv;

			// Create worker threads
			for (int i=0; i<S.getNumberOfThreads(); i++)
			{
				
				threads.emplace_back([this, i, &rayQueue, &queueMutex, &cv]()
					{
						T rt(S);
						while (true)
						{

							RayBase* ray;
							{
								std::unique_lock<std::mutex> lock(queueMutex);
								cv.wait(lock, [&rayQueue]() { return !rayQueue.empty(); });
								ray = rayQueue.front();
								rayQueue.pop();
								
							} // The lock is released here, allowing other threads to access the queue
							if (ray == nullptr) return; // Exit the thread if a nullptr is encountered
							int Reflexions = 0;
							int recursions = 0;
							rt.traceOneRay(ray, Reflexions, recursions);
							delete ray; // Clean up the ray after processing
						}
					});
			}
			


			for (int iLS = 0; iLS < S.getNumberOfLightSources(); iLS++)
			{
				int c = 0;
				do 
				{
					   
						RayBase* ray = getNextRay(iLS);
						if (ray == nullptr)
							break;

						{
							std::lock_guard<std::mutex> lock(queueMutex);
							rayQueue.push(ray);
						}
						c++;
						if (c > 100000)
						{
							std::cout << "Size of queue: " << rayQueue.size() << std::endl;
							c = 0;
						}
					cv.notify_one();
				} while (true);
			}
			for (int i = 0; i < S.getNumberOfThreads(); i++)
			{
				{
					std::lock_guard<std::mutex> lock(queueMutex);
					rayQueue.push(nullptr);
				}

				cv.notify_one();
			}
			for (auto& t : threads) // Wait for all threads to finish
				t.join();
		}
	}
}