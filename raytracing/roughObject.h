#pragma once
#include "objectshape.h"
#include "surface.h"
#include "ellipsoid.h"
#include "box.h"
#include <random>
namespace GOAT 
	{
		namespace raytracing 
		{
            class roughInterface
            {
            public:
                virtual double getSigma() const = 0;
                virtual void setSigma(double s) = 0;
            };

            template <class T>
            class roughObject : public T, public roughInterface
            {
            public:
                roughObject(const T& obj, double sigma)
                    : T(obj)
                {
                    setSigma(sigma);
                    this->rough = true;
                    dist = std::uniform_real_distribution<double>(-1.0, 1.0);
                }

                roughObject(T* obj, double sigma)
                    : T(*obj)
                {
                    setSigma(sigma);
                    this->rough = true;
                }

                roughObject(const roughObject& F)
                    : T(F),
                    dist(F.dist),
                    sigma(F.sigma),
                    thetaDist(F.thetaDist)
                {
                    this->type = OBJECTSHAPE_ROUGH_OBJECT;
                    this->rough = true;
                }

                maths::Vector<double> norm(const maths::Vector<double>& P) override
                {
                    
                    maths::Vector<double> n = T::norm(P);
                    if (sigma == 0.0) return n;
                    maths::Vector<double> a(dist(rng), dist(rng), dist(rng));

                    maths::Vector<double> t = a - n * (a * n);
                    t /= abs(t);

                    double theta = thetaDist(rng);
                    return n * cos(theta) + t * sin(theta);
                }

                void setSigma(double s)
                {
                    sigma = s;
                    if (s>0)
                    thetaDist = std::normal_distribution<double>(0.0, sigma);
                }

				double getSigma() const override
				{
					return sigma;
				}

            private:
                std::uniform_real_distribution<double> dist;
                double sigma = 0.0;
                std::normal_distribution<double> thetaDist;

                inline static thread_local std::mt19937 rng{ std::random_device{}() };
            };
            ObjectShape* makeRough(ObjectShape* obj, double sigma);
            
		}
}