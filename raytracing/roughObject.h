#pragma once
#include "objectshape.h"
#include "surface.h"
#include "ellipsoid.h"
#include "box.h"
#include <random>
#include <fstream>
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
                   // dist = std::uniform_real_distribution<double>(-1.0, 1.0);
					dist = std::uniform_real_distribution<double>(0.0, 1.0);
                    phiDist = std::uniform_real_distribution<double>(0.0, 2.0 * M_PI);
				//	logFile.open("roughObject.log");
                }

                roughObject(T* obj, double sigma)
                    : T(*obj)
                {
                    setSigma(sigma);
                    this->rough = true;
                  //  logFile.open("roughObject.log");
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
                    maths::Vector<double> t1 = fabs(n[2]) < 0.9 ? n % maths::ez : n % maths::ex;
                    t1 /= abs(t1);
                    maths::Vector<double> t2 = n % t1;

                    double m = 1.0;
				
					double u = dist(rng);

                    double cosAlphaMax = std::cos(sigma);
                    double p = 1.0 / (m + 1.0);

/*                    double cosAlpha = std::pow(1.0 - u * (1.0 - std::pow(cosAlphaMax, m + 1.0)), p);

                    double sinAlpha = std::sqrt(1.0 - cosAlpha * cosAlpha);*/
                //    double theta = asin(u);
                    double theta=thetaDist(rng);
                    double sinTheta = sin(theta);
                    double cosTheta = cos(theta);
                    double phi = phiDist(rng);

                    maths::Vector<double> nMicro =
                        sinTheta * std::cos(phi) * t1
                        + sinTheta * std::sin(phi) * t2
                        + cosTheta * n;
                    nMicro/=abs(nMicro);
					return nMicro;
                }
               /* maths::Vector<double> norm(const maths::Vector<double>& P) override
                {
					maths::Vector<double> n = T::norm(P);
                    maths::Vector<double> t1 = n[2] < 0.9 ? n % maths::ez : n % maths::ex;
					t1 /= abs(t1);
					maths::Vector<double> t2 = n % t1;
					double theta = thetaDist(rng);
					double phi = dist(rng) * M_PI;
					double cosTheta = cos(theta);
					double sinTheta = sin(theta);
					double cosPhi = cos(phi);
					double sinPhi = sin(phi);
					maths::Vector<double> newNormal = cosTheta * n + sinTheta * (cosPhi * t1 + sinPhi * t2);
                 //   logFile << P << "\t" << n << "\t" << newNormal << "\t" << theta << "\t" << phi << std::endl;
                    return newNormal;
                } */

                void setSigma(double s)
                {
		    phiDist = std::uniform_real_distribution<double>(0.0, 2.0 * M_PI);

                    sigma = s;
                    if (s>0)
                    thetaDist = std::normal_distribution<double>(0.0, sigma);
                }

				double getSigma() const override
				{
					return sigma;
				}

            private:
				std::ofstream logFile;
                std::uniform_real_distribution<double> dist;
                std::uniform_real_distribution<double>  phiDist;
                double sigma = 0.0;
                std::normal_distribution<double> thetaDist;

                inline static thread_local std::mt19937 rng{ std::random_device{}() };
            };
            ObjectShape* makeRough(ObjectShape* obj, double sigma);
            
		}
}
