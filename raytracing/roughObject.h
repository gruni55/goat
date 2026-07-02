#pragma once
#include "objectshape.h"
#include "surface.h"
#include "ellipsoid.h"
#include "box.h"
#include <random>
#include <fstream>
#include "roughObjectInterface.h"

namespace GOAT 
	{
		namespace raytracing 
		{
          

			

            /**
			 @brief roughObject is a template class that adds roughness to an object shape.
			 This class applies different types of roughness to object shapes. 
             The following distributions are available: Gaussian, normal, and cosine. 
             In each case, the surface normal is varied according to the selected distribution, 
             which results in a variation in the direction of reflection or transmission. 
             The angle range can be limited by specifying theta_max (values between 0 and pi/2). 
             For the normal distribution, the standard deviation sigma is also specified in radians. 
            */

            template <class T>
            class roughObject : public T, public roughInterface
            {
            public:
                roughObject(const T& obj, double sigma)
                    : T(obj)
                {
                    setSigma(sigma);
                    this->rough = true;
                    init();
                }
                
                roughObject(T* obj, double sigma)
                    : T(*obj)
                {
                    setSigma(sigma);
                    this->rough = true;
                    init();
                }

                roughObject(const roughObject& F)
                    : T(F),
                    uniformDist(F.uniformDist),
					normalDist(F.normalDist),
					phiDist(F.phiDist),
                    sigma(F.sigma)
                {
                    this->type = OBJECTSHAPE_ROUGH_OBJECT;
                    this->rough = true;
                }

            
                /**
				* @brief gives the normal vector of the rough object at a given point P.
                * The surface normal is determined randomly, depending on the parameters sigma 
                * (standard deviation for a normal distribution) and thetaMax (maximum scattering angle).
				* @param P The point on the surface of the object where the normal vector is to be calculated.
                */
                maths::Vector<double> norm(const maths::Vector<double>& P) override
                {
					switch (scatteringType)
					{
					case ScatteringType::Gaussian:
						return normGaussian(P);
					case ScatteringType::UniformCone:
						return normUniformCone(P);
					case ScatteringType::CosineCone:
						return normCosineCone(P);
					
					}
				}

                /** @brief Sets the standard deviation for the normal distribution.
                 *  @param s The standard deviation in radians.
                 */
                void setSigma(double s)
                {
                    sigma = std::max(0.0, s);
                }

                /** @brief Sets the maximum scattering angle.
                 *  @param t The maximum scattering angle in radians.
                 */
                void setThetaMax(double t)
                {
                    constexpr double eps = 1E-12;

                    if (t < eps)
                        thetaMax = 0.0;
                    else if (t > 0.5 * M_PI)
                        thetaMax = 0.5 * M_PI;
                    else
                        thetaMax = t;
                }

                
				/** @brief Returns the standard deviation for the normal distribution.
				 *  @return The standard deviation in radians.
				 */
				double getSigma() const override
				{
					return sigma;
				}

                /** @brief Returns the maximum scattering angle.
                 *  @return The maximum scattering angle in radians.
                 */
				double getThetaMax() const
				{
					return thetaMax;
				}
               
                /** @brief Sets the scattering type.
                 *  @param type The scattering type.
                 */
				void setScatteringType(ScatteringType type)
				{
					scatteringType = type;
				}

                /** @brief Returns the scattering type.
                 *  @return The scattering type.
                 */
				ScatteringType getScatteringType() const
				{
					return scatteringType;
				}


            private:    
                
                void init()
                {
                  uniformDist = std::uniform_real_distribution<double>(0.0, 1.0);
				  normalDist = std::normal_distribution<double>(0.0, 1.0);
                  phiDist = std::uniform_real_distribution<double>(0.0, 2.0 * M_PI);
                }


                maths::Vector<double> normGaussian(const maths::Vector<double>& P)
                {
                    maths::Vector<double> n = T::norm(P);

                    maths::Vector<double> t1 = (std::abs(n[2]) < 0.9) ? n % maths::ez : n % maths::ex;

                    t1 /= abs(t1);

                    maths::Vector<double> t2 = n % t1;

                    maths::Vector<double> nMicro;
                    double theta;

                    do
                    {
                        double sx = sigma * normalDist(rng);
                        double sy = sigma * normalDist(rng);

                        nMicro = n + sx * t1 + sy * t2;
                        nMicro /= abs(nMicro);

                        theta = std::acos(std::clamp(nMicro * n, -1.0, 1.0));
                    } 
                    while (theta > thetaMax);

                    return nMicro;
                }


                maths::Vector<double> normUniformCone(const maths::Vector<double>& P)
                {
                    maths::Vector<double> n = T::norm(P);

                    if (thetaMax <= 0.0) return n;

                    maths::Vector<double> t1 = (std::abs(n[2]) < 0.9) ? n % maths::ez : n % maths::ex;
                    t1 /= abs(t1);
                    maths::Vector<double> t2 = n % t1;

                    double u = uniformDist(rng);

                    double cosTheta = std::cos(thetaMax) + u * (1.0 - std::cos(thetaMax));
                    double sinTheta = std::sqrt(1.0 - cosTheta * cosTheta);
                    double phi = 2.0 * M_PI * uniformDist(rng);

                    maths::Vector<double> nMicro = sinTheta * std::cos(phi) * t1 + sinTheta * std::sin(phi) * t2 + cosTheta * n;
                    nMicro /= abs(nMicro);
                    return nMicro;
                }

                maths::Vector<double> normCosineCone(const maths::Vector<double>& P) 
                {
                    maths::Vector<double> n = T::norm(P);

                    if (thetaMax <= 0.0) return n;

                    maths::Vector<double> t1 =  (std::abs(n[2]) < 0.9) ? n % maths::ez : n % maths::ex;
                    t1 /= abs(t1);

                    maths::Vector<double> t2 = n % t1;

                    double u = uniformDist(rng);
                    double cosThetaMax = std::cos(thetaMax);
                    double cosTheta = std::sqrt(1.0 - u * (1.0 - cosThetaMax * cosThetaMax));

                    double sinTheta = std::sqrt(1.0 - cosTheta * cosTheta);

                    double phi = 2.0 * M_PI * uniformDist(rng);

                    maths::Vector<double> nMicro = sinTheta * std::cos(phi) * t1 + sinTheta * std::sin(phi) * t2 + cosTheta * n;
                    nMicro /= abs(nMicro);

                    return nMicro;
                }
                
               

                

				std::ofstream logFile;
                std::uniform_real_distribution<double> uniformDist;
				std::normal_distribution<double> normalDist;

                std::uniform_real_distribution<double>  phiDist;
                double sigma = 0.0;
                double thetaMax = M_PI / 2.0 ;
                inline static thread_local std::mt19937 rng{ std::random_device{}() };
				ScatteringType scatteringType = ScatteringType::Gaussian;
            };


            ObjectShape* makeRough(ObjectShape* obj, double sigma);
            
		}
}
