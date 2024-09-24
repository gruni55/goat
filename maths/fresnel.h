#ifndef FRESNEL_H
#define FRESNEL_H
#include <complex>

#include "vector.h"
namespace GOAT
{
	namespace maths
	{
#ifndef SENKRECHT
#define SENKRECHT     1    ///< used to determine the direction of the polarization 
#define PERPENDICULAR 1    ///< used to determine the direction of the polarization  
#endif

#ifndef PARALLEL
#define PARALLEL   0       ///< used to determine the direction of the polarization 
#endif 

		double abs(std::complex<double>  x); ///< absolute value of the complex valued variable x
		double sqr(double x); ///< square of x

		/**
		* @brief Calculates the Fresnel coefficient for transmission 
		* This function calculates the amplitude Fresnel coefficient for the transmitted part of the ray
		* @param pol Direction of the polarization (PARALLEL or SENKRECHT/PERPENDICULAR
		* @param k Directional vector for the incident ray
		* @param n Surface normal
		* @param n1 Refractive index of the medium on the side of the incident ray
		* @param n2 Refractive index of the medium on the side of the outgoing(transmitted) ray
		*/
		std::complex<double> Fresnel_trans(int pol, Vector<double> k, Vector<double> n, std::complex<double>  n1, std::complex<double>  n2);
		/**
		* @brief Calculates the Fresnel coefficient for reflection
		* This function calculates the amplitude Fresnel coefficient for the reflected part of the ray
		* @param pol Direction of the polarization (PARALLEL or SENKRECHT/PERPENDICULAR
		* @param k Directional vector for the incident ray
		* @param n Surface normal
		* @param n1 Refractive index of the medium on the side of the incident ray
		* @param n2 Refractive index of the medium on the side of the outgoing(transmitted) ray
		*/
		std::complex<double> Fresnel_reflect(int pol, Vector<double> k, Vector<double> n, std::complex<double>  n1, std::complex<double>  n2);
		/**
		* @brief Calculates the Fresnel coefficient for reflection
		* This function calculates the amplitude Fresnel coefficient for the reflected part of the ray
		* @param pol Direction of the polarization (PARALLEL or SENKRECHT/PERPENDICULAR
		* @param alpha Angle of incidence (in radiants)		
		* @param n1 Refractive index of the medium on the side of the incident ray
		* @param n2 Refractive index of the medium on the side of the outgoing(transmitted) ray
		*/
		std::complex<double>  freflect(int pol, double alpha, std::complex<double>  n1, std::complex<double>  n2);
		/**
		* @brief Calculates the Fresnel coefficient for transmission
		* This function calculates the amplitude Fresnel coefficient for the transmitted part of the ray
		* @param pol Direction of the polarization (PARALLEL or SENKRECHT/PERPENDICULAR
		* @param alpha Angle of incidence (in radiants)
		* @param n1 Refractive index of the medium on the side of the incident ray
		* @param n2 Refractive index of the medium on the side of the outgoing(transmitted) ray
		*/
		std::complex<double>  ftrans(int pol, double alpha, std::complex<double>  n1, std::complex<double>  n2);
	}
}
#endif
