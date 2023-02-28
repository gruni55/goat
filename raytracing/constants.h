#pragma once
namespace GOAT
{
	namespace raytracing
	{
#ifndef c_light
#define c_light 299792458.0   // Vakuumlichtgeschwindigkeit in m/s
#endif

#define Z0 376.730313461    // Impedanz des Vakuums Z0=sqrt (mu_0/epsilon_0)=mu_0*c_light in Ohm 
#define mu0 4.0*M_PI*1E-7   // µ0 in N/A^2
#define eps0 8.854187817E-12 // Epsilon_0
#define Planck_h 6.62606896E-34 // Plancksches Wirkungsquantum in Js
#define Planck_hquer 1.054571628E-34 // h/2PI in Js 
		constexpr double C_LIGHT = 299792458; ///< speed of light in m*s^-1
		constexpr double C_LIGHT_MU = 299792458E+6; ///< speed of light in mum*s^-1
	}
}