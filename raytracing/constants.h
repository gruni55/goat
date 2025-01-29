/** @file */
#pragma once
namespace GOAT
{
	namespace raytracing
	{
// #ifndef c_light
constexpr double c_light=299792458.0;  ///< speed of light in vacuum in m/s
// #endif

#define Z0 376.730313461   ///< impedance of the free space \f$Z_0=\sqrt{\frac{\mu_0}{\epsilon_0}}=\mu_0\cdot c_{light}\f$ in Ohms
// #define mu0 4.0*M_PI*1E-7   ///< vacuum permeability in \f$F\cdot A^{-2}\f$
constexpr double mu0=1.25663706127E-6; ///< vacuum permeability in \f$F\cdot A^{-2}\f$
constexpr double eps0=8.854187817E-12; ///< vacuum permittvity in \f$F\cdot m^{-1}\f$
#define Planck_h 6.62606896E-34 ///< Planck constant in Js
#define Planck_hquer 1.054571628E-34 ///<  reduced Planck constant in Js
        constexpr double C_LIGHT = 299792458; ///< speed of light in \f$ m\cdot s^{-1}\f$
        constexpr double C_LIGHT_MU_FS = 0.299792458; ///< speed of light in \f$ m*fs^{-1}\f$
        constexpr double C_LIGHT_MU = 299792458E+6; ///< speed of light in \f$\mu m*s^{-1} \f$

	}
}
