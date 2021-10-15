#ifndef FRESNEL_H
#define FRESNEL_H
#include <complex>

#include "vector.h"

#ifndef SENKRECHT
#define SENKRECHT  1
#endif

#ifndef PARALLEL
#define PARALLEL   0
#endif 

double abs (std::complex<double>  x); ///< absolute value of the complex valued variable x
double sqr (double x); ///< square of x

std::complex<double> Fresnel_trans (int pol, Vector<double> k, Vector<double> n, std::complex<double>  n1, std::complex<double>  n2);
std::complex<double> Fresnel_reflect (int pol, Vector<double> k, Vector<double> n, std::complex<double>  n1, std::complex<double>  n2);
std::complex<double>  freflect (int pol, double alpha, std::complex<double>  n1, std::complex<double>  n2);
std::complex<double>  ftrans (int pol, double alpha, std::complex<double>  n1, std::complex<double>  n2);
#endif
