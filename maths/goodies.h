
#ifndef GOODIES_H
#define GOODIES_H
#include <complex>
namespace GOAT
{
	namespace maths
	{

/*#ifndef min
#define min(a,b) a<b ? a : b
#endif

#ifndef max
#define max(a,b) a>b ? a : b
#endif*/
 
		bool testnan(double x);
		double sqr(double x);
		double abs(std::complex<double> x);
	}
}
#endif
