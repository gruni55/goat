#include "goodies.h"
#include <iostream>
#include <math.h>
#include <string>
#include <sstream>

using namespace std;
namespace GOAT
{
	namespace maths
	{
		bool testnan(double x)
		{
			stringstream ss;
			string S;
			ss << x;
			ss >> S;
			return S == "nan";
		}

		double sqr(double x) { return x * x; }

		double abs(complex<double> x)
		{
			return real(sqrt(x * conj(x)));
		}

		std::vector<double> movingAvg(std::vector<double> d, int n)
		{			
			
			std::size_t size = d.size();			
			if (n % 2 == 0) n = n + 1;
			std::vector<double> avg(n, 0);
			std::vector<double> erg;
			int mid = (n - 1) / 2;
			int hi;

			for (std::size_t i=0; i<size; i++)
			{
				for (int l = 0; l < n; l++)
				{
					hi = i + (l - mid);
					 avg[l] = (hi<0) || (hi>=size) ? 0 : d[hi];					 
				}
				erg.push_back(mean(avg));
			}
			return erg;
		}

		double mean(std::vector<double> d)
		{
			double erg = 0;
			int n = d.size();
			for (int i = 0; i < n; i++)
				erg += d[i];
			erg /= n;
			return erg;
		}

		double findmax(std::vector<double> d, std::size_t& index)
		{
			std::size_t size = d.size();
			double erg = d[0];
			index = 0;
			for (std::size_t i=0; i<size; i++)
			 if (d[i]>erg) 
			 {
				 index = i;
				 erg = d[i];
			 }
			return erg;
		}

		std::size_t FWHM(std::vector<double> d, std::size_t index)
		{
			std::size_t leftI = index;
			std::size_t rightI = index;
			
			bool found = false;
			double maxValue = d[index];
			std::size_t size = d.size();


			for (std::size_t i=index; (i>0) && !found; i--)
				if (d[i] / maxValue <= 0.5) { found = true; leftI = i; }

			found = false;
			for (std::size_t i = index; (i < (size - 1)) && !found; i++)
				if (d[i] / maxValue <= 0.5) { found = true; rightI = i; }

			std::size_t erg = rightI - leftI;
			return erg;
		}

#include <iostream>

		double newton_root(std::function<double(double)> f, std::function<double(double)> df, double x0, double eps)
		{
			double xold, xnew;
			xold = x0;
			double dx;
			do
			{
				std::cout << "xold=" << xold << std::endl;
				xnew = xold - f(xold) / df(xold);
				dx = fabs(xold - xnew);
				xold = xnew;
			} while (dx > eps);
			return xnew;
		}

	}
}

