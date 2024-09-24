
#ifndef GOODIES_H
#define GOODIES_H
#include <complex>
#include <vector>
#include <functional>
#include <limits>
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
 
		bool testnan(double x); ///< tests x on NaN
		double sqr(double x); ///< square of x (i.e. x^2)
		double abs(std::complex<double> x); ///< absolute value of x (sqrt (re(x)^2+im(x)+^2)

		/**
		* @brief calculates the moving average of the data 
		* Calculates the (unweighted) moving average over the data stored in d. The average will 
		* be calculated over n elements
		* \param d Vector with the data
		* \param n Number of elements over which the average is calculated
		* if n is even the n+1 elements were used for averaging
		*/
		std::vector<double> movingAvg(std::vector<double> d, int n);
		double mean(std::vector<double> d); ///< Calculates the mean of the data stored in d
		/**
		* @brief finds maximum in data
		* This function finds the maximal value inside a vector 
		* \param[in] d vector with the data
		* \return Returns the highest value in d
		* \param[out] Index of the maximum inside d 
		*/
		double findmax(std::vector<double> d, std::size_t &index);


		/**
		* @brief estimates the full width at half maximum of a peak around at the position inside a vector
		* The function gives back the full width at half maximum (FWHM) by finding the first point with less than (or equal to) the half of the maximal value left and right from position
		* index inside the vector d
		* \param d data
		* \param index position of the peak 
		* \return FWHM in indices
		*/
		std::size_t FWHM(std::vector<double> d, std::size_t index);
		double newton_root(std::function<double(double)> f, std::function<double(double)> df, double x0=0, double eps=std::numeric_limits<double>::min());		
	}
}
#endif
