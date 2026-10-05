#pragma once
#include "propagator.h"
#include <memory>
namespace GOAT
{
	namespace raytracing
	{
		/**
		* @brief This class implements the angular spectrum method for propagating optical fields.
		* \note Currently, only one source is included in the calculation. Due to the calculation method, the resolution of the source and destination must be the same, and both must be arranged parallel to each other. The resolution and position of the destination are adjusted automatically!
		* As sour
		*/
		class AngularSpectrum : public Propagator
		{
			public:
			/**
			* @brief constructor with rectangular grid definition
			* /param wvl wavelength
			* /param P position of the propagator plane (center)
			* /param e1 first direction vector of the plane (length of the plane in this direction is determined by the length of e1)
			* /param e2 second direction vector of the plane (length of the plane in this direction is determined by the length of e2)
			* /param n1 number of points in the first direction
			* /param n2 number of points in the second direction
			*/
			AngularSpectrum(double wvl, maths::Vector<double> P, maths::Vector<double> e1, maths::Vector<double> e2, int n1, int n2);
			/**
			* @brief constructor for square grid definition
			* /param wvl wavelength
			* /param P position of the propagator plane (center)
			* /param n normal vector of the plane
			* /param d width of the plane in both directions
			* /param N number of points in one direction (so the array is N x N)
			*/

			AngularSpectrum(const AngularSpectrum&) = delete;
			AngularSpectrum& operator=(const AngularSpectrum&) = delete;

			AngularSpectrum(AngularSpectrum&&) noexcept = default;
			AngularSpectrum& operator=(AngularSpectrum&&) noexcept = default;

			AngularSpectrum(double wvl, maths::Vector<double>P, maths::Vector<double> n, double d, int N);
			/**
			* @brief destructor
			*/
			~AngularSpectrum();
			
			/**
			* @brief Makes a calculation. clear determines if the content of the propagator should be cleared before the calculation.
			*/
			void calc(bool clear = true);
			/**
			* @brief Propagates the field from one source detector to the propagator plane. 
			* clear determines if the content of the propagator should be cleared before the calculation.
			*  /param det source detector
			* /param clear if true, the content of the propagator will be cleared before the calculation. If false, the result of the calculation will be added to the existing content.
			*/
			void calcOne(DetectorPlane* det, bool clear);
			

		private:
			double dx = 0.0;
			double dy = 0.0;

			
			void adaptTargetGridToSource(DetectorPlane* src);
			bool computeRoiOffset(DetectorPlane* src, int& ix0, int& iy0) const;
			
			struct Impl;              
			std::unique_ptr<Impl> impl;
			std::vector<std::vector<maths::Vector<std::complex<double>>>> uD; ///< this is the field on the propagator plane before smoothing	
			bool isSmoothed = false; ///< this is true, if the field on the propagator plane has been smoothed (with a Gaussian filter) to avoid aliasing effects. It is false, if the field on the propagator plane is not smoothed (so it is the raw result of the calculation).
		};
	}
}
