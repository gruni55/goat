#pragma once
#include "iray.h"
/** This class provides a ray which carries a special power*/
class Ray_pow :
	public IRay
{
public:
	Ray_pow();
	/**
	 * Constructor.
	 * 
	 * \param pow power of the ray (in W)
	 * \param p position of the ray
	 * \param Pol Direction vector for the electric field
	 * \param K Direction of the ray
	 * \param n0 refractive index of the surroundings
	 * \param r0 radius of the calculation sphere
	 * \param k0 wavenumber, i.e. \f$ k_0=\frac{2\pi}{\lambda_0}\f$ with the wavelength in the vacuum \f$\lambda_0\f$
	 * \param numObjs number of objects (optional, default numObjs=0)
	 * \param Objs list with objects (optional, default Objs=NULL)
	 */
	Ray_pow(double pow, const Vector<double> &p,
         const Vector<std::complex<double> > &Pol, const Vector<double> &K,
         std::complex<double>  n0, double r0, double k0,
         const int numObjs, ObjectShape **Einschluss);
	Ray_pow(const Ray_pow &r) : IRay(r)
	{
		this->Pow=r.Pow;
	}
	/**
	 * reflects the ray on a surface, returns the transmitted ray
	 * \param n surface normal
	 * \param n1 refractive index (incident side)
	 * \param n2 refractive index (side of the transmitted ray)
	 */
	Ray_pow reflect(Vector<double> n, std::complex<double> n1, std::complex<double> n2); 
	/**
	 * Reflects the ray and gives the transmitted ray back.
	 * 
	 * \param[out] tray transmitted ray
	 * \param[in] n surface normal
	 * \param[in] n1 refractive index (incident side)
	 * \param[in] n2 refractive index (transmitted side)
	 */
	void reflectRay(RayBase*& tray, Vector<double> n, std::complex<double> n1, std::complex<double> n2);
	/**
	 * refracts the ray with help of the Fresnel matrix (for internal use only).
	 * 
	 * \param FT Fresnel matrix 
	 * \param N surface normal
	 * \param n1 refractive index (incident side)
	 * \param n2 refractive index (transmitted side)
	 */
	void refract(Matrix<std::complex<double> > FT, Vector<double> N, std::complex<double> n1, std::complex<double> n2);
	~Ray_pow(void);
	double Pow;
	friend std::ostream& operator << (std::ostream &os,Ray_pow S);
};

