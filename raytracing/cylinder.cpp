#include "cylinder.h"
namespace GOAT
{
	namespace raytracing
	{
		Cylinder::Cylinder() : ObjectShape()
		{
			type = OBJECTSHAPE_CYLINDER;
		}

		Cylinder::Cylinder(ObjectShape& os) : ObjectShape (os)
		{
			type = OBJECTSHAPE_CYLINDER;
		}

		Cylinder::Cylinder(Cylinder& c) : ObjectShape(c)
		{
			type = OBJECTSHAPE_CYLINDER;
			r = c.r;
			h = c.h;
			initQuad();
		}

		Cylinder::Cylinder(const maths::Vector<double>& P, double r, double h, std::complex<double> n, double r0, const maths::Matrix<std::complex<double>> alpha, const maths::Vector<double>& Ex, const maths::Vector<double>& Ey, const maths::Vector<double>& Ez) 
			: ObjectShape(P, n, alpha, Ex, Ey, Ez, OBJECTSHAPE_CYLINDER)
		{
			this->r = r;
			this->h = h;
			initQuad();
		}



		void Cylinder::setRadius(double r)
		{
			this->r = r;
		}

		void Cylinder::setHeight(double h)
		{
			this->h = h;
		}

		maths::Vector<double> Cylinder::norm(const maths::Vector<double>& ps)
		{
			maths::Vector<double> p = H * (ps - P);
			if (p[2] < EPS) return -R * maths::ez;
			if (h - p[2] < EPS) return R * maths::ez;
			return R * maths::Vector<double>(p[0], p[1], 0);
		}

		bool Cylinder::isInside(const maths::Vector<double>& p)
		{
			// to be implemented !!!!!!!!!!!!!
			return false;
		}

		double Cylinder::volume()
		{
			return M_PI * r * r * h;
		}



		void Cylinder::initQuad()
			// must be implemented in a better way, in the moment only valid for the unrotated cylinder
		{
			pul = maths::Vector<double>(-r, -r, 0) + P;
			por = maths::Vector<double>(r, r, h) + P;
		}

		maths::Matrix<double> Cylinder::computeInertia()
		{
			GOAT::maths::Matrix<double> I;
			double r2 = r * r;
			double h2 = h * h;
			I(0, 0) = 1.0 / 12.0 * (3.0 * r2 + h2);
			I(1, 1) = 1.0 / 12.0 * (3.0 * r2 + h2);
			I(2, 2) = 1.0 / 12.0 * (6.0 * r2);
			return maths::Matrix<double>();
		}

		bool Cylinder::next(const maths::Vector<double>& Ps, const maths::Vector<double>& K, maths::Vector<double>& pout)
		{
			double A, B, C, D;
			double l1, l2, l3, l;
			GOAT::maths::Vector<double> n, k, p = Ps - P;

			// First, the lateral surface
			p = H * p;
			k = H * K;

			A = p[0] * k[0] + p[1] * k[1];
			B = k[0] * k[0] + k[1] * k[1];
			C = p[0] * p[0] + p[1] * p[1] - r * r;

			D = A * A - B * C;
			if (D >= 0)
			{
				if (D == 0) l1 = -A / B;
				else
				{
					double sd = sqrt(D);
					double la, lb;
					la = (-A + sd) / B;
					lb = (-A - sd) / B;
					if (la <= 0) l1 = lb;
					else l1 = la;
				}
			}
			else l1 = -1;
			
			// now bottom face
			l2 = -p[2] / k[2];

			// top face
			l3 = (h - p[2]) / k[2];

			if (((l1 > l2) && (l2>0)) || (l1 < 0)) l = l2;
			else l = l1;

			if ((l < 0) || ((l > l3) && (l3 > 0))) l = l3;
			if (l <= 0) return false;

			pout = Ps + l * K;
			return true;
		}

		Cylinder& Cylinder::operator=(Cylinder& f)
		{
			if (this == &f) return *this;
			H = f.H;
			P = f.P;
			R = f.R;
			alpha = f.alpha;
			for (int i = 0; i < 3; i++)
				e[i] = f.e[i];
			n = f.n;
			por = f.por;
			pul = f.pul;
			r = f.r;
			h = f.h;
			r0 = f.r0;
			type = f.type;
			Ealpha = f.Ealpha;
			Ebeta = f.Ebeta;
			Egamma = f.Egamma;
			return *this;
		}

		Cylinder& Cylinder::operator=(Cylinder f)
		{
			if (this == &f) return *this;
			H = f.H;
			P = f.P;
			R = f.R;
			alpha = f.alpha;
			for (int i = 0; i < 3; i++)
				e[i] = f.e[i];
			n = f.n;
			por = f.por;
			pul = f.pul;
			r = f.r;
			h = f.h;
			r0 = f.r0;
			type = f.type;
			Ealpha = f.Ealpha;
			Ebeta = f.Ebeta;
			Egamma = f.Egamma;
			return *this;
		}
		void Cylinder::scale(double sf)
		{
			r = sf * r;
			h = sf * h;
			initQuad();
		}
	}
}
