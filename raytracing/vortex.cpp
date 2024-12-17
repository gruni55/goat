#include "cylinder.h"
#include "vortex.h"
namespace GOAT
{
	namespace raytracing
	{
		VortexPlate::VortexPlate() : ObjectShape()
		{
			type = OBJECTSHAPE_CYLINDER;
		}

		VortexPlate::VortexPlate(ObjectShape& os) : ObjectShape (os)
		{
			type = OBJECTSHAPE_CYLINDER;
		}

		VortexPlate::VortexPlate(VortexPlate& c) : ObjectShape(c)
		{
			type = OBJECTSHAPE_CYLINDER;
			r = c.r;
			h = c.h;
			initQuad();
		}

		VortexPlate::VortexPlate(const maths::Vector<double>& P, double r, double h, double dh, int m, std::complex<double> n, double r0, const maths::Matrix<std::complex<double>> alpha, const maths::Vector<double>& Ex, const maths::Vector<double>& Ey, const maths::Vector<double>& Ez) 
			: ObjectShape(P, n, alpha, Ex, Ey, Ez, OBJECTSHAPE_CYLINDER)
		{
			this->r = r;
			this->h = h;
			this->dh = dh;
			this->m = m;
			initQuad();
		}



		void VortexPlate::setRadius(double r)
		{
			this->r = r;
		}

		void VortexPlate::setHeight(double h)
		{
			this->h = h;
		}

		void VortexPlate::setr0(double r0)
		{
			this->r0 = r0;
			initQuad();				
		}

		void VortexPlate::setm(int m)
		{
			this->m = m;
			initQuad();
		}

		void VortexPlate::binWrite(std::ofstream& os)
		{
			P.binWrite(os);
			H.binWrite(os);
			R.binWrite(os);
			os.write((char*)&n, (char)sizeof(n));
			alpha.binWrite(os);
			pul.binWrite(os);
			por.binWrite(os);
			for (int i = 0; i < 3; i++)
				e[i].binWrite(os);
			os.write((char*)&Ealpha, (char)sizeof(Ealpha));
			os.write((char*)&Ebeta, (char)sizeof(Ebeta));
			os.write((char*)&Egamma, (char)sizeof(Egamma));
			os.write((char*)&sf, (char)sizeof(sf));
			os.write((char*)&r0, (char)sizeof(r0));
			os.write((char*)&r, (char)sizeof(r));
			os.write((char*)&h, (char)sizeof(h));
		}

		void VortexPlate::binRead(std::ifstream& is)
		{
			type = OBJECTSHAPE_VORTEX_PLATE;
			P.binRead(is);
			H.binRead(is);
			R.binRead(is);
			is.read((char*)&n, (char)sizeof(n));
			alpha.binRead(is);
			pul.binRead(is);
			por.binRead(is);
			for (int i = 0; i < 3; i++)
				e[i].binRead(is);
			is.read((char*)&Ealpha, (char)sizeof(Ealpha));
			is.read((char*)&Ebeta, (char)sizeof(Ebeta));
			is.read((char*)&Egamma, (char)sizeof(Egamma));
			is.read((char*)&sf, (char)sizeof(sf));
			is.read((char*)&r0, (char)sizeof(r0));
			is.read((char*)&r, (char)sizeof(r));
			is.read((char*)&h, (char)sizeof(h));
		}



		maths::Vector<double> VortexPlate::norm(const maths::Vector<double>& ps)
		{
			maths::Vector<double> p = H * (ps - P);
			if (p[2] < EPS) return -R * maths::ez;
			// if (h - p[2] < EPS) return R * maths::ez;
			if (h - p[2] < dh + EPS)
			{
				double phi = atan2(p[1], p[0]);
				double r = sqrt(p[0] * p[0] + p[1] * p[1]);
				maths::Vector<double> n(-dh * (double)m / (2.0 * M_PI) * sin(phi), dh * (double)m / (2.0 * M_PI) * cos(phi), r);
				n /= -abs(n);
				return R * n;
			}
			return R * maths::Vector<double>(p[0], p[1], 0);
		}

		bool VortexPlate::isInside(const maths::Vector<double>& p)
		{
			// to be implemented !!!!!!!!!!!!!
			return false;
		}

		double VortexPlate::volume()
		{
			return M_PI * r * r * h;
		}



		void VortexPlate::initQuad()
			// must be implemented in a better way, in the moment only valid for the unrotated cylinder
		{
			pul = maths::Vector<double>(-r, -r, 0) + P;
			por = maths::Vector<double>(r, r, h+dh) + P;
		}

		maths::Matrix<double> VortexPlate::computeInertia()
		{
			GOAT::maths::Matrix<double> I;
			double r2 = r * r;
			double h2 = h * h;
			I(0, 0) = 1.0 / 12.0 * (3.0 * r2 + h2);
			I(1, 1) = 1.0 / 12.0 * (3.0 * r2 + h2);
			I(2, 2) = 1.0 / 12.0 * (6.0 * r2);
			return maths::Matrix<double>();
		}

		bool VortexPlate::next(const maths::Vector<double>& Ps, const maths::Vector<double>& K, maths::Vector<double>& pout)
		{
			double A, B, C, D;
			double l1, l2, l3, l;
			GOAT::maths::Vector<double> n, k, p = Ps - P;

			// First, the lateral surface
			p = H * p;
			k = H * K;

			A = k[0] * k[0] + k[1] * k[1];
			B = p[0] * k[0] + p[1] * k[1];			
			C = p[0] * p[0] + p[1] * p[1] - r * r;

			D = B * B - A * C;
			if (D >= 0)
			{
				if (D == 0) if (B == 0) l1 = -1; else l1 = -1 / B;
				else
				{
					double sd = sqrt(D);
					double la, lb;
					la = (-B + sd) / A;
					lb = (-B - sd) / A;
					if ((la < lb) && (la > EPS) || (lb < EPS)) l1 = la;
					else l1 = lb;					
				}
				pout = p + l1 * k;
				if ((pout[2] < 0) || (pout[2] > h)) l1 = -1;
			}
			else l1 = -1;
			
			// now bottom face
			l2 = -p[2] / k[2];
			if (l2 > EPS)
			{
				pout = p + l2 * k;
				if (pout[0] * pout[0] + pout[1] * pout[1] > r * r) l2 = -1;
			}

			// top face
			l3 = (h - p[2]) / k[2];
			if (l3 > EPS)
			{
				pout = p + l3 * k;
				if (pout[0] * pout[0] + pout[1] * pout[1] > r * r) l3 = -1;
				else
				{
					double phi = atan2(pout[1],pout[0]);
					double z = h + dh * (m * phi / (2.0 * M_PI) - floor(m * phi / (2.0 * M_PI)));
					l3 = (z - p[2]) / k[2];
				}
			}

			if (((l1 > l2) && (l2>EPS)) || (l1 < 0)) l = l2;
			else l = l1;

			if ((l < EPS) || ((l > l3) && (l3 > EPS))) l = l3;
			if (l <= EPS) return false;

			pout = Ps + l * K;
			return true;
		}

		VortexPlate& VortexPlate::operator=(VortexPlate& f)
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

		VortexPlate& VortexPlate::operator=(VortexPlate f)
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
		void VortexPlate::scale(double sf)
		{
			r = sf * r;
			h = sf * h;
			initQuad();
		}
	}
}
