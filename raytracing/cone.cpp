#include "cone.h"
namespace GOAT
{
	namespace raytracing
	{
		Cone::Cone(maths::Vector<double> Pos, double radius, double height, std::complex<double> n, double r0, const maths::Matrix<std::complex<double>> alpha, const maths::Vector<double>& Ex, const maths::Vector<double>& Ey, const maths::Vector<double>& Ez) 
			: ObjectShape (Pos,n,alpha,Ex,Ey,Ez, OBJECTSHAPE_CONE)
		{
			this->height = height;
			this->radius = radius;
			setr0(r0);
			init();
		}

		Cone::~Cone()
		{
		}

		void Cone::init()
		{
			v = -maths::ez;
			normv = -maths::ez;
			V = P - height * v;
			coneAngle = atan(radius / height);
			cosCA = cos(coneAngle);
			tan2CA = tan(coneAngle);
			tan2CA *= tan2CA;
			sideLen = height / cosCA;
		}

		bool Cone::next(const maths::Vector<double>& ps, const maths::Vector<double>& ks, maths::Vector<double>& pout)
		{
			// test with lateral surface 
			maths::Vector<double> poutC;

			maths::Vector<double> p = H * ps;
			maths::Vector<double> k = H * ks;
			double lC = nextCone(p, k, poutC);

			// std::cout << poutC << std::endl;

			// test with bottom surface 
			double lB = (P - p) * v / (k * v);
			double l;

			if (((lC < lB) && (lC > 0)) || (lB<=0))
			{
				if (lC <= 0) return false;
				pout = lC * ks + ps;
				return true;
			}
			
			pout = lB * ks + ps;
			if (abs(pout - P) > radius) return false;
			return true;
		}


		double Cone::nextCone(const maths::Vector<double>& p, const maths::Vector<double>& k, maths::Vector<double>& pout)
		{
			maths::Vector<double> n = GOAT::maths::norm(k % (V - p));
			maths::Vector<double> DV = V - p;
			double nv = n * v;

			if (nv > 0) { n = -n; nv = -nv; }
			double cos2theta = 1.0 - (nv * nv);
			double costheta = sqrt(cos2theta);

			if (costheta < cosCA) return -1;

			maths::Vector<double> u, w;
			u = GOAT::maths::norm(v % n);
			w = GOAT::maths::norm(u % v);

			/*double tan2theta = nv * nv / cos2theta;
			double tantheta = sqrt(tan2theta);*/
			double tantheta = std::acos(costheta);
			double tan2theta = tantheta * tantheta;
			maths::Vector<double> D = sqrt(tan2CA - tan2theta) * u;
			maths::Vector<double> d1 = v + tantheta * w + D;
			maths::Vector<double> d2 = v + tantheta * w - D;
			double r = ((V - p) % d1) * (k % d1) / abs2(k % d1);
			double s = ((V - p) % d2) * (k % d2) / abs2(k % d2);
	
			if ((r < s) && (r > 0))
			{
				pout = p + r * k;
				if (abs(pout - V) <= sideLen) return r;
				return -1;
			}
			
			if (s > 0) 
			{ 
				pout = p + s * k; 
				maths::Vector<double> h = pout - V;
				if ((abs(h) <= sideLen) && (h*v>=0)) return s; 
				return -1;
			}			
			return -1;
		}

		maths::Vector<double> Cone::norm(const maths::Vector<double>& p)
		{
			maths::Vector<double> Ps = H * p;	
			double dv = (Ps - P) * v;
			if (fabs((Ps - P) * v) > 0)
			{
				maths::Vector<double> res, h, VPs;
				VPs = Ps - V;
				h = VPs % v;
				res = VPs % h;
				res /= abs(res);
				return R * res;
			}
			return normv;
		}

		bool Cone::isInside(const maths::Vector<double>& p)
		{
			maths::Vector<double> Ps = H * p;
			maths::Vector<double> Vp = Ps - V;
			double h = Vp * normv;
			if ((h < 0) || (h > height)) return false;
			if (h / abs(Vp) < cosCA) return false;

			return true;
		}

		double Cone::volume()
		{
			return M_PI * radius * radius * height / 3.0;
		}

		void Cone::initQuad()
		{
			pul = maths::Vector<double>(P[0] - radius, P[1] - radius, P[2]);
			por = maths::Vector<double>(P[0] + radius, P[1] + radius, P[2] + height);
		}

		void Cone::setPos(maths::Vector<double> r)
		{
			P = r;
			init();
			initQuad();
		}

		void Cone::setPos(double x, double y, double z)
		{
			maths::Vector<double> r(x, y, z);
			setPos(r);
		}

		void Cone::setr0(double r0)
		{
			initQuad();
			this->r0 = r0;
		}

		maths::Vector<double> Cone::calcCoM()
		{			
			return maths::Vector<double>(P[0],P[1],3.0/4.0*height);
		}
		
		void Cone::binWrite(std::ofstream& os)
		{
			P.binWrite(os);
			H.binWrite(os);
			R.binWrite(os);
			os.write((char*)&height, (char)sizeof(height));
			os.write((char*)&radius, (char)sizeof(radius));

			
		}
		void Cone::binRead(std::ifstream& os)
		{
			P.binRead(os);
			H.binRead(os);
			R.binRead(os);
			os.read((char*)&height, (char)sizeof(height));
			os.read((char*)&radius, (char)sizeof(radius));
		}
		void Cone::scale(double sf)
		{
			height = height / this->sf * sf;
			radius = radius / this->sf * sf;
			v = v / this->sf * sf;
			init();
		}
                 
               void Cone::setConeAngle (double coneAngle)
               {
                height=radius * tan (coneAngle);
                init();
                initQuad();
               }
                
               double Cone::getConeAngle ()
               {
                return coneAngle;
               }
            
	}
}
