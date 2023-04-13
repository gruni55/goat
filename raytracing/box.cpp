#include "box.h"
#include <limits>

constexpr double FloatInf = std::numeric_limits<double>::infinity();


#ifndef EPS
#define EPS 1E-10
#endif
namespace GOAT 
  {
	namespace raytracing
	{
		Box::Box()
		{
			type = OBJECTSHAPE_BOX;
		}

		Box::Box(const ObjectShape& F) : ObjectShape(F)
		{
			type = OBJECTSHAPE_BOX;
			isOctree = ((Box)F).isOctree;
		}

		Box::Box(const Box& B) : ObjectShape(B)
		{
			d = B.d;
			bounds[0] = B.bounds[0];
			bounds[1] = B.bounds[1];
			calcDiag();
			type = OBJECTSHAPE_BOX;
			isOctree = B.isOctree;
		}

		Box::Box(const maths::Vector<double>& P,
			const maths::Vector<double>& d,
			std::complex<double> n,
			double r0,
			const maths::Matrix<std::complex<double> > alpha,
			const maths::Vector<double> Ex,
			const maths::Vector<double> Ey,
			const maths::Vector<double> Ez) : ObjectShape(P, n, alpha, Ex, Ey, Ez)
		{
		//	bounds[0] = P - d / 2.0;
		//	bounds[1] = P + d / 2.0;

			bounds[0] = -1.0 * d / 2.0;
			bounds[1] =  d / 2.0;
			this->d = d;
			this->r0 = r0;
			type = OBJECTSHAPE_BOX;
			calcDiag();
		}



		Box::~Box()
		{
		}

		void Box::binWrite(std::ofstream& os)
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
			d.binWrite(os);
			bounds[0].binWrite(os);
			bounds[1].binWrite(os);
		}

		void Box::binRead(std::ifstream& is)
		{
			type = OBJECTSHAPE_BOX;
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
			d.binRead(is);
			bounds[0].binRead(is);
			bounds[1].binRead(is);
			calcDiag();
		}


		void Box::scale(double sf)
		{
			double sfold = this->sf;
			this->sf = sf;
			d = d * sf / sfold;
			bounds[0] = P - d / 2.0;
			bounds[1] = P + d / 2.0;
			calcDiag();
			initQuad();
		}

		inline void swap(double& x, double& y)
		{
			double h = y;
			y = x;
			x = h;
		}

		double Box::distance(const maths::Vector<double>& ps, const maths::Vector<double>& K)
		{
			maths::Vector<double> p = ps - P;
			maths::Vector<double> k;
			p = H * p;
			k = H * K;
			double L, l[2], lm, lp;
			int found = 0;
			maths::Vector<double> h;
			maths::Vector<double> d2 = d / 2.0;
			// Erst mal mit den Seitenflächen x=const. anfangen 

			lp = (d[0] / 2.0 - p[0]) / k[0];
			lm = (-d[0] / 2.0 - p[0]) / k[0];
			if (lp > EPS)
			{
				h = p + lp * k;
				if ((fabs(h[1]) < d2[1]) && (fabs(h[2]) < d2[2])) { found++; l[found - 1] = lp; }
			}

			if (lm > EPS)
			{
				h = p + lm * k;
				if ((fabs(h[1]) < d2[1]) && (fabs(h[2]) < d2[2])) { found++; l[found - 1] = lm; }
			}

			if (found == 2) { L = (l[0] < l[1]) ? l[0] : l[1]; return L; }

			// jetzt das Ganze für y=const. ...
			lp = (d[1] / 2.0 - p[1]) / k[1];
			lm = (-d[1] / 2.0 - p[1]) / k[1];

			if (lp > EPS)
			{
				h = p + lp * k;
				if ((fabs(h[0]) < d2[0]) && (fabs(h[2]) < d2[2])) { found++; l[found - 1] = lp; }
			}

			if (lm > EPS)
			{
				h = p + lm * k;
				if ((fabs(h[0]) < d2[0]) && (fabs(h[2]) < d2[2])) { found++; l[found - 1] = lm; }
			}

			if (found == 2) { L = (l[0] < l[1]) ? l[0] : l[1]; return L; }

			// und schließlich für z=const.
			lp = (d[2] / 2.0 - p[2]) / k[2];
			lm = (-d[2] / 2.0 - p[2]) / k[2];
			if (lp > EPS)
			{
				h = p + lp * k;
				if ((fabs(h[1]) < d2[1]) && (fabs(h[0]) < d2[0])) { found++; l[found - 1] = lp; }
			}

			if (lm > EPS)
			{
				h = p + lm * k;
				if ((fabs(h[1]) < d2[1]) && (fabs(h[0]) < d2[0])) { found++; l[found - 1] = lm; }
			}

			if (found == 0) { return FloatInf; }
			if (found == 1) { return l[0]; }
			if (found == 2) { L = (l[0] < l[1]) ? l[0] : l[1]; return L; }
			return FloatInf;
		}

		double findmin(double t1, double t2)
		{
			double t;
			t = (t1 < t2) && (t1 > 0) ? t1 : t2;
			return (t < 0) ? INFINITY : t;
		}



		/*bool Box::next(const Vector<double> &ps, const Vector<double> &K,Vector<double> &pout, const int inside)
		{
		 Vector<double> p = ps - P;
		 Vector<double> k;
		 p = H*p;
		 k = H*K;
		 double L,l[2],lm, lp;
		 int found=0;
		 Vector<double> h;
		 Vector<double> d2=d/2.0;
		 // Erst mal mit den Seitenflächen x=const. anfangen

		  lp=(d[0]/2.0-p[0])/k[0];
		  lm=(-d[0]/2.0-p[0])/k[0];
		  if (lp>EPS)
					{
					  h=p+lp*k;
					  if ((fabs(h[1])<d2[1]) && (fabs(h[2])<d2[2])) { found++; l[found-1]=lp; }
					}

		  if (lm>EPS)
					{
					  h=p+lm*k;
					  if ((fabs(h[1])<d2[1]) && (fabs(h[2])<d2[2])) { found++; l[found-1]=lm; }
					}

		 if (found == 2)  {  L=(l[0]<l[1]) ? l[0] : l[1]; pout=ps+L*K; return true; }

		 // jetzt das Ganze für y=const. ...
		  lp=(d[1]/2.0-p[1])/k[1];
		  lm=(-d[1]/2.0-p[1])/k[1];

		  if (lp>EPS)
					{
					  h=p+lp*k;
					  if ((fabs(h[0])<d2[0]) && (fabs(h[2])<d2[2])) { found++; l[found-1]=lp; }
					}

		  if (lm>EPS)
					{
					  h=p+lm*k;
					  if ((fabs(h[0])<d2[0]) && (fabs(h[2])<d2[2])) { found++; l[found-1]=lm; }
					}

		 if (found == 2)  {  L=(l[0]<l[1]) ? l[0] : l[1]; pout=ps+L*K; return true; }

		// und schließlich für z=const.
		  lp=(d[2]/2.0-p[2])/k[2];
		  lm=(-d[2]/2.0-p[2])/k[2];
		  if (lp>EPS)
					{
					  h=p+lp*k;
					  if ((fabs(h[1])<d2[1]) && (fabs(h[0])<d2[0])) { found++; l[found-1]=lp; }
					}

		  if (lm>EPS)
					{
					  h=p+lm*k;
					  if ((fabs(h[1])<d2[1]) && (fabs(h[0])<d2[0])) { found++; l[found-1]=lm; }
					}

		 if (found == 0) { pout=ps; return false; }
		 if (found == 1) { pout=ps+l[0]*K; return true; }
		 if (found == 2)  {  L=(l[0]<l[1]) ? l[0] : l[1]; pout=ps+L*K; return true; }



		 return true;
		}
		*/

#define BOX_EPS 1E-10

		bool Box::next(const maths::Vector<double>& ps, const maths::Vector<double>& K, maths::Vector<double>& pout)
		{
			maths::Vector<double> k = H * K;
			maths::Vector<double> p = H * (ps - P);

			double tmin, tmax;
			if (k[0] >= 0)
			{
				tmin = (bounds[0][0] - p[0]) / k[0];
				tmax = (bounds[1][0] - p[0]) / k[0];
			}
			else
			{
				tmin = (bounds[1][0] - p[0]) / k[0];
				tmax = (bounds[0][0] - p[0]) / k[0];
			}

			double tymin, tymax;
			if (k[1] >= 0)
			{
				tymin = (bounds[0][1] - p[1]) / k[1];
				tymax = (bounds[1][1] - p[1]) / k[1];
			}
			else
			{
				tymin = (bounds[1][1] - p[1]) / k[1];
				tymax = (bounds[0][1] - p[1]) / k[1];
			}

			//        cout << ps << "    " << tmin << "  " << tymax << "   " << tymax << "   " << tymin << "   " << tmax << endl; 

			if ((tmin > tymax) || (tymin > tmax)) { pout = ps; return false; }

			if (tymin > tmin) tmin = tymin;
			if (tymax < tmax) tmax = tymax;

			double tzmin, tzmax;
			if (k[2] >= 0)
			{
				tzmin = (bounds[0][2] - p[2]) / k[2];
				tzmax = (bounds[1][2] - p[2]) / k[2];
			}
			else
			{
				tzmin = (bounds[1][2] - p[2]) / k[2];
				tzmax = (bounds[0][2] - p[2]) / k[2];
			}
			if ((tmin > tzmax) || (tzmin > tmax)) { pout = ps; return false; }
			if (tzmin > tmin) tmin = tzmin;
			if (tzmax < tmax) tmax = tzmax;
			if (tmin > BOX_EPS)
			{
				pout = ps + K * tmin;
				return true;
			}

			if (tmax > BOX_EPS)
			{
				pout = ps + K * tmax;
				return true;
			}

			pout = ps;
			return false;
		}

		/*
		 double tnear=-INFINITY;
		 double tfar = INFINITY;

		 double t1,t2;

		 Vector<double> p = ps - P;
		 Vector<double> k;
		 p = H*p;
		 k = H*K;
		 // cout << "p=" << p << "   k=" << k << endl;
		 t1=(bounds[0][0]-p[0])/k[0];
		 t2=(bounds[1][0]-p[0])/k[0];

		 if (t1 > t2) swap(t1,t2);
		 if (t1>tnear) tnear=t1;
		 if (t2<tfar)  tfar=t2;

		 if ((tfar < tnear) || (tfar<0)) return false;

		 t1=(bounds[0][1]-p[1])/k[1];
		 t2=(bounds[1][1]-p[1])/k[1];

		 if (t1 > t2) swap(t1,t2);
		 if (t1>tnear) tnear=t1;
		 if (t2<tfar)  tfar=t2;

		 if ((tfar < tnear) || (tfar<0)) return false;

		 t1=(bounds[0][2]-p[2])/k[2];
		 t2=(bounds[1][2]-p[2])/k[2];

		 if (t1 > t2) swap(t1,t2);
		 if (t1>tnear) tnear=t1;
		 if (t2<tfar)  tfar=t2;

		 if ((tfar < tnear) || (tfar < EPS)) return false;

		 pout=ps+tnear*K;
		 // cout << "tnear=" << tnear << endl;
		 return true;
		 */

		 // }


		 /*
		 bool Box::next(const Vector<double> &ps, const Vector<double> &K,Vector<double> &pout, const int inside)
		 {
			 double tmin, tmax, tymin, tymax, tzmin, tzmax;
			 double t0 = 0;
			 double t1 = 2.0*r0;
			 Vector<double> p = ps - P;
			 Vector<double> k;
			 p = H*p;
			 k = H*K;
		 cout << "Box::next ---------------------" << endl << "ps=" << ps << "    K=" << K << endl;
			 if (k[0] >= 0)
			 {
				 tmin = (bounds[0][0] - p[0]) / k[0];
				 tmax = (bounds[1][0] - p[0]) / k[0];
			 }
			 else
			 {
				 tmin = (bounds[1][0] - p[0]) / k[0];
				 tmax = (bounds[0][0] - p[0]) / k[0];
			 }
		 cout << "A : tmin=" << tmin << "    tmax=" << tmax << endl;
			 if (k[1] >= 0)
			 {
				 tymin = (bounds[0][1] - p[1]) / k[1];
				 tymax = (bounds[1][1] - p[1]) / k[1];
			 }
			 else
			 {
				 tymin = (bounds[1][1] - p[1]) / k[1];
				 tymax = (bounds[0][1] - p[1]) / k[1];
			 }
		 cout << "B : tymin=" << tymin << "    tymax=" << tymax << endl;

			 if ((tmin > tymax) || (tymin > tmax)) { cout << "not found! 1 " << endl; return false; }
			 if (tymin > tmin) tmin = tymin;
		 cout << "C : tymin=" << tmin << endl;
			 if ((tymax < tmax ) && (tymax > EPS) ) tmax = tymax;
		 cout << "D : tymax=" << tmax << endl;

			 if (k[2] >= 0)
			 {
				 tzmin = (bounds[0][2] - p[2]) / k[2];
				 tzmax = (bounds[1][2] - p[2]) / k[2];
			 }
			 else
			 {
				 tzmin = (bounds[1][2] - p[2]) / k[2];
				 tzmax = (bounds[0][2] - p[2]) / k[2];
			 }
		 cout << "B : tzmin=" << tzmin << "    tzmax=" << tzmax << endl;

			 if ((tmin>tzmax) || (tzmin > tmax)) { cout << "not found! A" << endl; return false; }
			 if (tzmin > tmin) tmin = tzmin;
			 if ((tmax > tzmax) && (tzmax > EPS)) tmax = tzmax;
			 double t = max(tmax, tmin);

				 if (t<EPS) { cout << "not found! B " << endl; return false; }
			 pout=ps + t*K;
				 cout << " t=" << t << "   pout=" << pout << endl;
			 return t>EPS;
		 }
		 */
maths::Vector<double> Box::norm(const maths::Vector<double>& ps)
		{
	maths::Vector<double>p = ps - P;
			p = H * p;
			/*	if (fabs(p[0] - d[0] / 2.0) < EPS) { cout << "n=" << ex << endl; return ex; }
				if (fabs(p[0] + d[0] / 2.0) < EPS) { cout << "n=" << -ex << endl; return -ex; }
				if (fabs(p[1] - d[1] / 2.0) < EPS) { cout << "n=" << ey << endl; return  ey; }
				if (fabs(p[1] + d[1] / 2.0) < EPS) { cout << "n=" << -ey << endl; return -ey; }
				if (fabs(p[2] - d[2] / 2.0) < EPS) { cout << "n=" << ez << endl; return ez; }
				if (fabs(p[2] + d[2] / 2.0) < EPS) { cout << "n=" << -ez << endl; return -ez;  } */

			if (fabs(p[0] - d[0] / 2.0) < EPS) return R * maths::ex;
			if (fabs(p[0] + d[0] / 2.0) < EPS) return -R * maths::ex;
			if (fabs(p[1] - d[1] / 2.0) < EPS) return  R * maths::ey;
			if (fabs(p[1] + d[1] / 2.0) < EPS) return -R * maths::ey;
			if (fabs(p[2] - d[2] / 2.0) < EPS) return R * maths::ez;
			if (fabs(p[2] + d[2] / 2.0) < EPS) return -R * maths::ez;
			//        cout << "n=" << zero << "   ps=" << ps << "   p=" << p << "   d=" << d << "    p[0] - d[0] / 2.0=" << p[0] - d[0] / 2.0 << endl;
			return maths::one;
		}

		bool Box::isInside(const maths::Vector<double>& Ps)
		{
			maths::Vector<double> p = P - Ps;
			p = H * p;
			bool xtest = fabs(p[0]) < d[0] / 2.0;
			bool ytest = fabs(p[1]) < d[1] / 2.0;
			bool ztest = fabs(p[2]) < d[2] / 2.0;
			return xtest && ytest && ztest;
		}

		void Box::initQuad()
		{
			maths::Vector<double> b0 = H * bounds[0];
			maths::Vector<double> b1 = H * bounds[1];
			por = maths::Vector<double>(std::max(b0[0], b1[0]), std::max(b0[1], b1[1]), std::max(b0[2], b1[2])) + P;
			pul = maths::Vector<double>(std::min(b0[0], b1[0]), std::min(b0[1], b1[1]), std::min(b0[2], b1[2])) + P;
		}

		void Box::setr0(double r0)
		{
			/*	d = d / this->r0 * r0;
				bounds[0] = P - d / 2.0;
				bounds[1] = P + d / 2.0;
				calcDiag();
				initQuad();*/
			this->r0 = r0;
		}

		std::ostream& operator<< (std::ostream& os, Box B)
		{
			os << B.P << "  " << B.d << std::endl;
			return os;
		}
	}
}
