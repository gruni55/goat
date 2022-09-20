#include "triangle.h"
namespace GOAT
{
	namespace raytracing
	{
		triangle::triangle()
		{
			//  P1=Vector<double>(0.0,0.0,0.0);
			//  P2=Vector<double>(0.0,0.0,0.0);
			//  P3=Vector<double>(0.0,0.0,0.0);
			//  u=0.0;
			//  v=0.0;
		}

		triangle::triangle(maths::Vector<double> ip1, maths::Vector<double> ip2, maths::Vector<double> ip3)
		{
			P[0] = ip1; P[1] = ip2; P[2] = ip3;
			calcSideVectors();
			setnorm();
		}

		triangle::triangle(const triangle& d)
		{
			this->n = d.n;
			for (int i = 0; i < 3; i++)
			{
				this->P[i] = d.P[i];
				this->f[i] = d.f[i];
			}

			this->u = d.u;
			this->v = d.v;
		}

		triangle& triangle::operator=(const triangle& dr)
		{
			for (int i = 0; i < 3; i++)
			{
				this->f[i] = dr.f[i];
				this->P[i] = dr.P[i];
			}
			this->n = dr.n;
			this->u = dr.u;
			this->v = dr.v;
			return *this;
		}

		void triangle::setnorm()
		{
			double hilf;
			n = (P[0] - P[2]) % (P[1] - P[2]);
			hilf = P[2] * n;
#ifdef __GNUC__
			n = n * copysign(1.0, hilf);
#else
			n = n * _copysign(1.0, hilf);
#endif
			n /= abs(n);
		}

		maths::Vector<double> triangle::getnorm(void)
		{
			// setnorm();
			return n;
		}

		triangle::triangle(maths::Vector<double> ip1, maths::Vector<double> ip2, maths::Vector<double> ip3, maths::Vector<double> P0)
		{
			P[0] = ip1 + P0; P[1] = ip2 + P0; P[2] = ip3 + P0;
			u = 0.0; v = 0.0;
			setnorm();
			calcSideVectors();
		}
		/*
		int triangle::calcIntersectionPoint(Vector<double> r, Vector<double> k, Vector<double>& p, double eps)
		{
			double t;
			Vector<double> e1 = P[1] - P[0];
			Vector<double> e2 = P[2] - P[0];
			double d = (e1 % k) * e2;
			if (fabs(d) < eps) return 0;
			Vector<double> pv = r - P[0];
			u = pv * e2 / d;
			if ( (u < eps) || (u > 1.0) ) return 0;
			v = -pv * e1 / d;
			if ( (v < eps) || (v > 1.0) ) return 0;
			p = P[0] + u * e1 + v * e2;
			t = (p - r) * k / abs(k);
			if (t < eps) return 0;
			return 1;
		}


		int triangle::calcIntersectionPoint(Vector<double> r, Vector<double> k, double& t, Vector<double>& p, double eps)
		{
			Vector<double> e1 = P[1] - P[0];
			Vector<double> e2 = P[2] - P[0];
			double d = (e1 % k) * e2;
			if (fabs(d) < eps) return 0;
			Vector<double> pv = r - P[0];
			u = pv * e2 / d;
			if ( (u < eps) || (u > 1.0) ) return 0;
			v = -pv * e1 / d;
			if ( (v < eps) || (v > 1.0) ) return 0;
			p = P[0] + u * e1 + v * e2;
			t = (p - r) * k / abs(k);
			if (t < eps) return 0;
			return 1;
		}

		*/
		int triangle::calcIntersectionPoint(maths::Vector<double> r, maths::Vector<double> k, maths::Vector<double>& p, double eps)
		{
			maths::Vector<double> edge1, edge2, pvec, qvec, tvec;
			double det, t, inv_det;

			edge1 = P[1] - P[0]; edge2 = P[2] - P[0];

			pvec = k % edge2;
			det = edge1 * pvec;

#ifdef TEST_CULL

			if (det < eps)
				return 0;

			tvec = r - P[0];
			u = tvec * pvec;

			if (u < 0.0 || u > det)
				return 0;

			qvec = tvec % edge1;
			v = k * qvec;

			if (v < 0.0 || u + v > det)
				return 0;

			t = edge2 * qvec;

			inv_det = 1.0 / det;

			t *= inv_det;
			u *= inv_det;
			v *= inv_det;

#else

			if (det > -eps && det < eps)
				return 0;

			inv_det = 1.0 / det;
			tvec = r - P[0];

			u = tvec * pvec * inv_det;

			if (u < 0.0 || u > 1.0)
				return 0;

			qvec = tvec % edge1;

			v = k * qvec * inv_det;

			if (v < 0.0 || u + v > 1.0)
				return 0;

			t = edge2 * qvec * inv_det;

#endif


			p = (1 - u - v) * P[0] + u * P[1] + v * P[2];

			if (t < eps)
				return 0;

			return 1;
		}

		double triangle::distance(maths::Vector<double> p, maths::Vector<double> k)
		{
			double t;
			maths::Vector<double> pout;
			if (calcIntersectionPoint(p, k, t, pout) == 0) return INFINITY;
			return t;
		}




		int triangle::calcIntersectionPoint(maths::Vector<double> r, maths::Vector<double> k,
			double& t, maths::Vector<double>& p, double eps)
		{
			maths::Vector<double> edge1, edge2, pvec, qvec, tvec;
			double det, inv_det;
			edge1 = P[1] - P[0]; edge2 = P[2] - P[0];

			pvec = k % edge2;
			det = edge1 * pvec;

#ifdef TEST_CULL

			if (det < eps)
				return 0;

			tvec = r - P[0];
			u = tvec * pvec;

			if (u < 0.0 || u > det)
				return 0;

			qvec = tvec % edge1;
			v = k * qvec;

			if (v < 0.0 || u + v > det)
				return 0;

			t = edge2 * qvec;

			inv_det = 1.0 / det;

			t *= inv_det;
			u *= inv_det;
			v *= inv_det;

#else

			if (det > -eps && det < eps)
				return 0;

			inv_det = 1.0 / det;
			tvec = r - P[0];

			u = tvec * pvec * inv_det;

			if (u < 0.0 || u > 1.0)
				return 0;

			qvec = tvec % edge1;

			v = k * qvec * inv_det;

			if (v < 0.0 || u + v > 1.0)
				return 0;

			t = edge2 * qvec * inv_det;

#endif


			p = (1 - u - v) * P[0] + u * P[1] + v * P[2];


			if (t < eps)
				return 0;

			return 1;
		}

		triangle::~triangle()
		{
		}

		maths::Vector<double>& triangle::operator[](int i)
		{
			return P[i];
		}

		const maths::Vector<double>& triangle::operator[](int i)
			const
		{
			return P[i];
		}

		void triangle::binWrite(std::ofstream& os)
		{
			for (int i = 0; i < 3; i++)
				P[i].binWrite(os);
		}

		void triangle::binRead(std::ifstream& is)
		{
			for (int i = 0; i < 3; i++)
				P[i].binRead(is);
		}



		std::ostream& operator << (std::ostream& os, const triangle& dr)
		{
			os << dr[0] << "  " << dr[1] << "  " << dr[2]; // << "  " << dr.n;
			return os;
		}

		triangle operator + (const triangle& dr, const maths::Vector<double>& v)
		{
			triangle h;
			h[0] = dr[0] + v;
			h[1] = dr[1] + v;
			h[2] = dr[2] + v;
			h.setnorm();

			return h;
		}

		triangle operator + (const maths::Vector<double>& v, const triangle& dr)
		{
			triangle h;
			h[0] = v + dr[0];
			h[1] = v + dr[1];
			h[2] = v + dr[2];
			h.setnorm();


			return h;
		}


		triangle operator - (const triangle& dr, const maths::Vector<double>& v)
		{
			triangle h;
			h[0] = dr[0] - v;
			h[1] = dr[1] - v;
			h[2] = dr[2] - v;
			h.setnorm();

			return h;
		}

		triangle operator - (const maths::Vector<double>& v, const triangle& dr)
		{
			triangle h;
			h[0] = v - dr[0];
			h[1] = v - dr[1];
			h[2] = v - dr[2];
			h.setnorm();

			return h;
		}

		triangle operator * (const maths::Matrix<double>& M, const triangle& dr)
		{
			triangle h;
			h[0] = M * dr[0];
			h[1] = M * dr[1];
			h[2] = M * dr[2];
			h.setnorm();

			return h;
		}

		triangle operator / (const triangle& dr, double a)
		{
			triangle h;
			h[0] = dr[0] / a;
			h[1] = dr[1] / a;
			h[2] = dr[2] / a;
			h.setnorm();
			return h;
		}


		triangle operator * (const triangle& dr, double a)
		{
			triangle h;
			h[0] = a * dr[0];
			h[1] = a * dr[1];
			h[2] = a * dr[2];
			h.setnorm();
			return h;
		}

		triangle operator * (double a, const triangle& dr)
		{
			triangle h;
			h[0] = a * dr[0];
			h[1] = a * dr[1];
			h[2] = a * dr[2];
			h.setnorm();
			return h;
		}
	}
}