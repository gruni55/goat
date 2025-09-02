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

#define EPS_T        1e-10   // nur t > EPS_T akzeptieren (Ray bias)
#define TAU          1e-10   // ~0 für Koeffizienten
#define TOL          1e-10   // Geometrische Toleranz
#define TOL_D        1e-16   // Diskriminanten-Toleranz
#define GRAZE_COS_MIN 1e-6   // Minimale |cos(theta)|, um Grazing zu akzeptieren (klein => schärfere Ablehnung)

		bool Cone::next(const maths::Vector<double>& ps,
			const maths::Vector<double>& ks,
			maths::Vector<double>& pout)
		{
			// Welt -> lokales Kegelsystem (Basiszentrum am Ursprung, Spitze bei z=h)
			maths::Vector<double> p = H * (ps - P);
			maths::Vector<double> k = H * ks;

			const double r = radius, h = height;
			const double s = r / h, s2 = s * s;

			auto in_z_range = [&](double t) {
				double z = p[2] + t * k[2];
				return z >= -TOL && z <= h + TOL;
				};

			// ---- Mantel-Quadratik: A t^2 + B t + C = 0
			double hpz = h - p[2];
			double A = (k[0] * k[0] + k[1] * k[1]) - s2 * (k[2] * k[2]);
			double B = 2.0 * (p[0] * k[0] + p[1] * k[1] + s2 * hpz * k[2]);
			double C = (p[0] * p[0] + p[1] * p[1]) - s2 * (hpz * hpz);

			// Kandidaten t's
			double t_lat = std::numeric_limits<double>::infinity();
			double t_bas = std::numeric_limits<double>::infinity();

			// --- Mantel: A!=0 und D>0 (Tangente & Sliding vermeiden)
			if (std::fabs(A) >= TAU) {
				double D = B * B - 4.0 * A * C;
				if (D > TOL_D) { // **Tangentialfall ausgeschlossen**
					double sD = std::sqrt(D);
					double t1 = (-B - sD) / (2.0 * A);
					double t2 = (-B + sD) / (2.0 * A);

					auto accept_mantel = [&](double t) {
						if (t <= EPS_T || !in_z_range(t)) return false;
						// **Sliding-Guard**: |n·k| / (|n||k|) >= GRAZE_COS_MIN
						double x = p[0] + t * k[0];
						double y = p[1] + t * k[1];
						double z = p[2] + t * k[2];
						maths::Vector<double> n(x, y, s2 * (h - z));  // ∇F/2
						double nk = std::fabs(n[0] * k[0] + n[1] * k[1] + n[2] * k[2]);
						double nn = std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
						double kk = std::sqrt(k[0] * k[0] + k[1] * k[1] + k[2] * k[2]);
						if (nn < TAU || kk < TAU) return false;
						double cosAbs = nk / (nn * kk);
						if (cosAbs < GRAZE_COS_MIN) return false; // **verwirf Strahlen entlang der Fläche**
						return true;
						};

					if (accept_mantel(t1)) t_lat = std::min(t_lat, t1);
					if (accept_mantel(t2)) t_lat = std::min(t_lat, t2);
				}
				// sonst: kein Manteltreffer (Tangente/keiner)
			}
			// A≈0 => Richtung „Mantel-generatrix“ -> **Sliding** ⇒ bewusst verwerfen

			// ---- Boden: z=0, nur echte Schnitte (keine Tangente/Parallelität)
			if (std::fabs(k[2]) > TAU) {
				double tB = -p[2] / k[2];
				if (tB > EPS_T) {
					double x = p[0] + tB * k[0];
					double y = p[1] + tB * k[1];
					double rr = x * x + y * y;
					// Kreis **ohne Rand** (Tangente ausschließen)
					if (rr < r * r - TOL) {
						// **Sliding-Guard am Boden**: |kz| groß genug (nicht nahezu parallel)
						double kk = std::sqrt(k[0] * k[0] + k[1] * k[1] + k[2] * k[2]);
						double cosAbs = std::fabs(k[2]) / (kk > TAU ? kk : 1.0);
						if (cosAbs >= GRAZE_COS_MIN) {
							t_bas = tB;
						}
					}
				}
			}
			else {
				// k nahezu parallel zur Bodenebene -> Sliding -> verwerfen
				// (Optional: 2D-Fall, wenn p[2]≈0, könntest du echte Durchschnitte prüfen, hier gewollt verwerfen.)
			}

			// ---- kleinsten gültigen Treffer wählen
			double t = std::min(t_lat, t_bas);
			if (!std::isfinite(t)) return false;

			// Weltkoordinate zurückgeben (immer ks benutzen, nicht lokales k!)
			pout = ps + t * ks;
			return true;
		}

		/*bool Cone::next(const maths::Vector<double>& ps, const maths::Vector<double>& ks, maths::Vector<double>& pout)
		{
			maths::Vector<double> n,n1, e1, e2;
			maths::Vector<double> p = H * (ps - P);
			maths::Vector<double> k = H * ks;
			double l1, l2;

			e2 = maths::ez;
			n = maths::norm(p % e2);
			e1 = maths::norm(n % e2);

			maths::Vector<double> Ph = e2 * height;
			maths::Vector<double> Psh = (Ph-p);
			maths::Vector<double> u, v;
			u =   -radius * e1 - Ph;
			v =   radius * e1 - Ph;
			double ku = k * u;
			double kv = k * v;
			double s2 = sideLen * sideLen;
			double lambda1, lambda2, lambdaB, lambda;
			lambda = -1;
			if (fabs(k[2]) > EPS)
			{
				l1 = (Psh * (e2 * ku - k[2] * u)) / (s2 * k[2] + height * ku);
				l2 = (Psh * (e2 * kv - k[2] * v)) / (s2 * k[2] + height * kv);
				// testing intersection with "left" side:
				lambda1 = -1;
				if ((EPS <= l1) && (l1 <= 1.0))
				{
					lambda1 = (Psh[2] - height * l1) / k[2];
	//				std::cout << "%-> " << Ph + l1 * u << "\t" << ps + lambda1 * ks << "\t" << k << std::endl;
				}
				// testing intersection with "right" side:
				lambda2 = -1;
				if ((EPS <= l2) && (l2 <= 1.0))
				{
					lambda2 = (Psh[2] - height * l2) / k[2];
			//		std::cout << "%->> " << Ph + l2 * v << "\t" << ps + lambda2 * ks <<  "\t" << k << std::endl;
				}
				
				// now testing intersection with bottom area
				lambdaB = -p[2] / k[2];

				maths::Vector<double> Pt = p + lambdaB * k;
				if (Pt[0] * Pt[0] + P[1] * Pt[1] > radius * radius) lambdaB = -1;

				
				if ((lambda1 < lambda2) && (lambda1 > EPS) || (lambda2<0)) lambda = lambda1;
				else if (lambda2 > EPS) lambda = lambda2;

				if ((lambdaB < lambda) && (lambdaB > EPS) || ((lambda < EPS) && (lambdaB > EPS)))
				{
					pout = ps + lambdaB * ks;
					return true;
				}
				else
					if (lambda > EPS)
					{
						pout = ps + lambda * ks;
 						return true;
					}				
			}
			else
			{
				l1 = Psh[2] / height;
				l2 = l1;
				lambda = -1;
				if ((EPS <= l1) && (l1 <= 1.0))
				{
					lambda1 = (Psh * u + s2 * l1) / ku;
					lambda2 = (Psh * v + s2 * l2) / kv;
		//			std::cout << "% " << Ph + l1 * u << "\t" << ps + lambda1 * ks << "\t" << k << std::endl;
					if ((lambda1 <= lambda2) && (lambda1 > 0)) lambda = lambda1;
					else  lambda = lambda2;

					if (lambda > EPS)
					{
						pout = ps + lambda * ks;
						return true;
					}

				}
			}
			pout = ps;
			return false;

		}
		*/
		/*
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

		*/
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

		maths::Vector<double> Cone::norm(const maths::Vector<double>& ps) 
		{
			maths::Vector<double> n;
			maths::Vector<double> p = H * (ps - P);
			
			if (fabs(p[2]) < 1E-10) return -R * maths::ez;
			// const double TOL = 1e-12;

			// Punkt ins lokale Kegelsystem
			maths::Vector<double> pL = H * (ps - P);

			// Spitze (Apex) im lokalen System
			maths::Vector<double> T(0.0, 0.0, height);

			// Achse (Einheitsvektor in z)
			maths::Vector<double> A(0.0, 0.0, 1.0);

			// Generator-Richtung (Apex -> Punkt)
			maths::Vector<double> v = pL - T;
			double vlen = std::sqrt(v * v);
			if (vlen < TOL) return maths::Vector<double>(0, 0, 0);
			maths::Vector<double> khat = v / vlen;

			// Azimutale Tangente (Achse × radial)
			maths::Vector<double> r = v - (v * A) * A;
			maths::Vector<double> ttheta = A % r;
			double tlen = std::sqrt(ttheta * ttheta);
			if (tlen < TOL) return maths::Vector<double>(0, 0, 0);
			ttheta /= tlen;

			// Normale (ttheta × khat)
			maths::Vector<double> nL = ttheta % khat;
			double nlen = std::sqrt(nL * nL);
			if (nlen < TOL) return maths::Vector<double>(0, 0, 0);
			nL /= nlen;

			// zurück in den Weltraum (H ist Rotation ⇒ inverse = transpose)
			maths::Vector<double> nW = R * nL;
			double nWlen = std::sqrt(nW * nW);
			if (nWlen > TOL) nW /= nWlen;

			return nW;
		
		}



		

		void Cone::setRadius(double radius)
		{
			this->radius = radius;
			coneAngle = atan(radius / height);
			initQuad();
		}

		double Cone::getRadius()
		{
			return radius;
		}

		void Cone::setHeight(double height)
		{
			this->height = height;
			coneAngle = atan(radius / height);
			initQuad();
		}

		double Cone::getHeight()
		{
			return height;
		}

	/*	maths::Vector<double> Cone::norm(const maths::Vector<double>& p)
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
		*/
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
                height=radius / tan (coneAngle);
                init();
                initQuad();
               }
                
               double Cone::getConeAngle ()
               {
                return coneAngle;
               }
            
	}
}
