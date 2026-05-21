#include "detector.h"
#include "matrix.h"
#include "constants.h"
#include <iostream>
#include <filesystem>
#include <numbers>

namespace GOAT
{
	namespace raytracing
	{

		Detector::Detector(void)
		{
			d1 = 0;
			d2 = 0;
			n1 = 0;
			n2 = 0;			
		}

		Detector::Detector(int n1, int n2)
		{
			init(n1, n2);
		}

		
		void Detector::init(int n1, int n2)
		{	
			D.clear();
			D.shrink_to_fit();	  			
			this->n1 = n1;
			this->n2 = n2;
			  D.resize(n1);
    				for (int i = 0; i < n1; ++i)
        			D[i].resize(n2);
		}

		Detector::Detector(const Detector& Det)
		{
			init (Det.n1, Det.n2);
				init(Det.n1, Det.n2);
				for (int i1 = 0; i1 < n1; i1++)
					for (int i2 = 0; i2 < n2; i2++)
						D[i1][i2] = Det.D[i1][i2];			
			n = Det.n;
			e1 = Det.e1;
			e2 = Det.e2;
			P = Det.P;
			type=Det.type;
		}

		Detector& Detector::operator = (const Detector& Det)
		{
			if (this != &Det)
			{
				init(Det.n1, Det.n2);
				for (int i1 = 0; i1 < n1; i1++)
					for (int i2 = 0; i2 < n2; i2++)
						D[i1][i2] = Det.D[i1][i2];
			
			n = Det.n;
			e1 = Det.e1;
			e2 = Det.e2;
			P = Det.P;
			type=Det.type;
			}
			return *this;
		}

		void Detector::mult(std::complex<double> fac)
		{
			if (n1 > 0)
			{
				for (int i1 = 0; i1 < n1; i1++)
					for (int i2 = 0; i2 < n2; i2++)
						D[i1][i2] *= fac;
			}
		}

		Detector::~Detector(void)
		{
			clear();
		}


		void Detector::clean()
		{
			if (n1 > 0)
				for (int i1 = 0; i1 < n1; i1++)
					for (int i2 = 0; i2 < n2; i2++)
						D[i1][i2] = GOAT::maths::czero;
		}

		void Detector::clear()
		{			
				D.clear();
				D.shrink_to_fit();
				n1 = 0;
				n2 = 0;			
		}

		int Detector::N1() { return n1; }
		int Detector::N2() { return n2; }
        void Detector::setN1(int n1)
        {
			clear();
			this->n1=n1;
			init(n1,n2);
        }

		void Detector::setN2(int n2)
        {
			clear();
			this->n2=n2;
			init(n1,n2);
        }

        void Detector::setN(int n1, int n2)
        {
			clear();
			this->n1=n1;
			this->n2=n2;
			init(n1,n2);
        }

        double Detector::D1() { return d1; }
        double Detector::D2() { return d2; }
		void Detector::setD(double d1, double d2) { this->d1 = d1; this->d2 = d2; }
		void Detector::setD1(double d1) { this->d1 = d1;}
		void Detector::setD2(double d2) { this->d2 = d2; }
		void Detector::setID(std::string ID) { this->ID = ID; }
		std::string GOAT::raytracing::Detector::getID() { return ID; }

		void Detector::save(const char* fn)
		{
			std::ofstream os;
			os.open(fn);
			os << "%n1 " << n1 << std::endl;
			os << "%n2 " << n2 << std::endl;
			for (int i1 = 0; i1 < n1; i1++)
			{
				for (int i2 = 0; i2 < n2; i2++)
				{

					os << D[i1][i2] << std::endl;
				}
			}
			os.close();
		}

		bool Detector::load(const char* fn)
		{
			// Datei existiert?
			if (!std::filesystem::exists(fn)) return false;


			// Vorherige Daten sicher löschen
			D.clear();
			D.shrink_to_fit();
			n1 = 0;
			n2 = 0;

			
			std::ifstream is(fn);
			if (!is.is_open()) return false;

			std::string str;
			int nl1 = 0, nl2 = 0;

			// %n1 lesen
			is >> str;
			if (str != "%n1" || !(is >> nl1) || nl1 <= 0) return false;

			// %n2 lesen
			is >> str;
			if (str != "%n2" || !(is >> nl2) || nl2 <= 0) return false;

			// Speicher allozieren
			D.resize(nl1);
			for (int i = 0; i < nl1; ++i)
				D[i].resize(nl2);

			// Daten einlesen
			for (int i1 = 0; i1 < nl1; ++i1)
			{
				for (int i2 = 0; i2 < nl2; ++i2)
				{
					if (!(is >> D[i1][i2]))
					{
						std::cerr << "Fehler beim Einlesen von D[" << i1 << "][" << i2 << "]" << std::endl;
						return false;
					}
				}
			}

			// Nur bei Erfolg interne Dimensionen setzen
			n1 = nl1;
			n2 = nl2;

			is.close();
			return true;
		}

		void Detector::saveabs(const char* fn)
		{
			std::ofstream os;
			double h;
			os.open(fn);
			for (int i1 = 0; i1 < n1; i1++)
			{
				for (int i2 = 0; i2 < n2; i2++)
				{
					if (type >= 3) h = abs(D[i1][i2]);
					else
						h = abs(D[i1][i2]);
					os << h << "   ";
				}
				os << std::endl;
			}
			os.close();
		}


		void Detector::savereal(const char* fn, int type)
		{
			std::ofstream os;
			double h;
			os.open(fn);
			for (int i1 = 0; i1 < n1; i1++)
			{
				for (int i2 = 0; i2 < n2; i2++)
				{
					if (type >= 3) h = abs(D[i1][i2]);
					else
						h = real(D[i1][i2][type]);
					os << h << "   ";
				}
				os << std::endl;
			}
			os.close();
		}


		void Detector::saveimag(const char* fn, int type)
		{
			std::ofstream os;
			double h;
			os.open(fn);
			for (int i1 = 0; i1 < n1; i1++)
			{
				for (int i2 = 0; i2 < n2; i2++)
				{
					if (type >= 3) h = abs(D[i1][i2]);
					else
						h = imag(D[i1][i2][type]);
					os << h << "   ";
				}
				os << std::endl;
			}
			os.close();
		}


		void Detector::savePhase(const char* fn, int type)
		{
			std::ofstream os;
			double h;
			os.open(fn);
			for (int i1 = 0; i1 < n1; i1++)
			{
				for (int i2 = 0; i2 < n2; i2++)
				{
					h = arg(D[i1][i2][type]);
					os << h << "   ";
				}
				os << std::endl;
			}
			os.close();
		}

		/*bool Detector::cross(GOAT::maths::Vector<double> P, GOAT::maths::Vector<double> k, int& i1, int& i2, double& l)
		{
		  switch (type)
		  {
		  case DETECTOR_PLANE : return ((DetectorPlane *)this)->cross(P,k,i1,i2,l); break;
		  }
		  return false;
		}
		*/

		DetectorPlane::DetectorPlane(void)
		{
			n1 = 0;
			n2 = 0;			
			type = DETECTOR_PLANE;
		}

		DetectorPlane::DetectorPlane(maths::Vector<double> P, maths::Vector<double> n, double d1, double d2, int n1, int n2)
		{
			init(n1, n2);
			type = DETECTOR_PLANE;
			this->n = n / abs(n);
			if (abs(this->n % GOAT::maths::ex) > 1E-5)
				e1 = GOAT::maths::ex - (GOAT::maths::ex * n) * n;
			else
				e1 = GOAT::maths::ey - (GOAT::maths::ey * n) * n;
			e1 = e1 / abs(e1);
			e2 = n % e1;
			e2 = e2 / abs(e2);
			this->d1 = d1;
			this->d2 = d2;
			this->P = P;
			this->n1 = n1;
			this->n2 = n2;
		}

		DetectorPlane::DetectorPlane(GOAT::maths::Vector<double> P, GOAT::maths::Vector<double> n, double d, int N)
		{			
			init(N, N);
			type = DETECTOR_PLANE;
			this->n = n / abs(n);
			if (abs(this->n % GOAT::maths::ex) > 1E-5)
				e1 = GOAT::maths::ex - (GOAT::maths::ex * n) * n;
			else
				e1 = GOAT::maths::ey - (GOAT::maths::ey * n) * n;
			e1 = e1 / abs(e1);
			e2 = n % e1;
			e2 = e2 / abs(e2);			
			d1 = d;
			d2 = d;
			this->P = P;
			n1 = N;
			n2 = N;
			
		}

		DetectorPlane::DetectorPlane(GOAT::maths::Vector<double> P, GOAT::maths::Vector<double> e1, GOAT::maths::Vector<double> e2, int n1, int n2)
		{
			type = DETECTOR_PLANE;
			init(n1, n2);
			this->P = P;
			
			d1 = abs(e1);
			d2 = abs(e2);
			this->e1 = e1 / abs(e1);
			this->e2 = e2 / abs(e2);
		/*for (int i = 0; i<3; i++)
			{
			 if (e1[i]!=0)
			 this->e1[i]=1.0/e1[i]*n1;
			 else this->e1[i]=0;

			 if (e2[i]!=0)
			 this->e2[i]=1.0/e2[i]*n2;
			 else this->e2[i]=0;
			}*/
           /*if (n1 == 1) this->e1 = e1;
			else this->e1 = e1 / (double)(n1 - 1);

			if (n2 == 1) this->e2 = e2;
			else this->e2 = e2 / (double)(n2 - 1);
			*/
			this->n = this->e1 % this->e2;
			this->n /= abs(this->n);			
		}

		void DetectorPlane::setNorm(maths::Vector<double> n)
		{
			this->n = n / abs(n);
			if (abs(this->n % GOAT::maths::ex) > 1E-5)
				e1 = GOAT::maths::ex - (GOAT::maths::ex * n) * n;
			else
				e1 = GOAT::maths::ey - (GOAT::maths::ey * n) * n;
			e1 = e1 / abs(e1);
			e2 = n % e1;
			e2 = e2 / abs(e2);						
		}

		/*bool DetectorPlane::cross(GOAT::maths::Vector<double> P, GOAT::maths::Vector<double> k, int& i1, int& i2, double& l)
		{
			/*
			*  Ebene: n*(PE-P)=0
			*  Strahl: P=k*l+PS
			*  in Ebenengleichung n*PE-n*(k*l+PS)=0 => n*(PE-PS)=n*k*l => l=n*(PE-PS)/(n*k)
			*/
		/*	GOAT::maths::Vector<double> dP = this->P - P;
			double kn = k * n;
			if (kn == 0) return false;
			l = (dP * n) / kn;

			if (l<0) return false;

			/* Index berechnen */
		/*	GOAT::maths::Vector<double> Ph = P + l * k;
			dP = Ph - this->P;
			i1 = dP * e1 * n1 / d1 + n1 / 2.0;
			i2 = dP * e2 * n2 / d2 + n2 / 2.0;
         //  std::cout << "n1=" << n1 << "\tn2=" << n2 << "\ti1=" << i1 << "\ti2=" << i2 << std::endl;
			if ((i1 < 0) || (i1 >= n1) || (i2 < 0) || (i2 >= n2)) return false;
			return true;
		}*/


		bool DetectorPlane::cross(
			GOAT::maths::Vector<double> P,
			GOAT::maths::Vector<double> k,
			int& i1,
			int& i2,
			double& l)
		{
			GOAT::maths::Vector<double> dP = this->P - P;

			double kn = k * n;

			constexpr double eps = 1E-12;
			if (std::abs(kn) < eps)
				return false;

			l = (dP * n) / kn;

			if (l < 0)
				return false;

			GOAT::maths::Vector<double> Ph = P + l * k;
			dP = Ph - this->P;

			double u = dP * e1;
			double v = dP * e2;

			int ii1 = static_cast<int>(std::floor((u / d1 + 0.5) * n1));
			int ii2 = static_cast<int>(std::floor((v / d2 + 0.5) * n2));

			if (ii1 < 0 || ii1 >= n1 || ii2 < 0 || ii2 >= n2)
				return false;

			i1 = ii1;
			i2 = ii2;

			return true;
		}


		void DetectorPlane::smooth(double sigma)
		{
			uD = D; // unsmoothed field speichern
			const int nx = static_cast<int>(D.size());
			const int ny = static_cast<int>(D[0].size());

			const int r = static_cast<int>(std::ceil(3.0 * sigma));
			for (int x = 0; x < nx; ++x)
			{
				for (int y = 0; y < ny; ++y)
				{
					GOAT::maths::Vector<std::complex<double>> sum(0, 0, 0);
					double norm = 0.0;

					for (int dx = -r; dx <= r; ++dx)
					{
						int xx = std::clamp(x + dx, 0, nx - 1);

						for (int dy = -r; dy <= r; ++dy)
						{
							int yy = std::clamp(y + dy, 0, ny - 1);

							double w = std::exp(-(dx * dx + dy * dy) / (2.0 * sigma * sigma));

							sum = sum + uD[xx][yy] * w;
							norm += w;
						}
					}

					D[x][y] = sum / norm;
				}
			}
			isSmoothed_ = true;
		}

		
		void DetectorPlane::holographicField(std::vector<std::vector<GOAT::maths::Vector<std::complex<double>>>> &result) const
		{
			result = D; // Kopie des gespeicherten Feldes
			size_t x0 = -0.5 * (n1 - 1) * d1;
			size_t y0 = -0.5 * (n2 - 1) * d2;
			for (size_t i2= 0; i2 < n2; ++i2)
			{
				for (size_t i1 = 0; i1< n1; ++i1)
				{
					double x = x0 + i1 * d1;
					double y = y0 + i2 * d2;

					double phi = 2.0 * std::numbers::pi * (fx * x + fy * y);
					auto phase = std::polar(1.0, phi);

					result[i1][i2][0] *= phase;
					result[i1][i2][1] *= phase;
					result[i1][i2][2] *= phase;
				}
			}
		}

		std::ostream& operator << (std::ostream& os, Detector& D)
		{
			for (int i1 = 0; i1 < D.n1; i1++)
			{
				for (int i2 = 0; i2 < D.n2; i2++)
					os << D.D[i1][i2] << "  ";
				os << std::endl;
			}
			return os;
		}
		
	}
}
