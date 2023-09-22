#include "lightsrc_mc.h"
#include <random>

namespace GOAT
{
    namespace raytracing
    {
        GOAT::maths::Vector<double> LightSrcGauss_mc::genStartingPos ()
        {
            std::random_device rd;
            std::mt19937 gen(rd());
            std::normal_distribution<double> nd (0,stddev);

            double x,y;

            do 
            {
               x=nd(gen);
            } while ((x<-D1/2.0) || (x>D1/2.0));

            do 
            {
               y=nd(gen);
            } while ((y<-D2/2.0) || (y>D2/2.0));

            GOAT::maths::Vector<double> P=Pos + x*e1 + y*e2;
            return P;
        }

        LightSrcGauss_mc::LightSrcGauss_mc (const LightSrcGauss_mc & L) : LightSrcGauss(L)
        {
            double z=abs(Pos-focuspos); 
            double w=calcw(z);
            stddev=w*M_SQRT1_2;
			D1 = L.D1;
			D2 = L.D2;
			type = LIGHTSRC_SRCTYPE_GAUSS_MC;
        }

        LightSrcGauss_mc::LightSrcGauss_mc(maths::Vector<double> Pos, int N, double wvl, double w0, maths::Vector<double> focuspos, double D, maths::Vector<std::complex<double> > Pol, int raytype, double r0) 
                                       : LightSrcGauss(Pos,N,wvl,w0,focuspos,D,Pol,raytype,r0)
        {
          stddev=w*M_SQRT1_2;
		  D1 = D;
		  D2 = D;
		  type = LIGHTSRC_SRCTYPE_GAUSS_MC;
        }

        int LightSrcGauss_mc::next(Ray_pow& S)
        {                    
			maths::Vector<double> fp, P = genStartingPos() + Pos;

			
			double x1, x2, x3; // Hilfsgrößen
			x1 = P * e1;
			x2 = P * e2;
			double r2 = x1 * x1 + x2 * x2;  // Quadrat des Abstands von der Laserstrahlachse

			double s2 = w0 * w0 / log(2.0);  // w0 : Strahltaille 
			double g;
			double L = 0.1;
			double R;

			maths::Vector<double> h, hk;
			maths::Matrix<double> DM;
			maths::Vector<std::complex<double> > E;
			std::complex<double> E0;
			double absh, gamma;
			double zeta;
			double absfp;
			
			fp = focuspos - P;  // Vector from the starting point to the focus
			hk = fp / abs(fp);  // Normalized vector pointing from start to the focus

			h = k % hk;
			absh = abs(h);

			x3 = fp * k; // This is the coordinate along the laser axis (i.e. in the laser coordinate system: z)

			R = x3 * (1.0 + z0 * z0 / (x3 * x3)); // R(z)
			if (z0 == 0) zeta = M_PI / 2.0;
			else zeta = atan(x3 / z0);


			if (absh == 0)
			{
				hk = k;
				gamma = 0;
			}
			else
			{
				h /= absh;
				gamma = std::acos(k * hk / (abs(k) * abs(hk)));
			}

			
			DM = rotMatrix(h, gamma);
			S = Ray_pow(1, P, Pol, hk, n0, r0, 2.0 * M_PI / wvl, numObjs, Obj);
			S.suppress_phase_progress = suppress_phase_progress;
			S.k = k;
			S.n = n0;
			absfp = fp * k;
			g = gaussw(-fp * k, 2.0 * M_PI / std::real(S.k0), w0); // w(z)  


			r2 = S.P[0] * S.P[0] + S.P[1] * S.P[1];

			// S.Pow = 1.0;
			S.Pow = P0;

		    S.E1 = Pol;
			S.E2 = Pol;

			std::complex<double> I(0.0,1.0); 
			std::complex<double> Phase=exp (I*(S.k0*S.n0*r2/(2.0*R)+S.k0*S.n0*x3-zeta)); // Phase factor
			maths::Vector<double> F = focuspos;
			S.k = fp / abs(fp);   // directional vector (normalized vector, that points in the direction of the focus)
			h = k % S.k;
			absh = abs(h);
			if (absh != 0)
			{
				h /= absh;
				gamma = std::acos(S.k * k);
				DM = rotMatrix(h, gamma);
				S.E1 = Phase * DM * S.E1;
				S.E2 = Phase * DM * S.E2;
			}
			
			S.n = n0;
			i1++;
			
			Pall += abs2(S.E2);
			rayCounter++;
			if ((rayCounter >= N) && (N > -1)) return LIGHTSRC_IS_LAST_RAY;
			return LIGHTSRC_NOT_LAST_RAY;
        }

        int LightSrcGauss_mc::next (IRay& S)
        {
            maths::Vector<double> fp, P = genStartingPos();
			double x1, x2;
			x1 = P * e1;
			x2 = P * e2;
			double r2 = x1 * x1 + x2 * x2;

			double s2 = w0 * w0 / log(2.0);
			double g;
			double L = 0.1;			

			maths::Vector<double> h, hk;
			maths::Matrix<double> DM;
			maths::Vector<std::complex<double> > E;
			std::complex<double> E0;
			double absh, gamma;
		
			fp = focuspos - P;  // Hier beginnt der Strahl
			hk = fp / abs(fp);  // Normierter Richtungsvektor vom Startpunkt zum Fokus


			h = k % hk;
			absh = abs(h);

			
			if (absh == 0)
			{
				hk = k;
				gamma = 0;
			}
			else
			{
				h /= absh;
				gamma = std::acos(k * hk / (abs(k) * abs(hk)));
			}

			

			DM = rotMatrix(h, gamma);
			S = IRay(P, Pol, hk, 1.0, r0, 2.0 * M_PI / wvl, numObjs, Obj);
			S.suppress_phase_progress = suppress_phase_progress;
			S.k = k;
			S.n = n0;
			
		
			 E0=sqrt(P0);
			maths::Vector<double> F = focuspos;

			S.k = focuspos - S.P;   // Richtungsvektor auf den Fokus gerichtet
			double z = abs(S.k);
			S.k /= z;
			z = -z;
			h = k % S.k;
			absh = abs(h);
			if (absh != 0)
			{
				h /= absh;
				gamma = std::acos(S.k * k);
				DM = rotMatrix(h, gamma);
				std::complex<double> fak = exp(-std::complex<double>(0, 1) * S.k0 * (z + r2 / (2.0 * z) - M_PI_2));
				S.E1 = E0 * DM * Pol * fak;
				S.E2 = E0 * DM * Pol2 * fak;
			}
			else
			{
				S.E1 = E0 * Pol;
				S.E2 = E0 * Pol2;
			}
			S.n = n0;

            
			rayCounter++;
            if ( (rayCounter >= N) && (N>-1)) return LIGHTSRC_IS_LAST_RAY;
			return LIGHTSRC_NOT_LAST_RAY;
        }

		int LightSrcGauss_mc::next(tubedRay& S)
		{

			maths::Vector<double> Ph = genStartingPos();
			maths::Vector<double> P = Pos + Ph;
			maths::Vector<double> k = focuspos - P;  // Richtung des Strahles
			k = k / abs(k);
			double E0 = sqrt(P0);
			S = tubedRay(P, density, density, Pol, k, 1.0, r0, 2.0 * M_PI / wvl, numObjs, Obj);
			S.suppress_phase_progress = suppress_phase_progress;
			S.setN0(n0);
			S.n = n0;

			maths::Vector<double> h;
			maths::Matrix<double> DM;
			double absh, gamma;

			for (int i = 0; i < 5; i++)
			{
				S.k[i] = focuspos - S.P[i];   // Richtungsvektor auf den Fokus gerichtet
				S.k[i] /= abs(S.k[i]);
				h = k % S.k[i];
				absh = abs(h);
				if (absh != 0)
				{
					h /= absh;
					gamma = std::acos(S.k[i] * k);
					DM = rotMatrix(h, gamma);
					S.E[i] = E0 * DM * Pol;					
				}
				else
				{
					S.E[i] = E0 * Pol;
				}
				S.n = n0;
			}
			rayCounter++;
			if ((rayCounter >= N) && (N > -1)) return LIGHTSRC_IS_LAST_RAY;
			return LIGHTSRC_NOT_LAST_RAY;
		}

        void LightSrcGauss_mc::reset() 
        {
            LightSrc::reset();
            rayCounter=0;
        }


		LightSrcPlane_mc::LightSrcPlane_mc(const LightSrcPlane_mc &L) : LightSrcPlane(L)
		{
			D1 = L.D1;
			D2 = L.D2;
			type = LIGHTSRC_SRCTYPE_PLANE_MC;
		}

		LightSrcPlane_mc::LightSrcPlane_mc(maths::Vector<double> Pos, int N, double wvl, double D, 
                                  maths::Vector<std::complex<double> > Pol, int raytype, double r0) : LightSrcPlane(Pos,N,wvl,D,Pol,raytype,r0)
		{
			D1 = D;
			D2 = D;
			type = LIGHTSRC_SRCTYPE_PLANE_MC;
		}

		void LightSrcPlane_mc::reset()
		{
			rayCounter=0;
		}

		int LightSrcPlane_mc::next (IRay &S)
		{
			Plane E;

			maths::Vector<double> P = genStartingPos();
			E.e1 = e1;
			E.e2 = e2;
			E.n = k;
			S = IRay(P, Pol, k, 1.0, r0, 2.0 * M_PI / wvl, numObjs, Obj);
			S.suppress_phase_progress = suppress_phase_progress;
			S.E1 = Pol / N;
			S.E2 = Pol2 / N;
			// S.init_Efeld(E,Pol);
			rayCounter++;
			if ((rayCounter >= N) && (N > -1)) return LIGHTSRC_IS_LAST_RAY;
			return LIGHTSRC_NOT_LAST_RAY;			
		}

		int LightSrcPlane_mc::next(Ray_pow& S)
		{

			Plane E;
			double Pow;

			maths::Vector<double> P = genStartingPos();
			E.e1 = e1;
			E.e2 = e2;
			E.n = k;
			// Pow = 1.0 / ((double)(N * N) * D * D);
			Pow = 1.0;
			S = Ray_pow(Pow, P, Pol, k, 1.0, r0, 2.0 * M_PI / wvl, numObjs, Obj);
			S.suppress_phase_progress = suppress_phase_progress;
			S.initElectricField(E, Pol);
			S.P = P;
			S.E1 = Pol;
			S.E2 = sqrt(Pow) * Pol / (double)(N * N);
			S.k = k;
			i1++;
			Pall += abs2(S.E2);
			rayCounter++;
			if ((rayCounter >= N) && (N > -1)) return LIGHTSRC_IS_LAST_RAY;
			return LIGHTSRC_NOT_LAST_RAY;			
		}


		int LightSrcPlane_mc::next(tubedRay& S)
		{
			double Pow = 1.0;
			maths::Vector<double> P = genStartingPos();
			S = tubedRay(P, density, density, sqrt(Pow) * Pol, k, 1.0, r0, 2.0 * M_PI / wvl, numObjs, Obj);			
			S.suppress_phase_progress = suppress_phase_progress;
			S.setN0(n0);
			i1++;
			rayCounter++;
			if ((rayCounter >= N) && (N > -1)) return LIGHTSRC_IS_LAST_RAY;
			return LIGHTSRC_NOT_LAST_RAY;			
		}



		GOAT::maths::Vector<double>  LightSrcPlane_mc::genStartingPos ()
		{
			std::random_device rd;
            std::mt19937_64 gen(rd());
			std::uniform_real_distribution<double> udx(-D1 / 2.0, D1 / 2.0);
			std::uniform_real_distribution<double> udy(-D2 / 2.0, D2 / 2.0);

            double x,y;

               x=udx(gen);
			   y=udy(gen);            
            GOAT::maths::Vector<double> P=Pos + x*e1 + y*e2;
            return P;
		}

		LightSrcRing_mc::LightSrcRing_mc(const LightSrcRing_mc& L) : LightSrcPlane(L)
		{
			rmin = L.rmin;
			rmax = L.rmax;
			type = LIGHTSRC_SRCTYPE_RING_MC;
		}

		LightSrcRing_mc::LightSrcRing_mc( maths::Vector<double> Pos, int N, double wvl,double rmin, double rmax,
			maths::Vector<std::complex<double> > Pol , int raytype, double r0) : LightSrcPlane (Pos,N,wvl,rmax,Pol,raytype,r0)
		{
			this->rmin = rmin;
			this->rmax = rmax;
			type = LIGHTSRC_SRCTYPE_RING_MC;
		}

		GOAT::maths::Vector<double> LightSrcRing_mc::genStartingPos()
		{
			std::random_device rd;
			std::mt19937_64 gen(rd());
			std::uniform_real_distribution<double> uphi(0, 2.0 * M_PI);
			std::uniform_real_distribution<double> ur((rmin*rmin)/(rmax*rmax), 1.0);
			double r = std::sqrt(ur(gen));
			double phi = uphi(gen);
			double x, y;
			x = r * cos(phi);
			y = r * sin(phi);
			GOAT::maths::Vector<double> P = Pos + x * e1 + y * e2;
			return P;
		}

		int LightSrcRing_mc::next(IRay& S)
		{
			Plane E;

			maths::Vector<double> P = genStartingPos();
			E.e1 = e1;
			E.e2 = e2;
			E.n = k;
			S = IRay(P, Pol, k, 1.0, r0, 2.0 * M_PI / wvl, numObjs, Obj);
			S.suppress_phase_progress = suppress_phase_progress;
			S.E1 = Pol / N;
			S.E2 = Pol2 / N;
			// S.init_Efeld(E,Pol);
			rayCounter++;
			if ((rayCounter >= N) && (N > -1)) return LIGHTSRC_IS_LAST_RAY;
			return LIGHTSRC_NOT_LAST_RAY;
		}

		int LightSrcRing_mc::next(Ray_pow& S)
		{

			Plane E;
			double Pow;

			maths::Vector<double> P = genStartingPos();
			E.e1 = e1;
			E.e2 = e2;
			E.n = k;
			// Pow = 1.0 / ((double)(N * N) * D * D);
			Pow = 1.0;
			S = Ray_pow(Pow, P, Pol, k, 1.0, r0, 2.0 * M_PI / wvl, numObjs, Obj);
			S.suppress_phase_progress = suppress_phase_progress;
			S.initElectricField(E, Pol);
			S.P = P;
			S.E1 = Pol;
			S.E2 = sqrt(Pow) * Pol / (double)(N * N);
			S.k = k;
			i1++;
			Pall += abs2(S.E2);
			rayCounter++;
			if ((rayCounter >= N) && (N > -1)) return LIGHTSRC_IS_LAST_RAY;
			return LIGHTSRC_NOT_LAST_RAY;
		}


		int LightSrcRing_mc::next(tubedRay& S)
		{
			double Pow = 1.0;
			maths::Vector<double> P = genStartingPos();
			S = tubedRay(P, density, density, sqrt(Pow) * Pol, k, 1.0, r0, 2.0 * M_PI / wvl, numObjs, Obj);
			S.suppress_phase_progress = suppress_phase_progress;
			S.setN0(n0);
			i1++;
			rayCounter++;
			if ((rayCounter >= N) && (N > -1)) return LIGHTSRC_IS_LAST_RAY;
			return LIGHTSRC_NOT_LAST_RAY;
		}
    }
}