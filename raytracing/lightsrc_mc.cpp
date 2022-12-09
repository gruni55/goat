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
            } while ((x>=-D/2.0) && (x<=D/2.0));

            do 
            {
               y=nd(gen);
            } while ((y>=-D/2.0) &&(y<=D/2.0));

            GOAT::maths::Vector<double> P=Pos + x*e1 + y*e2;
            return P;
        }

        LightSrcGauss_mc::LightSrcGauss_mc (const LightSrcGauss_mc & L) : LightSrcGauss(L)
        {
            double z=abs(Pos-focuspos); 
            double w=calcw(z);
            stddev=w*M_SQRT1_2;
        }

        LightSrcGauss_mc::LightSrcGauss_mc(maths::Vector<double> Pos, int N, double wvl, double w0, maths::Vector<double> focuspos, double D, maths::Vector<std::complex<double> > Pol, int raytype, double r0) 
                                       : LightSrcGauss(Pos,N,wvl,w0,focuspos,D,Pol,raytype,r0)
        {
          stddev=w*M_SQRT1_2;
        }

        int LightSrcGauss_mc::next(Ray_pow& S)
        {                    
			maths::Vector<double> fp, P;
			//  maths::Vector<double> fp,P=Pos+(i1*density-D/2.0)*e1; // nur zu TESTZWECKEN !!!
				// P : Startort, density=Anzahl Strahlen/Längeneinheit, D: Breite der Lichtquelle
				//     i1=momentaner Index in e1-Richtung (normalerweise x-Richtung), i2=momentaner Index in e2-Richtung (normalerweise y-Richtung)

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
			
			fp = focuspos - P;  // Hier beginnt der Strahl
			hk = fp / abs(fp);  // Normierter Richtungsvektor vom Startpunkt zum Fokus

			h = k % hk;
			absh = abs(h);

			x3 = fp * k;
			R = x3 * (1.0 + z0 * z0 / (x3 * x3));
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

			// Pol=maths::Vector<std::complex<double> >(0,1,0);  // y-Polarisation

			DM = rotMatrix(h, gamma);
			S = Ray_pow(1, P, Pol, hk, n0, r0, 2.0 * M_PI / wvl, numObjs, Obj);
			S.k = k;
			S.n = n0;
			absfp = fp * k;
			g = gaussw(-fp * k, 2.0 * M_PI / std::real(S.k0), w0); // w(z)  


			r2 = S.P[0] * S.P[0] + S.P[1] * S.P[1];

			S.Pow = P0 * k * hk * exp(-2.0 * r2 / g / g) / ((double)N * N) * D * D / g;//*sqrt(8.0/M_PI);// *real(Normfak);

			// cout << "%w=" << g << "   Pow=" << S.Pow << "   P0=" << P0 << endl;
#ifdef MIT_NORMIERUNG
			S.Pow = k * hk * exp(-2.0 * r2 / g / g) / ((double)N * N) * D * D / g;
#endif
			S.E1 = Pol;
			S.E2 = sqrt(S.Pow) * Pol;

			maths::Vector<double> F = focuspos;
			S.k = focuspos - S.P;   // Richtungsvektor auf den Fokus gerichtet
			S.k /= abs(S.k);
			h = k % S.k;
			absh = abs(h);
			if (absh != 0)
			{
				h /= absh;
				gamma = std::acos(S.k * k);
				DM = rotMatrix(h, gamma);
				S.E1 = DM * S.E1;
				S.E2 = DM * S.E2;
			}
			// else S.E1=E0*Pol;
			S.n = n0;
			i1++;
			if (type == LIGHTSRC_SRCTYPE_TOPHAT)
			{
				S.E2 = S.E2 / abs(S.E2);
				S.Pow = 1.0;
			}
			Pall += abs2(S.E2);
			//  if (i1*density>D) { return LIGHTSRC_IS_LAST_RAY; } // nur zu TESTZWECKEN !!!!
			if (i1 * density > D) { i1 = 0; i2++; }
			if (i2 * density > D) return LIGHTSRC_IS_LAST_RAY;
			return LIGHTSRC_NOT_LAST_RAY;
        }

        int LightSrcGauss_mc::next (IRay& S)
        {
            maths::Vector<double> fp, P = genStartingPos();
			// maths::Vector<double> P=Pos+(i1*density-D/2.0)*e1; // NUR ZU TESTZWECKEN !!!!!!!!!! 
			double x1, x2, x3;
			x1 = P * e1;
			x2 = P * e2;
			double r2 = x1 * x1 + x2 * x2;

			double s2 = w0 * w0 / log(2.0);
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
			S.k = k;
			S.n = n0;
			
		
			 E0=1.0;
			maths::Vector<double> F = focuspos;

			S.k = focuspos - S.P;   // Richtungsvektor auf den Fokus gerichtet
			S.k /= abs(S.k);
			h = k % S.k;
			absh = abs(h);
			if (absh != 0)
			{
				h /= absh;
				gamma = std::acos(S.k * k);
				DM = rotMatrix(h, gamma);
				S.E1 = E0 * DM * Pol;
				S.E2 = E0 * DM * Pol2;
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

        void LightSrcGauss_mc::reset() 
        {
            LightSrc::reset();
            rayCounter=0;
        }
    }
}