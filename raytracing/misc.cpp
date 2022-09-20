/***************************************************************************
                          misc.cpp  -  description
                             -------------------
    begin                : Mit Mär 12 2003
    copyright            : (C) 2003 by Thomas Weigel
    email                : weigel@lat.ruhr-uni-bochum.de
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "misc.h"
#include <iostream>


// extern long int deletedItems, createdItems;
namespace GOAT
{
	namespace raytracing
	{
		void initInc(ObjectShape* E)
		{
			switch (E->type)
			{
			case OBJECTSHAPE_ELLIPSOID: ((Ellipsoid*)E)->initQuad(); break;
			case OBJECTSHAPE_SURFACE: ((surface*)E)->initQuad(); break;
			}
		}

		void setR0(ObjectShape* E, double r0)
		{
			switch (E->type)
			{
			case OBJECTSHAPE_ELLIPSOID: ((Ellipsoid*)E)->setr0(r0); break;
			case OBJECTSHAPE_SURFACE: ((surface*)E)->setr0(r0); break;
			}
		}

		void copyInc(ObjectShape*& d, ObjectShape* s)
		{
			// createdItems++;
			switch (s->type)
			{
			case OBJECTSHAPE_ELLIPSOID: d = new Ellipsoid(*((Ellipsoid*)s)); break;
			case OBJECTSHAPE_SURFACE: d = new surface(*((surface*)s)); break;
			case OBJECTSHAPE_BOX: d = new Box(*((Box*)s)); break;
			}
		}

		void copyFormList(ObjectShape**& d, ObjectShape** s, int anz)
		{
			/*
			  Kopiert die Einschluss-Liste s in die neue Liste d
			*/
			//	cout << "COPY OBJECT LIST" << endl;
			if ((anz > 0) && (s != 0))
			{
				d = (ObjectShape**)malloc(sizeof(ObjectShape*) * anz);
				for (int i = 0; i < anz; i++)
				{
					switch (s[i]->type)
					{
					case OBJECTSHAPE_ELLIPSOID: d[i] = new Ellipsoid(*((Ellipsoid*)s[i])); break;
					case OBJECTSHAPE_SURFACE: d[i] = new surface(*((surface*)s[i])); break;
					case OBJECTSHAPE_BOX: d[i] = new Box(*((Box*)s[i])); break;
					}
					d[i]->initQuad();
					//   sprintf (d[i]->Beschreibung,"%s",s[i]->Beschreibung); 
				}
			}
		}

		void binWriteIncList(std::ofstream& os, ObjectShape** E, int anz)
		{
			os.write((char*)&anz, (char)sizeof(anz));
			if (anz > 0)
				for (int i = 0; i < anz; i++)
					binWriteInc(os, E[i]);
		}

		void binWriteInc(std::ofstream& os, ObjectShape* E)
		{
			os.write((char*)&E->type, (char)sizeof(E->type));
			switch (E->type)
			{
			case OBJECTSHAPE_ELLIPSOID: ((Ellipsoid*)E)->binWrite(os); break;
			case OBJECTSHAPE_SURFACE: ((surface*)E)->binWrite(os); break;
			case OBJECTSHAPE_BOX: ((Box*)E)->binWrite(os); break;
			}
		}

		void binReadIncList(std::ifstream& is, ObjectShape**& E, int anz)
		{
			int Anz;
			is.read((char*)&Anz, (char)sizeof(Anz));
//			std::cout << "ANZ=" << Anz << std::endl;
			if (anz > 0)
			{
				E = (ObjectShape**)malloc(sizeof(ObjectShape*) * anz);
				for (int i = 0; i < anz; i++)
				{
					binReadInc(is, E[i], true);
				}
			}
			else E = 0;
		}

		void binReadInc(std::ifstream& is, ObjectShape*& E, bool isNew)
		{
			int type;
			is.read((char*)&type, (char)sizeof(type));
			// createdItems++;
			switch (type)
			{
			case OBJECTSHAPE_ELLIPSOID: if (isNew) E = new Ellipsoid; ((Ellipsoid*)E)->binRead(is); break;
			case OBJECTSHAPE_SURFACE: if (isNew) E = new surface; ((surface*)E)->binRead(is); break;
			case OBJECTSHAPE_BOX: if (isNew) E = new Box; ((Box*)E)->binRead(is); break;
			}
		}

		void deleteInc(ObjectShape* E)
		{
			// deletedItems++;
			/*
			 switch (E->type)
			 {
			  case FUNSURF : delete ((Funsurf *) E); break;
			  case OBJECTSHAPE_ELLIPSOID : delete ((Ellipsoid *) E); break;
			  case SUPERELLIPSOID_D : delete ((Superellipsoid_D *)E); break;
			  case SUPERELLIPSOID : delete ((Superellipsoid *)E); break;
			  case OBJECTSHAPE_SURFACE : delete ((surface *)E); break;
			  case ZYLINDER : delete ((Zylinder *)E); break;
			  case COMPOUND : delete ((Compound *)E); break;
			 }*/
			delete E;
			E = 0;
		}
		double NA2w0(double lambda, double NA, std::complex<double> n)
		{
			double nr = real(n);
			double w0 = lambda / M_PI * sqrt(nr * nr - NA * NA) / NA;
			return w0;
		}

		maths::Vector<double> force(maths::Vector<double> norm, tubedRay Se, tubedRay Sr, tubedRay St, double df)
			// hier wird mu_x als 1 angenommen !
			// norm zeigt entgegen der Einfallsrichtung des Strahls !!!!
		{
			maths::Vector<double> F, fe, fr, ft, n;
			// double dF=Se.flaeche();
			double dF = df;
			double ne = std::real(Se.n);
			double nr = std::real(Sr.n);
			double nt = real(St.n);
			// double h=Se.k0/Z0; 
			double h = 1.0 / c_light;
			if (!Se.inObject) fe = ne * dF * abs2(Se.E[4]) * Se.k[4] * fabs(norm * Se.k[4]);
			if (!Sr.inObject) fr = nr * dF * abs2(Sr.E[4]) * Sr.k[4] * fabs(norm * Sr.k[4]);
			if (!St.inObject) ft = nt * dF * abs2(St.E[4]) * St.k[4] * fabs(norm * St.k[4]);

			F = h * (fe - fr - ft);

			//  cout << Se.P[4] << "    " << F << endl;
			return F;
		}

		/*Vector<double> force (Vector<double> norm, Strahl Se, Strahl Sr, Strahl St, double df)
		{
			double dA=Se.flaeche();
			double phi_r,phi_t;
			Vector<double> F;
			phi_r=0.5*c_light*real(Sr.n)*eps0*dA*abs2(Sr.E[4]);
			phi_t=0.5*c_light*real(St.n)*eps0*dA*abs2(St.E[4]);
			F=real(Se.n)/c_light*((Sr.k[4]-Se.k[4])*phi_r+(St.k[4]-Se.k[4])*phi_t);
			return F;
		}*/

		double gaussw(double z, double wvl, double w0)
		{
			double z0 = M_PI * w0 * w0 / wvl;
			double w = w0 * sqrt(1 + z * z / (z0 * z0));
			return w;
		}

		std::complex<double> gaussphase(maths::Vector<double> P, maths::Vector<double> F, maths::Vector<double> k, double w0, double k0)
			// Achtung k wird als normiert angenommen !
		{
			maths::Vector<double> h = F - P;
			double z = k * h;
			double r = abs(h - z * k);
			double z0 = k0 * w0 / 2.0;
			double rz = z / z0;
			double zeta = atan(rz);
			double R = z * (1 + rz * rz);
			return exp(-I * k0 * r * r / 2.0 / R) * exp(I * (zeta - k0 * z));
		}

		float readLE_float32(std::istream& is)
		{
			char* d;
			char h;
			float f;
			d = (char*)&f;
			is.read(d, 4);
			/*h=d[3];
			d[3]=d[0];
			d[0]=h;
			h=d[1];
			d[1]=d[2];
			d[2]=h;*/
			return f;
		}

		int readLE_int32(std::istream& is)
		{
			unsigned char d[4];
			is.read((char*)d, 4);
			int i = d[0] + d[1] * 256 + d[2] * 65536 + d[3] * 16777216;
			return i;
		}
	}
}