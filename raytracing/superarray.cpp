/***************************************************************************
                          supergitter.cpp  -  description
                             -------------------
    begin                : Thu Dec 14 2000
    copyright            : (C) 2000 by Thomas Weigel
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

#include "superarray.h"

namespace GOAT
{
  namespace raytracing 
  {
      void saveExPhase(SuperArray<maths::Vector<std::complex<double> > > &S, char* FName, int i)
      {
          maths::Vector<int> Pi;
          maths::Vector<std::complex<double> > E;
          std::ofstream os;
          os.open(FName);
          std::complex<double> phase;
          if (S.type == IN_OBJECT)
          {
              os << "%Dimensionen " << S.n[i][0] << "  x  " << S.n[i][1] << "  x  " << S.n[i][2] << std::endl;

              for (int ix = 0; ix < S.n[i][0]; ix++)
                  for (int iy = 0; iy < S.n[i][1]; iy++)
                      for (int iz = 0; iz < S.n[i][2]; iz++)
                      {
                          os << atan2(imag(S.G[i][ix][iy][iz][0]), real(S.G[i][ix][iy][iz][0])) << std::endl;
                      }
          }
          else
          {
              for (int ix = 0; ix < S.nges[0]; ix++)
                  for (int iy = 0; iy < S.nges[1]; iy++)
                      for (int iz = 0; iz < S.nges[2]; iz++)
                      {
                          Pi = S.kugelindex(maths::Vector<int>(ix, iy, iz));
                          if (S.Error != NO_ERRORS) os << 0.0 << std::endl;
                          else
                          {
                              E = S.K[Pi[0]][Pi[1]][Pi[2]];
                              os << atan2(imag(E[0]), real(E[0])) << std::endl;
                          }
                      }
          }
          os.close();
      }

      void saveEyPhase(SuperArray<maths::Vector<std::complex<double> > > &S, char* FName, int i)
      {
          maths::Vector<int> Pi;
          maths::Vector<std::complex<double> > E;
          std::ofstream os;
          os.open(FName);
          std::complex<double> phase;
          if (S.type == IN_OBJECT)
          {
              os << "%Dimensionen " << S.n[i][0] << "  x  " << S.n[i][1] << "  x  " << S.n[i][2] << std::endl;

              for (int ix = 0; ix < S.n[i][0]; ix++)
                  for (int iy = 0; iy < S.n[i][1]; iy++)
                      for (int iz = 0; iz < S.n[i][2]; iz++)
                      {
                          os << atan2(imag(S.G[i][ix][iy][iz][1]), real(S.G[i][ix][iy][iz][1])) << std::endl;
                      }
          }
          else
          {
              for (int ix = 0; ix < S.nges[0]; ix++)
                  for (int iy = 0; iy < S.nges[1]; iy++)
                      for (int iz = 0; iz < S.nges[2]; iz++)
                      {
                          Pi = S.kugelindex(maths::Vector<int>(ix, iy, iz));
                          if (S.Error != NO_ERRORS) os << 0.0 << std::endl;
                          else
                          {
                              E = S.K[Pi[0]][Pi[1]][Pi[2]];
                              os << atan2(imag(E[1]), real(E[1])) << std::endl;
                          }
                      }
          }
          os.close();
      }

      void saveEzPhase(SuperArray<maths::Vector<std::complex<double> > > &S, char* FName, int i)
      {
          maths::Vector<int> Pi;
          maths::Vector<std::complex<double> > E;
          std::ofstream os;
          os.open(FName);
          std::complex<double> phase;
          if (S.type == IN_OBJECT)
          {
              os << "%Dimensionen " << S.n[i][0] << "  x  " << S.n[i][1] << "  x  " << S.n[i][2] << std::endl;

              for (int ix = 0; ix < S.n[i][0]; ix++)
                  for (int iy = 0; iy < S.n[i][1]; iy++)
                      for (int iz = 0; iz < S.n[i][2]; iz++)
                      {
                          os << atan2(imag(S.G[i][ix][iy][iz][2]), real(S.G[i][ix][iy][iz][2])) << std::endl;
                      }
          }
          else
          {
              for (int ix = 0; ix < S.nges[0]; ix++)
                  for (int iy = 0; iy < S.nges[1]; iy++)
                      for (int iz = 0; iz < S.nges[2]; iz++)
                      {
                          Pi = S.kugelindex(maths::Vector<int>(ix, iy, iz));
                          if (S.Error != NO_ERRORS) os << 0.0 << std::endl;
                          else
                          {
                              E = S.K[Pi[0]][Pi[1]][Pi[2]];
                              os << atan2(imag(E[2]), real(E[2])) << std::endl;
                          }
                      }
          }
          os.close();
      }


      void saveExPol(SuperArray < maths::Vector < std::complex<double> > > &S, char* FName, int i)
      {
          maths::Vector<int> Pi;
          maths::Vector<std::complex<double> > E;
          std::ofstream os;
          os.open(FName);
          if (S.type == IN_OBJECT)
          {
              os << "%Dimensionen " << S.n[i][0] << "  x  " << S.n[i][1] << "  x  " << S.n[i][2] << std::endl;

              for (int ix = 0; ix < S.n[i][0]; ix++)
                  for (int iy = 0; iy < S.n[i][1]; iy++)
                      for (int iz = 0; iz < S.n[i][2]; iz++)
                          os << abs(S.G[i][ix][iy][iz][0]) << std::endl;
          }
          else
          {
              for (int ix = 0; ix < S.nges[0]; ix++)
                  for (int iy = 0; iy < S.nges[1]; iy++)
                      for (int iz = 0; iz < S.nges[2]; iz++)
                      {
                          Pi = S.kugelindex(maths::Vector<int>(ix, iy, iz));
                          if (S.Error != NO_ERRORS) os << 0.0 << std::endl;
                          else
                          {
                              E = S.K[Pi[0]][Pi[1]][Pi[2]];
                              os << abs(E[0]) << std::endl;
                          }
                      }
          }
          os.close();
      }

      void saveEyPol(SuperArray < maths::Vector < std::complex<double> > > &S, char* FName, int i)
      {
          maths::Vector<int> Pi;
          maths::Vector<std::complex<double> > E;
          std::ofstream os;
          os.open(FName);
          if (S.type == IN_OBJECT)
          {
              os << "%Dimension " << S.n[i][0] << "  x  " << S.n[i][1] << "  x  " << S.n[i][2] << std::endl;
              for (int ix = 0; ix < S.n[i][0]; ix++)
                  for (int iy = 0; iy < S.n[i][1]; iy++)
                      for (int iz = 0; iz < S.n[i][2]; iz++)
                          os << abs(S.G[i][ix][iy][iz][1]) << std::endl;
          }
          else
          {
              for (int ix = 0; ix < S.nges[0]; ix++)
                  for (int iy = 0; iy < S.nges[1]; iy++)
                      for (int iz = 0; iz < S.nges[2]; iz++)
                      {
                          Pi = S.kugelindex(maths::Vector<int>(ix, iy, iz));
                          if (S.Error != NO_ERRORS) os << 0.0 << std::endl;
                          else
                          {
                              E = S.K[Pi[0]][Pi[1]][Pi[2]];
                              os << abs(E[1]) << std::endl;
                          }
                      }


          }


          os.close();
      }

      void saveEzPol(SuperArray < maths::Vector < std::complex<double> > > &S, char* FName, int i)
      {
          maths::Vector<int> Pi;
          maths::Vector<std::complex<double> > E;
          std::ofstream os;
          os.open(FName);
          if (S.type == IN_OBJECT)
          {
              os << "%Dimension " << S.n[i][0] << "  x  " << S.n[i][1] << "  x  " << S.n[i][2] << std::endl;


              for (int ix = 0; ix < S.n[i][0]; ix++)
                  for (int iy = 0; iy < S.n[i][1]; iy++)
                      for (int iz = 0; iz < S.n[i][2]; iz++)
                          //   os << abs(G[i][ix][iy][iz][2])/abs(G[i][ix][iy][iz]) << std::endl;
                          // os << real(G[i][ix][iy][iz][2]) << std::endl;
                          os << abs(S.G[i][ix][iy][iz][2]) << std::endl;
          }
          else
          {
              for (int ix = 0; ix < S.nges[0]; ix++)
                  for (int iy = 0; iy < S.nges[1]; iy++)
                      for (int iz = 0; iz < S.nges[2]; iz++)
                      {
                          Pi = S.kugelindex(maths::Vector<int>(ix, iy, iz));
                          if (S.Error != NO_ERRORS) os << 0.0 << std::endl;
                          else
                          {
                              E = S.K[Pi[0]][Pi[1]][Pi[2]];
                              os << abs(E[2]) << std::endl;
                          }
                      }

          }

          os.close();
      }

      void saveabsE(SuperArray < maths::Vector < std::complex<double> > > &S,  std::string FName, int i)
      {
          maths::Vector<std::complex<double> > E;
          std::ofstream os;
          os.open(FName);
          maths::Vector<int> Pi;
          double x;
          if (S.type == IN_OBJECT)
          {
              os << "%Dimension " << S.n[i][0] << "  x  " << S.n[i][1] << "  x  " << S.n[i][2] << std::endl;


              for (int ix = 0; ix < S.n[i][0]; ix++)
                  for (int iy = 0; iy < S.n[i][1]; iy++)
                      for (int iz = 0; iz < S.n[i][2]; iz++)
                      {
                          x = real(S.G[i][ix][iy][iz][0] * conj(S.G[i][ix][iy][iz][0]));
                          x += real(S.G[i][ix][iy][iz][1] * conj(S.G[i][ix][iy][iz][1]));
                          x += real(S.G[i][ix][iy][iz][2] * conj(S.G[i][ix][iy][iz][2]));
                          os << x << std::endl;
                      }
          }
          else
          {
              for (int ix = 0; ix < S.nges[0]; ix++)
                  for (int iy = 0; iy < S.nges[1]; iy++)
                      for (int iz = 0; iz < S.nges[2]; iz++)
                      {
                          Pi = S.kugelindex(maths::Vector<int>(ix, iy, iz));
                          if (S.Error != NO_ERRORS) os << 0.0 << std::endl;
                          else
                          {
                              E = S.K[Pi[0]][Pi[1]][Pi[2]];
                              os << real(E * conj(E)) << std::endl;
                          }
                      }
          }

          os.close();
      }

      void saveabsEbin(SuperArray < maths::Vector < std::complex<double> > >& S, std::string FName, int i)
      {
          maths::Vector<std::complex<double> > E;
          std::ofstream os;
          os.open(FName,std::ios_base::binary);
          maths::Vector<int> Pi;
          double x;
          double erg;
          if (S.type == IN_OBJECT)
          {              
               for (int ix = 0; ix < S.n[i][0]; ix++)
                  for (int iy = 0; iy < S.n[i][1]; iy++)
                      for (int iz = 0; iz < S.n[i][2]; iz++)
                      {
                          x = real(S.G[i][ix][iy][iz][0] * conj(S.G[i][ix][iy][iz][0]));
                          x += real(S.G[i][ix][iy][iz][1] * conj(S.G[i][ix][iy][iz][1]));
                          x += real(S.G[i][ix][iy][iz][2] * conj(S.G[i][ix][iy][iz][2]));
                          os.write((char*)&x, sizeof(x));                          
                      }
          }
          else
          {
              for (int ix = 0; ix < S.nges[0]; ix++)
                  for (int iy = 0; iy < S.nges[1]; iy++)
                      for (int iz = 0; iz < S.nges[2]; iz++)
                      {
                          Pi = S.kugelindex(maths::Vector<int>(ix, iy, iz));
                          if (S.Error != NO_ERRORS) os << 0.0 << std::endl;
                          else
                          {
                              E = S.K[Pi[0]][Pi[1]][Pi[2]];
                              erg = real(E * conj(E));
                              os.write((char*)&erg, sizeof(erg));                              
                          }
                      }
          }

          os.close();
      }

      void saveFullE(SuperArray < maths::Vector < std::complex<double> > > &S, std::string FName, int i)
      {
          maths::Vector<std::complex<double> > E;
          std::ofstream os;
          os.open(FName);
          maths::Vector<int> Pi;
          
          if (S.type == IN_OBJECT)
          {
              os << "%Dimension " << S.n[i][0] << "  x  " << S.n[i][1] << "  x  " << S.n[i][2] << std::endl;


              for (int ix = 0; ix < S.n[i][0]; ix++)
                  for (int iy = 0; iy < S.n[i][1]; iy++)
                      for (int iz = 0; iz < S.n[i][2]; iz++)
                      {
                          os << real(S.G[i][ix][iy][iz][0]) << "\t" << imag(S.G[i][ix][iy][iz][0]) << "\t";
                          os << real(S.G[i][ix][iy][iz][1]) << "\t" << imag(S.G[i][ix][iy][iz][1]) << "\t";
                          os << real(S.G[i][ix][iy][iz][2]) << "\t" << imag(S.G[i][ix][iy][iz][2]) << "\t" << std::endl;
                      }
          }
          else
          {
              for (int ix = 0; ix < S.nges[0]; ix++)
                  for (int iy = 0; iy < S.nges[1]; iy++)
                      for (int iz = 0; iz < S.nges[2]; iz++)
                      {
                          Pi = S.kugelindex(maths::Vector<int>(ix, iy, iz));
                          if (S.Error != NO_ERRORS) os << 0.0 << std::endl;
                          else
                          {
                              E = S.K[Pi[0]][Pi[1]][Pi[2]];
                              os << real(E[0]) << "\t" << imag(E[0]) << "\t";
                              os << real(E[1]) << "\t" << imag(E[1]) << "\t";
                              os << real(E[2]) << "\t" << imag(E[2]) << "\t" << std::endl;
                          }
                      }
          }

          os.close();
      }




      double sumabs(const SuperArray<maths::Vector<std::complex<double> > >& S)
      {
          double Erg, h;
          int anzx2 = S.nges[0] / 2;
          Erg = 0.0;
          if (S.type == IN_OBJECT)
          {
              for (int i = 0; i < S.numObjs; i++)
                  for (int ix = 0; ix < S.n[i][0]; ix++)
                      for (int iy = 0; iy < S.n[i][1]; iy++)
                          for (int iz = 0; iz < S.n[i][2]; iz++)
                          {
                              h = abs(S.G[i][ix][iy][iz]);
                              Erg += h;
                          }
          }
          else
          {
              for (int k = 0; k < anzx2; k++)
              {
                  for (int l = 0; l < S.ywerte[k]; l++)
                  {

                      for (int m = 0; m < S.zwerte[k][l]; m++)
                      {
                          h = abs(S.K[anzx2 - 1 - k][S.ywerte[k] - 1 - l][m]);
                          Erg += h;
                          h = abs(S.K[anzx2 + k][S.ywerte[k] - 1 - l][m]);
                          Erg += h;
                          h = abs(S.K[anzx2 - 1 - k][S.ywerte[k] + l][m]);
                          Erg += h;
                          h = abs(S.K[anzx2 + k][S.ywerte[k] + l][m]);
                          Erg += h;
                          h = abs(S.K[anzx2 - 1 - k][S.ywerte[k] - 1 - l][2 * S.zwerte[k][l] - 1 - m]);
                          Erg += h;
                          h = abs(S.K[anzx2 + k][S.ywerte[k] - 1 - l][2 * S.zwerte[k][l] - 1 - m]);
                          Erg += h;
                          h = abs(S.K[anzx2 - 1 - k][S.ywerte[k] + l][2 * S.zwerte[k][l] - 1 - m]);
                          Erg += h;
                          h = abs(S.K[anzx2 + k][S.ywerte[k] + l][2 * S.zwerte[k][l] - 1 - m]);
                          Erg += h;
                      }
                  }
              }

          }

          return Erg;
      }

      double sumabs2(const SuperArray<maths::Vector<std::complex<double> > >& S)
      {
          double Erg, h;
          int anzx2 = S.nges[0] / 2;
          Erg = 0.0;
          if (S.type == IN_OBJECT)
          {
              for (int i = 0; i < S.numObjs; i++)
                  for (int ix = 0; ix < S.n[i][0]; ix++)
                      for (int iy = 0; iy < S.n[i][1]; iy++)
                          for (int iz = 0; iz < S.n[i][2]; iz++)
                          {
                              h = abs(S.G[i][ix][iy][iz]);
                              Erg += h * h;
                          }
          }
          else
          {
              for (int k = 0; k < anzx2; k++)
              {
                  for (int l = 0; l < S.ywerte[k]; l++)
                  {

                      for (int m = 0; m < S.zwerte[k][l]; m++)
                      {
                          h = abs(S.K[anzx2 - 1 - k][S.ywerte[k] - 1 - l][m]);
                          Erg += h * h;
                          h = abs(S.K[anzx2 + k][S.ywerte[k] - 1 - l][m]);
                          Erg += h * h;
                          h = abs(S.K[anzx2 - 1 - k][S.ywerte[k] + l][m]);
                          Erg += h * h;
                          h = abs(S.K[anzx2 + k][S.ywerte[k] + l][m]);
                          Erg += h * h;
                          h = abs(S.K[anzx2 - 1 - k][S.ywerte[k] - 1 - l][2 * S.zwerte[k][l] - 1 - m]);
                          Erg += h * h;
                          h = abs(S.K[anzx2 + k][S.ywerte[k] - 1 - l][2 * S.zwerte[k][l] - 1 - m]);
                          Erg += h * h;
                          h = abs(S.K[anzx2 - 1 - k][S.ywerte[k] + l][2 * S.zwerte[k][l] - 1 - m]);
                          Erg += h * h;
                          h = abs(S.K[anzx2 + k][S.ywerte[k] + l][2 * S.zwerte[k][l] - 1 - m]);
                          Erg += h * h;
                      }
                  }
              }

          }
          return Erg;
      }

      double abs2sum(const SuperArray<maths::Vector<std::complex<double> > >& S)
      {
          maths::Vector<std::complex<double> > h;
          double Erg;
          int anzx2 = S.nges[0] / 2;
          if (S.type == IN_OBJECT)
          {
              for (int i = 0; i < S.numObjs; i++)
                  for (int ix = 0; ix < S.n[i][0]; ix++)
                      for (int iy = 0; iy < S.n[i][1]; iy++)
                          for (int iz = 0; iz < S.n[i][2]; iz++)
                          {
                              h += S.G[i][ix][iy][iz];
                          }
          }
          else
          {
              for (int k = 0; k < anzx2; k++)
              {
                  for (int l = 0; l < S.ywerte[k]; l++)
                  {

                      for (int m = 0; m < S.zwerte[k][l]; m++)
                      {
                          h += S.K[anzx2 - 1 - k][S.ywerte[k] - 1 - l][m];
                          h += S.K[anzx2 + k][S.ywerte[k] - 1 - l][m];
                          h += S.K[anzx2 - 1 - k][S.ywerte[k] + l][m];
                          h += S.K[anzx2 + k][S.ywerte[k] + l][m];
                          h += S.K[anzx2 - 1 - k][S.ywerte[k] - 1 - l][2 * S.zwerte[k][l] - 1 - m];
                          h += S.K[anzx2 + k][S.ywerte[k] - 1 - l][2 * S.zwerte[k][l] - 1 - m];
                          h += S.K[anzx2 - 1 - k][S.ywerte[k] + l][2 * S.zwerte[k][l] - 1 - m];
                          h += S.K[anzx2 + k][S.ywerte[k] + l][2 * S.zwerte[k][l] - 1 - m];
                      }
                  }
              }

          }
          Erg = abs(h);
          return Erg * Erg;
      }

      void save(SuperArray<GOAT::raytracing::gridEntry > S, std::string FName)
      {
          std::ofstream os(FName);
          os << "% Superarray: " << S.numObjs << " objects" << std::endl;
          int vecnumObj = S.G.size();
          int nx, ny, nz;
          for (int i = 0; i < vecnumObj; i++)
          {
              os << "% ------------ object no. " << i << " ----------------" << std::endl;              
              nx = S.G[i].size();              
              for (int ix = 0; ix < nx; ix++)
              {
                  ny = S.G[i][ix].size();
                  for (int iy = 0; iy < ny; iy++)
                  {
                      nz = S.G[i][ix][iy].size();
                      os << "% " << nx << " x " << ny << " x " << nz << std::endl;
                      for (int iz = 0; iz < nz; iz++)
                      {
                          os << "% " << S.G[i][ix][iy][iz].E << "  ";
                          os << "% steps: " << S.G[i][ix][iy][iz].step.size() << std::endl;
                          for (auto se : S.G[i][ix][iy][iz].step)
                              os << ix << "\t" << iy << "\t" << iz << "\t" << se.l << "\t" << se.matIndex << std::endl;                          
                      }
                  }

              }

          }
          os.close();
      }

      void save(SuperArray<std::vector<GOAT::raytracing::gridEntry > > S, std::string FName)
      {
          std::ofstream os(FName);
          os << "% Superarray: " << S.numObjs << " objects" << std::endl;
          int vecnumObj = S.G.size();
          std::cout << "vecnumObj=" << vecnumObj << std::endl; 
          int nx, ny, nz;
          for (int i = 0; i < vecnumObj; i++)
          {
              os << "% ------------ object no. " << i << " ----------------" << std::endl;
              nx = S.G[i].size();
              for (int ix = 0; ix < nx; ix++)
              {
                  ny = S.G[i][ix].size();
                  for (int iy = 0; iy < ny; iy++)
                  {
                      nz = S.G[i][ix][iy].size();                      
                      for (int iz = 0; iz < nz; iz++)
                      {
                          os << "% Array entry " << ix << "," << iy << "," << iz << std::endl;
                          for (auto ge : S.G[i][ix][iy][iz])
                          {
                              os << "% E=" << ge.E << "\t" << ge.step.size() << " steps." << std::endl;
                              for (auto se : ge.step)
                                  os << se.l << "\t" << se.matIndex << std::endl;
                              
                          }
                      }
                  }
              }
              os.close();
          }

      }
      template<> bool SuperArray<GOAT::maths::Vector<std::complex<double> > >::write(std::string fname)
      {
          std::ofstream os;
          os.open(fname, std::ios::out | std::ios::binary);
          if (os.good())
          {
              os.write((char*)&numObjs, sizeof(numObjs)); // number of objects (objects itself were not saved yet !)
              os.write((char*)&type, sizeof(type)); 
              
              if (type == IN_OBJECT)
              {
                  bool active;
                  for (int i = 0; i < numObjs; i++)
                  {
                      active = Obj[i]->isActive();
                      os.write((char*)&active, sizeof(active));
                      if (active)
                      {
                          os.write((char*)&n[i][0], sizeof(n[i][0]));
                          os.write((char*)&n[i][1], sizeof(n[i][1]));
                          os.write((char*)&n[i][2], sizeof(n[i][2]));
                          for (int ix = 0; ix < n[i][0]; ix++)
                              for (int iy = 0; iy < n[i][1]; iy++)
                                  for (int iz = 0; iz < n[i][2]; iz++)
                                      os.write((char*)&G[i][ix][iy][iz], sizeof(G[i][ix][iy][iz]));
                      }

                  }

              }
              else
              {
                  os.write((char*)&nges[0], sizeof(nges[0]));
                  os.write((char*)&nges[1], sizeof(nges[1]));
                  os.write((char*)&nges[2], sizeof(nges[2]));

                  int anzx, anzx2;
                  anzx = nges[0];
                  anzx2 = nges[0] / 2;
                  for (int i = 0; i < ywerte.size(); i++)
                      os.write((char*)&ywerte[i], sizeof(ywerte[i]));

                  for (int i = 0; i < zwerte.size(); i++)
                      for (int j = 0; j < zwerte[0].size(); j++)
                          os.write((char*)&zwerte[i][j], sizeof(zwerte[i][j]));


                  for (int k = 0; k < anzx2; k++)
                      for (int l = 0; l < ywerte[k]; l++)
                          for (int m = 0; m < zwerte[k][l]; m++)
                          {
                              os.write((char*)&K[anzx2 - 1 - k][ywerte[k] - 1 - l][zwerte[k][l] - 1 - m], sizeof(K[anzx2 - 1 - k][ywerte[k] - 1 - l][zwerte[k][l] - 1 - m]));
                              os.write((char*)&K[anzx2 + k][ywerte[k] - 1 - l][zwerte[k][l] - 1 - m], sizeof(K[anzx2 + k][ywerte[k] - 1 - l][zwerte[k][l] - 1 - m]));
                              os.write((char*)&K[anzx2 - 1 - k][ywerte[k] + l][zwerte[k][l] - 1 - m], sizeof(K[anzx2 - 1 - k][ywerte[k] + l][zwerte[k][l] - 1 - m]));
                              os.write((char*)&K[anzx2 + k][ywerte[k] + l][zwerte[k][l] - 1 - m], sizeof(K[anzx2 + k][ywerte[k] + l][zwerte[k][l] - 1 - m]));
                              os.write((char*)&K[anzx2 - 1 - k][ywerte[k] - 1 - l][zwerte[k][l] + m], sizeof(K[anzx2 - 1 - k][ywerte[k] - 1 - l][zwerte[k][l] + m]));
                              os.write((char*)&K[anzx2 + k][ywerte[k] - 1 - l][zwerte[k][l] + m], sizeof(K[anzx2 + k][ywerte[k] - 1 - l][zwerte[k][l] + m]));
                              os.write((char*)&K[anzx2 - 1 - k][ywerte[k] + l][zwerte[k][l] + m], sizeof(K[anzx2 - 1 - k][ywerte[k] + l][zwerte[k][l] + m]));
                              os.write((char*)&K[anzx2 + k][ywerte[k] + l][zwerte[k][l] + m], sizeof(K[anzx2 + k][ywerte[k] + l][zwerte[k][l] + m]));
                          }
              }
              return true;
          }
          return false;
      }

      template<> bool SuperArray<GOAT::maths::Vector<std::complex<double> > >::read(std::string fname)
      {
          std::ifstream is;
          is.open(fname, std::ios::out | std::ios::binary);
          if (is.good())
          {
              is.read((char*)&numObjs, sizeof(numObjs)); // number of objects (objects itself were not saved yet !)
              is.read((char*)&type, sizeof(type));

              if (type == IN_OBJECT)
              {
                  GOAT::maths::Vector<int> hiv;
                  bool active;
                  int hi;
                  for (int i = 0; i < numObjs; i++)
                  {
                      is.read((char*)&active, sizeof(active));
                      is.read((char*)&hi, sizeof(hi));
                      hiv[0] = hi;
                      is.read((char*)&hi, sizeof(hi));
                      hiv[1] = hi;
                      is.read((char*)&hi, sizeof(hi));
                      hiv[2] = hi;
                      n.push_back(hiv);
                       
                      std::complex<double> cx, cy, cz;
                      G.resize(numObjs);
                      
                      for (int ix = 0; ix < n[i][0]; ix++)
                          for (int iy = 0; iy < n[i][1]; iy++)
                              for (int iz = 0; iz < n[i][2]; iz++)
                              {
                                  is.read((char*)&cx, sizeof(cx));
                                  is.read((char*)&cy, sizeof(cy));
                                  is.read((char*)&cz, sizeof(cz));

                                  is.read((char*)&G[i][ix][iy][iz], sizeof(G[i][ix][iy][iz]));
                              }
                  }

              }
              else
              {
                  int intx, inty, intz;
                  is.read((char*)&intx, sizeof(intx));
                  is.read((char*)&inty, sizeof(inty));
                  is.read((char*)&intz, sizeof(intz));
                  nges = GOAT::maths::Vector<int>(intx, inty, intz);

                  int anzx, anzx2;
                  anzx = nges[0];
                  anzx2 = nges[0] / 2;
                  for (int i = 0; i < ywerte.size(); i++)
                      is.read((char*)&ywerte[i], sizeof(ywerte[i]));

                  for (int i = 0; i < zwerte.size(); i++)
                      for (int j = 0; j < zwerte[0].size(); j++)
                          is.read((char*)&zwerte[i][j], sizeof(zwerte[i][j]));

                  for (int k = 0; k < anzx2; k++)
                      for (int l = 0; l < ywerte[k]; l++)
                          for (int m = 0; m < zwerte[k][l]; m++)
                          {
                              is.read((char*)&K[anzx2 - 1 - k][ywerte[k] - 1 - l][zwerte[k][l] - 1 - m], sizeof(K[anzx2 - 1 - k][ywerte[k] - 1 - l][zwerte[k][l] - 1 - m]));
                              is.read((char*)&K[anzx2 + k][ywerte[k] - 1 - l][zwerte[k][l] - 1 - m], sizeof(K[anzx2 + k][ywerte[k] - 1 - l][zwerte[k][l] - 1 - m]));
                              is.read((char*)&K[anzx2 - 1 - k][ywerte[k] + l][zwerte[k][l] - 1 - m], sizeof(K[anzx2 - 1 - k][ywerte[k] + l][zwerte[k][l] - 1 - m]));
                              is.read((char*)&K[anzx2 + k][ywerte[k] + l][zwerte[k][l] - 1 - m], sizeof(K[anzx2 + k][ywerte[k] + l][zwerte[k][l] - 1 - m]));
                              is.read((char*)&K[anzx2 - 1 - k][ywerte[k] - 1 - l][zwerte[k][l] + m], sizeof(K[anzx2 - 1 - k][ywerte[k] - 1 - l][zwerte[k][l] + m]));
                              is.read((char*)&K[anzx2 + k][ywerte[k] - 1 - l][zwerte[k][l] + m], sizeof(K[anzx2 + k][ywerte[k] - 1 - l][zwerte[k][l] + m]));
                              is.read((char*)&K[anzx2 - 1 - k][ywerte[k] + l][zwerte[k][l] + m], sizeof(K[anzx2 - 1 - k][ywerte[k] + l][zwerte[k][l] + m]));
                              is.read((char*)&K[anzx2 + k][ywerte[k] + l][zwerte[k][l] + m], sizeof(K[anzx2 + k][ywerte[k] + l][zwerte[k][l] + m]));
                          }
              }
              return true;
          }
          return false;
      }
  }
}
