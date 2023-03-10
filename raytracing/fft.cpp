#include "fft.h"
namespace GOAT
{
    namespace raytracing
    {
        Trafo::Trafo()
        {
            tp.nR = INEL_MAX_NREFLEX;
        }


        Trafo::Trafo(TrafoParms tp)
        {
            this->tp = tp;
            omegastart = 2.0 * M_PI * C_LIGHT_MU / tp.lstop;
            omegastop = 2.0 * M_PI * C_LIGHT_MU / tp.lstart;
            twoSigma2 = tp.dt * tp.dt / (4.0 * M_LN2);
        }

        void Trafo::initResult(SuperArray<maths::Vector<std::complex<double> >>& SA)
        {
            SAres.clear();
            SAres = SuperArray<maths::Vector<std::complex<double> > >(SA.r0, SA.nges[0], SA.nges[1], SA.nges[2], SA.Obj, SA.numObjs);
        }

        void Trafo::initResult(double r0, int nx, int ny, int nz, ObjectShape** Obj, int numObjs)
        {
            SAres.clear();
            SAres = SuperArray<maths::Vector<std::complex<double> > >(r0, nx, ny, nz, Obj, numObjs);
        }


      
        void Trafo::calc(std::vector<SuperArray <std::vector<gridEntry> > >& SA, double t)
        {
            maths::Vector<std::complex<double> > E;
            std::complex<double> phase;
            int nsteps;
            double dwvl; // step size for the integration inside the subrange in wavelength
            double Dwvl; // width of one subrange
            double domega; // step size for the integration inside the subrange in angular frequency
            double omega; // angular frequency
            double wvl;

            initResult(SA[0].r0, SA[0].nges[0], SA[0].nges[1], SA[0].nges[2], SA[0].Obj, SA[0].numObjs);

            SAres.fill(maths::czero); // empty the whole result array
            Dwvl = (tp.lstop - tp.lstart) / (double)tp.nI;
            /* -----  Loops over positions ------ */
            for (int iR = 0; iR < tp.nR; iR++) // loop over reflection order
                for (int i = 0; i < SA[iR].numObjs; i++) // loop over object number (i.e. over Sub-Array in SuperArray)
                    for (int ix = 0; ix < SA[iR].n[i][0]; ix++) // loops over x-,y- and z- indices
                        for (int iy = 0; iy < SA[iR].n[i][1]; iy++)
                            for (int iz = 0; iz < SA[iR].n[i][2]; iz++)
                            {
                                SAres.G[i][ix][iy][iz] = integrate(t, SA[iR].G[i][ix][iy][iz], omegastart, omegastop, tp.nS);
                            }


        }

        void Trafo::setRefractiveIndexFunctions(std::vector<std::function<std::complex<double>(double) > > nList)
        {
            tp.nList = nList;
        }

        void Trafo::createLTexpo()
        {
            double dl = (tp.lstop - tp.lstart) / (double)(tp.nS - tp.nI - 1);
            int N = tp.nS * tp.nI;
            double l;
            for (int i = 0; i < N; i++)
            {
                l = tp.lstart + i * dl;
                freq.push_back(2.0 * M_PI * C_LIGHT_MU / l);
            }
        }

        std::complex<double> Trafo::calcPhase(std::vector<stepEntry> steps, double k0)
        {
            std::complex<double> sum = 0;
            for (stepEntry se : steps)
                sum += k0 * tp.nList[se.matIndex](2.0 * M_PI / k0) * se.l;
            return sum;
        }



        GOAT::maths::Vector<std::complex<double> >  Trafo::integrate(double t, std::vector<gridEntry> vge, double omegastart, double omegastop, int nsteps)
        {
            double k0;
        GOAT:maths::Vector<std::complex<double> > E;
            std::complex<double> phase;
            double omega0 = 2.0 * M_PI * C_LIGHT_MU / tp.wvl;
            double omega;
            double domega = (omegastop - omegastart) / ((double)nsteps - 1.0);
            double weight;
            for (auto ge : vge) // Loop over all rays which hit the cell
            {
                for (int iomega = 0; iomega < nsteps; iomega++)
                {
                    omega = omegastart + iomega * domega;
                    k0 = omega / C_LIGHT_MU;
                    phase = exp(I * calcPhase(ge.step, k0));
                    phase *= exp(-I * omega * t);
                    weight = exp(-twoSigma2 * (omega - omega0) * (omega - omega0));
                    E += weight * ge.E * phase;
                }
            }
            return E;
        }




        /*

        unsigned int bitReverse(unsigned int x, int log2n)  ///< function used for FFT algorithm
        {
            int n = 0;
            int mask = 0x1;

            for (int i = 0; i < log2n; i++) {
                n <<= 1;
                n |= (x & 1);
                x >>= 1;
            }
            return n;
        }

        std::vector<unsigned int> bitReverseLT(int log2n)
        {
            std::vector<unsigned int> LT;
            int n = 1 << log2n;
            for (int x = 0; x < n; x++)
            {
                LT.push_back(bitReverse(x, log2n));
            }
            return LT;
        }


        template<class Iter_T>
        void fft(Iter_T a, Iter_T b, int log2n)
        {
            typedef typename iterator_traits<Iter_T>::value_type complex;
            const complex J(0, 1);
            int n = 1 << log2n;
            for (unsigned int i = 0; i < n; ++i) {
                b[bitReverse(i, log2n)] = a[i];
            }
            for (int s = 1; s <= log2n; ++s) {
                int m = 1 << s;
                int m2 = m >> 1;
                complex w(1, 0);
                complex wm = exp(-J * (M_PI / m2));
                for (int j = 0; j < m2; ++j) {
                    for (int k = j; k < n; k += m) {
                        complex t = w * b[k + m2];
                        complex u = b[k];
                        b[k] = u + t;
                        b[k + m2] = u - t;
                    }
                    w *= wm;
                }
            }
        }

        template<class Iter_T>
        void fft(Iter_T a, Iter_T b, int log2n, std::vector<unsigned int> LT)
        {
            typedef typename iterator_traits<Iter_T>::value_type complex;
            const complex J(0, 1);
            int n = 1 << log2n;
            for (unsigned int i = 0; i < n; ++i) {
                b[LT[i]] = a[i];
            }
            for (int s = 1; s <= log2n; ++s) {
                int m = 1 << s;
                int m2 = m >> 1;
                complex w(1, 0);
                complex wm = exp(-J * (M_PI / m2));
                for (int j = 0; j < m2; ++j) {
                    for (int k = j; k < n; k += m) {
                        complex t = w * b[k + m2];
                        complex u = b[k];
                        b[k] = u + t;
                        b[k + m2] = u - t;
                    }
                    w *= wm;
                }
            }
        }
        */
    }
}