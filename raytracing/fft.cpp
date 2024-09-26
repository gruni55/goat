#include "fft.h"
#include <chrono>
#include <math.h>
#include <omp.h>

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
            twoSigma2 = tp.dt * tp.dt / (4.0 * M_LN2);
            
            sigma2 = tp.dt * tp.dt / (4.0 * M_LN2);
            prefactor = 1.0 / sqrt(twoSigma2 * 2.0 * M_PI);            
        }

        void Trafo::clear()
        {
            SAres.clear();
            D.clear();
        }

        void Trafo::empty()
        {

            SAres.fill(maths::czero); // empty the whole result array
        }

        void Trafo::initResult(SuperArray<maths::Vector<std::complex<double> >>& SA)
        {           
            clear();
            SAres = SuperArray<maths::Vector<std::complex<double> > >(SA.r0, SA.nges[0], SA.nges[1], SA.nges[2], SA.Obj, SA.numObjs);             
            for (int i=0; i<numObjs; i++)
              D.push_back(0);
            this->numObjs=numObjs; 
        }

        void Trafo::initResult(double r0, INDEX_TYPE nx, INDEX_TYPE ny, INDEX_TYPE nz, ObjectShape** Obj, int numObjs)
        {
            clear();
            SAres = SuperArray<maths::Vector<std::complex<double> > >(r0, nx, ny, nz, Obj, numObjs);
            for (int i=0; i<numObjs; i++)
              D.push_back(0);
            this->numObjs=numObjs;  
        }

       

        GOAT::maths::Vector<std::complex<double> >  Trafo::integrate(double t, std::vector<gridEntry> vge, double omegastart, double omegastop)
        {
        GOAT:maths::Vector<std::complex<double> > E(0,0,0);
            if (vge.size() > 0)
            {
                double k0;
                double dw;
                std::complex<double> phase;
                std::complex<double> phase1;
                std::complex<double> phase2;

                double omega;
                double domega = 2.0 * M_PI / fabs(tref - t);
                double Domega = omegastop - omegastart;
// //                 double dw;
                int nsteps = tp.nS;
                    domega = Domega / ((double)(nsteps));
                double weight;

                for (int iomega = 0; iomega < nsteps; iomega++)
                {
                    omega = omegastart + iomega * domega;
                    
                    dw = (omega - tp.omega0);
                    double ws = dw * dw * sigma2 / 2.0;
                    weight = exp(-ws); // this weighting is due to the Fourier transform of the gaussian (temporal) pulse shape
                    k0 = omega / C_LIGHT_MU_FS;
                    setCurrNList(2.0 * M_PI / k0);
                    for (auto ge : vge) // Loop over all rays which hit the cell
                    {
                        phase1 = exp(I * calcPhase(ge.step, k0));
                        phase = phase1 * exp(I * (-omega * t));
                        E += weight * phase * ge.E;
                    }
                }
            }
            return E;
        }

      
        void Trafo::calc(std::vector < std::vector<SuperArray <std::vector<gridEntry> > > >& SA, double t)
        {
            
            maths::Vector<std::complex<double> > E; 
            std::complex<double> phase;
            int nsteps;
            double dwvl;   // step size for the integration inside the subrange in wavelength
            double Dwvl;   // width of one subrange
            double domega; // step size for the integration inside the subrange in angular frequency
            double omega;  // angular frequency
            double omegaStart, omegaEnd;
            double wvl;
            double lambda;
            double omega0 = 2.0 * M_PI * C_LIGHT_MU_FS / tp.wvl;
            initResult(SA[0][0].r0, SA[0][0].nges[0], SA[0][0].nges[1], SA[0][0].nges[2], SA[0][0].Obj, SA[0][0].numObjs);
            double Sigma = 2.3548 / tp.dt;
            double Domega = 4.0 * Sigma;            

            domega = Domega / (double)tp.nI;
                        
             for (int iOmega = 0; iOmega < tp.nI; iOmega++) // loop over the spectral ranges
            {               
                 auto start = std::chrono::high_resolution_clock::now();
                 omega = tp.omegaStart + iOmega * domega - domega / 2.0;
                omegaStart = omega - domega / 2.0;
                omegaEnd = omega + domega / 2.0;
                
                for (int iR = 0; iR < tp.nR; iR++)   // loop over reflection order
                    for (int i = 0; i < SA[iOmega][iR].numObjs; i++)        // loop over object number (i.e. over Sub-Array in SuperArray)
                        if (SAres.Obj[i]->isActive())
                        {
// #if defined(_OPENMP)
                           //omp_set_num_threads(tp.number_of_threads);
                            // std::cout << "number of threads used:" << omp_get_num_threads() << std::endl;
                             #pragma omp parallel for


// #endif
                            // #pragma omp parallel for num_threads (7)
                            for (int ix = 0; ix < SA[iOmega][iR].n[i][0]; ix++) // loops over x-,y- and z- indices
                            {
                                std::cout << ix << std::endl << std::flush;
                                for (int iy = 0; iy < SA[iOmega][iR].n[i][1]; iy++)
                                    for (int iz = 0; iz < SA[iOmega][iR].n[i][2]; iz++)
                                    {
                                        /*int ix = 0;
                                        int iy = 0;
                                        int iz = 0;*/

                                        SAres.G[i][ix][iy][iz] += integrate(t, SA[iOmega][iR].G[i][ix][iy][iz], omegaStart, omegaEnd);
                                       
                                    }
                                                    std::cout  << ix << "  " << GOAT::maths::abs2(SAres.G[i][ix][2][2]) << std::endl;
                            }
                        }
                auto end = std::chrono::high_resolution_clock::now();
                std::cout << "%integration time: " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000000 << " s" << std::endl;
            }

        }

        void Trafo::calc(std::vector<SuperArray <std::vector<gridEntry> > > & SA, double omegaStart, double omegaEnd, double t, bool do_clear)
        {
              //  auto start = std::chrono::high_resolution_clock::now();
            int counter=0;
  //          omp_set_num_threads(tp.number_of_threads);

//            std::cout << "number of threads used:" << omp_get_num_threads() << "(" << tp.number_of_threads << ")" << std::endl;
            maths::Vector<std::complex<double> > hint;
                for (int iR = 0; iR < tp.nR; iR++)   // loop over reflection order
                    for (int i = 0; i < SA[iR].numObjs; i++)        // loop over object number (i.e. over Sub-Array in SuperArray)
                        if (SAres.Obj[i]->isActive())
                        {
                         D[i]=0;  
//  #pragma omp parallel for num_threads(6)
#pragma omp parallel for
                            for (INDEX_TYPE ix = 0; ix < SA[iR].n[i][0]; ix++) // loops over x-,y- and z- indices
                            {                    
                                // std::cout << ix << std::endl << std::flush;
                                for (INDEX_TYPE iy = 0; iy < SA[iR].n[i][1]; iy++)
                                    for (INDEX_TYPE iz= 0; iz < SA[iR].n[i][2]; iz++)
                                    {
                                      //  std::cout << ix << "\t" << iy << "\t" << iz << std::endl << std::flush;
                                        if (SA[iR].G[i][ix][iy][iz].size() > 0)
                                        {
                                            maths::Vector<std::complex<double>> hint=integrate(t, SA[iR].G[i][ix][iy][iz], omegaStart, omegaEnd);
                                            // double absold = abs(SAres.G[i][ix][iy][iz]);
                                            //if (absold > 1E-20) { D += abs(hint) / absold; }
                                            SAres.G[i][ix][iy][iz] += hint;
                                            D[i]+=abs2(SAres.G[i][ix][iy][iz]);
                                            counter++;
                                            if (do_clear) SA[iR].G[i][ix][iy][iz].clear();
                                        }
                                    }
                             }
                          D[i]/=counter; 
                        }
                 // D/=counter;
              //  auto end = std::chrono::high_resolution_clock::now();
              //  std::cout << "%integration time: " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000000 << " s" << std::endl;
        }

        void Trafo::setRefractiveIndexFunctions(std::vector<std::function<std::complex<double>(double) > > nList)
        {
            tp.nList = nList;
            this->nList = nList;
            currNList.resize(nList.size());
        }


        void Trafo::setTrafoParms(TrafoParms trafoparms)
        {
            tp = trafoparms;
            
            //double sigma = 2.0 * sqrt(M_LN2) / tp.dt;
            sigma2 = tp.dt*tp.dt/(4.0*M_LN2);            
        }
/*
        void Trafo::createLTexpo()
        {
            double dl = (tp.lstop - tp.lstart) / (double)(tp.nS - tp.nI - 1);
            int N = tp.nS * tp.nI;
            double l;s
            for (int i = 0; i < N; i++)
            {
                l = tp.lstart + i * dl;
                freq.push_back(2.0 * M_PI * C_LIGHT_MU_FS / l);
            }
        }
        */
       inline std::complex<double> Trafo::calcPhase(std::vector<stepEntry> steps, double k0)
        {
            std::complex<double> sum = 0;
            for (stepEntry se : steps)
               sum += k0 * currNList[se.matIndex] * se.l;
            return sum;
        }





        void Trafo::setReferenceTime(double tref) { this->tref = tref; }

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
