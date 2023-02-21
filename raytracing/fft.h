#pragma once
#include <functional>
#include <vector>
#include "superarray.h"
#include "raytrace_usp.h"

namespace GOAT
{
	namespace raytracing
	{
#ifndef I 
#define I std::complex<double> (0.0,1.0)
#endif
        constexpr int N_INTEGRAL = 100;
        /**
        * Structure, which acts as a container for all informations needed to process the calculation
        * @var TrafoParms::lstart 
        * defines the start of the integration range
        * @var TrafoParms::lstop 
        *  defines end of the integration range
        *  @var TrafoParms::nI 
        *  number of subranges 
        *  @var TrafoParms::nS
        * number of steps to consider within a subrange
        * @var TrafoParms::nR
        * number of reflections considered in each SuperArray
        */
        typedef struct TrafoParms
        {
            int nI;
            int nS=1000;
            int nR= INEL_MAX_NREFLEX;
            double lstart;
            double lstop;
            std::vector<std::function<std::complex<double>(double) > > nList;
        };

        constexpr double C_LIGHT = 299792458; ///< speed of light in m*s^-1
        /**
        * @brief This class calculates the time dependence of a field which calculated before
        * The calculation of the time dependence is made in two steps. In the first step, the fields are calculate 
        * for different wavelengths. In this class, the second step will be made. Here the time dependence is calculated
        * with a kind of Fourier transformation. The class requires an array of SuperArray elements which contain the 
        * fields for the different wavelengths. The process assumes, that the wavelength range can be subdivided into 
        * subranges where the paths of the rays remain nearly the same. Since the phase is much more sensitive to changes
        * of the refractive index, only the phase change within this subranges is considered. The structure TrafoParms 
        * carries all information needed for the calculation
        */
        template <class T> class Trafo
        {
        public:
            Trafo();
            Trafo(TrafoParms);  ///< Contructor for initialization          
          //  ~Trafo();
            /**
            * Perform a calculation at the time t
            * @param SA Array of SuperArray elements. Contains the fields for the different wavelengths  
            * defined by the parameters in the corresponding TrafoParms structure, given in the constructor
            */
            void calc(std::vector<SuperArray <std::vector<gridEntry> > >& SA, double t); 
            SuperArray<maths::Vector<std::complex<double> > >SAres; ///< Container for the last result       


        private:
            void initResult(SuperArray<T>& SA);
            void initResult(double r0, int nx, int ny, int nz, ObjectShape** Obj, int numObjs);
            void createLTexpo();
            std::complex<double> calcPhase(std::vector<stepEntry> steps, double k0);
            GOAT::maths::Vector<std::complex<double> > calcOne(std::vector<stepEntry> steps, double t);
            GOAT::maths::Vector<std::complex<double> > integrate(double t, std::vector<gridEntry> ge, double omegastart, double omegastop,  int nsteps=100);

            std::vector<double> freq;
            int nI; ///< Number of wavelength intervalls
            int nS; ///< Number of steps per intervall
            double lstart; ///< Beginning of the spectral range in wavelengths
            double lstop; ///< End of the spectral range in wavelengths
            double omegastart; ///< Beginning of the spectral range in frequencies
            double omegastop; ///< End of the spectral range in frequencies
            int  nR= INEL_MAX_NREFLEX; ///< Number of reflections
            std::vector<std::function<std::complex<double>(double) > > nList; ///< List of the refractive index functions             
        };

/* ------------------ IMPLEMENTATION ------------------------ */

        template<class T> Trafo<T>::Trafo()
        {
            nR = INEL_MAX_NREFLEX;
        }

        template<class T> Trafo<T>::Trafo(TrafoParms tp)
        {
            nI = tp.nI;
            nS = tp.nS;
            nR = tp.nR;
            lstart = tp.lstart;
            lstop = tp.lstop;
            nList = tp.nList;
            omegastart = 2.0 * M_PI * C_LIGHT / lstop;
            omegastop = 2.0 * M_PI * C_LIGHT / lstart;
        }

        /*template <class T> Trafo<T>::~Trafo()
        {
        }        
        */
        template<class T> void Trafo<T>::initResult(SuperArray<T>& SA)
        {
            SAres = SuperArray<maths::Vector<std::complex<double> > >(SA.r0, SA.nx, SA.ny, SA.nz, SA.Obj, SA.numObjs);
        }

       template<class T> void Trafo<T>::initResult(double r0, int nx, int ny, int nz, ObjectShape** Obj, int numObjs)
        {
           SAres.clear();
           SAres = SuperArray<maths::Vector<std::complex<double> > >(r0, nx, ny, nz, Obj, numObjs);
        }
        
       

        template<class T> GOAT::maths::Vector<std::complex<double> > Trafo<T>::calcOne(std::vector<stepEntry> steps, double t)
        {
            
        } 
        

        template<class T> void Trafo<T>::calc(std::vector<SuperArray <std::vector<gridEntry> > >& SA, double t)
        {
            maths::Vector<std::complex<double> > E;
            std::complex<double> phase;
            int nsteps;
            double dwvl; // step size for the integration inside the subrange in wavelength
            double Dwvl; // width of one subrange
            double domega; // step size for the integration inside the subrange in angular frequency
            double omega; // angular frequency
            double wvl;
            
            initResult(SA[0].r0,SA[0].nges[0], SA[0].nges[1], SA[0].nges[2],SA[0].Obj,SA[0].numObjs);
            
            SAres.fill(maths::czero); // empty the whole result array
            Dwvl = (lstop - lstart) / (double)nI; 
            /* -----  Loops over positions ------ */
             for (int iR = 0; iR < nR; iR++) // loop over reflection order
            for (int i = 0; i < SA[iR].numObjs; i++) // loop over object number (i.e. over Sub-Array in SuperArray)
                for (int ix = 0; ix < SA[iR].n[i][0]; ix++) // loops over x-,y- and z- indices
                    for (int iy = 0; iy < SA[iR].n[i][1]; iy++)
                        for (int iz = 0; iz < SA[iR].n[i][2]; iz++)
                        {
                            SAres.G[i][ix][iy][iz] = integrate(t, SA[iR].G[i][ix][iy][iz], omegastart, omegastop);
                        }
 

        }

        template<class T> void Trafo<T>::createLTexpo()
        {
            double dl = (lstop - lstart) / (double) (nS-nI-1);
            int N = nS * nI;
            double l;
            for (int i = 0; i < N; i++)
            {
                l = lstart + i * dl;
                freq.push_back(2.0 * M_PI * C_LIGHT / l);                
            }                
        }

        template<class T> std::complex<double> Trafo<T>::calcPhase(std::vector<stepEntry> steps, double k0)
        {
            std::complex<double> sum = 0;            
            for (stepEntry se : steps)
                sum += k0 * nList[se.matIndex](2.0*M_PI/k0) * se.l;
            return sum;
        }

        

        template<class T> GOAT::maths::Vector<std::complex<double> >  Trafo<T>::integrate(double t, std::vector<gridEntry> vge, double omegastart, double omegastop, int nsteps)
        {
            GOAT:maths::Vector<std::complex<double> > E;
            std::complex<double> phase;
            double k0;
            double omega;
            double domega = (omegastop - omegastart) / ((double)nsteps - 1.0);
            for (auto ge : vge) // Loop over all rays which hit the cell
            {
                for (int iomega = 0; iomega < nsteps; iomega++)
                {                    
                    omega = omegastart + iomega * domega;
                    k0 = omega / C_LIGHT;
                    phase = exp(I * calcPhase(ge.step, k0));
                    phase *= exp(-I * omega * t);
                    E += ge.E * phase;
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
