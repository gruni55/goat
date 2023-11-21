#pragma once
#include <stdlib.h>
#include <fstream>
#include "vector.h"
#include "resutil.h"
#include "objectshape.h"
#include "ellipsoid.h"
#include "error.h"
#include <vector>
#include "gridentry.h"

// #define INDEX_TYPE long long int
namespace GOAT
{
   namespace raytracing 
   {
       typedef long long int INDEX_TYPE;
      /** 
       * This template class provides a virtual 3D-grid of complex vectors which can be used to store e.g. the electric field in a volume. 
       *  It circumscribes a sphere with radius r0 and virtual means here, that only  the parts which are needed were allocated, 
       *  but for the programmer it seems like a whole array or grid. 
       *  By adding an object to the grid, the memory for the grid around the object is allocated, so one can store the electric field 
       * inside the object. The cells of the grid can be addressed with help of the bracket()-operator with different arguments 
       **/ 
      template <class T> class SuperArray
      {
         public :
         SuperArray ();
         /**
          * @brief Constructor for SuperArray
          * @param r0 Radius of the circumscribed sphere
          * @p nx, @p ny, @p nz number of cells in the x-, y-, z- direction
          * @param type describes wether the whole space of the grid is allocated (IN_HOST) or in the objects only (IN_OBJECT) 
         */
         SuperArray(double r0, INDEX_TYPE nx, INDEX_TYPE ny, INDEX_TYPE nz, const int typ=IN_OBJECT);
         SuperArray(double r0, INDEX_TYPE nx, INDEX_TYPE ny, INDEX_TYPE nz,ObjectShape **Obj, int numObjs, const bool isAbsolute=false, const int typ=IN_OBJECT);
         SuperArray(const SuperArray& S)
         {
             Error = NO_ERRORS;
             
             type = S.type;
             ywerte = S.ywerte;
             zwerte = S.zwerte;
             numObjs = S.numObjs;
             Obj = S.Obj;
             G = S.G;
             K = S.K;
             Pul = S.Pul;
             n = S.n;
             nges = S.nges;
             d = S.d;
             r0 = S.r0;
             isequal = S.isequal;
             iscleared = S.iscleared;
             pc = S.pc;
             H = S.H;
             R = S.R;
         }
         ~SuperArray () { clear();}
         
         /**
          * @brief Copies SuperArray object
          * Here, \p S will be copied into the existing SuperArray and all elements will be overridden
         */
         void copy (const SuperArray &S);
         /**
          * 
          *  @brief Adds another SuperArray object. 
          * In this function, all array elements from \p S are added to the existing array elements.
          *           *
          */ 
         void add (const SuperArray& S);
         /**
          *
          *  @brief Subtract another SuperArray object.
          * In this function, all array elements from \p S are subtracted from the existing array elements.
          *
          */
         void sub (const SuperArray& S);
      //   friend double sumabs (const SuperArray &S); ///< The absolute values of the cell contents are summed up
      //   friend double sumabs2 (const SuperArray &S); ///< The squared absolute values of the cell contents are summed up
         /**
         * @brief Summing up SuperArray
         * Here, all cell contents are first added up. The absolute value is taken from the result and squared.
         */
      //   friend double abs2sum (const SuperArray &S);
         /**
         * @brief Add objects to SuperArray
         * Adds a list of objects to the SuperArray.
         * An object is added to the SuperArray by allocating memory that corresponds to the volume of the circumscribing cuboid.
         * @param Obj The object list
         * @param numObj Number of objects to add
         * @param isAbsolute true: Location/size information in absolute coordinates
         * 
         */
         bool addInc(ObjectShape **Obj,int numObj, const bool isAbsolute=false); // numObjs Einschluesse hinzufuegen
                                           // isAbsolute=true => Orts-/Groessenangaben in absoluten Koordinaten 
         /**
         * @brief Add object to SuperArray
         * Add one objects to the SuperArray.
         * An object is added to the SuperArray by allocating memory that corresponds to the volume of the circumscribing cuboid.
         * @param E The object 
         * @param isAbsolute true: Location/size information in absolute coordinates
         *
         */
         bool addInc (ObjectShape *E,const bool isAbsolute=false);  // Einschluss hinzufuegen (isAbsolute s.o.) 
         maths::Vector<INDEX_TYPE> gitterpunkt (maths::Vector<double> P);
         bool inObject (maths::Vector<double> P, int i); ///< checks if \p P is inside the i-th object (p in real coordinates)
         bool inObject (maths::Vector<INDEX_TYPE> Pi, int i); ///< checks if \p Pi is inside the i-th object (pi in indices)
         bool inObject (INDEX_TYPE ix, INDEX_TYPE iy, INDEX_TYPE iz, int i); ///< checks if a point indicated by its indices (\p ix, \p iy, \p iz) is inside the \p i -th object
         T& operator () (INDEX_TYPE ix, INDEX_TYPE iy, INDEX_TYPE iz); ///< gives the content of the cell[ix][iy][iz]
         /**
          * @brief returns the content of the i-th object, with the grid-coordinates ix,iy,iz
          * @param i number of the object
          * @param ix 
         */
         T& operator () (int i, INDEX_TYPE ix, INDEX_TYPE iy, INDEX_TYPE iz, bool isObjectCoordinate=true); ///< gives back the contents of the cell from the i-th object with the indices (ix,iy,iz) (faster)
         T& operator () (maths::Vector<INDEX_TYPE> Pi); ///< gives back the contents of the cell with indices stored in Pi
         T& operator () (int i,maths::Vector<INDEX_TYPE> Pi); ///< gives back the contents of the cell with indices stored in Pi from the i-th object (faster)
         T& operator () (maths::Vector<double> P); ///< gives back the contents of the cell at P 
         T& operator () (int i, maths::Vector<double> P); ///< gives back the contents of the cell at P from the i-th object (faster)
         
          void makeReal ();
         void fill(const T &x); ///< Fill the whole SuperArray with value \p x
         SuperArray& operator = (const SuperArray &S); ///< Assignment operator
         void clear(); ///< Clear the SuperArray and release the allocated memory
         void allockugel(); 

         /**
         * @brief Checks if the point \p Pi (indicated by the indices) is inside the calculation sphere
         * The function returns the Vector Pi, if it is inside the calculation sphere otherwise a the vector (-1,-1,-1) is returned 
         * (for internal use only)
         */
         maths::Vector<INDEX_TYPE> kugelindex(maths::Vector<INDEX_TYPE> Pi); ///< for internal use only
         T kugelwert(maths::Vector<int> Pi); ///< for internal use only
         T kugelwert(int ix, int iy, int iz); ///< for internal use only
 
         int Error;  ///< Holds an error number 
         ObjectShape** Obj=NULL; ///< here are the objects
         int numObjs; ///< Number of objects
         int type; ///< Mainly used for inelastic scattering. type=IN_HOST means the grid is stored in the whole volume, type=IN_OBJECT means grid is only used in the (active) objects
         std::vector<int> ywerte;
         std::vector<std::vector<int> > zwerte;
         std::vector <std::vector <std::vector <std::vector <T> > > > G; ///< Here, the data is stored. G[i][ix][iy][iz], whereas i: index of the object, ix,iy,iz: indices of the grid around object i
         std::vector < std::vector < std::vector <T> > >  K;
         T dummy;
         std::vector<maths::Vector<INDEX_TYPE>> Pul;
         std::vector<maths::Vector<INDEX_TYPE>> n; ///< n[i]: Dimensions, i.e. number of subdivision in x-, y- and z-direction
         maths::Vector<INDEX_TYPE> nges; ///< Vector which contains the number of subdivisions in x-, y- and z-direction for the whole (virtual) array

         maths::Vector<double> d;  ///<  Edge length of one cell in x-, y- and z-direction
         double r0; ///< Radius of the calculation sphere
         bool isequal;
         bool iscleared = true; ///< Array is cleared (needed by the clear() function)
         T pc;        
         maths::Matrix<double> H, R; ///< Transformation matrices in the local array coordinate system and backwards
         bool write(std::string fname);
         bool read(std::string fname);
    };

    /* Specialization of Superarray for use together with std::vector<gridEntry> */
    template<> void SuperArray<std::vector<GOAT::raytracing::gridEntry> >::clear();
   
    /*------------------------- IMPLEMENTATION --------------------------------------*/

    template<> bool SuperArray<GOAT::maths::Vector<std::complex<double> > >::write(std::string fname);

    template <class T> SuperArray<T>::SuperArray()
    {
        H = maths::unity();
        R = maths::unity();
        Error = NO_ERRORS;
        numObjs = 0;
        isequal = false;

        type = IN_HOST;        
    }

    template <class T> SuperArray<T>::SuperArray(double r0, INDEX_TYPE nx, INDEX_TYPE ny, INDEX_TYPE nz, const int typ)
    {
        Error = NO_ERRORS;
        double b = 2.0 * r0;
        isequal = false;
        numObjs = 0;
        Obj = 0;
        this->r0 = r0;



        nges = maths::Vector<INDEX_TYPE>(nx, ny, nz);
        d = maths::Vector<double>(b / (double)nx, b / (double)ny, b / (double)nz);
        iscleared = true;
        type = typ;
        if (type & IN_HOST)
            allockugel();
    }

    template <class T> SuperArray<T>::SuperArray(double r0, INDEX_TYPE nx, INDEX_TYPE ny, INDEX_TYPE nz, ObjectShape** Obj, int numObjs, const bool isAbsolute, const int typ)
    {
        Error = NO_ERRORS;
        double b = 2.0 * r0;
        isequal = isAbsolute;
        this->numObjs = 0;
        this->r0 = r0;
        nges = maths::Vector<INDEX_TYPE>(nx, ny, nz);
        type = typ;
        d = maths::Vector<double>(b / nx, b / ny, b / nz);
        iscleared = true;
        if (type & IN_HOST)
        {
            allockugel();
        }
        addInc(Obj, numObjs, isAbsolute);
    }


    template <class T> bool SuperArray<T>::addInc(ObjectShape** Obj, int numObjs, const bool isAbsolute)
    {
        bool ok = true;
        for (int i = 0; (i < numObjs)/* && (canCalc)*/ && (ok); i++) ok = addInc(Obj[i], isAbsolute);
        return ok;
    }

    template <class T> maths::Vector<INDEX_TYPE> SuperArray<T>::gitterpunkt(maths::Vector<double> P)
    {
        maths::Vector<INDEX_TYPE> pi, ph;
        maths::Vector<double> h;
        h = H * P + maths::Vector<double>(r0, r0, r0);
        pi = maths::Vector<INDEX_TYPE>((INDEX_TYPE)floor(h[0] / d[0]), (INDEX_TYPE)floor(h[1] / d[1]), (INDEX_TYPE)floor(h[2] / d[2]));

        /*
        if (type & IN_HOST) //????????
        {
            pi = ph;
        }
        else
        {
            pi = ph;
        }*/
        return pi;
    }

    template <class T> bool SuperArray<T>::addInc(ObjectShape* E, const bool isAbsolute)
    {
       
        //SysMemInfo smi;
        //MemInfo mi;
        long int allocMem;
        maths::Vector<INDEX_TYPE> pul, por, N, hn;
        maths::Vector<double> h, O;
        O = maths::Vector<double>(-r0, -r0, -r0);

       /* smi = sysmem();
        mi = memstat();
        */
        h = ceil(ediv(E->por, d)) - floor(ediv(E->pul, d));

        hn = maths::Vector<INDEX_TYPE>((INDEX_TYPE)h[0]+1, (INDEX_TYPE)h[1]+1, (INDEX_TYPE)h[2]+1); // Gr��e des 3D-Gitters in die drei Koordinatenrichtungen

        /* Berechne den tats�chlichen Bedarf */
        allocMem = sizeof(T***) + 2 * sizeof(maths::Vector<int>) + sizeof(ObjectShape*)
            + hn[0] * sizeof(T**)
            + hn[0] * hn[1] * sizeof(T*)
            + hn[0] * hn[1] * hn[2] * sizeof(T);

        {
            if (numObjs < 1)  // Es ist der erste Einschluss der hinzugef�gt wird
            {   
                G.resize(1);                        
                Obj = (ObjectShape**)malloc(sizeof(ObjectShape*));
                if (Obj == NULL) { error(MALLOC_ERR, "SuperArray::addInc Obj=.."); return false; }
            }
            else
            {
                G.resize(numObjs + 1);
                Obj = (ObjectShape**)realloc(Obj, (numObjs + 1) * sizeof(ObjectShape*));
                if (Obj == NULL) { error(REALLOC_ERR, "SuperArray::addInc Obj=.."); return false; }
            }
        }
      
        n.push_back(hn);
        Obj[numObjs] = E;
        h = floor(ediv(Obj[numObjs]->pul + maths::Vector<double>(r0, r0, r0), d));
        Pul.push_back(maths::Vector<INDEX_TYPE>((INDEX_TYPE)h[0], (INDEX_TYPE)h[1], (INDEX_TYPE)h[2]));
        
        if (E->isActive())  // Ist der Einschluss �berhaupt inelastisch aktiv ? 
        {
//            std::cout << "n[" << numObjs << "]=" << n[numObjs] << std::endl;
            G[numObjs].resize(n[numObjs][0] + 1);
            
            for (INDEX_TYPE ix = 0; ix < n[numObjs][0] + 1; ix++)
            {
                G[numObjs][ix].resize(n[numObjs][1] + 1);                
                
                for (INDEX_TYPE iy = 0; iy < n[numObjs][1] + 1; iy++)
                {
                    G[numObjs][ix][iy].resize(n[numObjs][2] + 1) ;                                       
                }
            }
        }
        

        int i = 0;
        numObjs++;
        iscleared = false;
        return true;
    }

    template <class T> bool SuperArray<T>::inObject(maths::Vector<double> P, int i) // da muss noch was gemacht werden
    {
        return Obj[i]->isInside(P);
    }

    template <class T> bool SuperArray<T>::inObject(maths::Vector<INDEX_TYPE> Pi, int i)
    {
        return (Obj[i]->isInside(emult(Pi, d)));
    }

    template <class T> bool SuperArray<T>::inObject(INDEX_TYPE ix, INDEX_TYPE iy, INDEX_TYPE iz, int i)
    {
        return (Obj[i]->isInside(emult(maths::Vector<double>(ix, iy, iz), d)));
    }

    template <class T> T& SuperArray<T>::operator () (INDEX_TYPE ix, INDEX_TYPE iy, INDEX_TYPE iz)
    {
        maths::Vector<int> Pi = maths::Vector<int>(ix, iy, iz);
        int i = 0;
        bool found = false;
        // dummy = maths::Vector<std::complex<double> >(0.0, 0.0, 0.0);
        if (type == IN_OBJECT)
        {
            do
            {
                found = inObject(ix, iy, iz, i);
                i++;
            } while ((i < numObjs) && (!found));
            i--;
            if (!found)
            {
                Error = NOT_FOUND;
                return dummy;
            }
            else
            {
                Pi = Pi - Pul[i];
                Pi = kugelindex(Pi);
                if (Error == NO_ERRORS)
                    return G[i][Pi[0]][Pi[1]][Pi[2]];
                else
                {
                    Error = NOT_FOUND;
                    return dummy;
                }
            }
        }
        else
        {
            Pi = kugelindex(Pi);
            if (Error == NO_ERRORS)
                return K[Pi[0]][Pi[1]][Pi[2]];
            else
            {
                Error = NOT_FOUND;
                return dummy;
            }
        }
    }

    template <class T> T& SuperArray<T>::operator () (maths::Vector<INDEX_TYPE> Pi)
    {
        // dummy = maths::Vector<std::complex<double> >(0.0, 0.0, 0.0);
       
        if (type == IN_OBJECT)
        {
            int i = 0;
            bool found = false;
            do
            {
                found = inObject(Pi, i);
                i++;
            } while ((i < numObjs) && (!found));
            i--;
            
            if (!found)
            {
                Error = NOT_FOUND;
                return dummy;
            }
            else
            {
                Pi = Pi - Pul[i];
                Pi = kugelindex(Pi);                
                if (Error == NO_ERRORS)
                    return G[i][Pi[0]][Pi[1]][Pi[2]];
                else
                {
                    Error = NOT_FOUND;
                    return dummy;
                }
            }
        }
        else
        {
            Pi = kugelindex(Pi);
            if (Error == NO_ERRORS)
                return K[Pi[0]][Pi[1]][Pi[2]];
            else
            {
                Error = NOT_FOUND;
                return dummy;
            }
        }
    }


    template <class T> T& SuperArray<T>::operator () (int i, INDEX_TYPE ix, INDEX_TYPE iy, INDEX_TYPE iz, bool isEinKoord)
    {
        T dummy;
        maths::Vector<int> Pi = maths::Vector<int>(ix, iy, iz);
        if (type == IN_OBJECT)
        {
            if (G[i].size() == 0)
            {
                Error = NOT_FOUND;
                return dummy;
            }
            if (!isEinKoord)  Pi = Pi - Pul[i];
            return G[i][Pi[0]][Pi[1]][Pi[2]];
        }
        else
        {
            Pi = kugelindex(Pi);
            if (Error == NO_ERRORS)
                return K[Pi[0]][Pi[1]][Pi[2]];
            else
            {
                Error = NOT_FOUND;
                return dummy;
            }
        }
    }

    template <class T> T& SuperArray<T>::operator () (int i, maths::Vector<INDEX_TYPE> Pi)
    {
        T dummy;
        if (type == IN_OBJECT)
        {
            Pi = Pi - Pul[i];
            if (G[i].size() == 0) { Error = NOT_FOUND; return dummy; }
            if (Pi[0] < 0) { Error = NOT_FOUND; return dummy; }//maths::Vector<std::complex<double> > (0,0,0);
            if (Pi[1] < 0) { Error = NOT_FOUND; return dummy; } //maths::Vector<std::complex<double> > (0,0,0);
            if (Pi[2] < 0) { Error = NOT_FOUND; return dummy; } // maths::Vector<std::complex<double> > (0,0,0);
            if (Pi[0] >= n[i][0])
            {
                Error = SUPERGITTER;
                 return dummy; //maths::Vector<std::complex<double> > (0,0,0);
            }

            if (Pi[1] >= n[i][1])
            {
                Error = SUPERGITTER;
                 return dummy; //maths::Vector<std::complex<double> > (0,0,0);
            }

            if (Pi[2] >= n[i][2])
            {
                Error = SUPERGITTER;
                 return dummy;
            }

            Error = NO_ERRORS;
           return G[i][Pi[0]][Pi[1]][Pi[2]];
        }
        else
        {
            Pi = kugelindex(Pi);
            if (Error == NO_ERRORS)
                return K[Pi[0]][Pi[1]][Pi[2]];
            else
            {
                Error = NOT_FOUND; return dummy;
            }
        }
    }

    template <class T> T& SuperArray<T>::operator () (maths::Vector<double> P)
    {
        int i;
        maths::Vector<int> Pi;
        maths::Vector<double> h;
        h = H * P - maths::Vector<double>(-r0, -r0, -r0);
        dummy = maths::Vector<std::complex<double> >(0.0, 0.0, 0.0);
        for (i = 0; i < 3; i++)
            Pi[i] = h[i] / d[i];

        if (type == IN_OBJECT)
        {
            bool found = false;
            i = 0;
            do
            {
                found = inObject(Pi, i);
                i++;
            } while ((i < numObjs) && (!found));
            i--;
            // if (!found) return maths::Vector<std::complex<double> > (INF,INF,INF);
            if (!found)
            {
                Error = NOT_FOUND;
                return dummy;
            }
            else
            {
                if (G[i] == NULL)
                {
                    Error = NOT_FOUND;
                    return dummy;
                }
                Pi = Pi - Pul[i];
                std::cout << "PI:" << Pi << std::endl;
                return G[i][Pi[0]][Pi[1]][Pi[2]];
            }
        }
        else
        {
            Pi = kugelindex(Pi);
            if (Error == NO_ERRORS)
                return K[Pi[0]][Pi[1]][Pi[2]];
            else
            {
                Error = NOT_FOUND;
                return dummy;
            }
        }
    }

    template <class T> T& SuperArray<T>::operator () (int i, maths::Vector<double> P)
    {
        maths::Vector<int> Pi;
        maths::Vector<double> h;
        h = H * P - maths::Vector<double>(-r0, -r0, -r0);
        T dummy;
        // dummy = maths::Vector<std::complex<double> >(0.0, 0.0, 0.0);

        for (int j = 0; j < 3; j++)
            Pi[j] = h[j] / d[j];

        if (type == IN_OBJECT)
        {
            if (G[i] == NULL)
            {
                Error = NOT_FOUND;
                return dummy;
            }
            Pi = Pi - Pul[i];
            return G[i][Pi[0]][Pi[1]][Pi[2]];
        }
        else
        {
            Pi = kugelindex(Pi);
            if (Error == NO_ERRORS)
                return K[Pi[0]][Pi[1]][Pi[2]];
            else
            {
                Error = NOT_FOUND;
                return dummy;
            }
        }
    }


    template <class T> void SuperArray<T>::clear()
    {
        int anzx, anzx2;
        anzx = nges[0];
        anzx2 = nges[0] / 2;
        if (!iscleared)
        if (type == IN_OBJECT)
        {
            if (numObjs > 0)
            {
                for (int i = numObjs - 1; i >= 0; i--)
                {
                    if (G[i].size()>0)
                    {
                        for (int ix = n[i][0]; ix >= 0; ix--)
                        {
                            for (int iy = n[i][1]; iy >= 0; iy--)
                                G[i][ix][iy].clear();
                            G[i][ix].clear();                            
                        } // for ix    
                        G[i].clear();
                    } // if (G[i]!=NULL)
                } // for i
                n.clear(); 
                Pul.clear();
              /*  if (Obj != NULL)
                free(Obj);
                numObjs = 0;
                Obj = NULL;*/
                iscleared = true;
                
            }
        }
        else
        {
            if (K.size()>0)
            {
                for (int k = 0; k < anzx2; k++)
                {
                    /*for (int l = 0; l < ywerte[k]; l++)
                    {
                        delete[] K[anzx2 - 1 - k][ywerte[k] + l];
                        delete[] K[anzx2 + k][ywerte[k] + l];
                        delete[] K[anzx2 - 1 - k][ywerte[k] - 1 - l];
                        delete[] K[anzx2 + k][ywerte[k] - 1 - l];
                    }
                    */
                    //delete[] K[anzx2 - 1 - k];
                    //delete[] K[anzx2 + k];
                    zwerte[k].clear();                    
                }

                ywerte.clear();
                zwerte.clear();

               /* if (numObjs > 0)
                    free(Obj);*/
            }
            
            // numObjs = 0;
            iscleared = true;
        }
    }

    template <class T> void SuperArray<T>::copy(const SuperArray<T>& S)  // MUSS DRINGEND GE�NDERT WERDEN !
    {
        clear();
        addInc(S.Obj, S.numObjs);
        for (int i = 0; i < numObjs; i++)
            if (Obj[i]->isActive())
                for (int ix = 0; ix < n[i][0]; ix++)
                    for (int iy = 0; iy < n[i][1]; iy++)
                        for (int iz = 0; iz < n[i][2]; iz++)
                            G[i][ix][iy][iz] = S.G[i][ix][iy][iz];
        K = S.K;
        Pul = S.Pul;
        n = S.n;
        nges = S.nges;
        d = S.d;
        r0 = S.r0;
        isequal = S.isequal;
        iscleared = S.iscleared;
        pc = S.pc;
        H = S.H;
        R = S.R;
        
    }

    template <class T> SuperArray<T>& SuperArray<T>::operator = (const SuperArray<T>& S)
    {
        int i = 0;
        int anzx, anzx2;


        if (this == &S) return *this;
        clear();
        if (this->numObjs > 0)
        {
        /*  for (int i = 0; i < this->numObjs; i++)
                delete this->Obj[i];*/
            delete[] this->Obj;
            this->numObjs = 0;
        }
        r0 = S.r0;
        d = S.d;
        nges = S.nges;
        type = S.type;

        anzx = nges[0]; anzx2 = nges[0] / 2;

        if (type == IN_OBJECT)
        {
            isequal = true; addInc(S.Obj, S.numObjs, true);
            isequal = true;
            G = S.G;
          /*  if (S.numObjs != 0)
                do
                {
                    for (int ix = 0; ix < n[i][0]; ix++)
                        for (int iy = 0; iy < n[i][1]; iy++)
                            for (int iz = 0; iz < n[i][2]; iz++)
                                G[i][ix][iy][iz] = S.G[i][ix][iy][iz];
                    i++;
                } while (i < S.numObjs);*/
                isequal = false;
        }
        else
        {
            allockugel();
            isequal = true; addInc(S.Obj, S.numObjs, true);
            isequal = true;
            
            for (int k = 0; k < anzx2; k++)
            {
                ywerte[k] = S.ywerte[k];                            
                for (int l = 0; l < S.ywerte[k]; l++)
                {
                    zwerte[k][l] = S.zwerte[k][l];
                    for (int m = 0; m < S.zwerte[k][l]; m++)
                    {
                        K[anzx2 - 1 - k][ywerte[k] - 1 - l][m] = S.K[anzx2 - 1 - k][ywerte[k] - 1 - l][m];
                        K[anzx2 + k][ywerte[k] - 1 - l][m] = S.K[anzx2 + k][ywerte[k] - 1 - l][m];
                        K[anzx2 - 1 - k][ywerte[k] + l][m] = S.K[anzx2 - 1 - k][ywerte[k] + l][m];
                        K[anzx2 + k][ywerte[k] + l][m] = S.K[anzx2 + k][ywerte[k] + l][m];
                        K[anzx2 - 1 - k][ywerte[k] - 1 - l][2 * zwerte[k][l] - 1 - m] = S.K[anzx2 - 1 - k][ywerte[k] - 1 - l][2 * zwerte[k][l] - 1 - m];
                        K[anzx2 + k][ywerte[k] - 1 - l][2 * zwerte[k][l] - 1 - m] = S.K[anzx2 + k][ywerte[k] - 1 - l][2 * zwerte[k][l] - 1 - m];
                        K[anzx2 - 1 - k][ywerte[k] + l][2 * zwerte[k][l] - 1 - m] = S.K[anzx2 - 1 - k][ywerte[k] + l][2 * zwerte[k][l] - 1 - m];
                        K[anzx2 + k][ywerte[k] + l][2 * zwerte[k][l] - 1 - m] = S.K[anzx2 + k][ywerte[k] + l][2 * zwerte[k][l] - 1 - m];
                    }
                }
            }
        }
        return *this;
    }



    template <class T> void SuperArray<T>::add(const SuperArray<T>& S)
    {
        int anzx2 = nges[0] / 2;
        isequal = false;
        H = S.H;
        R = S.R;
        if (type == IN_OBJECT)
        {
            for (int i = 0; i < S.numObjs; i++)
                for (int ix = 0; ix < S.n[i][0]; ix++)
                    for (int iy = 0; iy < S.n[i][1]; iy++)
                        for (int iz = 0; iz < S.n[i][2]; iz++)
                            G[i][ix][iy][iz] += S.G[i][ix][iy][iz];
        }
        else
        {
            for (int k = 0; k < anzx2; k++)
            {
                for (int l = 0; l < ywerte[k]; l++)
                {

                    for (int m = 0; m < zwerte[k][l]; m++)
                    {
                        K[anzx2 - 1 - k][ywerte[k] - 1 - l][m] += S.K[anzx2 - 1 - k][ywerte[k] - 1 - l][m];
                        K[anzx2 + k][ywerte[k] - 1 - l][m] += S.K[anzx2 + k][ywerte[k] - 1 - l][m];
                        K[anzx2 - 1 - k][ywerte[k] + l][m] += S.K[anzx2 - 1 - k][ywerte[k] + l][m];
                        K[anzx2 + k][ywerte[k] + l][m] += S.K[anzx2 + k][ywerte[k] + l][m];
                        K[anzx2 - 1 - k][ywerte[k] - 1 - l][2 * zwerte[k][l] - 1 - m] += S.K[anzx2 - 1 - k][ywerte[k] - 1 - l][2 * zwerte[k][l] - 1 - m];
                        K[anzx2 + k][ywerte[k] - 1 - l][2 * zwerte[k][l] - 1 - m] += S.K[anzx2 + k][ywerte[k] - 1 - l][2 * zwerte[k][l] - 1 - m];
                        K[anzx2 - 1 - k][ywerte[k] + l][2 * zwerte[k][l] - 1 - m] += S.K[anzx2 - 1 - k][ywerte[k] + l][2 * zwerte[k][l] - 1 - m];
                        K[anzx2 + k][ywerte[k] + l][2 * zwerte[k][l] - 1 - m] += S.K[anzx2 + k][ywerte[k] + l][2 * zwerte[k][l] - 1 - m];
                    }
                }
            }
        }
    }

    template <class T> void SuperArray<T>::sub(const SuperArray<T>& S)
    {
        int anzx2 = nges[0] / 2;

        isequal = false;
        if (type == IN_OBJECT)
        {
            for (int i = 0; i < S.numObjs; i++)
                for (int ix = 0; ix < S.n[i][0]; ix++)
                    for (int iy = 0; iy < S.n[i][1]; iy++)
                        for (int iz = 0; iz < S.n[i][2]; iz++)
                            G[i][ix][iy][iz] -= S.G[i][ix][iy][iz];
        }
        else
        {
            for (int k = 0; k < anzx2; k++)
            {
                for (int l = 0; l < ywerte[k]; l++)
                {

                    for (int m = 0; m < zwerte[k][l]; m++)
                    {
                        K[anzx2 - 1 - k][ywerte[k] - 1 - l][m] -= S.K[anzx2 - 1 - k][ywerte[k] - 1 - l][m];
                        K[anzx2 + k][ywerte[k] - 1 - l][m] -= S.K[anzx2 + k][ywerte[k] - 1 - l][m];
                        K[anzx2 - 1 - k][ywerte[k] + l][m] -= S.K[anzx2 - 1 - k][ywerte[k] + l][m];
                        K[anzx2 + k][ywerte[k] + l][m] -= S.K[anzx2 + k][ywerte[k] + l][m];
                        K[anzx2 - 1 - k][ywerte[k] - 1 - l][2 * zwerte[k][l] - 1 - m] -= S.K[anzx2 - 1 - k][ywerte[k] - 1 - l][2 * zwerte[k][l] - 1 - m];
                        K[anzx2 + k][ywerte[k] - 1 - l][2 * zwerte[k][l] - 1 - m] -= S.K[anzx2 + k][ywerte[k] - 1 - l][2 * zwerte[k][l] - 1 - m];
                        K[anzx2 - 1 - k][ywerte[k] + l][2 * zwerte[k][l] - 1 - m] -= S.K[anzx2 - 1 - k][ywerte[k] + l][2 * zwerte[k][l] - 1 - m];
                        K[anzx2 + k][ywerte[k] + l][2 * zwerte[k][l] - 1 - m] -= S.K[anzx2 + k][ywerte[k] + l][2 * zwerte[k][l] - 1 - m];
                    }
                }
            }
        }
    }

    template <class T> std::ostream& operator << (std::ostream& os, const SuperArray<T>& S)
    {
        os << "r0=" << S.r0 << std::endl;
        os << "Ausdehnung:    nx=" << S.nges[0] << "  ny=" << S.nges[1] << "  nz=" << S.nges[2] << std::endl;
        os << S.numObjs << " Einschluesse" << std::endl;
        if (S.numObjs > 0)
            for (int i = 0; i < S.numObjs; i++)
            {
                os << "====================== Einschluss Nr. " << i << " ======================" << std::endl;
                for (int ix = 0; ix < S.n[i][0]; ix++)
                    for (int iy = 0; iy < S.n[i][1]; iy++)
                        for (int iz = 0; iz < S.n[i][2]; iz++)
                            os << abs(S.G[i][ix][iy][iz]) << std::endl;
            }
        return os;
    } 
  
    std::ostream& operator << (std::ostream& os, const SuperArray<std::vector<gridEntry> >& S);
    

    template <class T> void SuperArray<T>::fill(const T& x)
    {
        int anzx, anzx2;
        anzx = nges[0]; anzx2 = nges[0] / 2;

        if (type == IN_OBJECT)
        {
            for (int i = 0; i < numObjs; i++)
                if (Obj[i]->isActive())
                for (int ix = 0; ix < n[i][0]; ix++)
                    for (int iy = 0; iy < n[i][1]; iy++)
                        for (int iz = 0; iz < n[i][2]; iz++)
                            G[i][ix][iy][iz] = x;
        }
        else
        {
            for (int k = 0; k < anzx2; k++)
            {
                for (int l = 0; l < ywerte[k]; l++)
                {
                    for (int m = 0; m < zwerte[k][l]; m++)
                    {
                        K[anzx2 - 1 - k][ywerte[k] - 1 - l][zwerte[k][l] - 1 - m] = x;
                        K[anzx2 + k][ywerte[k] - 1 - l][zwerte[k][l] - 1 - m] = x;
                        K[anzx2 - 1 - k][ywerte[k] + l][zwerte[k][l] - 1 - m] = x;
                        K[anzx2 + k][ywerte[k] + l][zwerte[k][l] - 1 - m] = x;
                        K[anzx2 - 1 - k][ywerte[k] - 1 - l][zwerte[k][l] + m] = x;
                        K[anzx2 + k][ywerte[k] - 1 - l][zwerte[k][l] + m] = x;
                        K[anzx2 - 1 - k][ywerte[k] + l][zwerte[k][l] + m] = x;
                        K[anzx2 + k][ywerte[k] + l][zwerte[k][l] + m] = x;
                    }
                }
            }
        }
    }


    template <class T> void SuperArray<T>::makeReal()
    {
        for (int i = 0; i < numObjs; i++)
            for (int ix = 0; ix < n[i][0]; ix++)
                for (int iy = 0; iy < n[i][1]; iy++)
                    for (int iz = 0; iz < n[i][2]; iz++)
                        for (int j = 0; j < 3; j++)
                            G[i][ix][iy][iz][j] = abs(G[i][ix][iy][iz][j]);
    }


    

    template <class T> void SuperArray<T>::allockugel()
    {
        int anzx, anzx2;
        double dx, dy, dz, hilf;
        anzx = nges[0];
        anzx2 = nges[0] / 2;
        // cout << "allockugel" << std::endl;
      
        /*dx=d[0];
        dy=d[1];
        dz=d[2];*/
        dx = 2.0 / nges[0];
        dy = 2.0 / nges[1];
        dz = 2.0 / nges[2];

        zwerte.resize(anzx2);
        for (int i = 0; i < anzx2; i++)
        {
            ywerte.push_back(ceil(sqrt(1.0 - (i * dx) * (i * dx)) / dy));
            for (int j = 0; j < anzx2; j++)
            {
                hilf = 1.0 - (i * dx) * (i * dx) - (j * dy) * (j * dy);
                //      cout << "hilf:" << hilf << std::endl;
                if (hilf >= 0)
                    //        zwerte[i][j]=int(ceil(sqrt(1.0-(i*dx)*(i*dx)-(j*dy)*(j*dy))*nges[2]/2.0));
                    zwerte[i].push_back(ceil(sqrt(hilf) / dz));
                else
                    zwerte[i].push_back(1);

                //      cout << "zwerte[" << i << "][" << j << "]:" <<  zwerte[i][j] << std::endl;
            }
        }

        K.resize(anzx);
        for (int k = 0; k < anzx2; k++)
        {
            K[anzx2 - 1 - k].resize(2 * ywerte[k]);
            K[anzx2 + k].resize(2 * ywerte[k]);
            for (int l = 0; l < ywerte[k]; l++)
            {
                K[anzx2 - 1 - k][ywerte[k] + l].resize(2 * zwerte[k][l]);                
                K[anzx2 + k][ywerte[k] + l].resize(2 * zwerte[k][l]);
                K[anzx2 - 1 - k][ywerte[k] - 1 - l].resize(2 * zwerte[k][l]);
                K[anzx2 + k][ywerte[k] - 1 - l].resize(2 * zwerte[k][l]);                
            }
        }
    }

    template <class T> maths::Vector<INDEX_TYPE> SuperArray<T>::kugelindex(maths::Vector<INDEX_TYPE> Pi)
    {
        int anzx, jx, jy, jz, ix, iy, iz, jxh, jyh;
        maths::Vector<INDEX_TYPE> pk;
        anzx = nges[0];
        ix = Pi[0]; iy = Pi[1]; iz = Pi[2];
        Error = NO_ERRORS;
        jx = ix;
        jxh = abs(int(jx - (anzx - 1) / 2.0));
        jy = iy - (anzx / 2 - ywerte[jxh]);
        if ((jy < 0) || (jy > (2 * ywerte[jxh] - 1)))   // Punkt ausserhalb 
        {
            Error = SUPERGITTER;
            //    cout << "ix:" << ix << ", iy:" << iy << ", iz:" << iz <<  ", jx:" << jx <<
            //             ", jxh:" << jxh << ", jy:" << jy << std::endl;
            //    cout << "ywerte[" << jxh << "]:" << ywerte[jxh] << std::endl;
            return maths::Vector<INDEX_TYPE>(-1, -1, -1);
        }
        else
        {
            jyh = abs(INDEX_TYPE(jy - (2 * ywerte[jxh] - 1) / 2.0));
            jz = iz - (anzx / 2 - zwerte[jxh][jyh]);

            if ((jz < 0) || (jz > (2 * zwerte[jxh][jyh]) - 1)) // Punkt ausserhalb 
            {
                Error = SUPERGITTER;
                //    cout << "ix:" << ix << ", iy:" << iy << ", iz:" << iz << ", jx:" << jx << ", jy:" << jy << 
                //             ", jxh:" << jxh << ", jyh:" << jyh << ", jz:" << jz << std::endl;
                //    cout << "ywerte[" << jxh << "]:" << ywerte[jxh] << ", zwerte[" << jxh << "][" << jyh << "]:" <<
                //    zwerte[jxh][jyh] << std::endl;
                return maths::Vector<INDEX_TYPE>(-1, -1, -1);
                return maths::Vector<INDEX_TYPE>(-1, -1, -1);
            }
            else
            {
                return  maths::Vector<INDEX_TYPE>(jx, jy, jz);
            }
        }
        return pk;
    }

    template<class T> T SuperArray<T>::kugelwert(maths::Vector<int> Pi)
    {
        int anzx, jx, jy, jz, ix, iy, iz, jxh, jyh;
        // maths::Vector<std::complex<double> > pc;
        anzx = nges[0];
        ix = Pi[0]; iy = Pi[1]; iz = Pi[2];
        jx = ix;
        jxh = abs(int(jx - (anzx - 1) / 2.0));
        jy = iy - (anzx / 2 - ywerte[jxh]);
        pc = maths::Vector<std::complex<double> >(0.0, 0.0, 0.0);
        if ((jy < 0) || (jy > (2 * ywerte[jxh] - 1)))
        {
            std::cout << "Fehler!! Index iy falsch" << std::endl;
            return pc;
        }
        else
        {
            jyh = abs(int(jy - (2 * ywerte[jxh] - 1) / 2.0));
            jz = iz - (anzx / 2 - zwerte[jxh][jyh]);

            if ((jz < 0) || (jz > (2 * zwerte[jxh][jyh]) - 1))
            {
                return pc;
            }
            return K[jx][jy][jz];
        }
        return pc;
    }

    template<class T> T SuperArray<T>::kugelwert(int ix, int iy, int iz)
    {
        int anzx, jx, jy, jz, jxh, jyh;
        T pc;
        anzx = nges[0];


        jx = ix;
        jxh = abs(int(jx - (anzx - 1) / 2.0));
        jy = iy - (anzx / 2 - ywerte[jxh]);
        if ((jy < 0) || (jy > (2 * ywerte[jxh] - 1)))
        {
            std::cout << "Fehler!! Index iy falsch" << std::endl;
           // pc = maths::Vector<std::complex<double> >(0.0, 0.0, 0.0);
            return pc;
        }
        else
        {
            jyh = abs(int(jy - (2 * ywerte[jxh] - 1) / 2.0));
            jz = iz - (anzx / 2 - zwerte[jxh][jyh]);
            if ((jz < 0) || (jz > (2 * zwerte[jxh][jyh]) - 1))
            {
                std::cout << "Fehler!! Index iz falsch" << std::endl;
               // pc = maths::Vector<std::complex<double> >(0.0, 0.0, 0.0);
                return pc;
            }
            pc = K[jx][jy][jz];
        }
        return pc;
    }

    

    bool saveExPhase(SuperArray<maths::Vector<std::complex<double> > > &S, char* FName, int i = 0);
    bool saveEyPhase(SuperArray<maths::Vector<std::complex<double> > > &S, char* FName, int i = 0);
    bool saveEzPhase(SuperArray<maths::Vector<std::complex<double> > > &S, char* FName, int i = 0);
    bool saveExPol(SuperArray < maths::Vector < std::complex<double> > > &S, char* FName, int i = 0);
    bool saveEyPol(SuperArray < maths::Vector < std::complex<double> > > &S, char* FName, int i = 0);
    bool saveEzPol(SuperArray < maths::Vector < std::complex<double> > > &S, char* FName, int i = 0);
    bool saveabsE(SuperArray < maths::Vector < std::complex<double> > > &S, std::string FName, int i = 0);
    bool saveabsEbin(SuperArray < maths::Vector < std::complex<double> > >& S, std::string FName, int i = 0);
    bool saveFullE(SuperArray < maths::Vector < std::complex<double> > > &S, std::string FName, int i = 0);
    double sumabs(const SuperArray<maths::Vector<std::complex<double> > >& S);
    double sumabs2(const SuperArray<maths::Vector<std::complex<double> > >& S);
    double abs2sum(const SuperArray<maths::Vector<std::complex<double> > >& S);
    bool save(SuperArray<GOAT::raytracing::gridEntry > S, std::string FName);
    bool save(SuperArray<std::vector<GOAT::raytracing::gridEntry > > S, std::string FName);
      }
}
