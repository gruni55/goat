
#pragma once

#include <stdlib.h>
#include <fstream>
#include "vector.h"
#include "resutil.h"
#include "objectshape.h"
#include <vector>

namespace GOAT
{
    namespace raytracing
    {
        constexpr int SUPERGRID_NOERRORS = -1;
        /**
         * This class provides a virtual 3D-grid of complex vectors which can be used to store e.g. the electric field in a volume.
         *  It circumscribes a sphere with radius r0 and virtual means here, that only  the parts which are needed were allocated,
         *  but for the programmer it seems like a whole array or grid.
         *  By adding an object to the grid, the memory for the grid around the object is allocated, so one can store the electric field
         * inside the object. The cells of the grid can be addressed with help of the bracket()-operator with different arguments
         **/
       template<class T> class SuperGrid
        {
        public:
            SuperGrid();
            /**
             * @brief Constructor for SuperArray
             * @param r0 Radius of the circumscribed sphere
             * @p nx, @p ny, @p nz number of cells in the x-, y-, z- direction
             * @param type describes wether the whole space of the grid is allocated (IN_HOST) or in the objects only (IN_OBJECT)
            */
            SuperGrid(double r0, int nx, int ny, int nz, const int typ = IN_HOST);
            SuperGrid(double r0, int nx, int ny, int nz, ObjectShape** Obj, int numObjs, const bool isAbsolute = false, const int typ = 0);
            ~SuperGrid() { clear(); };

            /**
             * @brief Copies SuperArray object
             * Here, \p S will be copied into the existing SuperArray and all elements will be overridden
            */
            void copy(const SuperGrid& S);
            /**
             *
             *  @brief Adds another SuperGrid object.
             * In this function, all array elements from \p S are added to the existing array elements.
             *
             */
            void add(const SuperGrid& S);
            /**
             *
             *  @brief Subtract another SuperArray object.
             * In this function, all array elements from \p S are subtracted from the existing array elements.
             *
             */
            void sub(const SuperGrid& S);
            
            bool addPlane(maths::Vector<double> P, maths::Vector<double> n, double d1, double d2, maths::Vector<double> cellSize);
            /**
            * @brief Add objects to SuperGrid
            * Adds a list of objects to the SuperGrid.
            * An object is added to the SuperGrid by allocating memory that corresponds to the volume of the circumscribing cuboid.
            * @param Obj The object list
            * @param numObj Number of objects to add
            * @param isAbsolute true: Location/size information in absolute coordinates
            *
            */
            bool addObject(ObjectShape** Obj, int numObj, const bool isAbsolute = false); // numObjs Einschluesse hinzufuegen
            // isAbsolute=true => Orts-/Groessenangaben in absoluten Koordinaten 
/**
* @brief Add object to SuperGrid
* Add one objects to the SuperGrid.
* An object is added to the SuperGrid by allocating memory that corresponds to the volume of the circumscribing cuboid.
* @param E The object
* @param isAbsolute true: Location/size information in absolute coordinates
*
*/
            bool addObject(ObjectShape* E, const bool isAbsolute = false);  ///< add new object to the grid 
            maths::Vector<int> gitterpunkt(maths::Vector<double> P);
            bool inObject(maths::Vector<double> P, int i); ///< checks if \p P is inside the i-th object (p in real coordinates)
            bool inObject(maths::Vector<int> Pi, int i); ///< checks if \p Pi is inside the i-th object (pi in indices)
            bool inObject(int ix, int iy, int iz, int i); ///< checks if a point indicated by its indices (\p ix, \p iy, \p iz) is inside the \p i -th object
            T & operator () (int ix, int iy, int iz); ///< gives the content of the cell[ix][iy][iz]
            /**
             * @brief returns the content of the i-th object, with the grid-coordinates ix,iy,iz
             * @param i number of the object
             * @param ix
            */
            T& operator () (int i, int ix, int iy, int iz, bool isObjectCoordinate = true);
            T& operator () (maths::Vector<int> Pi);
            T& operator () (int i, maths::Vector<int> Pi);
            T& operator () (maths::Vector<double> P);
            T& operator () (int i, maths::Vector<double> P);
            
            void fill(const T& x); ///< Fill the whole SuperArray with value \p x
            SuperGrid& operator = (const SuperGrid& S); ///< Assignment operator
            void clear(); ///< Clear the SuperGrid and release the allocated memory
            void allockugel();

            /**
            * @brief Checks if the point \p Pi (indicated by the indices) is inside the calculation sphere
            * The function returns the Vector Pi, if it is inside the calculation sphere otherwise a the vector (-1,-1,-1) is returned
            * (for internal use only)
            */
            maths::Vector<int> kugelindex(maths::Vector<int> Pi);
            T kugelwert(maths::Vector<int> Pi);
            T kugelwert(int ix, int iy, int iz);
            friend std::ostream& operator << (std::ostream& os, const SuperArray& S);

            /*
               Hier kommen die eigentlichen Daten:
               Ein : Array mit den Einschluessen (Orte und Ausdehnungen der Einschluesse werden absolut angegeben)
               numObjs : Anzahl Einschluesse
               type: Verteilung der aktiven Molek�le
                     0: nur im Einschluss
                     1: nur im Host
                     2: �berall
               ywerte,
               zwerte : Arrays zur Bestimmung der Dimensionen von K
               r0  : Radius des Host-Partikels=halbe Breite des Uebergitters
               G   : Array mit den Untergittern, wird ben�tigt falls type 0 ist; ObjectShape : G[Einschlussindex][ix][iy][iz]
               K   : Kugelgitter, wird ben�tigt falls type 1 oder 2 ist
               Pul : Array mit den Basispunkten der Untergitter
               n   : Array mit den Ausdehnungen der    "   n[Einschlussindex][nx][ny][nz] (nx,ny,nz sind als Vektor zusammengefasst)
               nges : Ausdehnung des Uebergitters als Vektor (nx,ny,nz) zusammengefasst
               d   : Gitterkonstanten als Vektor (dx,dy,dz) zusammengefasst
            */


            std::vector<ObjectShape*> Obj;      ///< List of the objects      
            int numObjs; ///< Number of objects
            int type;   ///< Type of grid (0: Field is stored only in the objects, 1: field is stored in the surroundings only, 2: field is stored in the whole volume                 
            int* ywerte;                           
            int** zwerte;
           // *std::vector<std:vector<std::vector<T> > > grid;
            std::vector<std::vector<std::vector<std::vector<T> > > > G;
            std::vector<std::vector<std::vector<T> > > K;
            T dummy;
            std::vector<GOAT::maths::Vector<int> >Pul;
            std::vector<GOAT::maths::Vector<int> >n;
            GOAT::maths::Vector<int> nges;            
            maths::Vector<double> d;
            double r0;
            bool isequal;
            bool iscleared;
            maths::Vector<std::complex<double> > pc;
            int Error;
            maths::Matrix<double> H, R;
        };

/* -------------------- IMPLEMENTATION  ---------------------- */




        template <class T>
        SuperGrid<T>::SuperGrid()
        {
            H = maths::unity();
            R = maths::unity();
            Error = SUPERGRID_NOERRORS;
            numObjs = 0;
            isequal = false;

            type = IN_HOST;
            ywerte = 0;
            zwerte = 0;
        }


        template <class T> SuperGrid<T>::SuperGrid(double r0, int nx, int ny, int nz, const int typ)
        {
            Error = SUPERGRID_NOERRORS;
            double b = 2.0 * r0;
           // G = 0;
            isequal = false;
            numObjs = 0;
      //      Obj = 0;
            ywerte = 0;
            zwerte = 0;
           // K = 0;
            this->r0 = r0;



            nges = maths::Vector<int>(nx, ny, nz);
            d = maths::Vector<double>(b / nx, b / ny, b / nz);
            iscleared = true;
            type = typ;
            if (type & IN_HOST)
                allockugel();
        }

        template <class T> SuperGrid<T>::SuperGrid(double r0, int nx, int ny, int nz, ObjectShape** Obj, int numObjs, const bool isAbsolute, const int typ)
        {
            Error = SUPERGRID_NOERRORS;
            double b = 2.0 * r0;
            isequal = isAbsolute;
            numObjs = 0;
            G = 0;
            this->r0 = r0;
            nges = maths::Vector<int>(nx, ny, nz);
            type = typ;
            d = maths::Vector<double>(b / nx, b / ny, b / nz);
            iscleared = true;
            if (type & IN_HOST)
            {
                allockugel();
            }
            addObject(Obj, numObjs, isAbsolute);
        }


        template <class T> bool SuperGrid<T>::addObject(ObjectShape** Obj, int numObjs, const bool isAbsolute)
        {
            bool ok = true;
            for (int i = 0; (i < numObjs) && (ok); i++) ok = addObject(Obj[i], isAbsolute);
            return ok;
        }

        template <class T> maths::Vector<int>  SuperGrid<T>::gitterpunkt(maths::Vector<double> P)
        {
            int ix, iy, iz, jx, jy, jz, jxh, jyh, anzx;
            maths::Vector<int> pi, ph;
            maths::Vector<double> h;
            h = H * P + maths::Vector<double>(r0, r0, r0);
            ph = maths::Vector<int>(floor(h[0] / d[0]), floor(h[1] / d[1]), floor(h[2] / d[2]));

            if (type & IN_HOST)
            {
                pi = ph;
            }
            else
            {
                pi = ph;
            }
            return pi;
        }

        template<class T> bool SuperGrid<T>::addPlane(maths::Vector<double> P, maths::Vector<double> norm, double d1, double d2, maths::Vector<double> cellSize)
        {
            GOAT::maths::Vector<double> e1, e2;
            if (abs(norm % GOAT::maths::ex) > 1E-5)
                e1 = GOAT::maths::ex - (GOAT::maths::ex * norm) * norm;
            else
                e1 = GOAT::maths::ey - (GOAT::maths::ey * norm) * norm;
            e1 = e1 / abs(e1);
            e2 = norm % e1;
            e2 = e2 / abs(e2);

            maths::trafo(e1, e2, norm, H, R);


            if (numObjs < 1)  // Es ist der erste Einschluss der hinzugef�gt wird
            {
                G = std::vector < std::vector < std:vector <std::vector> > >(1);
                Pul = std::vector < GOAT::maths::Vector<int> >(1);
                n = std::vector < GOAT::maths::Vector<int> >(1);
                Obj = std::vector <ObjectShape*>(1);
            }
            else
            {
                G.push_back(std::vector<std::vector<std::vector<T> > >(1));
                Pul.push_back(GOAT::maths::Vector<int>  (1));
                n.push_back(GOAT::maths::Vector<int>  (1));
                Obj.push_back(ObjectShape*);                
            }
            n[numObjs] = GOAT::maths::Vector<int>(d1/cellSize[0], d2/cellSize[1], 1);

            GOAT::maths::Vector<double> h;
            GOAT::maths::Vector<int> hn;
            h = ceil(ediv(E->por, d)) - floor(ediv(E->pul, d));
            hn = maths::Vector<int>((int)h[0], (int)h[1], (int)h[2]);
            Pul[numObjs] = maths::Vector<int>((int)h[0], (int)h[1], (int)h[2]);

            G[numObjs] = std::vector<std::vector<std::vector<T> > >(n[numObjs][0]);
            for (int ix = 0; ix < n[numObjs][0] + 1; ix++)
            {
                G[numObjs][ix] = std::vector<std::vector<T> >(n[numObjs][1]);
                for (int iy = 0; iy < n[numObjs][1] + 1; iy++)
                {
                    G[numObjs][ix][iy] = std::vector<T>(n[numObjs][2]);
                    /*for (int iz = 0; iz < n[numObjs][2] + 1; iz++)
                        G[numObjs][ix][iy][iz] = maths::Vector<std::complex<double> >(0.0, 0.0, 0.0);*/
                    G[numObjs][ix][iy][0] = maths::Vector<std::complex<double> >(0.0, 0.0, 0.0);
                }
            }
            numObjs++;

        }

        template<class T> bool SuperGrid<T>::addObject(ObjectShape* E, const bool isAbsolute)
        {
            int hmem, hmem2;
            char str[500];
            SysMemInfo smi;
            MemInfo mi;
            long int allocMem;
            maths::Vector<int> pul, por, N, hn;
            maths::Vector<double> h, O;
            O = maths::Vector<double>(-r0, -r0, -r0);

            smi = sysmem();
            mi = memstat();
            h = ceil(ediv(E->por, d)) - floor(ediv(E->pul, d));

            hn = maths::Vector<int>((int)h[0], (int)h[1], (int)h[2]); // Gr��e des 3D-Gitters in die drei Koordinatenrichtungen
            
                if (numObjs < 1)  // Es ist der erste Einschluss der hinzugef�gt wird
                {
                    G = std::vector < std::vector < std:vector <std::vector> > >(1);
                    Pul = std::vector < GOAT::maths::Vector<int> >(1);
                    n = std::vector < GOAT::maths::Vector<int> >(1);                    
                    Obj = std::vector <ObjectShape*>(1);
                }
                else
                {
                    G.push_back(std::vector<std::vector<std::vector<T> > >(1));
                    Pul.push_back(GOAT::maths::Vector<int>  (1));
                    n.push_back(GOAT::maths::Vector<int>  (1));
                    Obj.push_back(ObjectShape *);                    
                }
            
            n[numObjs] = hn;

            Obj[numObjs] = E;
            h = floor(ediv(Obj[numObjs]->pul + maths::Vector<double>(r0, r0, r0), d));
            Pul[numObjs] = maths::Vector<int>((int)h[0], (int)h[1], (int)h[2]);  // ???????????
            Ellipsoid* H = (Ellipsoid*)E;
            if (E->isActive())  // Ist der Einschluss �berhaupt inelastisch aktiv ? 
            {
                G[numObjs]=std::vector<std::vector<std::vector<T> > > (n[numObjs][0]);
                for (int ix = 0; ix < n[numObjs][0] + 1; ix++)
                {
                    G[numObjs][ix]=std::vector<std::vector<T> > (n[numObjs][1]);                    
                    for (int iy = 0; iy < n[numObjs][1] + 1; iy++)
                    {
                        G[numObjs][ix][iy] = std::vector<T>(n[numObjs][2]);                        
                        for (int iz = 0; iz < n[numObjs][2] + 1; iz++)
                            G[numObjs][ix][iy][iz] = maths::Vector<std::complex<double> >(0.0, 0.0, 0.0);
                    }
                }
            }
            // else G[numObjs] = NULL;

            int i = 0;
            numObjs++;
            iscleared = false;
            return true;
        }

        template <class T> bool SuperGrid<T>::inObject(maths::Vector<double> P, int i) // da muss noch was gemacht werden
        {
            return Obj[i]->isInside(P);
        }

        template <class T> bool SuperGrid<T>::inObject(maths::Vector<int> Pi, int i)
        {
            return (Obj[i]->isInside(emult(Pi, d)));
        }

        template <class T> bool SuperGrid<T>::inObject(int ix, int iy, int iz, int i)
        {
            return (Obj[i]->isInside(emult(maths::Vector<double>(ix, iy, iz), d)));
        }

        template <class T> T& SuperGrid<T>::operator () (int ix, int iy, int iz)
        {
            maths::Vector<int> Pi = maths::Vector<int>(ix, iy, iz);
            int i = 0;
            bool found = false;
            T dummy;
            if (type == IN_OBJECT)
            {
                do
                {
                    found = inObject(ix, iy, iz, i);
                    i++;
                } while ((i < numObjs) && (!found));
                i--;
                if (!found) return dummy;
                else
                {
                    Pi = Pi - Pul[i];
                    Pi = kugelindex(Pi);
                    if (Error == SUPERGRID_NOERRORS)
                        return G[i][Pi[0]][Pi[1]][Pi[2]];
                    else return dummy;
                }
            }
            else
            {
                Pi = kugelindex(Pi);
                if (Error == SUPERGRID_NOERRORS)
                    return K[Pi[0]][Pi[1]][Pi[2]];
                else return dummy;
            }
        }

        template<class T> T& SuperGrid<T>::operator () (maths::Vector<int> Pi)
        {
            T dummy;
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

                if (!found) return dummy;
                else
                {
                    Pi = Pi - Pul[i];
                    Pi = kugelindex(Pi);
                    if (Error == SUPERGRID_NOERRORS)
                        return G[i][Pi[0]][Pi[1]][Pi[2]];
                    else
                        return dummy;
                }
            }
            else
            {
                Pi = kugelindex(Pi);
                if (Error == SUPERGRID_NOERRORS)
                    return K[Pi[0]][Pi[1]][Pi[2]];
                else
                    return dummy;
            }
        }

        template <class T> T& SuperGrid<T>::operator () (int i, int ix, int iy, int iz, bool isEinKoord)
        {
            T dummy;
            maths::Vector<int> Pi = maths::Vector<int>(ix, iy, iz);
            if (type == IN_OBJECT)
            {
                if (G[i] == NULL) return dummy;
                if (!isEinKoord)  Pi = Pi - Pul[i];
                return G[i][Pi[0]][Pi[1]][Pi[2]];
            }
            else
            {
                Pi = kugelindex(Pi);
                if (Error == SUPERGRID_NOERRORS)
                    return K[Pi[0]][Pi[1]][Pi[2]];
                else
                    return dummy;
            }
        }

        template <class T> T& SuperGrid<T>::operator () (int i, maths::Vector<int> Pi)
        {
            dummy = maths::Vector<std::complex<double> >(0.0, 0.0, 0.0);

            if (type == IN_OBJECT)
            {
                Pi = Pi - Pul[i];
                if (G[i] == NULL) return dummy;
                if (Pi[0] < 0) return dummy; //maths::Vector<std::complex<double> > (0,0,0);
                if (Pi[1] < 0) return dummy; //maths::Vector<std::complex<double> > (0,0,0);
                if (Pi[2] < 0) return dummy; // maths::Vector<std::complex<double> > (0,0,0);
                if (Pi[0] > n[i][0])
                {
                    Error = SUPERGITTER;
                    return dummy; //maths::Vector<std::complex<double> > (0,0,0);
                }

                if (Pi[1] > n[i][1])
                {
                    Error = SUPERGITTER;
                    return dummy; //maths::Vector<std::complex<double> > (0,0,0);
                }

                if (Pi[2] > n[i][2])
                {
                    Error = SUPERGITTER;
                }

                Error = SUPERGRID_NOERRORS;
                return G[i][Pi[0]][Pi[1]][Pi[2]];
            }
            else
            {
                Pi = kugelindex(Pi);
                if (Error == SUPERGRID_NOERRORS)
                    return K[Pi[0]][Pi[1]][Pi[2]];
                else
                    return dummy;
            }
        }

        template<class T> T& SuperGrid<T>::operator () (maths::Vector<double> P)
        {
            int i;
            maths::Vector<int> Pi;
            maths::Vector<double> h;
            h = H * P - maths::Vector<double>(-r0, -r0, -r0);
            T dummy;
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
                
                if (!found) return dummy;
                else
                {
                    if (G[i] == NULL) return dummy;
                    Pi = Pi - Pul[i];
                    return G[i][Pi[0]][Pi[1]][Pi[2]];
                }
            }
            else
            {
                Pi = kugelindex(Pi);
                if (Error == SUPERGRID_NOERRORS)
                    return K[Pi[0]][Pi[1]][Pi[2]];
                else return dummy;
            }
        }

        template <class T> T& SuperGrid<T>::operator () (int i, maths::Vector<double> P)
        {
            maths::Vector<int> Pi;
            maths::Vector<double> h;
            h = H * P - maths::Vector<double>(-r0, -r0, -r0);
            T dummy;

            for (int j = 0; j < 3; j++)
                Pi[j] = h[j] / d[j];

            if (type == IN_OBJECT)
            {
                if (G[i] == NULL) return dummy;
                Pi = Pi - Pul[i];
                return G[i][Pi[0]][Pi[1]][Pi[2]];
            }
            else
            {
                Pi = kugelindex(Pi);
                if (Error == SUPERGRID_NOERRORS)
                    return K[Pi[0]][Pi[1]][Pi[2]];
                else return dummy;
            }
        }

        template <class T> void SuperGrid<T>::clear()
        {
         /*   int anzx, anzx2;
            anzx = nges[0];
            anzx2 = nges[0] / 2;

            if (type == IN_OBJECT)
            {
                if (numObjs > 0)
                {
                    for (int i = numObjs - 1; i >= 0; i--)
                    {
                        if (G[i] != NULL)
                        {
                            for (int ix = n[i][0]; ix >= 0; ix--)
                            {
                                for (int iy = n[i][1]; iy >= 0; iy--)
                                    if (G[i][ix][iy] != 0)
                                        free(G[i][ix][iy]);
                                if (G[i][ix] != 0)
                                {
                                    free(G[i][ix]);
                                } // if
                            } // for ix    
                            free(G[i]);
                        } // if (G[i]!=NULL)
                    } // for i
                    free(n);
                    free(Pul);
                    free(Obj);
                    numObjs = 0;
                    iscleared = true;
                    G = 0;
                }
            }
            else
            {
                if (K != 0)
                {
                    for (int k = 0; k < anzx2; k++)
                    {
                        for (int l = 0; l < ywerte[k]; l++)
                        {
                            delete[] K[anzx2 - 1 - k][ywerte[k] + l];
                            delete[] K[anzx2 + k][ywerte[k] + l];
                            delete[] K[anzx2 - 1 - k][ywerte[k] - 1 - l];
                            delete[] K[anzx2 + k][ywerte[k] - 1 - l];
                        }
                        delete[] K[anzx2 - 1 - k];
                        delete[] K[anzx2 + k];
                        delete[] zwerte[k];
                    }

                    delete[] ywerte;
                    delete[] zwerte;

                    if (numObjs > 0)
                        free(Obj);
                }
                K = 0;
                zwerte = 0;
                ywerte = 0;
                numObjs = 0;
                iscleared = true;
            }*/
        }

        template <class T> void SuperGrid<T>::copy(const SuperGrid<T>& S)  // MUSS DRINGEND GE�NDERT WERDEN !
        {
            for (int i = 0; i < numObjs; i++)
                if (S.G[i].size()>0)
                {
                    for (int ix = 0; ix < n[i][0]; ix++)
                        for (int iy = 0; iy < n[i][1]; iy++)
                            for (int iz = 0; iz < n[i][2]; iz++)
                                G[i][ix][iy][iz] = S.G[i][ix][iy][iz];
                }
        }

        template<class T> SuperGrid<T>& SuperGrid<T>::operator = (const SuperGrid<T>& S)
        {
            int i = 0;
            int anzx, anzx2;


            if (this == &S) return *this;
            clear();
            r0 = S.r0;
            d = S.d;
            nges = S.nges;
            type = S.type;

            anzx = nges[0]; anzx2 = nges[0] / 2;

            if (type == IN_OBJECT)
            {
                isequal = true; addObject(S.Obj, S.numObjs, true);
                isequal = true;
                if (S.numObjs != 0)
                    do
                    {
                        for (int ix = 0; ix < n[i][0]; ix++)
                            for (int iy = 0; iy < n[i][1]; iy++)
                                for (int iz = 0; iz < n[i][2]; iz++)
                                    G[i][ix][iy][iz] = S.G[i][ix][iy][iz];
                        i++;
                    } while (i < S.numObjs);
                    isequal = false;
            }
            else
            {
                allockugel();
                isequal = true; addObject(S.Obj, S.numObjs, true);
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



        template <class T> void SuperGrid<T>::add(const SuperGrid<T>& S)
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

        template <class T> void SuperGrid<T>::sub(const SuperGrid<T>& S)
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

        template<class T> std::ostream& operator << (std::ostream& os, const SuperGrid<T> & S)
        {
            /*os << "r0=" << S.r0 << std::endl;
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
                }*/
            return os;
        }

        template <class T> void SuperGrid<T>::fill(const T& x)
        {
            int anzx, anzx2;
            anzx = nges[0]; anzx2 = nges[0] / 2;

            if (type == IN_OBJECT)
            {
                for (int i = 0; i < numObjs; i++)
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


        template <class T> void SuperGrid<T>::allockugel()
        {
            int anzx, anzx2;
            double dx, dy, dz, hilf;
            anzx = nges[0];
            anzx2 = nges[0] / 2;
            // cout << "allockugel" << std::endl;
            ywerte = new int[anzx2];
            zwerte = new int* [anzx2];

            /*dx=d[0];
            dy=d[1];
            dz=d[2];*/
            dx = 2.0 / nges[0];
            dy = 2.0 / nges[1];
            dz = 2.0 / nges[2];

            // cout << "dx=" << dx << std::endl; 
            for (int i = 0; i < anzx2; i++)
            {
                ywerte[i] = int(ceil(sqrt(1.0 - (i * dx) * (i * dx)) / dy));
                zwerte[i] = new int[anzx2];
                for (int j = 0; j < anzx2; j++)
                {
                    hilf = 1.0 - (i * dx) * (i * dx) - (j * dy) * (j * dy);
                    //      cout << "hilf:" << hilf << std::endl;
                    if (hilf >= 0)
                        //        zwerte[i][j]=int(ceil(sqrt(1.0-(i*dx)*(i*dx)-(j*dy)*(j*dy))*nges[2]/2.0));
                        zwerte[i][j] = ceil(sqrt(hilf) / dz);
                    else
                        zwerte[i][j] = 1;

                    //      cout << "zwerte[" << i << "][" << j << "]:" <<  zwerte[i][j] << std::endl;
                }
            }
        }



      /*  template<class T> void SuperGrid<T>::makeReal()
        {
            for (int i = 0; i < numObjs; i++)
                for (int ix = 0; ix < n[i][0]; ix++)
                    for (int iy = 0; iy < n[i][1]; iy++)
                        for (int iz = 0; iz < n[i][2]; iz++)
                            for (int j = 0; j < 3; j++)
                                G[i][ix][iy][iz][j] = abs(G[i][ix][iy][iz][j]);
        }
        */
    }
}