#ifndef SURFACE_H
#define SURFACE_H

#include <complex>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <math.h>
#include "vector.h"
#include "objectshape.h"
#include <string.h>
#include "triangle.h"
#ifdef WITH_OCTREE
#include "octree.h"
#endif


namespace GOAT
{
    namespace raytracing
    {
#ifndef SQRT3
#define SQRT3 1.73205080756887729352744634150587236694280525381038062805580  /// SQRT(3)
#endif



        /* Klasse, die eine ObjectShape repräsentiert, die aus Dreiecken zusammengebaut ist. Die Dreiecke werden durch die Klasse \ref refdreieck "triangle"
        beschrieben */
        /**
        @brief This class represents objects those surface is described by triangles

        The surface class represents a volume object, whose surface is described by triangles. This method is often used by CAD programs, therefore a function
        is provided to import binary stl files. The surface class has also a simple proprietary file format:
        The first entry is the number of triangles in the file. Than, the information about the triangles follows. Each triangle is described by the components
        of the corner vectors. Therefore the file looks like:
        \f$
        \newline
        \begin{array}{ccccccccc}
        \text{<number of triangles>} & & & & & & & & \\
        \text{<number of triangles>} & & & & & & & & \\
        x_{00} & y_{00} & z_{00} & x_{01} & y_{01} & z_{01} &  x_{02} &  y_{02} & z_{02} \\
        x_{10} & y_{10} & z_{10} & x_{11} & y_{11} & z_{11} &  x_{12} &  y_{12} & z_{12} \\
        \vdots & \vdots  & & & & & & &
        \end{array}
        \f$
        The intersection points are calculated with help of an octree procedure to enhance the calculation speed for the intersection point finding.
        */
        class surface : public ObjectShape
        {
        public:

            // Konstruktoren und Destruktor
            //   surface(int num); ///< Constructor for a surface object with num empty triangles
            surface();
            surface(maths::Vector<double> Oh); ///< Empty surface object with position Oh
            // surface(const ObjectShape &); 
            surface(const surface& Su);
            /**
             * Constructor for an empty surface object.
             *
             * \param P position of the object
             * \param n refractive index
             * \param alpha polarisability
             * \param Ex direction of the first axis of the local coordinate system (usually ex)
             * \param Ey direction of the second axis of the local coordinate system (usually ey)
             * \param Ez direction of the third axis of the local coordinate system (usually ez)
             */
            surface(const maths::Vector<double>& P,
                std::complex<double>  n,
                const maths::Matrix<std::complex<double> > alpha = maths::CUNITY,
                const maths::Vector<double>& Ex = maths::ex,
                const maths::Vector<double>& Ey = maths::ey,
                const maths::Vector<double>& Ez = maths::ez
            );

            /**
             * Constructor for an object, whose surface is described by a list of triangles
             *
             * \param P position of the object
             * \param n refractive index
             * \param num number of triangles
             * \param list list of triangles
             * \param alpha polarisability
              * \param Ex direction of the first axis of the local coordinate system (usually ex)
             * \param Ey direction of the second axis of the local coordinate system (usually ey)
             * \param Ez direction of the third axis of the local coordinate system (usually ez)
             */
            surface(const maths::Vector<double>& P,
                std::complex<double>  n,
                int num, triangle* list,
                const maths::Matrix<std::complex<double> > alpha = maths::CUNITY,
                const maths::Vector<double>& Ex = maths::ex,
                const maths::Vector<double>& Ey = maths::ey,
                const maths::Vector<double>& Ez = maths::ez);

            // ~surface();
            maths::Vector<double> calcCoM(); /// calculates center of mass
            /**
             * Returns the smallest and largest edge length of the triangles stored in the surface object
             * @param min smallest edge length
             * @param max largest edge length
             */
            void getMinMax(double& min, double& max);
            triangle* S = 0; ///< list of all triangles
            int numTriangles = 0; ///< number of triangles 

            // Erzeugen der Dreiecksliste
            // 1. interaktiv
            int createsurface(); ///< Creates a list of triangles, which describes the surface, from standard input
            // Skaliere die Dreiecke
            /**
             * @brief scale the lengths of all triangles
             * @param sf scaling factor
             */
            void scale(double sf);

            // 2. Einlesen aus Datei FName
            int createsurface(std::string FName); ///< reads the triangle list from a file (proprietary file format, SRF)
            int importBinSTL(std::string FName); ///< import triangle list from binary STL-file
            void exportSRF(std::string FName);   ///< export triangle list to proprietary file format SRF
            // Andere Hilfsfunktionen
            // 1. Rückgabe einer Klasse mit leerem S
            surface nosurface(); ///< returns copy of the object without triangles
            // 2. Hinzufügen einer Liste S 
            /**
             * @brief Adds triangles to the internal list of triangles.
             * @param  S list of triangles
             * @param  num number of triangles to add
             */
            void addTriangle(triangle* S, int num);

            void clearS(); ///< clears the internal list of triangles
            //  string getFName() {return FName; }
            std::string getFName() { return FName; } ///< returns the current file name
            void setFilename(std::string FName) ///< set the current file name
            {
                this->FName = FName;
                //if (this->FName!=0) delete[] this->FName;
                // this->FName=new char[strlen(FName)+1];
                // strcpy(this->FName,FName);
            }
            // Operatoren 
            surface& operator=(const surface&);

            // Die virtuellen Funktionen der ObjectShape-Klasse
              /**
               * @brief Calculates the intersection point between a straight line and the surface object.
               * The straight line is represented by a point on the line and its direction
               * @param r point on the line
               * @param direction of the line
               * @param p calculated intersection point (zero, if there is no intersection at all)
               * @return true if there is an intersection point otherwise false
               */
            bool next(const maths::Vector<double>& r, const maths::Vector<double>& k, maths::Vector<double>& p);
            maths::Vector<double> norm(const maths::Vector<double>& P); ///< Surface normal at point P (There is no check, if P is on the surface. This function is used for internal purposes)
            bool isInside(const maths::Vector<double>& p);
            void setr0(double r0); ///< sets the radius of the calculation sphere (necessary only if this class is used outside Scene)
            void initQuad(); ///< Boundaries of the circumferent cuboid is calculated. The cuboid is represented by the upper right corner (por) and the lower left corner (pul). Used only for inelastic scattering
            maths::Matrix<double> computeInertia(); ///< Calculates the inertia of the object
            void setPos(maths::Vector<double> r); ///< Set position of the object to r.
            void setPos(double x, double y, double z) { setPos(maths::Vector<double>(x, y, z)); } ///< Set position of the object to a position vector represented by its components x,y and z
            //double isInHost(void);
            void binWrite(std::ofstream& os); ///< Writes object to binary file represented by ostream os
            void binRead(std::ifstream& is); ///< Reads object from binary file represented by istream is
            // double volume(); 
            double volume();///< Calculates the volume of the object
            void setCenter(maths::Vector<double> P);
            void setCenter2CoM(); ///< Sets the Position to the center of mass (CoM)
            int getCurrentIndex() { return currentIndex; } ///< Returns the index of the triangle that was last hit.
            triangle& getTriangle(int i) { return S[i]; } ///< Returns the i-th triangle in the internal triangle list.

            // protected:
            maths::Vector<double> currentnorm;  ///< Normal vector of the triangle that was last hit
            int currentIndex; ///< Index of the triangle that was last hit.
            std::string FName; ///< File name (used to save the internal triangle list)
            // char *FName ;
            void initBounds(maths::Vector<double>& pul, maths::Vector<double>& por); ///< Calculates the bounding box (=circumferent cuboid) but without rotation. The cuboid is represented by the upper right corner (por) and the lower left corner (pul).

#ifdef WITH_OCTREE
            // #include "Octree.hpp"
            //

            Octree<triangle> Tree;
            // Octree<triangle>
          //::Section cube;
#endif
        };

        std::ostream& operator << (std::ostream& os, const surface& su);

        /** @name operators on surface
         * @brief Operators, which act directly on the surface class
         */
         ///{
        surface operator + (const surface& s, const maths::Vector<double>& v); ///< adds vector v to all corners of the triangle list of the surface object s
        surface operator - (const surface& s, const maths::Vector<double>& v); ///< subtracts vector v from all corners of the triangle list of the surface object s
        surface operator + (const maths::Vector<double>& v, const surface& s); ///< adds vector v to all corners of the triangle list of the surface object s
        surface operator - (const maths::Vector<double>& v, const surface& s); ///< subtracts vector v from all corners of the triangle list of the surface object s
        surface operator * (const maths::Matrix<double>& M, const surface& s); ///< multiply every corner of the triangles from surface s by the matrix M
        ///}

        /**
         * @brief Generates a cylinder with hexagonal cross section
         * @param a side length
         * @param h height
         * @param M rotation matrix
         */
        surface generateHexagonCylinder(double a, double h, maths::Matrix<double> M = maths::UNITY);

        /**
         * @brief generates a cylinder with hexagonal cross section. At the upper and lower faces, there are pyramidal dips.
         * @param a side length
         * @param h height
         * @param t depth of the dip
         * @param M rotation matrix
         */
        surface generateHollowHexagonalCylinder(double a, double h, double t, maths::Matrix<double> M = maths::UNITY);

        /**
         * @brief generates a cylinder with hexagonal cross section with a pyramid-shaped tip on the upper surface.
         * @param a side length
         * @param h total height
         * @param t height of the tip
         * @param M rotation matrix
         */
        surface& generateBullet(double a, double h, double t, maths::Matrix<double> M = maths::UNITY); /// Generiert einen hexagonalen Zylinder mit Seitenlänge a und Höhe h und einer pyramidalen Spitze 

        /**
         * @brief generates a triangulated ellipsoid
         * @param a first half axis
         * @param b second half axis
         * @param N resolution
         */
        surface generateEllipsoid(double a, double b, int N, double r0 = 1.0, maths::Matrix<double> M = maths::UNITY); // generiert einen Rotationsellipsoiden mit Halbachsen a (z-Achse) und b

        /**
         * @brief Generates a cylinder with ellipsoidal caps
         * @param  a hei
         */
        surface generatePill(double a, double b, double h, int N, double r0 = 1.0, maths::Matrix<double> M = maths::UNITY);

        surface generatePolyBullet(int n, double a, double h, double t, maths::Matrix<double> M = maths::UNITY);

        /**
         * @brief calculates the longest and shortest side length
         * @param numTriangles number of triangles
         * @param triangle List of the triangles
         * @return min shortest side length
         * @return max longest side length
         */
        void getMinMax(int numTriangles, triangle *S, double& min, double& max);
    }
}
#endif

