#pragma once
#include "objectshape.h"
namespace GOAT
{
    namespace raytracing
    {
        typedef struct
        {
            double k; ///< conic constant
            double R; ///< Radius
            std::vector<double> A = { 0,0,0,0,0,0,0 };
        } asphericLensSide;

        typedef struct
        {
            asphericLensSide left,right;
            double height;
            double radius;
        } asphericLensParms;

        class asphericLens :
            public ObjectShape
        {
         public:
             asphericLens(asphericLensParms lensParms);
            void binWrite();
            void binRead();
            void scale(double sf);
            bool next(const maths::Vector<double>& p, const maths::Vector<double>& k,
                 maths::Vector<double>& pout);
            maths::Vector<double> norm(const maths::Vector<double>& P);
            bool isInside(const maths::Vector<double>& p);
            double volume();

        private:
            asphericLensParms lensParms;
        };
    }
}
