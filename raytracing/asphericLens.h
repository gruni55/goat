#pragma once
#include "objectshape.h"
namespace GOAT
{
    namespace raytracing
    {
        /**
        * @brief Structure to describe one side of an aspheric lens
        * The surface is described by 
        * \f$ z(r)=\frac{r^2}{R\left( 1+\sqrt{1-(1+\kappa)\frac{r^2}{R^2}}\right)}+\sum_l\alpha_lr^l\f$
        */
        typedef struct
        {
            double k; ///< conic constant
            double R; ///< Radius
            std::vector<double> A;
            double shift=0; ///< how much must the curve be shifted (with consideration of the offset)
            bool isPlano = false;///< true, if the corresponding side is plano
        } asphericLensSide;


        /**
        * @Structure, which holds the full information about the aspheric lens
        */
        typedef struct
        {
            asphericLensSide left,right;
            double offset;
            double radius;
        } asphericLensParms;
        
   

        /**
        * @brief Representation of aspheric lens
        * The left and the right side is described by a formula \f$ z(r)=\frac{r^2}{R\left( 1+\sqrt{1-(1+\kappa)\frac{r^2}{R^2}}\right)}+\sum_l \alpha_l r^l \f$
        * or the side can also be plano. 
        */

        class asphericLens :
            public ObjectShape
        {
         public:
             asphericLens(const maths::Vector<double>& P,
                 std::complex<double>  n, asphericLensParms lensParms, 
                 const maths::Matrix<std::complex<double> > alpha = maths::CUNITY,
                 const maths::Vector<double>& Ex = maths::ex,
                 const maths::Vector<double>& Ey = maths::ey,
                 const maths::Vector<double>& Ez = maths::ez);
            void binWrite();
            void binRead();
            void scale(double sf);
            bool next(const maths::Vector<double>& p, const maths::Vector<double>& k,
                 maths::Vector<double>& pout);
            maths::Vector<double> norm(const maths::Vector<double>& P);
            bool isInside(const maths::Vector<double>& p);
            double volume();///< not yet implemented
            void binWrite(std::ofstream& os) { } ///< not yet implemented
            void binRead(std::ifstream& os) { } ///< not yet implemented
            void initQuad() { } ///< not yet implemented
            void setr0(double r0) { } ///< not yet implemented
            void setPos(maths::Vector<double> r) { } ///< not yet implemented
            void setPos(double x, double y, double z) { } ///< not yet implemented
            maths::Vector<double> calcCoM() {
                return maths::Vector<double> (0, 0, 0);
            } ///< not yet implemented

         
            enum Side {left,right,lateral};
            asphericLensParms lensParms;
            double zleft(maths::Vector<double> P, double lambda);
            double zright(maths::Vector<double> P, double lambda);
            double z(double r, asphericLensSide side);

            double dzdrleft(maths::Vector<double> P, double lambda);
            double dzdrright(maths::Vector<double> P, double lambda);
         
/*
            double dzright(double lambda);
            double dzleft(double lambda);
            */
            double drdlambda;
            maths::Vector<double> P, k;
           maths::Vector<double> currentnorm;
        };
    }
}
