#pragma once

namespace GOAT
{
    namespace raytracing
    {
        enum ScatteringType
        {
            None = -1,
            Gaussian = 0,
            UniformCone = 1,
            CosineCone = 2
        };
        /** @brief Interface for rough surfaces.
        * This interface class is used to allow access to the parameters
        * of `roughObject` without having to cast them.
        */
        class roughInterface
        {
        public:
            virtual double getSigma() const = 0;
            virtual void setSigma(double s) = 0;
            virtual double getThetaMax() const = 0;
            virtual void setThetaMax(double t) = 0;
            virtual void setScatteringType(ScatteringType type) = 0;
            virtual ScatteringType getScatteringType() const = 0;
			virtual maths::Vector<double> norm(const maths::Vector<double>& P) = 0;
        };
    }
}