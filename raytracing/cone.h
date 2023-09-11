#pragma once
#include "objectshape.h"

namespace GOAT
{
	namespace raytracing
	{
		/**
		 * @brief This class represents a cone.
		 * The cone is defined by the reference point (position), the height and the radius
		 * \image html cone.png
		 * \image latex cone.png "Cone" width=\textwidth
		 * 
		 */
		class Cone : public ObjectShape
		{
		public:
			Cone();
			/**
			 * @brief Constructor which defines the cone by the the reference point, radius and height			 
			 * @param Pos: Reference point in the center of the base (acts as position of the cone)
			 * @param radius: Radius of the base area
			 * @param height: Height of the cone
			 * The cone is oriented in z-direction.
			 */
			Cone(
				maths::Vector<double> Pos, 
				double radius, 
				double height,
				std::complex<double>  n,
				double r0 = 1.0,
				const maths::Matrix<std::complex<double> > alpha = maths::CUNITY,
				const maths::Vector<double>& Ex = maths::ex,
				const maths::Vector<double>& Ey = maths::ey,
				const maths::Vector<double>& Ez = maths::ez
			);
			~Cone();

			void binWrite(std::ofstream& os);
			void binRead(std::ifstream& os);
			void scale(double sf);
			/**
			 * @brief Calculates next crossing point with ray.
			 * This method calculates the next crossing point with the entire cone. 
			 * \param[in] p Reference point of the ray 
			 * \param[in] k Directional vector of the ray (normalized)
			 * \param[out] pout crossing point (if found)
			 * \return True, if a crossing point was found
			 */
			bool next(const maths::Vector< double >& p, const maths::Vector< double >& k, maths::Vector< double >& pout);

			/**
			 * @brief Calculates the next crossing point with the lateral surface of the cone.
			 * The method calculates the next crossing point with the lateral surface of the cone according to the algorithm described
			 * in the article: Ching-Kuang Shene, "V.1 - Computing the Intersection of a Line and a Cone" in Graphics Gems V, p. 227--231, Academic Press, 1995.
			 * \param[in] p Reference point of the ray 
			 * \param[in] k Directional vector of the ray (normalized)
			 * \param[out] pout crossing point (if found)
			 * \return distance of p and the crossing point if found, otherwise -1 
			 */
			double nextCone(const maths::Vector< double >& p, const maths::Vector< double >& k, maths::Vector< double >& pout); 
			maths::Vector< double > norm(const maths::Vector< double >& P);
			
			bool isInside(const maths::Vector< double >& p);
			double volume();
			void initQuad();
			void setPos(maths::Vector< double > r);
			void setPos(double x, double y, double z);
			void setr0(double r0);
			maths::Vector<double> calcCoM();

		private:
			maths::Vector<double> V; ///< Position of the tip of the cone
			maths::Vector<double> v; ///< Vector from the tip cone to the reference point Pos
			maths::Vector<double> normv; ///< normalized Vector from the tip cone to the reference point Pos
			double coneAngle; ///< cone angle
			double cosCA; ///< cosine of the cone angle
			double tan2CA; ///< square of the tangens of the cone angle
			double sideLen; ///< length of the side of the cone
			double radius; ///< radius of the base area
			double height; ///< height of the cone

			void init(); ///< Initialisation of all internal variables 
		};

		
	}
}