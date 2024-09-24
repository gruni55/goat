#include "xmltoken.h"
#include "xml.h"
#include<algorithm>
#include "refractive_index_functions.h"
#include "xmltoken.h"

namespace GOAT
{
	namespace XML
	{
		int mapString2LightSourceToken(std::string tstr)
		{
			std::string str = str_tolower(tstr);
			bool found = false;
			int i;
			for (i = 0; (i < numLightSourceToken) && (!found); i++)
				found = (str.compare(lightSourceToken[i])==0);
			if (found) return i-1;
			return TOKEN_NOT_FOUND;
		}
		
		int mapString2ObjectToken(std::string tstr)
		{
			std::string str = str_tolower(tstr);
			bool found = false;
			int i;
			for (i = 0; (i < numObjectToken) && (!found); i++)
				found = (str.compare(objectToken[i])==0);
			if (found) return (i-1) + 100;
			return TOKEN_NOT_FOUND;
		}

		int mapString2DetectorToken(std::string tstr)
		{
			std::string str = str_tolower(tstr);
			bool found = false;
			int j;
			for (j = 0; (j < numDetectorToken) && (!found); j++)
				found = (str.compare(detectorToken[j]) == 0);
			if (found) return (j - 1) + 150;
			return TOKEN_NOT_FOUND;
		}

		int mapString2CalculationToken(std::string tstr)
		{
			std::string str = str_tolower(tstr);
			bool found = false;
			int i;
			for (i = 0; (i < numCalculationToken) && (!found); i++)
				found = (str.compare(calculationToken[i]) == 0);
			if (found) return (i - 1) + 200;
			return TOKEN_NOT_FOUND;
		}

		int mapString2RefractiveIndexToken(std::string tstr)
		{
			std::string str = str_tolower(tstr);
			bool found = false;
			int i;
			for (i = 0; (i < numRefractiveIndexToken) && (!found); i++)
				found = (str.compare(refractiveIndexToken[i]) == 0);
			if (found) return (i - 1) + 300;
			return TOKEN_NOT_FOUND;
		}

		std::string str_tolower(std::string s)
		{
			std::transform(s.begin(), s.end(), s.begin(),				
				[](unsigned char c) { return std::tolower(c); } // correct
			);
			return s;
		}

		bool addFunction2IndexList(std::vector<std::function<std::complex<double>(double)>>& nList, int refIndexToken)
		{
			switch (refIndexToken)
			{
				case TOKEN_REFRACTIVE_INDEX_AIR: nList.push_back(GOAT::raytracing::n_Air); break;
				case TOKEN_REFRACTIVE_INDEX_VACUUM: nList.push_back(GOAT::raytracing::n_Vacuum); break;
				case TOKEN_REFRACTIVE_INDEX_BK7: nList.push_back(GOAT::raytracing::n_BK7); break;
				case TOKEN_REFRACTIVE_INDEX_GLASS: nList.push_back(GOAT::raytracing::n_Glass); break;
				case TOKEN_REFRACTIVE_INDEX_LASF55: nList.push_back(GOAT::raytracing::n_LASF55); break;									
				default: std::cerr << "Pulse calculation: Wrong refractive index function name! Calculation stopped!" << std::endl; return false;
			}
			return true;
		}
	}
}
