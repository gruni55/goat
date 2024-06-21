#include "xml.h"
#include<algorithm>
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
			for (i = 0; (i < numLightSourceToken) && (!found); i++)
				found = (str.compare(objectToken[i])==0);
			if (found) return (i-1) + 100;
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

		std::string str_tolower(std::string s)
		{
			std::transform(s.begin(), s.end(), s.begin(),				
				[](unsigned char c) { return std::tolower(c); } // correct
			);
			return s;
		}
	}
}