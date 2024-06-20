#include "xml.h"
namespace GOAT
{
	namespace XML
	{
		int mapString2LightSourceToken(std::string str)
		{
			bool found = false;
			int i;
			for (i = 0; (i < numLightSourceToken) && (!found); i++)
				found = (str.compare(lightSourceToken[i])==0);
			if (found) return i-1;
			return TOKEN_NOT_FOUND;
		}

		int mapString2ObjectToken(std::string str)
		{
			bool found = false;
			int i;
			for (i = 0; (i < numLightSourceToken) && (!found); i++)
				found = (str.compare(objectToken[i])==0);
			if (found) return (i-1) + 100;
			return TOKEN_NOT_FOUND;
		}
	}
}