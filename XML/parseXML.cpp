#include "parseXML.h"

namespace GOAT
{
	namespace XML
	{
		xmlParser::xmlParser()
		{}

		xmlParser::xmlParser(const raytrace::Scene &S)
		{
			this->S = S;
		}

		void xmlParser::read(std::string str)
		{
			XMLDocument doc;
			doc.LoadFile(str);
			XMLElement* pRootElement = doc.RootElement();
			if (pRootElemnt != NULL)
			{
				XMLElement* pLightSrcs = pRootElement->FirstChildElement("Light sources");
				while (pLightSrcs)
				{
					XMLElement* pLightSrc = pLightSrcs->FirstChildElement("Gaussian");
					if (pLightSrc != NULL)
					{

					}
				}
			}
 			
			
		}

		void xmlParser::write(std::string str)
		{
			document.LoadFile(str);
			*tinyxml2::XMLElement rootElement = document.RootElement();
			document.NewElement()
				
		}
	}
}
