#pragma once
#include "parseXML.h"
#include "raytrace.h"
#include "tinyxml2.h"

namespace GOAT
{
	namespace XML
	{
		class xmlParser
		{
			public:
				xmlParser();
				xmlParser(const raytrace::Scene& S);
				void read(std::string str);
				void write(std::string str);
			private:
				Scene S;
				tinyxml2::XMLDocument document;
		};
	}
}