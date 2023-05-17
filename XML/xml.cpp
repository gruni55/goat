
#include "ticpp.h"
#include "xml.h"


namespace GOAT
{
	namespace XML
	{
		void xmlReader::readXML(std::string fname, GOAT::raytracing::Scene& S)
		{
			ticpp::Document doc(fname);
			try
			{
				doc.LoadFile();

			}
			catch (std::string err)
			{
				std::cerr << err << std::endl;				
			}
				rootElement = doc.FirstChildElement();
				ticpp::Element* ell = rootElement->FirstChildElement();;
				std::string str;
				do
				{
					
					if (ell != 0)
					{
						str = ell->Value();
						switch (checkRootElements(ell->Value()))
						{
						    case XML_LIGHTSOURCES: readLightSrc(ell); break;
							case XML_OBJECTS: break;
							case XML_DETECTORS: break;
						}
					}
					ell = ell->NextSiblingElement();
				} while (ell != 0);


		}

		int xmlReader::checkRootElements(std::string str)
		{
			for (int i = 0; i < numXMLRootElements; i++)
				if (rootXMLElements[i].compare(str) == 0) return i;
			return XML_NONE;
		}
		
		int xmlReader::checkLSElements(std::string str)
		{
			for (int i = 0; i < numXML_LS_TYPES; i++)
				if (LSXMLTYPES[i].compare(str) == 0) return i;
			return XML_NONE;
		}

		GOAT::maths::Vector<double> xmlReader::readDVector(ticpp::Element* ell)
		{
			std::string str;
			GOAT::maths::Vector<double> Pos;
			if (ell != 0)
			{
				double x, y, z;

				str = ell->GetAttribute("x");
				if (str.empty()) x = 0.0;
				else x = std::stoi(str);

				str = ell->GetAttribute("y");
				if (str.empty()) y = 0.0;
				else y = std::stoi(str);

				str = ell->GetAttribute("z");
				if (str.empty()) z = 0.0;
				else z = std::stoi(str);

				Pos = GOAT::maths::Vector<double>(x, y, z);
			}
			return Pos;
		}

		void xmlReader::readLightSrcDefaults(ticpp::Element* ell, double& size, int& numRays, double& wvl, GOAT::maths::Vector<double>& Pos)
		{
			std::string str;
			str = ell->GetAttribute("Wavelength");
			if (str.empty()) wvl = 1.0;
			else wvl = std::stod(str);

			str = ell->GetAttribute("Size");
			if (str.empty()) size = 1.0;
			else size = std::stod(str);

			str = ell->GetAttribute("NumRays");
			if (str.empty()) numRays = 100;
			else numRays = std::stoi(str);
			
			ticpp::Element *ellPos = ell->NextSiblingElement("Position");
			Pos = readDVector(ellPos);
		}

		int xmlReader::readLightSrc(ticpp::Element *ell)
		{
			GOAT::raytracing::LightSrc *LS;
			GOAT::maths::Vector<double> focuspos;
			GOAT::maths::Vector<double> Pos;
			double wavelength;
			int numRays;
			double size;
			readLightSrcDefaults(ell, size, numRays, wavelength, Pos);
			std::string str;
			int result = checkLSElements(ell->GetAttribute("Type"));
			
			switch (result)
			{
			case XML_LS_TYPE_PLANE: LS = new GOAT::raytracing::LightSrcPlane(Pos, numRays, wavelength, size); S.addLightSource(LS);  break;
				case XML_LS_TYPE_GAUSS: 
										double w0;
										str = ell->GetAttribute("w0");
										if (str.empty()) w0 = 1.0;
										else w0 = std::stod(str);
										
										focuspos = readDVector(ell);

					                    LS = new GOAT::raytracing::LightSrcGauss(Pos, numRays, wavelength, 1.0, focuspos, size);

										str = ell->GetAttribute("NA");
										if (!str.empty()) ((GOAT::raytracing::LightSrcGauss_mc*)LS)->setNA(std::stod(str));
										S.addLightSource(LS);										
					                     break;

				case XML_LS_TYPE_PLANE_MC: break;
				case XML_LS_TYPE_GAUSS_MC:break;
			}
			return result;
		}
	}
}