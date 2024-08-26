#include <array>
#include <string>
#include "lightsrc.h"
#pragma once

namespace GOAT
{
	namespace XML
	{
#define numLightSourceToken 8
#define numObjectToken	    6 
#define numDetectorToken    1
#define numCalculationToken 5
#define numRefractiveIndexToken 5

#define TOKEN_NOT_FOUND					-1
#define TOKEN_LIGHTSOURCE_PLANE			0
#define TOKEN_LIGHTSOURCE_GAUSSIAN		1
#define TOKEN_LIGHTSOURCE_RING			2
#define TOKEN_LIGHTSOURCE_TOPHAT		3
#define TOKEN_LIGHTSOURCE_PLANE_MC		4
#define TOKEN_LIGHTSOURCE_GAUSSIAN_MC	5
#define TOKEN_LIGHTSOURCE_RING_MC		6
#define TOKEN_LIGHTSOURCE_GAUSSIAN_RING_MC 7

#define TOKEN_OBJECT_ELLIPSOID		  100
#define TOKEN_OBJECT_SURFACE		  101
#define TOKEN_OBJECT_CONE			  102
#define TOKEN_OBJECT_ASPHERIC_LENS	  103
#define TOKEN_OBJECT_SPHERIC_LENS	  104
#define TOKEN_OBJECT_BOX			  105

#define TOKEN_DETECTOR_PLANE 		  150

#define TOKEN_CALCULATION_PURE		  200
#define TOKEN_CALCULATION_PATH		  201
#define TOKEN_CALCULATION_PULSE		  202
#define TOKEN_CALCULATION_PULSE_FIELD 203
#define TOKEN_CALCULATION_INELASTIC   204

#define TOKEN_REFRACTIVE_INDEX_AIR	  300
#define TOKEN_REFRACTIVE_INDEX_GLASS  301
#define TOKEN_REFRACTIVE_INDEX_BK7	  302
#define TOKEN_REFRACTIVE_INDEX_LASF55 303
#define TOKEN_REFRACTIVE_INDEX_VACUUM 304




        const std::vector<std::string> lightSourceToken = { "plane","gaussian","ring","tophat","plane_mc","gaussian_mc","ring_mc","gaussian_ring_mc" };
		const std::vector<std::string> objectToken = { "ellipsoid","surface","cone","aspheric_lens","spheric_lens","box"};
		const std::vector<std::string> detectorToken = {"plane"};
        const std::vector<std::string> calculationToken = { "pure","path","pulse","pulse_field","inelastic"};
		const std::vector<std::string> refractiveIndexToken = { "air","glass","bk7","lasf55","vacuum" };

		int mapString2LightSourceToken(std::string str);
		int mapString2ObjectToken(std::string str);
		int mapString2DetectorToken(std::string str);
		int mapString2CalculationToken(std::string str);
		int mapString2RefractiveIndexToken(std::string str);
		std::string str_tolower(std::string s); ///< converts all letters in s into lower case (taken from https://en.cppreference.com/w/cpp/string/byte/tolower)
		bool addFunction2IndexList(std::vector< std::function< std::complex< double >(double) > >& nList, int refIndexToken);
	}
}

