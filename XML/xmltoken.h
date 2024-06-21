#include <array>
#include <string>
#include "lightsrc.h"
#pragma once

namespace GOAT
{
	namespace XML
	{
#define numLightSourceToken 7
#define numObjectToken	    5 
#define numCalculationToken 3
#define TOKEN_NOT_FOUND					-1
#define TOKEN_LIGHTSOURCE_PLANE			0
#define TOKEN_LIGHTSOURCE_GAUSSIAN		1
#define TOKEN_LIGHTSOURCE_RING			2
#define TOKEN_LIGHTSOURCE_TOPHAT		3
#define TOKEN_LIGHTSOURCE_PLANE_MC		4
#define TOKEN_LIGHTSOURCE_GAUSSIAN_MC	5
#define TOKEN_LIGHTSOURCE_RING_MC		6

#define TOKEN_OBJECT_ELLISPOID		  100
#define TOKEN_OBJECT_SURFACE		  101
#define TOKEN_OBJECT_CONE			  102
#define TOKEN_OBJECT_ASPHERIC_LENS	  103
#define TOKEN_OBJECT_SPHERIC_LENS	  104
#define TOKEN_OBJECT_BOX			  105

#define TOKEN_CALCULATION_PURE		  200
#define TOKEN_CALCULATION_PATH		  201
#define TOKEN_CALCULATION_PULSE		  202

		const std::vector<std::string> lightSourceToken = { "plane","gaussian","ring","tophat","plane_mc","gaussian_mc","ring_mc" };
		const std::vector<std::string> objectToken = { "ellipsoid","surface","cone","aspheric_lens","spheric_lens" };	
		const std::vector<std::string> calculationToken = { "pure","path","pulse" };

		int mapString2LightSourceToken(std::string str);
		int mapString2ObjectToken(std::string str);
		int mapString2CalculationToken(std::string str);
		std::string str_tolower(std::string s); ///< converts all letters in s into lower case (taken from https://en.cppreference.com/w/cpp/string/byte/tolower)
	}
}

