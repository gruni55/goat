#include "vector.h"
#include "box.h"
#include "triangle.h"
#pragma once
int triBoxOverlap(const Vector<double> boxcenter, const Vector<double> boxhalfsize, const Vector<double> triverts[3]);
int triBoxOverlap(double boxcenter[3], double boxhalfsize[3], double triverts[3][3]);
bool triangleAABBIntersection(const Box &B, const triangle &d);
inline double min3(double x, double y, double z)
{
	double erg=x;
	if (erg > y) erg = y; 
	if (erg > z) erg = z;
	return erg;
}

inline double max3(double x, double y, double z)
{
	double erg = x;
	if (erg < y) erg = y;
	if (erg < z) erg = z;
	return erg;
}