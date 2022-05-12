#include "vector.h"
#include "box.h"
#include "triangle.h"
#include "triangle_box_intersection.h"

/********************************************************/
/* AABB-triangle overlap test code                      */
/* by Tomas Akenine-M�ller                              */
/* Function: int triBoxOverlap(double boxcenter[3],      */
/*          double boxhalfsize[3],double triverts[3][3]); */
/* History:                                             */
/*   2001-03-05: released the code in its first version */
/*   2001-06-18: changed the order of the tests, faster */
/*                                                      */
/* Acknowledgement: Many thanks to Pierre Terdiman for  */
/* suggestions and discussions on how to optimize code. */
/* Thanks to David Hunt for finding a ">="-bug!         */
/********************************************************/

#include <math.h>
#include <stdio.h>

#define X 0
#define Y 1
#define Z 2

#define CROSS(dest,v1,v2) \
dest[0] = v1[1] * v2[2] - v1[2] * v2[1]; \
dest[1] = v1[2] * v2[0] - v1[0] * v2[2]; \
dest[2] = v1[0] * v2[1] - v1[1] * v2[0];



#define DOT(v1,v2) (v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])

#define SUB(dest,v1,v2) \
dest[0] = v1[0] - v2[0]; \
dest[1] = v1[1] - v2[1]; \
dest[2] = v1[2] - v2[2];



#define FINDMINMAX(x0,x1,x2,min,max) \
min = max = x0;   \
if (x1<min) min = x1; \
	if (x1>max) max = x1; \
		if (x2<min) min = x2; \
			if (x2>max) max = x2;



int planeBoxOverlap(double normal[3], double vert[3], double maxbox[3])	// -NJMP-
{
	int q;
	double vmin[3], vmax[3], v;
	for (q = X; q <= Z; q++)
	{
		v = vert[q];					// -NJMP-
		if (normal[q]>0.0f)
		{
			vmin[q] = -maxbox[q] - v;	// -NJMP-
			vmax[q] = maxbox[q] - v;	// -NJMP-
		}
		else
		{
			vmin[q] = maxbox[q] - v;	// -NJMP-
			vmax[q] = -maxbox[q] - v;	// -NJMP-
		}
	}
	if (DOT(normal, vmin)>0.0f) return 0;	// -NJMP-
	if (DOT(normal, vmax) >= 0.0f) return 1;	// -NJMP-
	return 0;
}





/*======================== X-tests ========================*/

#define AXISTEST_X01(a, b, fa, fb)			   \
p0 = a*v0[Y] - b*v0[Z];			       	   \
p2 = a*v2[Y] - b*v2[Z];			       	   \
if (p0<p2) { min = p0; max = p2; } \
else { min = p2; max = p0; } \
rad = fa * boxhalfsize[Y] + fb * boxhalfsize[Z];   \
if (min>rad || max<-rad) return 0;



#define AXISTEST_X2(a, b, fa, fb)			   \
p0 = a*v0[Y] - b*v0[Z];			           \
p1 = a*v1[Y] - b*v1[Z];			       	   \
if (p0<p1) { min = p0; max = p1; } \
else { min = p1; max = p0; } \
rad = fa * boxhalfsize[Y] + fb * boxhalfsize[Z];   \
if (min>rad || max<-rad) return 0;



/*======================== Y-tests ========================*/

#define AXISTEST_Y02(a, b, fa, fb)			   \
p0 = -a*v0[X] + b*v0[Z];		      	   \
p2 = -a*v2[X] + b*v2[Z];	       	       	   \
if (p0<p2) { min = p0; max = p2; }\
else { min = p2; max = p0; } \
rad = fa * boxhalfsize[X] + fb * boxhalfsize[Z];   \
if (min>rad || max<-rad) return 0;



#define AXISTEST_Y1(a, b, fa, fb)			   \
p0 = -a*v0[X] + b*v0[Z];		      	   \
p1 = -a*v1[X] + b*v1[Z];	     	       	   \
if (p0<p1) { min = p0; max = p1; } \
else { min = p1; max = p0; } \
rad = fa * boxhalfsize[X] + fb * boxhalfsize[Z];   \
if (min>rad || max<-rad) return 0;



/*======================== Z-tests ========================*/



#define AXISTEST_Z12(a, b, fa, fb)			   \
p1 = a*v1[X] - b*v1[Y];			           \
p2 = a*v2[X] - b*v2[Y];			       	   \
if (p2<p1) { min = p2; max = p1; } \
else { min = p1; max = p2; } \
rad = fa * boxhalfsize[X] + fb * boxhalfsize[Y];   \
if (min>rad || max<-rad) return 0;



#define AXISTEST_Z0(a, b, fa, fb)			   \
p0 = a*v0[X] - b*v0[Y];				   \
p1 = a*v1[X] - b*v1[Y];			           \
if (p0<p1) { min = p0; max = p1; } \
else { min = p1; max = p0; } \
rad = fa * boxhalfsize[X] + fb * boxhalfsize[Y];   \
if (min>rad || max<-rad) return 0;



int triBoxOverlap(double boxcenter[3], double boxhalfsize[3], double triverts[3][3])
{

	/*    use separating axis theorem to test overlap between triangle and box */
	/*    need to test for overlap in these directions: */
	/*    1) the {x,y,z}-directions (actually, since we use the AABB of the triangle */
	/*       we do not even need to test these) */
	/*    2) normal of the triangle */
	/*    3) crossproduct(edge from tri, {x,y,z}-directin) */
	/*       this gives 3x3=9 more tests */
	double v0[3], v1[3], v2[3];
	//   double axis[3];
	double min, max, p0, p1, p2, rad, fex, fey, fez;		// -NJMP- "d" local variable removed
	double normal[3], e0[3], e1[3], e2[3];

	/* This is the fastest branch on Sun */
	/* move everything so that the boxcenter is in (0,0,0) */

	SUB(v0, triverts[0], boxcenter);
	SUB(v1, triverts[1], boxcenter);
	SUB(v2, triverts[2], boxcenter);



	/* compute triangle edges */

	SUB(e0, v1, v0);      /* tri edge 0 */
	SUB(e1, v2, v1);      /* tri edge 1 */
	SUB(e2, v0, v2);      /* tri edge 2 */

						  /* Bullet 3:  */
						  /*  test the 9 tests first (this was faster) */
	fex = fabsf(e0[X]);
	fey = fabsf(e0[Y]);
	fez = fabsf(e0[Z]);
	AXISTEST_X01(e0[Z], e0[Y], fez, fey);
	AXISTEST_Y02(e0[Z], e0[X], fez, fex);
	AXISTEST_Z12(e0[Y], e0[X], fey, fex);

	fex = fabsf(e1[X]);
	fey = fabsf(e1[Y]);
	fez = fabsf(e1[Z]);
	AXISTEST_X01(e1[Z], e1[Y], fez, fey);
	AXISTEST_Y02(e1[Z], e1[X], fez, fex);
	AXISTEST_Z0(e1[Y], e1[X], fey, fex);

	fex = fabsf(e2[X]);
	fey = fabsf(e2[Y]);
	fez = fabsf(e2[Z]);
	AXISTEST_X2(e2[Z], e2[Y], fez, fey);
	AXISTEST_Y1(e2[Z], e2[X], fez, fex);
	AXISTEST_Z12(e2[Y], e2[X], fey, fex);


	/* Bullet 1: */
	/*  first test overlap in the {x,y,z}-directions */
	/*  find min, max of the triangle each direction, and test for overlap in */
	/*  that direction -- this is equivalent to testing a minimal AABB around */
	/*  the triangle against the AABB */

	/* test in X-direction */
	FINDMINMAX(v0[X], v1[X], v2[X], min, max);
	if (min>boxhalfsize[X] || max<-boxhalfsize[X]) return 0;

	/* test in Y-direction */
	FINDMINMAX(v0[Y], v1[Y], v2[Y], min, max);
	if (min>boxhalfsize[Y] || max<-boxhalfsize[Y]) return 0;

	/* test in Z-direction */
	FINDMINMAX(v0[Z], v1[Z], v2[Z], min, max);
	if (min>boxhalfsize[Z] || max<-boxhalfsize[Z]) return 0;

	/* Bullet 2: */
	/*  test if the box intersects the plane of the triangle */
	/*  compute plane equation of triangle: normal*x+d=0 */

	CROSS(normal, e0, e1);

	// -NJMP- (line removed here)

	if (!planeBoxOverlap(normal, v0, boxhalfsize)) return 0;	// -NJMP-

	return 1;   /* box and triangle overlaps */
}



int triBoxOverlap(const Vector<double> boxcenter, const Vector<double> boxhalfsize, const Vector<double> triverts[3])
{
	double bc[3] = { boxcenter[0],boxcenter[1],boxcenter[2] };
	double bhs[3] = { boxhalfsize[0], boxhalfsize[1],boxhalfsize[2] };
	double tv[3][3];
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			tv[i][j] = triverts[i][j];
	return triBoxOverlap(bc, bhs, tv);
}

#ifndef min
#define min(x,y)  x<y ? x : y
#endif

#ifndef max
#define max(x,y)  x>y ? x : y
#endif
#ifndef EPS_D 
#define EPS_D 1E-10 
#endif

#define OBBTEST(p,R) \
 if ((min3(p[0],p[1],p[2])>R) || (max3(p[0],p[1],p[2])<-R)) return false;

bool triangleAABBIntersection(const Box &B, const triangle &T)
{
	
	double R;
	Vector<double> E[3],N;
	Vector<double> D = T.P[0] - B.P;
	Vector<double> A[3] = { ex,ey,ez };
	
	int i,j;
	int l, m;

	E[0] = T.P[1] - T.P[0];
	E[1] = T.P[2] - T.P[0];
	E[2] = E[1] - E[0];
	N = E[0] % E[1];
	double a[3] = { B.d[0] / 2.0,B.d[1] / 2.0,B.d[2] / 2.0 };
	
	// Achse N Test
	R = a[0] * fabs(N*A[0]) + a[1] * fabs(N*A[1]) + a[2] * fabs(N*A[2]);
	if (abs(N*D) > R) return false;
//	cout << "% N Test" << endl;

	// test Ai % Ej
	double p;
	double d;
	for (int i = 0; i < 9; i++)
	{
		switch (i)
		{
		case 0: p = (A[0] % E[0])*D; d =  A[0] * N;  R = a[1] * abs(A[2] * E[0]) + a[2] * abs(A[1] * E[0]); break;
		case 1: p = (A[0] % E[1])*D; d = -A[0] * N;  R = a[1] * abs(A[2] * E[1]) + a[2] * abs(A[1] * E[1]); break;
		case 2: p = (A[0] % E[2])*D; d=  -A[0] * N;  R = a[1] * abs(A[2] * E[2]) + a[2] * abs(A[1] * E[2]); break;
		case 3: p = (A[1] % E[0])*D; d =  A[1] * N;  R = a[0] * abs(A[2] * E[0]) + a[2] * abs(A[0] * E[0]); break;
		case 4: p = (A[1] % E[1])*D; d = -A[1] * N;  R = a[0] * abs(A[2] * E[1]) + a[2] * abs(A[0] * E[1]); break;
		case 5: p = (A[1] % E[2])*D; d = -A[1] * N;  R = a[0] * abs(A[2] * E[2]) + a[2] * abs(A[0] * E[2]); break;
		case 6: p = (A[2] % E[0])*D; d =  A[2] * N;  R = a[0] * abs(A[1] * E[0]) + a[1] * abs(A[0] * E[0]); break;
		case 7: p = (A[2] % E[1])*D; d = -A[2] * N;  R = a[0] * abs(A[1] * E[1]) + a[1] * abs(A[0] * E[1]); break;
		case 8: p = (A[2] % E[2])*D; d = -A[2] * N;  R = a[0] * abs(A[1] * E[2]) + a[1] * abs(A[0] * E[2]); break;
		}

		if (p > R)
		{
			if (d >= 0) return false;
			if (p + d > R) return false;
		}
		else if (p < -R)
		{
			if (d <= 0) return false;
			if (p + d < -R) return false;
		}
	}

	// test Ak
	double d0, d1;
	for (int i = 0; i < 3; i++)
	{
		p = A[i] * D;
		d0 = A[i] * E[0];
		d1 = A[i] * E[1];
		R = a[i];
		if (p > R)
		{
			if (d0 >= 0)
			{
				if (d1 >= 0) return false;
				if (p + d1 > R) return false;
			}
			else if (d1 <= d0)
			{
				if (p + d1 > R) return false;
			}
			else
			{
				if (p + d0 > R) return false;
			}
		}
		else if (p < -R)
		{
			if (d0 <= 0)
			{
				if (d1 <= 0) return false;
				if (p + d1 < -R) return false;
			}
			else if (d1 >= d0)
			{
				if (p + d1 < -R) return false;
			}
			else
			{
				if (p + d0 < -R) return false;
			}
		}
	}
	return true;
}


