/*****************************************************************//**
 * \file   vector.h
 * \brief  This file contains the Vector template class and some useful functions around this class
 * 
 * 
 * \author Thomas Weigel
 * \date   July 2021
 *********************************************************************/

#ifndef VEKTOR_H
#define VEKTOR_H
#include <ostream>
#include <stdlib.h>
#include <stdio.h>
#include <complex>
#include <fstream>
#ifdef _WIN32
#define _USE_MATH_DEFINES
#endif
#include<math.h> 
#include<iostream>

double arctan (double y, double x);
std::complex<double>  tan(std::complex<double>  z);

inline std::complex<double> operator + (std::complex<double> x, int y) { return x+(double)y; }
inline std::complex<double> operator + (int y, std::complex<double> x) { return (double)y+x; }
inline std::complex<double> operator - (std::complex<double> x, int y) { return x-(double)y; }
inline std::complex<double> operator - (int y, std::complex<double> x) { return (double)y-x; }

inline std::complex<double> operator * (std::complex<double> x, int y) { return x*(double)y; }
inline std::complex<double> operator * (int y, std::complex<double> x) { return (double)y*x; }
inline std::complex<double> operator / (std::complex<double> x, int y) { return x/(double)y; }
inline std::complex<double> operator / (int y, std::complex<double> x) { return (double)y/x; }

/** @brief Template class for threedimensional vectors
* 
* This class provides standard operations and functions for threedimensional vectors. It is intended for use with int, double and complex<double> based vectors. 
* Also, operators with mixed types are given (e.g. operators for double and complex vectors).
*/
template <class T>
class Vector          
{
 public :

 
 void binWrite (std::ofstream &os) ///< writes vector to a binary file, represented by the output stream os
 {
  for (int i=0; i<3; i++)
  os.write ((char *) &data[i], (char) sizeof (data[i]));
 }


 void binRead (std::ifstream &is) ///< reads vector from binary file, represented by the input stream is
 {
  for (int i=0; i<3; i++)
  {
   is.read ((char *) &data[i], (char) sizeof (data[i]));   
  }  
 }
 
 /** @brief  Standard constructor (all components are initialized with 0)*/
 Vector () 
  { 
   data[0]=0;
   data[1]=0;
   data[2]=0;
  }

 /**
  * @brief Constructor with the three components as parameters
  * 
  * \param x first component
  * \param y second component
  * \param z third component
  */
 Vector (T x, T y, T z) 
  { 
   data[0]=x; 
   data[1]=y;
   data[2]=z;
  }
 
/** Copy constructor */ 
inline Vector (const Vector& r)
 { 
  for (int i=0; i<3; i++)
  data[i]=r.data[i];
 }
  
/**
 * @brief Addition of two vectors.
 * 
 * \param r1 first summand
 * \param r2 second summand
 * \return sum of r1 and r2 
 */
inline  friend Vector operator + (const Vector& r1, const Vector& r2)
 {
 Vector<T> h; 
  h.data[0]=r1.data[0] + r2.data[0];
  h.data[1]=r1.data[1] + r2.data[1];
  h.data[2]=r1.data[2] + r2.data[2];
 return h;
 }

/** sign minus */
inline Vector operator - ()
 {
  Vector Erg;
  for (int i=0; i<3; i++)
  Erg.data[i]=-data[i];
  return Erg;
 } 
 
/**
 * @brief Subtraction of two vectors.
 * 
 * \param r1 minuend
 * \param r2 subtrahend
 * \return 
 */
inline friend Vector operator - (const Vector& r1, const Vector& r2)
 {
  Vector<T> h;
  h.data[0]=r1.data[0] - r2.data[0];
  h.data[1]=r1.data[1] - r2.data[1];
  h.data[2]=r1.data[2] - r2.data[2];
  return h; 
 }
 

/**
 * @brief cross product between two vectors.
 * 
 * \param r1 first vector
 * \param r2 second vector
 * \return cross product
 */
inline friend Vector operator % (const Vector &r1,const Vector& r2)
 {
 Vector<T> h;
 h.data[0]=r1.data[1]*r2.data[2] - r1.data[2]*r2.data[1];
 h.data[1]=r1.data[2]*r2.data[0] - r1.data[0]*r2.data[2];
 h.data[2]=r1.data[0]*r2.data[1] - r1.data[1]*r2.data[0];
 return h;
 }

/**
 * @brief calculates the square root of a vector: result: sqrt((x,y,z)) \f$ \rightarrow (\sqrt{x},\sqrt{y},\sqrt{z}) \f$
 * 
 */
inline friend Vector sqrt (const Vector& r)
 {
  return Vector<T> (sqrt(r.data[0]),sqrt(r.data[1]),sqrt(r.data[2]));
 }

/**
 * @brief exponential function of a vector (x,y,z), defined by exp((x,y,z)) \f$\rightarrow (e^x,e^y,e^z)\f$
 */
inline  friend Vector exp(const Vector& r)
 {
  return Vector<T> (exp(r.data[0]),exp(r.data[1]),exp(r.data[2]));
 }

/** @brief Component-wise multiplication of two vectors */ 
inline friend Vector emult (const Vector& r1, const Vector& r2)
 {
  Vector<T> Erg;
  for (int i=0; i<3; i++)  
    Erg.data[i]=r1.data[i]*r2.data[i];
  return Erg;  
 } 

/** @brief Component-wise division of two vectors */
inline friend Vector ediv (const Vector& r1, const Vector& r2)
 {
  Vector<T> Erg;
  for (int i=0; i<3; i++)  
    Erg.data[i]=r1.data[i]/r2.data[i];
  return Erg;  
 } 

/**
 * @brief Component-wise division, scalar (s) by vector (r) 
 * 
 * \param s Scalar
 * \param r Vector
 * \return Vector with components: (s/r[0],s/r[1],s/r[2])
 */
inline friend Vector ediv (T s, const Vector& r)
{
	Vector<T> Erg;
	for (int i = 0; i<3; i++)
		Erg.data[i] = s / r.data[i];
	return Erg;
}

 
/** @brief Assignment operator */ 
 Vector& operator = (const Vector& r)
 {
	 // if (&r == NULL) return *this;
  if (this == &r) return *this;
  for (int i=0; i<3; i++)
  data[i]=r.data[i];
  return *this;
 }

/** @brief Assignment operator (with an ordinary threedimensional array */ 
 Vector& operator = (T r[3])
 {
  for (int i=0; i<3; i++)
  data[i]=r[i];
  return *this;
 }

/** @brief Comparison operators: two vectors are equal if their components are (exactly) the same */
  bool operator == (const Vector& r)
 {
  bool Erg;
  Erg=((data[0]==r.data[0]) && (data[1]==r.data[1]) && (data[2]==r.data[2]));
  return Erg;
 } 
 
/** @brief Inequality operator: Two vectors are inequal if at least one component is not equal */
 bool operator != (const Vector& r)
 {
  bool Erg;
  Erg=((data[0]!=r.data[0]) || (data[1]!=r.data[1]) || (data[2]!=r.data[2]));
  return Erg;
 }
 
/** @brief Addition operator*/
 Vector& operator += (const Vector& r)
 {
  for (int i=0; i<3; i++)
  data[i]=data[i]+r.data[i];
  return *this;
 }
 
/** @brief Subtraction operator*/
 Vector& operator -= (const Vector& r)
 {
  for (int i=0; i<3; i++)
  data[i]=data[i]-r.data[i];
  return *this;
 }
 
/** @brief Multiplication with scalar */
 Vector& operator *= (T r)
 {
  for (int i=0; i<3; i++)
  data[i]=data[i]*r;
  return *this;
 }
 
/** @brief Division by scalar */
 Vector& operator /= (T r) 
 {
   for (int i=0; i<3; i++)
   data[i]=data[i]/r;
   return *this;
 }
 
inline friend bool operator == (const Vector& A, const Vector& B) ///< Comparison operator. Two vectors are equal if all components are the same
 {
  return ((A[0]==B[0])&&(A[1]==B[1])&&(A[2]==B[2])); 
 }
 
inline friend bool operator != (const Vector& A, const Vector& B) ///< Inequality operator. Two vectors are unequal if they differ in at least one component.
 {
 return !((A[0]==B[0])&&(A[1]==B[1])&&(A[2]==B[2])); 
 }



inline friend T operator * (const Vector& a, const Vector& b) ///< Scalar product between two vectors
 {
  T h;
  h=0;
  for (int i=0; i<3; i++)
   h+=a[i]*b[i];
  return h;
 }

/**
 * @brief Multiplication of a scalar with a vector
 * 
 * \param x scalar
 * \param r vector
 * 
 */
 friend Vector operator * (T x, const Vector& r)
 {
  Vector<T> h;
  for (int i=0; i<3; i++)
  h.data[i]=r.data[i]*x;
  return h; 
 }

/**
 * @brief Multiplication of a vector with a scalar.
 * 
 * \param r vector
 * \param x scalar
 * 
 */
 friend Vector operator * ( const Vector& r, T x)
 {
  Vector<T> h;
  for (int i=0; i<3; i++)
  h.data[i]=r.data[i]*x;
  return h; 
 }
 

/**
 * @brief Division of a vector by a scalar.
 * 
 * \param r vector 
 * \param x scalar
 * 
 */
 friend Vector operator / (const Vector& r,T x)
 {
  Vector<T> Erg;
  for (int i=0; i<3; i++)
   Erg.data[i]=r.data[i] / x;
  return Erg; 
 }
 

/** @brief Bracket operator, gives access to the i-th component of the vector (as for a usual array).
* 
 * \param i Index of the component
 *  \return content of the i-th component
 */
 T& operator[] (int i) {return data[i]; };

 /** @brief Bracket operator, gives access to the i-th component of the vector (as for a usual array).
*
 * \param i Index of the component
 *  \return content of the i-th component
 */
 const T& operator[ ] (int i) const  { return data[i]; };
 T data[3]; 
 //int Dim;
};


/**
 * @brief Transformation from cartesian coordinates to spherical coordinates
 * 
 * \param x Vector in cartesian coordinates
 * \return Vector in spherical coordinates \f$(r,\vartheta,\varphi)\f$
 */
template <class T> Vector<T> cart2sphere(const Vector<T>& x)
{
 Vector<T> y,a;

  a[0] = x[0];
  a[1] = x[1];
  a[2] = x[2];
  y[1] = sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
  if (a[2] == 0.0)
    if (y[0] == 0.0)
      y[1] = 0;
    else
      y[1] = M_PI / 2.0;
    else
      {
        y[1] = atan(sqrt(a[0] * a[0] + a[1] * a[1]) / a[2]);
        if (a[2] < 0.0)
          y[1] = M_PI + y[1];
      }
  if (a[0] == 0.0)
    if (a[1] == 0.0)
      y[2] = 0.0;
    else
      if (a[1] > 0.0)
        y[2] = M_PI / 2.0;
      else
        y[2] = 3.0 * M_PI / 2.0;
  else
    {
      y[2] = atan(a[1] / a[0]);
      if (a[0] < 0.0)
        y[2] = M_PI + y[2];
    }
  return y;
}



/**
 * @brief Transforms vector from spherical coordinates into cartesian coordinates.
 * 
 * \param P Vector in spherical coordinates \f$(r,\vartheta,\varphi)\f$
 * \return Vector in cartesian coordinates
 */
template <class T> Vector<T> sphere2cart (const Vector<T>& P)
{


 Vector<T> Erg;
 Erg[0]=P[0]*sin(P[1])*cos(P[2]);
 Erg[1]=P[0]*sin(P[1])*sin(P[2]);
 Erg[2]=P[0]*cos(P[1]); 
 return Erg; 
}


std::istream& operator >> (std::istream &is, Vector<std::complex<double> > &r);
std::istream& operator >> (std::istream &is, Vector<double> &r);
std::istream& operator >> (std::istream &is, Vector<int> &r);

std::ostream& operator << (std::ostream &is, const Vector<std::complex<double> > &r);
std::ostream& operator << (std::ostream &is, const Vector<double> &r);
std::ostream& operator << (std::ostream &is, const Vector<int> &r);

void errmsg(char *S); //! Fehlermeldung ausgeben
double abs (const Vector<double> &r); ///< Calculates the absolute value of a double vector
double abs (const Vector<std::complex<double> > &r); ///< Calculates the absolute value of a complex vector \$f(abs(\vec{r})=\sqrt{\vec{r} \cdot \vec{r}*\$f})
double abs (const Vector<int> &r); ///< Calculates the absolute value of an integer vector 
double abs2 (const Vector<double> &r); ///< Calculates the squared absolute value of double vector
double abs2 (const Vector<std::complex<double> > &r);///< Calculates the squared absolute value of double vector \$f(abs2(\vec{r})=\vec{r} \cdot \vec{r}*\$f)
double  abs2 (std::complex<double>  x); ///< Absolute value of the complex number x \$f(abs(x)=x \cdot x*\$f)
double abs2 (double x); ///< Squared absolute value of x

/** @name Standard operators for vectors
 * @brief Standard operators for vectors with mixed types 
 */
 ///@{
Vector<double> operator - (const Vector<double>&r1, const Vector<int> &r2);
Vector<double> operator - (const Vector<int>&r1, const Vector<double> &r2);
Vector<std::complex<double> > operator - (const Vector<double> &r1,const Vector<std::complex<double> > &r2);
Vector<std::complex<double> > operator - (const Vector<std::complex<double> > &r1,const Vector<double> &r2);
Vector<std::complex<double> > operator - (const Vector<std::complex<double> >& r1, const Vector<int> &r2);
Vector<std::complex<double> > operator - (const Vector<int>& r1, const Vector<std::complex<double> > &r2);

Vector<std::complex<double> > operator + (const Vector<double> &r1,const Vector<std::complex<double> > &r2);
Vector<std::complex<double> > operator + (const Vector<std::complex<double> > &r1,const Vector<double> &r2);
Vector<double> operator + (const Vector<int>& r1, const Vector<double> &r2);
Vector<double> operator + (const Vector<double>& r1, const Vector<int> &r2);
Vector<std::complex<double> > operator + (const Vector<std::complex<double> >& r1, const Vector<int> &r2);
Vector<std::complex<double> > operator + (const Vector<int>& r1, const Vector<std::complex<double> > &r2);
Vector<double> operator * (int x, const Vector<double>& r);
Vector<double> operator * (const Vector<double>& r, int x);
Vector<double> operator * (double x, const Vector<int>& r);
Vector<double> operator * (const Vector<int>& r, double x);
Vector<std::complex<double> > operator * (std::complex<double>  x, const Vector<double>& r);
Vector<std::complex<double> > operator * (const Vector<double>& r, std::complex<double>  x);
Vector<std::complex<double> > operator * (std::complex<double>  x, const Vector<int>& r);
Vector<std::complex<double> > operator * (const Vector<int>& r, std::complex<double>  x);

Vector<std::complex<double> > operator / (const Vector<double>& r, std::complex<double>  x);
Vector<double> operator / (const Vector<double>& r, int x);
Vector<std::complex<double> > operator / (const Vector<std::complex<double> >& r, double x);
Vector<std::complex<double> > operator / (const Vector<std::complex<double> >& r, int x);
Vector<double> operator / (const Vector<int>& r, double x);
Vector<std::complex<double> > operator / (const Vector<int>& r, std::complex<double>  x);


std::complex<double>  operator * (const Vector<double>& a, const Vector<std::complex<double> >& b);
std::complex<double>  operator * (const Vector<std::complex<double> >& a, const Vector<double>& b);
double operator * (const Vector<double>& a, const Vector<int>& b);
double operator * (const Vector<int>& a, const Vector<double>& b);
std::complex<double>  operator * (const Vector<std::complex<double> >& a, const Vector<int>& b);
std::complex<double>  operator * (const Vector<int>& a, const Vector<std::complex<double> >& b);

Vector<std::complex<double> > operator % (const Vector<std::complex<double> >& a,
    const Vector<double>& b);
Vector<std::complex<double> > operator % (const Vector<double>& a,
    const Vector<std::complex<double> >& b);
Vector<std::complex<double> > operator % (const Vector<std::complex<double> >& a, const Vector<int>& b);
Vector<std::complex<double> > operator % (const Vector<int>& a, const Vector<std::complex<double> >& b);
Vector<std::complex<double> > operator % (const Vector<double>& a, const Vector<int>& b);
Vector<std::complex<double> > operator % (const Vector<int>& a, const Vector<double>& b);

///@}


Vector<std::complex<double> > conj (const Vector<std::complex<double> > &r); ///< Conjugate complex of the vector r 

/**
 * @brief This function returns an array with all components of the Vector-array r.
 * @param numV number of elements in vector list r
 * @param r List of vectors to convert
 * @return Array with all components of r: erg[0]=real(r[0]), erg[1]=imag(r[0]),erg[2]=real(r[1]), erg[2]=imag(r[1]),...
 */
double *conv2double (int numV, Vector<std::complex<double> > *r);

/**
 * @name element-wise multiplication
 * @brief Element-wise multiplication of two vectors with different types
 */
///@{
Vector<double> emult (const Vector<double> &r1, const Vector<int> &r2);
Vector<double> emult (const Vector<int> &r1, const Vector<double> &r2);
///@}
Vector<double> floor (const Vector<double> &r); ///< rounding up, Component-by-component 
Vector<int> ifloor (const Vector<double> &r);///< rounding up, component-by-component rounding
inline Vector <double> ceil (const Vector<double>& r) ///< rounding down, component - by - component rounding
{
 return Vector<double> ((int) ceil (r[0]),(int) ceil (r[1]), (int) ceil (r[2]));
}

inline Vector <int> iceil (const Vector<double>& r) ///< rounding down, component - by - component rounding
{
 return Vector<int> ((int) ceil (r[0]),(int) ceil (r[1]), (int) ceil (r[2]));
}

inline double sign(double x) ///< returns -1 if x<0, 0 if x equals to 0 and 1 if x>0
{
 double y=(x>0) ? 1.0 : 0.0 + (x<0) ? -1.0 : 0.0;
 return y;
}

/**
 * @brief This function converts a list of doubles into a list of complex vectors.
 * @param num Number of complex vectors to convert
 * @param r List of double values
 * @return List of complex vectors (erg): erg[0][0]=r[0]+i*r[1],erg[0][1]=r[2]+i*r[3],erg[0][2]=r[4]+i*r[5], erg[1][0]=r[6]+i*r[7],...
 */
Vector<std::complex<double> > *conv2vector (int num, double *r);
Vector<std::complex<double> > grad2d (double (*f(Vector<std::complex<double> >)),Vector<std::complex<double> > x, double dx);
Vector<std::complex<double> > convd2c (const Vector<double> &r); 
Vector<double> real(Vector <std::complex<double> > r); ///< returns a vector with the real parts of the components of the complex vector r
Vector<double> imag(Vector <std::complex<double> > r); ///< returns a vector with the imaginary parts of the components of the complex vector r
double sDreieck (Vector<double>,Vector<double>,Vector<double>);
double sViereck (Vector<double>,Vector<double>,Vector<double>,Vector<double>); 
void rotate (Vector<double> &r, double phi);
Vector<double> rotate (const Vector<double> &r, double dtheta, double dphi);
Vector<std::complex<double> > vdc (const Vector<double> &r);
Vector<std::complex<double> > makeReal (const Vector<std::complex<double> > &r);
/**
 * @brief Converts a vector, described by its cartesian components, into a vector in spherical coordinates (\f$(r, \vartheta, \varphi) \f$)
 * @param x x-component
 * @param y y-component
 * @param z z-component
 */
Vector<double> cart2sphere(double x, double y, double z);
/**
 * @brief Converts a vector, described by its spherical components (\f$(r, \vartheta, \varphi) \f$), into a vector in cartesian coordinates 
 * @param r r-component
 * @param theta \f$(\vartheta\f$)-component
 * @param phi \f$(\varphi\f$)-component
 */
Vector<double> sphere2cart (double r, double theta, double phi);

/**
 * @brief Calculates the three direction vectors of a coordinate system formed by the surface normal n and the propagation vector k 
 * @param n Surface normal
 * @param k propagation vector
 * @param e0 normalized direction vector for the first coordinate axis, equal to n 
 * @param e1 normalized direction vector for the first coordinate axis, \f$\vec{e_1}=\vec{e_2}\times\vec{e_0}\f$
 * @param e2 normalized direction vector for the third coordinate axis: \f$\vec{e_2}=\vec{n}\times\vec{k}\f$. If \f$\vec{e_2}=0\f$, \f$\vec{e_2}= \frac{\vec{n}\times\vec{e_x}}{\left|\vec{n}\times\vec{e_x}\right|\f$
 */ 
void getKSystem (const Vector<double> &n, const Vector<double> &k,
                Vector<double> &e1, Vector<double> &e2, Vector<double> &e3);

Vector<double> arg (Vector<std::complex<double> > &r);
std::complex<double>  asin (std::complex<double>  z);
std::complex<double>  tan  (std::complex<double>  z);
std::complex<double>  ihoch (int l);
void sphereunitv (Vector<double> &P, Vector<double> &er, Vector<double> &etheta, Vector<double> &ephi); 
Vector <double> emax (Vector<double> &P1, Vector<double> &P2);

#define ex Vector<double> (1.0,0.0,0.0)
#define ey Vector<double> (0.0,1.0,0.0)
#define ez Vector<double> (0.0,0.0,1.0)
#define dzero Vector<double> (0.0,0.0,0.0)
#define one Vector<double>  (1.0,1.0,1.0)
#define czero Vector<std::complex<double> > (0.0,0.0,0.0)
#define cone Vector<std::complex<double> >  (1.0,1.0,1.0)


Vector<double> er(double theta, double phi);
Vector<double> etheta(double theta, double phi);
Vector<double> ephi(double theta, double phi);

Vector<std::complex<double> > operator * (std::complex<double>  x,const Vector<double>& r);
Vector<std::complex<double> > operator * (const Vector<double>& r, std::complex<double>  x);
double *conv2double (int Anz, Vector<std::complex<double> > *r);
Vector<std::complex<double> > *conv2vector (int Anz, double *r);
Vector<std::complex<double> > grad2d (double (*f(Vector<std::complex<double> >)),Vector<std::complex<double> > x, double dx);
Vector<std::complex<double> > convd2c (const Vector<std::complex<double> > &r); 
void rotate (Vector<double> &r, double phi);
std::complex<double>  asin (std::complex<double>  z);
std::complex<double>  tan  (std::complex<double>  z);
Vector<std::complex<double> > conj (const Vector<std::complex<double> > &r);
Vector<std::complex<double> > norm (const Vector<std::complex<double> > &r);  // Normiert Vektor r
Vector<double> norm (const Vector<double> &r);  // Normiert Vektor r (r=r/|r|)




/**
 * @brief Converts a double vector into a string
 * @param S string in which the components of the vector P will be written. The components will be written in one row, separated by spaces. 
 * @param P Vector 
 */
char *toString (char *S,Vector<double> P);
/**
 * @brief Converts a complex vector into a string
 * @param S string in which the components of the vector P will be written. The components (first real part, then imaginary part) will be written in one row, separated by spaces.
 * @param P Vector
 */
char* toString (char *S, Vector<std::complex<double> > P);
 Vector<double> nanV (const char *tagp);   

 bool isnan(Vector<std::complex<double> > v); ///< returns true if one of the components of v is NaN
 bool isnan (Vector<double> v); ///< returns true if one of the components of v is NaN

 #endif

#ifndef complex_I
#define complex_I
const std::complex<double>  I=std::complex<double>  (0.0,1.0); ///< imaginary unit 
#endif

