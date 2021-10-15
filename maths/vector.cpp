#ifndef VEKTOR_CC
#define VEKTOR_CC
#ifdef _WIN32

#ifndef min
#define min(a,b) a<b ? a : b
#endif

#ifndef max
#define max(a,b) a>b ? a : b
#endif

#endif

#include <stdlib.h>
#include <iostream>
#include "vector.h"
#include <complex>
#ifdef _WIN32
#define _USE_MATH_DEFINES
#endif
#include <math.h>

using namespace std;

double abs (const Vector<double> &r)
{
 double Erg=sqrt (r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
 return Erg;
}

 double abs (const Vector<int> &r)
{
 double Erg=sqrt ((double)r[0]*r[0]+(double)r[1]*r[1]+(double)r[2]*r[2]);
 return Erg;
}

 double abs (const Vector<complex<double> > &r)
{
 double re,im;
 double Erg=0.0;
 for (int i=0; i<3; i++)
 {
  re=real(r[i]);
  im=imag(r[i]);
  Erg+=re*re+im*im;
 }
 return sqrt(Erg);
}

 double abs2(double x)
{
 return x*x;
}

 double  abs2 (complex<double>  x)
{
 return real(x*conj(x)); 
}

 double abs2 (const Vector<double> &r)
{
 return r[0]*r[0]+r[1]*r[1]+r[2]*r[2];
}

 double abs2 (const Vector<complex<double> > &r)
{
 return real(r*conj(r));
}

Vector<double> operator - (const Vector<double>&r1, const Vector<int> &r2)
{
 return Vector<double> (r1[0]-r2[0],r1[1]-r2[1],r1[2]-r2[2]);
}

Vector<double> operator - (const Vector<int>&r1, const Vector<double> &r2)
{
 return Vector<double> (r1[0]-r2[0],r1[1]-r2[1],r1[2]-r2[2]);
}

Vector<complex<double> > operator - (const Vector<complex<double> >& r1, const Vector<double> &r2)
{
 return Vector<complex<double> > (r1[0]-r2[0],r1[1]-r2[1],r1[2]-r2[2]);

}

Vector<complex<double> > operator - (const Vector<double>& r1, const Vector<complex<double> > &r2)
{
 return Vector<complex<double> > (r1[0]-r2[0],r1[1]-r2[1],r1[2]-r2[2]);

}

Vector<complex<double> > operator - (const Vector<complex<double> >& r1, const Vector<int> &r2)
{
 return Vector<complex<double> > (r1[0]-r2[0],r1[1]-r2[1],r1[2]-r2[2]);
}

Vector<complex<double> > operator - (const Vector<int>& r1, const Vector<complex<double> > &r2)
{
 return Vector<complex<double> > (r1[0]-r2[0],r1[1]-r2[1],r1[2]-r2[2]);
}

Vector<complex<double> > operator + (const Vector<complex<double> >& r1, const Vector<double> &r2)
{
 return Vector<complex<double> > (r1[0]+r2[0],r1[1]+r2[1],r1[2]+r2[2]);
}

Vector<complex<double> > operator + (const Vector<double>& r1, const Vector<complex<double> > &r2)
{
 return Vector<complex<double> > (r1[0]+r2[0],r1[1]+r2[1],r1[2]+r2[2]);
}

Vector<double> operator + (const Vector<int>& r1, const Vector<double> &r2)
{
 return Vector<double> (r1[0]+r2[0],r1[1]+r2[1],r1[2]+r2[2]);
}

Vector<double> operator + (const Vector<double>& r1, const Vector<int> &r2)
{
 return Vector<double> (r1[0]+r2[0],r1[1]+r2[1],r1[2]+r2[2]);
}

Vector<complex<double> > operator + (const Vector<complex<double> >& r1, const Vector<int> &r2)
{
 return Vector<complex<double> > (r1[0]+r2[0],r1[1]+r2[1],r1[2]+r2[2]);
}

Vector<complex<double> > operator + (const Vector<int>& r1, const Vector<complex<double> > &r2)
{
 return Vector<complex<double> > (r1[0]+r2[0],r1[1]+r2[1],r1[2]+r2[2]);
}

complex<double>  operator * (const Vector<double>& a, const Vector<complex<double> > &b)
{
 return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

complex<double>  operator * (const Vector<complex<double> >& a, const Vector<double> &b)
{
 return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

double operator * (const Vector<double>& a, const Vector<int> &b)
{
 return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

double operator * (const Vector<int>& a, const Vector<double> &b)
{
 return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

complex<double>  operator * (const Vector<complex<double> >& a, const Vector<int> &b)
{
 return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

complex<double>  operator * (const Vector<int>& a, const Vector<complex<double> > &b)
{
 return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

Vector<double> operator * (int x, const Vector<double> &r)
{
 Vector<double> h;
 for (int i=0; i<3; i++)
  h.data[i]=r.data[i]*x;
 return h; 
}

Vector<double> operator * (const Vector<double> &r,int x)
{
 Vector<double> h;
 for (int i=0; i<3; i++)
  h.data[i]=r.data[i]*x;
 return h; 
}

Vector<complex<double> > operator * (double x,const Vector<complex<double> >& r)
{
 Vector<complex<double> > h;
 for (int i=0; i<3; i++)
  h.data[i]=r.data[i]*x;
 return h; 
}

Vector<complex<double> > operator * (const Vector<complex<double> >& r,double x)
{
 Vector<complex<double> > h;
 for (int i=0; i<3; i++)
  h.data[i]=r.data[i]*x;
 return h; 
}

Vector<complex<double> > operator * (int x,const Vector<complex<double> >& r)
{
 Vector<complex<double> > h;
 for (int i=0; i<3; i++)
  h.data[i]=r.data[i]*x;
 return h; 
}

Vector<complex<double> > operator * (const Vector<complex<double> >& r,int x)
{
 Vector<complex<double> > h;
 for (int i=0; i<3; i++)
  h.data[i]=r.data[i]*x;
 return h; 
}

Vector<complex<double> > operator * (complex<double>  x,const Vector<double>& r)
{
 Vector<complex<double> > h;
 for (int i=0; i<3; i++)
  h.data[i]=r.data[i]*x;
 return h; 
}

Vector<complex<double> > operator * (const Vector<double>& r, complex<double>  x)
{
 Vector<complex<double> > h;
 for (int i=0; i<3; i++)
  h.data[i]=r.data[i]*x;
 return h; 
}

Vector<double> operator * (double x,const Vector<int>& r)
{
 Vector<double> h;
 for (int i=0; i<3; i++)
  h.data[i]=r.data[i]*x;
 return h; 
}

Vector<double> operator * (const Vector<int>& r, double x)
{
 Vector<double> h;
 for (int i=0; i<3; i++)
  h.data[i]=r.data[i]*x;
 return h; 
}

Vector<complex<double> > operator * (complex<double>  x,const Vector<int>& r)
{
 Vector<complex<double> > h;
 for (int i=0; i<3; i++)
  h.data[i]=r.data[i]*x;
 return h; 
}

Vector<complex<double> > operator * (const Vector<int>& r, complex<double>  x)
{
 Vector<complex<double> > h;
 for (int i=0; i<3; i++)
  h.data[i]=r.data[i]*x;
 return h; 
}

Vector<complex<double> > operator / (const Vector<double>& r, complex<double>  x)
{
 return Vector<complex<double> > (r.data[0]/x,r.data[1]/x,r.data[2]/x); 
}

Vector<double> operator / (const Vector<double>& r, int x)
{
 return Vector<double> (r.data[0]/(double)x,r.data[1]/(double)x,r.data[2]/(double)x); 
}

Vector<complex<double> > operator / (const Vector<complex<double> >& r, double x)
{
 return Vector<complex<double> > (r.data[0]/x,r.data[1]/x,r.data[2]/x); 
}

Vector<complex<double> > operator / (const Vector<complex<double> >& r, int x)
{
 return Vector<complex<double> > (r.data[0]/(double)x,r.data[1]/(double)x,r.data[2]/(double)x); 
}

Vector<double> operator / (const Vector<int>& r, double x)
{
 return Vector<double> (r.data[0]/x,r.data[1]/x,r.data[2]/x); 
}

Vector<complex<double> > operator / (const Vector<int>& r, complex<double>  x)
{
 return Vector<complex<double> > (r.data[0]/x,r.data[1]/x,r.data[2]/x); 
}

double *conv2double (int numV, Vector<complex<double> > *r)
{
 double *Erg;
 int j=1;
 
 Erg=new double[numV*2*3+1];
 for (int i=0; i<numV; i++)
 {
  Erg[j]=real(r[i][0]);
  Erg[j+1]=imag(r[i][0]);
  Erg[j+2]=real(r[i][1]);
  Erg[j+3]=imag(r[i][1]);
  Erg[j+4]=real(r[i][2]);
  Erg[j+5]=imag(r[i][2]);
  j+=6; 
 }
 return Erg;
} 

Vector<complex<double> > *conv2vector (int num, double *r)
{
 Vector<complex<double> > *Erg;
 int j=1;
 Erg=new Vector<complex<double> >[num];
 for (int i=0; i<num/3; i++)
 {
  Erg[i][1]=complex<double> (r[j+2],r[j+3]);
  Erg[i][2]=complex<double> (r[j+4],r[j+5]);
  j+=6;
 } 
 return Erg;
}


 Vector<complex<double> > grad2d (double (*f(Vector<complex<double> >)),Vector<complex<double> > x, double dx)  
{
 // Berechnet der Gradienten der Funktion f(x) in 2D 
 // dx : Genauigkeit (-> grad(f(x))=f(x+ex*dx)/dx+f(x+ey*dx)/dx...
 Vector<complex<double> > Erg,r2;
 double f1,f2;
 f1=*f(x);
 for (int i=0; i<2; i++)
 {
  r2=x;
  r2[i]=r2[i]+dx;
  f2=*f(r2);
  Erg[i]=(f2-f1)/dx;
 }
 return Erg; 
}

Vector<complex<double> > convd2c (const Vector<double> &r)
{
 Vector<complex<double> > Erg; 
 for (int i=0; i<3; i++)
  Erg[i]=r[i]; 
 return Erg;
}

Vector<double> real(Vector <complex<double> > r)
{
 Vector<double> Erg;
 for (int i=0; i<3; i++)
  Erg[i]=real(r[i]);
 return Erg;
}

Vector<double> imag(Vector <complex<double> > r)
{
 Vector<double> Erg;
 for (int i=0; i<3; i++)
  Erg[i]=imag(r[i]);
 return Erg;
}

void rotate (Vector<double> &r, double phi)
{
 // Drehung um die z-Achse
 // dreht den Vektor r um den Winkel phi 
 // phi < 0 -> Drehung nach links
 // phi > 0 -> Drehung nach rechts

 Vector<double> Erg;
 Erg[0]=  cos(phi) * r[0] + sin(phi) * r[1];
 Erg[1]= -sin(phi) * r[0] + cos(phi) * r[1];
 Erg[2]=r[2];
 r=Erg;
}

Vector<double> rotate (const Vector<double> &r, double dtheta, double dphi)
{
 Vector<double> Erg;
 Erg=cart2sphere(r);
 Erg[1]+=dtheta;
 Erg[2]+=dphi;
 return Erg;
}

complex<double>  asin (complex<double>  z)
{
 return (-I*log(I*z+sqrt(1.0-z*z)));
}

complex<double>  tan  (complex<double>  z)
{
 return sin(z)/cos(z);
}

Vector<complex<double> > conj (const Vector<complex<double> > &r)
 {
  Vector<complex<double> > Erg;
  Erg[0]=conj(r[0]);
  Erg[1]=conj(r[1]);
  Erg[2]=conj(r[2]);
 return Erg;  
 }

complex<double>  ihoch (int l)
{
 if (l>=0)
 {
  switch (l%4) 
  {
    case 0 : return 1.0; 
    case 1 : return I;
    case 2 : return -1.0;
    case 3 : return -I;   
  }
 }
 else
 {
  switch ((-l)%4)
  {
    case 0 : return 1.0;
    case 1 : return -I;
    case 2 : return -1.0;
    case 3 : return I;
  }
 }
}



istream& operator >> (istream &is, Vector<double> &r)
{
 is >> r[0] >> r[1] >> r[2];
 return is;
}

istream& operator >> (istream &is, Vector<int> &r)
{
 is >> r[0] >> r[1] >> r[2];
 return is;
}

istream& operator >> (istream &is, Vector<complex<double> > &r)
{
 double re,im;
 for (int i=0; i<3; i++)
 {
  is >> re >> im;
  r[i]=complex<double> (re,im);
 }
 return is; 
}

ostream& operator << (ostream& os, const Vector<double> &r)
{
 os << r[0] << "  " << r[1] << "  " << r[2];
 return os;
}

ostream& operator << (ostream& os, const Vector<int> &r)
{
 os << r[0] << "  " << r[1] << "  " << r[2];
 return os;
}

ostream& operator << (ostream& os, const Vector<complex<double> > &r)
{
 for (int i=0; i<3; i++)
  os << real(r[i]) << "  " << imag(r[i]) << "  ";
 return os; 
}


double sDreieck (Vector<double> A,Vector<double> B,Vector<double> C) 
{ 
 double a,b,c;
 double rA,rB,rC,r2,d,h; 
 double Erg; 

 rA=abs(A); rB=abs(B); rC=abs(C); r2=rA*rA;
 if ((rA==rB)&&(rB==rC)) 
 {
  d=abs(B-C)/2.0;
  h=sqrt(r2-d*d);
  a=2.0*atan2(d,h);
 
  d=abs(A-C)/2.0;
  h=sqrt(r2-d*d);
  b=2.0*atan2(d,h);

  d=abs(A-B)/2.0;
  h=sqrt(r2-d*d);
  c=2.0*atan2(d,h);

  Erg=(a+b+c-rA*M_PI)*rA; 
  return Erg; 
 }
 else errmsg ("sTriangle : Punkte liegen nicht auf einer Kugel !");
 return 0.0; 
}

double sViereck (Vector<double> A,Vector<double> B,Vector<double> C,Vector<double> D)
{
 return sDreieck(A,B,C)+sDreieck(A,C,D); 
}

Vector<complex<double> > vdc (const Vector<double> &r)
{
 Vector<complex<double> > Erg=Vector<complex<double> >(r[0],r[1],r[2]);
 return Erg;
}


Vector<double> cart2sphere(double X, double Y,  double Z)
{
 Vector<double> y,a;

  a[0] = X;
  a[1] = Y;
  a[2] = Z;
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




Vector<double> sphere2cart (double r, double theta, double phi)
{
 Vector<double> Erg;
 Erg[0]=r*sin(theta)*cos(phi);
 Erg[1]=r*sin(theta)*sin(phi);
 Erg[2]=r*cos(theta);
 return Erg;
}


void getKSystem (const Vector<double> &n, const Vector<double> &k,
                   Vector<double> &e0, Vector<double> &e1, Vector<double> &e2)
/*
  Berechnet aus der Flächennormalen n und des Ausbreitungsvektors k die Einheitsvektoren e1,e2,e3 des
  entsprechenden Koordinatensystems
*/
 {
  e0=n/abs(n);
  e2=n%k;
  if (abs(e2)==0) 
  {
   e2=n%ex;
   if (abs(e2)==0) e2=n%ey; 
  }
  e2/=abs(e2);
  e1=e2%e0;
  e1/=abs(e1);
 }


double arctan (double y, double x)
{
 double Erg;
 double ax=fabs(x);
 double ay=fabs(y);
 Erg=atan(ay/ax); 
 if (x>=0.0) 
 {
  if (y<0.0) Erg=2.0*M_PI-Erg;  
 }    
 else
 {
  if (y>=0.0) Erg=M_PI-Erg; 
  else Erg=M_PI+Erg;
 }
 return Erg; 
}

char* toString (char *S, Vector<double> P)
{  
 sprintf (S,"%f   %f   %f",P[0],P[1],P[2]);
 return S; 
}

char* toString (char *S, Vector<complex<double> > P)
{  
 sprintf (S,"%f   %f   %f   %f   %f   %f"
,real(P[0]),imag(P[0]),real(P[1]),imag(P[1]),real(P[2]),imag(P[2]));
 return S; 
}


Vector<double> emult (const Vector<double> &r1, const Vector<int> &r2)
{
 Vector<double> Erg;
 for (int i=0; i<3; i++)
 Erg[i]=r1[i]*r2[i];
 return Erg;  
}

Vector<double> emult (const Vector<int> &r1, const Vector<double> &r2)
{
 Vector<double> Erg;
 for (int i=0; i<3; i++)
 Erg[i]=r1[i]*r2[i];
 return Erg;  
}
 
Vector<double> floor (const Vector<double> &r)
{
 return Vector<double> (floor(r[0]),floor(r[1]),floor(r[2]));
}

Vector<int> ifloor (const Vector<double> &r)
{
 return Vector<int> ((int)floor(r[0]),(int)floor(r[1]),(int)floor(r[2]));
}

Vector<complex<double> > makeReal (const Vector<complex<double> > &r)
{
 Vector<complex<double> > Erg;
 for (int i=0; i<3; i++)
   Erg[i]=abs(r[i]);
 return Erg;
}

Vector<complex<double> > operator % (const Vector<complex<double> > &a, const Vector<double> &b)
{
 Vector <complex<double> > Erg;
 Erg[0]=a[1]*b[2]-a[2]*b[1];
 Erg[1]=a[2]*b[0]-a[0]*b[2];
 Erg[2]=a[0]*b[1]-a[1]*b[0];
 return Erg;
}

Vector<complex<double> > operator % (const Vector<double> &a, const Vector<complex<double> > &b)
{
 Vector <complex<double> > Erg;
 Erg[0]=a[1]*b[2]-a[2]*b[1];
 Erg[1]=a[2]*b[0]-a[0]*b[2];
 Erg[2]=a[0]*b[1]-a[1]*b[0];
 return Erg;
}

Vector<complex<double> > operator % (const Vector<complex<double> > &a, const Vector<int> &b)
{
 Vector <complex<double> > Erg;
 Erg[0]=a[1]*b[2]-a[2]*b[1];
 Erg[1]=a[2]*b[0]-a[0]*b[2];
 Erg[2]=a[0]*b[1]-a[1]*b[0];
 return Erg;
}

Vector<complex<double> > operator % (const Vector<int> &a, const Vector<complex<double> > &b)
{
 Vector <complex<double> > Erg;
 Erg[0]=a[1]*b[2]-a[2]*b[1];
 Erg[1]=a[2]*b[0]-a[0]*b[2];
 Erg[2]=a[0]*b[1]-a[1]*b[0];
 return Erg;
}

Vector<complex<double> > operator % (const Vector<double> &a, const Vector<int> &b)
{
 Vector <complex<double> > Erg;
 Erg[0]=a[1]*b[2]-a[2]*b[1];
 Erg[1]=a[2]*b[0]-a[0]*b[2];
 Erg[2]=a[0]*b[1]-a[1]*b[0];
 return Erg;
}

Vector<complex<double> > operator % (const Vector<int> &a, const Vector<double> &b)
{
 Vector <complex<double> > Erg;
 Erg[0]=a[1]*b[2]-a[2]*b[1];
 Erg[1]=a[2]*b[0]-a[0]*b[2];
 Erg[2]=a[0]*b[1]-a[1]*b[0];
 return Erg;
}



Vector<double> arg (Vector<complex<double> > &r)
{
//  return Vector<double> (arg(r[0]),arg(r[1]),arg(r[2]));
 complex <double> x,y,z;
 x=r[0];
 if (fabs(real(x))<1E-15) x=complex<double> (0.0,imag(x));
 if (fabs(imag(x))<1E-15) x=complex<double> (real(x),0.0);
 
 y=r[1];
 if (fabs(real(y))<1E-15) y=complex<double> (0.0,imag(y));
 if (fabs(imag(y))<1E-15) y=complex<double> (real(y),0.0);
 
 z=r[2]; 
 if (fabs(real(z))<1E-15) z=complex<double> (0.0,imag(z));
 if (fabs(imag(z))<1E-15) z=complex<double> (real(z),0.0);
 return Vector<double> (arg(x),arg(y),arg(z));
}

void sphereunitv (Vector<double> &P, Vector<double> &er, Vector<double> &etheta, Vector<double> &ephi)
{
 Vector <double> h=cart2sphere(P);
 er=P/abs(P);
 etheta=Vector<double> (cos(h[2])*cos(h[1]),sin(h[2])*cos(h[1]),-sin(h[1]));
 ephi=Vector<double> (-sin(h[2]),cos(h[2]),0.0);
}

Vector <double> emax (Vector<double> &P1, Vector<double> &P2)
{
 return Vector<double> (max (P1[0],P2[0]),max (P1[1],P2[1]),max (P1[2],P2[2]));
}
 

Vector<double> er(double theta, double phi)
{
 return Vector<double> (sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta));
}

Vector<double> etheta(double theta, double phi)
{
 return Vector<double> (cos(phi)*cos(theta),sin(phi)*cos(theta),-sin(theta));
}

Vector<double> ephi(double theta, double phi)
{
 return Vector<double> (-sin(phi),cos(phi),0);
}



void errmsg(char *S) 
{
 cerr << endl;
 cerr << "Fehler bei Vektorrechnung: " << S << endl;
 exit(1);
}

Vector<complex<double> > norm (const Vector<complex<double> > &r)
{
 double absr=abs(r);
 if (absr>0.0) return r/absr;
 else return czero;
}

Vector<double> norm (const Vector<double> &r)
{
 double absr=abs(r);
 if (absr>0.0) return r/absr;
 else return dzero;
}

Vector<double> nanV (const char *tagp)   
{ 
  double NaN=nan(tagp);
  return Vector<double> (NaN,NaN,NaN);
} 

bool isnan(Vector<std::complex<double> > v)
{
    bool Erg = std::isnan(real(v[0])) || std::isnan(real(v[1])) || std::isnan(real(v[2]));
    Erg = Erg || std::isnan(imag(v[0])) || std::isnan(imag(v[1])) || std::isnan(imag(v[2]));
    return Erg;
}

 bool isnan (Vector<double> v) 
{
 return std::isnan(v[0]) || std::isnan(v[1]) || std::isnan(v[2]);
}


#endif
