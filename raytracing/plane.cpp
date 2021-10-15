/***************************************************************************
                          ebene.cpp  -  description
                             -------------------
    begin                : Sat Feb 19 2000
    copyright            : (C) 2000 by Thomas Weigel
    email                : weigel@lat.ruhr-uni-bochum.de
 ***************************************************************************/


#include "plane.h"

Plane::Plane(){

}

Plane::Plane(const Vector<double>& P, const Vector<double>& e1, const Vector<double>& e2)
{
 this->P=P;
 this->e1=e1/abs(e1);
 this->e2=e2/abs(e2);
 n=e1%e2;
 n/=abs(n);
}

void Plane::norm()
{
 e1/=abs(e1);
 e2/=abs(e2);
 //n=e1%e2;
 n/=abs(n);
}

void Plane::Normalenform()
{
 Vector<double> n=e1%e2;
 Vector<double> Ps;
 Matrix<double> M1,M2;

 for (int i=0; i<3; i++)
 {
  M1(i,0)=n[i];
  M1(i,1)=-e1[i];
  M1(i,2)=-e2[i];
  }
  M2=invert(M1);
  Ps=M2*P;
  P=Ps[0]*P;
  norm();
}

void Plane::intersectSphere (Vector<double> O, double r)
{
 double a;
 double Pe1=P*e1;
  a=(-P*e1+sqrt(Pe1*Pe1-(P*P)*(e1*e1)+r*r*(e1*e1)))/(e1*e1);
}

Plane::~Plane(){
}

double Plane::distance (Vector<double> R)
{
  // Voraussetzung |n|=1 !

  double d,Erg;
  d=(P*n)/abs(n);
  Erg=fabs(n*R-d);
 return Erg;
}

void Plane::toString(char *S)
{
 sprintf (S," P=%f   %f   %f      n=%f  %f  %f\ne1=%f  %f  %f   e2=%f  %f  %f",
          P[0],P[1],P[2],n[0],n[1],n[2],e1[0],e1[1],e1[2],e2[0],e2[1],e2[2]);
}

void Plane::rotate(double dr, double dtheta,double dphi)
{
 Vector<double> P1,P2;
 Vector<double> Ps,P1s,P2s;
 P1=P+e1;
 P2=P+e2;

 Ps=cart2sphere(P);
 Ps[1]+=dtheta;
 Ps[2]+=dphi;
 P=sphere2cart(Ps);

 P1s=cart2sphere(P1);
 P1s[1]+=dtheta;
 P1s[2]+=dphi;
 P1=sphere2cart(P1s);

 P2s=cart2sphere(P2);
 P2s[1]+=dtheta;
 P2s[2]+=dphi;
 P2=sphere2cart(P2s);

 e1=P1-P;
 e2=P2-P;
 n=e1%e2;
 n/=abs(n);
 P+=dr*n;
}

void Plane::binWrite (std::ofstream &os)
{
 char svd=sizeof (Vector<double>);
 os.write((char *) &P, svd);
 os.write((char *) &e1, svd);
 os.write((char *) &e2, svd);
 os.write((char *) &n, svd);
}

void Plane::binRead (std::ifstream &is)
{
 char svd=sizeof (Vector<double>);
 is.read ((char *) &P, svd);
 is.read ((char *) &e1, svd);
 is.read ((char *) &e2, svd);
 is.read ((char *) &n, svd);
}
