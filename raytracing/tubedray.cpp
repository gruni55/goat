/***************************************************************************
                          strahl.cpp  -  description
                             -------------------
    begin                : Fri Oct 15 1999
    copyright            : (C) 1999 by Thomas Weigel
    email                : weigel@lat.ruhr-uni-bochum.de
 ***************************************************************************/

#include "tubedray.h"
#include "ellipsoid.h"
#include <fstream> 


namespace GOAT
{
    namespace raytracing
    { 
tubedRay::tubedRay(){
 g.isGauss=false;
 iR=0;
 objIndex=-1;
 OK= GOAT::maths::dzero;
 isValid=true;
 getunnelt = false;
 // KORR=1E-8;
}

tubedRay::tubedRay(const tubedRay& ray)
{
    r0 = ray.r0;
    OK=ray.OK;
    for (int i = 0; i < 5; i++)
    {
        this->k[i] = ray.k[i];
        this->P[i] = ray.P[i];
        this->E[i] = ray.E[i];
        this->phi[i] = ray.phi[i];
    }
    n = ray.n;
    n0 = ray.n0;
    k0 = ray.k0;
    r0 = ray.r0;
    
    
    numObj = ray.numObj;
    Obj = ray.Obj;
    objIndex = ray.objIndex;
    iR = ray.iR;
    pol = ray.pol;
    inObject = ray.inObject;
    isValid = ray.isValid;   
    getunnelt = false;
}

tubedRay::~tubedRay(){
}

void tubedRay::checkObjectIntersection(int Index[5], maths::Vector<double> *Pmin)
{
    maths::Vector<double> P1,n;
 double amin,a;
 bool found;
 bool inside=(objIndex>-1);
 for (int j=0; j<5; j++)
 {
  amin=100.0*r0;
  Pmin[j]=P[j];
  Index[j]=-1;
  for (int i=0; i<numObj; i++)
  {
   found=Obj[i]->next(P[j],k[j],P1);
  // cout << P[j] << "  /  " << P1 << "   k=" << k[j] << endl;
   if (found)
   {
    a=abs(P1-P[j]);
    if ((a<amin) && (a>0))
    {
      amin=a;
      Pmin[j]=P1;
      if (i!=objIndex) Index[j]=i;
    }
   }
  }
 } 
 
   found=true;
   for (int j=0; (j<=3) && found; j++)
    for (int i=j+1; (i<=4) && found; i++)
     found=found && (Index[j]==Index[i]);
   isValid=found;  
 } 


tubedRay::tubedRay(Plane E, const maths::Vector<double> &p, double dy, double dz,
         const maths::Vector<std::complex<double> > &Pol,
         std::complex<double>  n0, double r0, double k0,
         const int Anzein, ObjectShape **Einschluss,bool logRay)
{
 g.isGauss=false;
 GOAT::maths::Matrix<double> D;
 double phi=(double)rand()/(double)RAND_MAX*M_PI;
 this->logRay=logRay;
  ka=E.n;
 P[4]=p;
 D=rotMatrix(E.n,phi);
 P[0]=p+D*(-dy/2.0*E.e1+dz/2.0*E.e2); // links oben
 P[1]=p+D*(-dy/2.0*E.e1-dz/2.0*E.e2); // links unten
 P[2]=p+D*( dy/2.0*E.e1-dz/2.0*E.e2); // rechts unten
 P[3]=p+D*( dy/2.0*E.e1+dz/2.0*E.e2); // rechts oben

 
 for (int i=0; i<5; i++) k[i]=E.n;
 Obj=Einschluss;
 numObj=Anzein;
 n=n0;
 for (int i=0; i<5; i++)
   this->E[i]=Pol;
 this->r0=r0;
 this->k0=k0;
 inObject=false;
 objIndex=-1;
 OK= maths::dzero;
 rc=abs(P[4])/real(n);
 getunnelt=false;
 isValid=true;
}


tubedRay::tubedRay(const maths::Vector<double> &p, double dy, double dz,
         const maths::Vector<std::complex<double> > &Pol, const maths::Vector<double> &K,
         std::complex<double>  n0, double r0, double k0,
         const int Anzein, ObjectShape **Einschluss,bool logRay)
{
 g.isGauss=false;
 this->logRay=logRay;
  ka=K;
 
 P[4]=p;
 P[0]=p-dy/2.0* maths::ey+dz/2.0* maths::ex; // links oben
 P[1]=p-dy/2.0* maths::ey-dz/2.0* maths::ex; // links unten
 P[2]=p+dy/2.0* maths::ey-dz/2.0* maths::ex; // rechts unten
 P[3]=p+dy/2.0* maths::ey+dz/2.0* maths::ex; // rechts oben
 for (int i=0; i<5; i++) k[i]=K;
 Obj=Einschluss;
 numObj=Anzein;
 // KORR=1E-8;
 n=n0;
 for (int i=0; i<5; i++)
   E[i]=Pol;
 this->r0=r0;
 this->k0=k0;
 //init_Efeld(Pol);
 inObject=false;
 objIndex=-1;
 OK= maths::dzero;
 rc=abs(P[4])/real(n);
 getunnelt=false;
 isValid=true;
}

tubedRay::tubedRay(tubedRayBuffer &B)
{
 g.isGauss=false;
 iR=0;
 for (int i=0; i<5; i++)
 {
  P[i]=B.P[i];
  k[i]=B.k[i];
  E[i]=B.E[i];
 }
 ka=B.k[4];
 n=B.n;
 Obj=B.Obj;
 r0=B.r0;
 inObject=false;
 getunnelt=false;
 OK= maths::dzero;
 numObj=0;
 isValid=true;
}

bool tubedRay::next()
{
 int Index[5];
 maths::Vector<double> R[5];
 bool found2,found=true;
 bool Einschluss=true;
 if (!inObject)
 {
  if (numObj==0) Index[4]=-1;
  else 
  {    checkObjectIntersection(Index,R);
     if (!isValid) return false;		 
  } 

  if(/*Index[0]==-1 || Index[1]==-1 || Index[2]==-1 || Index[3]==-1 ||*/ Index[4]==-1) // Es wurde kein Einschluss gefunden !
  {
   
   for (int i=0; i<5; i++)
   {
      R[i]=nextP(P[i],k[i], maths::dzero,r0,found2);
	  found=found && found2;
   }
   if (!found) { objIndex=-1;  return found;}
  }
  
  objIndex=Index[4];
  
  for (int i = 0; i < 5; i++) E[i] = E[i]  *exp(I * k0 * n * abs(R[i] - P[i]));
  
 }
 else
 {
  for (int i=0;i<5;i++)
  {
   Obj[objIndex]->next(P[i],k[i],R[i]);
   E[i]=E[i]*exp(I*k0*Obj[objIndex]->n*abs(R[i]-P[i]));  
  }
 
 }


 
 if (abs(P[4]-R[4])>1E-15) for (int i=0;i<5;i++) this->P[i]=R[i];
 else objIndex=-1;

return found;
}




void tubedRay::refract(maths::Vector<double> *N, std::complex<double>  n1, std::complex<double>  n2)
{
 double s,det;
 maths::Matrix<double> H,R;
 maths::Matrix<std::complex<double> > T;
 maths::Matrix<double> D;
 maths::Vector <double> n,e0,e1,e2;
 double alpha, gamma;
 double  beta;
 isValid=true;
   for (int i=0; i<5; i++)
 {
  //int i=4;
  n=N[i]/abs(N[i]);
 
  getKSystem (n,k[i],e0,e1,e2);
  det=n*k[i]/(abs(n)*abs(k[i]));
  if (det<-1) det=-1;
  else
	  if (det>1) det=1;
  alpha= std::acos(det);
  if (alpha>M_PI/2.0) { alpha=M_PI-alpha; e2=-e2; }

  beta=real(asin((std::complex<double>)real(n1)/real(n2)*sin(alpha)));
  isValid=isValid ; //&& (fabs(beta)>EPS_WINKEL);
//  if (real(n1/n2*sin(alpha))>=1.0) cout << "P=" << P[4] << "    " << n1/n2*sin(alpha) << endl;
  gamma=beta-alpha;
  //if (real(n2)<real(n1)) { e2=-e2;}
  s=1.0;
  trafo(e0,e1,e2,H,R);

  D=rotMatrixA(e2,k[i],gamma);
  T=Fresnel_trans(alpha,beta,n1,n2);
  k[i]=D*k[i];
  E[i]=H*E[i];
  E[i]=T*E[i];
  E[i]=R*E[i];
  E[i]=D*E[i];
 }
}



tubedRay tubedRay::reflect(maths::Vector<double> *n, std::complex<double>  n1, std::complex<double>  n2)
/* Strahl wird an einer Oberflaeche reflektiert. Wird an einer Einschlussoberflaeche
   reflektiert (objIndex >-1 ), dann wird der transmittierte Strahl zurueckgegeben */
{
 tubedRay Erg;
 maths::Matrix <double> D,H,R;
 maths::Matrix<std::complex<double> > F;
 maths::Vector <double> Ph,e0,e1,e2;
 double det,alpha,gamma;
  if (inObject)
  {
  // Strahl will aus dem Einschluss raus
  Erg=*this;
  Erg.OK= maths::dzero;
  Erg.inObject=false;
  Erg.n=n2;
  Erg.r0=r0;
  Erg.objIndex=-1;
  Erg.refract(n,n1,n2);
  Erg.isValid=true;
  /*if (Ein[objIndex]->type==OBJECTSHAPE_ELLIPSOID)
  { */
/*   for (int i=0; i<5; i++)
   {
    Erg.P[i]=P[i]+Ein[objIndex]->norm (P[i])*KORR;
    P[i]=P[i]-Ein[objIndex]->norm (P[i])*KORR;
   }*/
//   }
  } // if inObject
  else
  {
   if (objIndex>-1)
   {
    // Einschluss wurde getroffen
    Erg=*this;
    Erg.OK=Obj[objIndex]->P;
    Erg.inObject=true;
    Erg.n=Obj[objIndex]->n;
    Erg.objIndex=objIndex;
    Erg.refract(n,n1,n2);
    Erg.getunnelt=false;
  /*  if (Ein[objIndex]->type==OBJECTSHAPE_ELLIPSOID)
    {*/
/*     for (int i=0; i<5; i++)
     {
       Erg.P[i]=P[i]-Ein[objIndex]->norm (P[i])*KORR;
       P[i]=P[i]+Ein[objIndex]->norm(P[i])*KORR;
     }*/
//    }
    objIndex=-1;
   } // if objIndex!=-1
   else
   {
    Erg=*this;
    Erg.OK= maths::dzero;
    Erg.inObject=false;
    Erg.n=n2;
    Erg.objIndex=-1;
    Erg.refract(n,n1,n2);
/*    for (int i=0; i<5; i++) 
     P[i]=P[i]/abs(P[i])*(1.0-KORR)*r0;*/
   }
  } // else if inObject
 // cout << "vorher: P=" << P[4] << "  k=" << k[4] << endl; 
  for (int i=0; i<5; i++)
 {
   n[i]=-n[i];
    getKSystem(n[i],k[i],e0,e1,e2);
  det=n[i]*k[i]/(abs(n[i])*abs(k[i]));
  if (det<-1) det=-1;
  else
	  if (det>1) det=1;
  alpha=std::acos(det);
  // isValid=isValid && (fabs(alpha)>EPS_WINKEL);
  gamma=M_PI-2.0*alpha;
  if (alpha>M_PI/2.0) { gamma=-gamma; alpha=M_PI-alpha;}
  if (objIndex>-1) e0=-e0;
  trafo(e0,e1,e2,H,R);
  F=Fresnel_reflect(alpha,n1,n2);
  D=rotMatrixA(e2,k[i],gamma);
  k[i]=D*k[i];
  E[i]=H*E[i];
  E[i]=F*E[i];
  E[i]=R*E[i];
  E[i]=D*E[i];
 } // for i
 /*cout << "gamma=" << gamma/M_PI*180.0 << endl;
 cout << "D=" << D << endl;
 cout << "nachher : k=" << k[4] << endl;
 cout << "------------------------------------------------------------" << endl;
*/
  iR++;
   
 return Erg;
}



void tubedRay::refract(maths::Vector<double> N, std::complex<double>  n1, std::complex<double>  n2)
{
 double s;
 maths::Matrix<double> H,R;
 maths::Matrix<std::complex<double> > T;
 maths::Matrix<double> D;
 maths::Vector <double> n,e0,e1,e2;
 double ha,alpha, gamma;
 double  beta;

  n=N/abs(N);
  isValid=true;
 for (int i=0; i<5; i++)
 {
  //int i=4;
  getKSystem (n,k[i],e0,e1,e2);
//  cout << "e0=" << e0 << "   e1=" << e1 << "   e2=" << e2 << endl;
  ha=n*k[i]/abs(k[i]); 
  if (fabs(ha)>=1) alpha=0;  // <- um Ungenauigkeiten auszugleichen (sollte real eigentlich nicht vorkommen... ;-) )
  else
  {
	  alpha=std::acos(ha);
       if (alpha>M_PI/2.0) { alpha=M_PI-alpha; e2=-e2; }
  }

  beta=real(asin((std::complex<double>)real(n1)/real(n2)*sin(alpha)));
 //  isValid=isValid && fabs(beta)>EPS_WINKEL;  // HÄH ????
  gamma=beta-alpha;
  s=1.0;
  trafo(e0,e1,e2,H,R);

  D=rotMatrixA(e2,k[i],gamma);
  T=Fresnel_trans(alpha,beta,n1,n2);
  
  k[i]=D*k[i];
  E[i]=H*E[i];
  E[i]=T*E[i];
  E[i]=R*E[i];
  E[i]=D*E[i];
 }
}



void tubedRay::reflectRay(RayBase* &tray, maths::Vector<double> n, std::complex<double> n1, std::complex<double> n2)
{

    tubedRay Tray= reflect(n, n1, n2);
    *(tubedRay *)tray = tubedRay(Tray);
}

tubedRay tubedRay::reflect(maths::Vector<double> n, std::complex<double>  n1, std::complex<double>  n2)
/* Strahl wird an einer Oberflaeche reflektiert. Wird an einer Einschlussoberflaeche
   reflektiert (objIndex >-1 ), dann wird der transmittierte Strahl zurueckgegeben */
{
 tubedRay Erg;
 maths::Matrix <double> D,H,R;
 maths::Matrix<std::complex<double> > F;
 maths::Vector <double> Ph,e0,e1,e2;
 double ha;
 double alpha,gamma;
 isValid=true;
  if (inObject)
  {
  // Strahl will aus dem Einschluss raus
  Erg=*this;
  Erg.OK= maths::dzero;
  Erg.inObject=false;
  Erg.n=n2;
  Erg.r0=r0;
  Erg.objIndex=-1;
  Erg.refract(n,n1,n2);
  /*if (Ein[objIndex]->type==OBJECTSHAPE_ELLIPSOID)
  { */
/*   for (int i=0; i<5; i++)
   {
    Erg.P[i]=P[i]+Ein[objIndex]->norm (P[i])*KORR;
    P[i]=P[i]-Ein[objIndex]->norm (P[i])*KORR;
   }*/
//   }
  } // if inObject
  else
  {
   if (objIndex!=-1)
   {
    // Einschluss wurde getroffen
    Erg=*this;
    Erg.OK=Obj[objIndex]->P;
    Erg.inObject=true;
    Erg.n=Obj[objIndex]->n;
    Erg.objIndex=objIndex;
    Erg.refract(n,n1,n2);
    Erg.getunnelt=false;
  /*  if (Ein[objIndex]->type==OBJECTSHAPE_ELLIPSOID)
    {*/
/*     for (int i=0; i<5; i++)
     {
       Erg.P[i]=P[i]-Ein[objIndex]->norm (P[i])*KORR;
       P[i]=P[i]+Ein[objIndex]->norm(P[i])*KORR;
     }*/
//    }
    objIndex=-1;
   } // if objIndex!=-1
   else
   {
    Erg=*this;
    Erg.OK= maths::dzero;
    Erg.inObject=false;
    Erg.n=n2;
    Erg.objIndex=-1;
    Erg.refract(n,n1,n2);
/*    for (int i=0; i<5; i++) 
     P[i]=P[i]/abs(P[i])*(1.0-KORR)*r0;*/
   }
  } // else if inObject

 n=-n;

 maths::Vector<std::complex<double> >hE=E[4];
 for (int i=0; i<5; i++)
 {
    getKSystem(n,k[i],e0,e1,e2);
	ha=n*k[i]/(abs(n)*abs(k[i]));
	if (fabs(ha)>1) alpha=0;
  else alpha=std::acos(ha);
  isValid=isValid && fabs(alpha)>EPS_WINKEL;
  gamma=M_PI-2.0*alpha;
  if (alpha>M_PI/2.0) { gamma=-gamma; alpha=M_PI-alpha;}
  if (objIndex>-1) e0=-e0;
  trafo(e0,e1,e2,H,R);
  F=Fresnel_reflect(alpha,n1,n2);
  D=rotMatrixA(e2,k[i],gamma);
    k[i]=D*k[i];
  E[i]=H*E[i];
  E[i]=F*E[i];
  E[i]=R*E[i];
  E[i]=D*E[i];
  
 } // for i


  iR++;
 return Erg;
}

GOAT::maths::Matrix<std::complex<double> > tubedRay::Fresnel_reflect (double alpha, std::complex<double>  n1, std::complex<double>  n2)
{
 std::complex<double>  rs,rp;
 GOAT::maths::Matrix<std::complex<double> > R;
 double  n12;
 double  beta;
 n12=real(n2)/real(n1);
 int l;
 double x;
 beta=real(asin((std::complex<double> ) sin(alpha) / n12));
  
 
  if (alpha==0) { rp=(n2-n1)/(n2+n1); rs=-rp; }
  else
  {
   rp=tan(alpha-beta)/tan(alpha+beta);
   rs=-sin(alpha-beta)/sin(alpha+beta);
  }
 
 
 R(0,0)=rp;
 R(1,1)=rp;
 R(2,2)=rs;
 return R;
}

GOAT::maths::Matrix<std::complex<double> > tubedRay::Fresnel_trans (double alpha,std::complex<double>  beta, std::complex<double>  n1, std::complex<double>  n2)
{
 double rbeta=real(beta);
 std::complex<double>  ts,tp;
 GOAT::maths::Matrix<std::complex<double> > T;
 
 if (alpha==0.0)
 {
  ts=2.0*n1/(n1+n2);
  tp=ts; 
 }
 else 
 {
   ts=2.0 * sin(rbeta) * cos(alpha) / sin(alpha+rbeta);
   tp=2.0 * sin(rbeta) * cos(alpha) / (sin(alpha+rbeta)*cos(alpha-rbeta));

    
 }

 T(0,0)=tp;
 T(1,1)=tp;
 T(2,2)=ts;
 return T;
}

void tubedRay::initElectricFieldGauss (int i, const Plane& Eb,
                                   const maths::Vector<std::complex<double> >& Pol)
{ 
 k[i]=g.F-P[i];
 k[i]=k[i]/abs(k[i]);
 double theta;
 double z=abs(P[i]-g.F);
 double r;
 double R;
 double l0=2.0*M_PI/k0;
 maths::Vector<double> hr,h;
  std::complex<double> E0;
  double w,G; 
  double z0=M_PI*g.w0*g.w0/l0;
  theta=atan2(g.w0,z0);
  r=sin(theta/2.0)*z;
  w=g.w0*sqrt(1.0+z*z/(z0*z0));
  R=z*sqrt(1.0+z*z/(z0*z0));
  G=atan2(z,z0);
  E0=g.w0/w*exp((-r*r/w*w)-I*k0*r*r/(2*R)-I*k0*z+I*G);  
  double l;
  double Pk=P[i]*k[i];
  l=-Pk-sqrt(Pk*Pk-P[i]*P[i]+r0*r0);
  P[i]=P[i]+l*k[i];  
  E[i]=Pol*E0*exp(I*k0*l);    
}


void tubedRay::initElectricField (const maths::Vector<std::complex<double> >& Pol,const int
AnzRays, double dx, Plane Eb)
{
  double dx2=dx/2.0;
  char Str[255];
  double w2,R,b,w02,z;
  double l,r,r2;
  maths::Vector<double> h,hr;
  bool found;
  isValid=true;
  double x1,x2,xn,p;
  maths::Vector<double> rh;
  int i;
  if (!getunnelt)
  { 
   for (i=0; i<5; i++)
   {
      k[i]/=abs(k[i]);
      if ((logRay)&&(i==4))
      {
      toString (Str,P[i]);
      }  // if logRay
    if (g.isGauss)
    { // ----------- Gaussstrahl ---------------
      for (i=0; i<5; i++) initElectricFieldGauss( i,Eb,Pol); 
      } // if Gaussstrahl
      else
      {
//       cout << "jau" << endl;
       // ----------- einfallende Plane Welle -------------------

	  hr=P[i];	   
	  h=P[i]-(P[i]*k[i])*k[i];
      P[i]=h-sqrt(r0*r0-abs(h)*abs(h))*k[i];
      E[i]=Pol*exp(I*abs(hr-P[i])*k0); 
      }


      if ((logRay) && (!getunnelt) && (i==4))
      {
       toString (Str,P[i]);
      }
   } // for
  }
  else // getunnelt!
  {
    /*Vector<double> kall, nall, Pall;   
    Vector< std::complex<double> > Eall;
    kall = k[4]/abs(k[4]);

       hr=P[4];
       h=P[4]-(P[4]*kall)*kall;
       Pall=h-sqrt(r0*r0-abs(h)*abs(h))*kall;
       Eall=Pol*exp(I*abs(hr-Pall)*k0);

    for (int i=0; i<5; i++)
    {
    E[i]=Eall;
    P[i]=Pall;
    k[i]=kall;*/

    double alpha,b;
    maths::Vector<double> er;
    for (int i=0; i<5; i++)
    {
     E[i]=Pol*exp(I*k0*r0);
     P[i]=(Eb.e1*P[i])*Eb.e1+(Eb.e2*P[i])*Eb.e2;
     b=abs(P[i]);
     alpha=std::acos(b/r0/real(n));
     er=P[i]/abs(P[i]); 
     k[i]=cos(alpha)*Eb.n-sin(alpha)*er;
     P[i]=r0*er;
    }
  }
}


void tubedRay::initElectricField (const maths::Vector<std::complex<double> >& Pol,int AnzRays)
{
  char Str[255];
  //QString HStr;
  double absh,w2,R,b,w02,z;
  double r,r2;
  maths::Vector<double> h,hr;
  bool found;
double x1,x2,xn,p;
maths::Vector<double> rh;
  isValid=true;
//   cout << "P Anfang:" << P[0] << P[1] << P[2] << P[3] << P[4] << endl;
  if (!getunnelt)
  {
   for (int i=0; i<5; i++)
   {
      k[i]/=abs(k[i]);
      if ((logRay)&&(i==4))
      {
      toString (Str,P[i]);
     /* HStr=Str;
      HStr+="   ";*/
      }  // if logRay



    if (g.isGauss)
    { // ----------- Gaussstrahl ---------------
      double zi0,h2,alpha,cosa,absh,phi,det,l,R;

      maths::Matrix<double> D;
      maths::Vector<double> N,Pkt,Ps,d,hk;

      hk=k[i];
      d=P[i]-g.F*r0;
      r=sqrt(d[0]*d[0]+d[1]*d[1]);
      w02=g.w0*g.w0*r0*r0;
      zi0=k0*w02/2.0;
      h2=k0*w02/(2*d[2]);
      R=d[2]*(1+h2*h2);

      // Normale zur Phasenfläche
      N[0]=k0*d[0]/R;
      N[1]=k0*d[1]/R;
      N[2]=-1/(zi0*(1+d[2]*d[2]))-k0+k0*r*r/(R*R)*(1-k0*k0*w02*w02/(2*d[2]*d[2]));
      k[i]=N/abs(N);

      // Nächsten Punkt auf der Kugel berechnen
      det=(P[i]*k[i])*(P[i]*k[i])-(abs2(P[i])-r0*r0);
      if (det>0) // Schnittpunkt gefunden
      {
       l=(-P[i]*k[i]-sqrt(det));
       P[i]=P[i]+l*k[i];
       h2=2*d[2]/(k0*w02);
       w2=w02*(1+h2*h2);
      alpha=abs(maths::ez*k[i]);
      phi=atan(d[2]/zi0)-k0*(d[2]+r*r/(2*R));
      E[i]=Pol*g.w0*r0/w2*exp(-r*r/w2)*exp(I*phi);
      D=rotMatrix(k[i]%hk,alpha);
      E[i]=D*E[i];
      }
      else isValid=false; // P[i]=Vector<double>(0,0,0);

      //else P[i]=0.0;
      } // if Gaussstrahl
      else
      {
//       cout << "jau" << endl;
       // ----------- einfallende Plane Welle -------------------
     /* x1=P[i]*ex;
      x2=P[i]*ey;
      xn=P[i]*ez;
      rh=P[i];
      p=x1*x1+x2*x2;
      P[i]=x1*ex+x2*ey-sqrt(fabs(r0*r0-p))*ez;
      E[i]=Pol*exp(I*abs(rh-P[i])*k0);
*/


      hr=P[i];	   
	  h=P[i]-(P[i]*k[i])*k[i];
      P[i]=h-sqrt(r0*r0-abs(h)*abs(h))*k[i];
      E[i]=Pol*exp(I*abs(hr-P[i])*k0); 
      }


      if ((logRay) && (!getunnelt) && (i==4))
      {
       toString (Str,P[i]);
       //HStr+=Str; 
       //ausgabe->print(HStr);
      }
   } // for
  }
  else // getunnelt!
  {
    double lambda, alpha,gamma;
    maths::Vector<double> N;
    maths::Matrix<double> D;
    lambda=-P[4]*k[4];
    P[4]=P[4]+lambda*k[4];
    alpha=std::acos(abs(P[4])/(r0*real(n)));
     N=P[4]/abs(P[4]);
    D=rotMatrixA(N,k[4],alpha);
    P[4]=P[4]/abs(P[4])*r0; 
    k[4]=D*k[4];
    for (int i=0; i<4; i++)
    {
     P[i]=P[4];
     k[i]=k[4];
     E[i]=D*Pol;
     E[i]*=exp(I*k0*n*lambda);
    }
  }

/*  else
  {
   // ------------- getunnelter Strahl ------------------
   Vector<double> e0,e1,e2,n;
   double alpha,d,gamma;
   std::complex<double>  ts,tp,rs,rp;
   double FR;
   int l;
   double m;
   double x;

   GOAT::maths::Matrix<double> D,H,R;
   GOAT::maths::Matrix<std::complex<double> > F;

   for (int i=0; i<5; i++)
   {
    E[i]=Pol*exp(I*k0*r0)/sqrt(AnzRays);
    P[i]=P[i]+k[i]*r0;
    d=abs(P[i]);
    P[i]=P[i]/d*r0*(1-KORR);
    n=P[i]/r0;
    getKSystem(n,k[i],e0,e1,e2);
   alpha=acos(d/(r0*real(Parms.n0)));
  if (objIndex>-1) e0=-e0;
    trafo(e0,e1,e2,H,R);

  D=rotMatrixA(e2,k[i],alpha);
    x=r0*k0;
//  l=rint(sin(alpha)*x*real(Parms.n0)-0.5);
   l=rint(d*Parms.k0-0.5);
  ts=teneu(l,SENKRECHT,1.0,Parms.n0,x);
  tp=teneu(l,PARALLEL,1.0,Parms.n0,x);

  F(0,0)=tp;
  F(1,1)=tp;
  F(2,2)=ts;


//    D=Dy(alpha);
  k[i]=D*k[i];
  E[i]=H*E[i];
  E[i]=F*E[i];
  E[i]=R*E[i];
  E[i]=D*E[i];
  rp=conj(ri(l,PARALLEL,1.0,Parms.n0,x));
  rs=conj(ri(l,SENKRECHT,1.0,Parms.n0,x));
  if (abs(rp)>abs(rs)) FR=abs(rp); else FR=abs(rs);
  m=0.5*log(Grenze)/log(FR);
  E[i]*=pow(FR,Parms.AnzReflex-m);
   }
  } */
}
/*void Strahl::init_Efeld (const Vector<std::complex<double> >& Pol,const int AnzRays)
{
  char Str[255];
  //QString HStr;
  double w2,R,b,w02,z;
  double r,r2;
  Vector<double> h,hr;
  bool found;
  isValid=true;
  if (!getunnelt)
  {
   for (int i=0; i<5; i++)
   {
      k[i]/=abs(k[i]);
      if ((logRay)&&(i==4))
      {
      toString (Str,P[i]);
    }  // if logRay



    if (g.isGauss)
    { // ----------- Gaussstrahl ---------------
      double zi0,h2,alpha,cosa,absh,phi,det,l,R;

      GOAT::maths::Matrix<double> D;
      Vector<double> n,Pkt,Ps,d,hk;

      hk=k[i];
      d=P[i]-g.F*r0;
      r=sqrt(d[0]*d[0]+d[1]*d[1]);
      w02=g.w0*g.w0*r0*r0;
      zi0=k0*w02/2.0;
      h2=k0*w02/(2*d[2]);
      R=d[2]*(1+h2*h2);

      // Normale zur Phasenflï¿½he
      n[0]=k0*d[0]/R;
      n[1]=k0*d[1]/R;
      n[2]=-1/(zi0*(1+d[2]*d[2]))-k0+k0*r*r/(R*R)*(1-k0*k0*w02*w02/(2*d[2]*d[2]));
      k[i]=n/abs(n);

      // Nï¿½hsten Punkt auf der Kugel berechnen
      det=(P[i]*k[i])*(P[i]*k[i])-(abs2(P[i])-r0*r0);
      if (det>0) // Schnittpunkt gefunden
      {
       l=(-P[i]*k[i]-sqrt(det));
       P[i]=P[i]+l*k[i];
       h2=2*d[2]/(k0*w02);
       w2=w02*(1+h2*h2);
      alpha=abs(ez*k[i]);
      phi=atan(d[2]/zi0)-k0*(d[2]+r*r/(2*R));
      E[i]=Pol*g.w0*r0/w2*exp(-r*r/w2)*exp(I*phi);
      D=rotMatrix(k[i]%hk,alpha);
      E[i]=D*E[i];
      P[i]=P[i]*(1.0-KORR);
      }
      else isValid=false; // P[i]=Vector<double>(0,0,0);

      //else P[i]=0.0;
      } // if Gaussstrahl
      else
      {
       // ----------- einfallende Plane Welle -------------------
       hr=P[i];
       h=P[i]-(P[i]*k[i])*k[i];
       P[i]=h-sqrt(r0*r0-abs(h)*abs(h))*k[i];
       E[i]=Pol*exp(I*abs(hr-P[i])*k0); 
      }


      if ((logRay) && (!getunnelt) && (i==4))
      {
       toString (Str,P[i]);
       //HStr+=Str;
       //ausgabe->print(HStr);
      }
   } // for
  }
  else
  {
   // ------------- getunnelter Strahl ------------------
   Vector<double> e0,e1,e2,n;
   double alpha,d,gamma;
   std::complex<double>  ts,tp,rs,rp;
   double FR;
   int l;
   double m;
   double x;

   GOAT::maths::Matrix<double> D,H,R;
   GOAT::maths::Matrix<std::complex<double> > F;

   for (int i=0; i<5; i++)
   {
    E[i]=Pol*exp(I*k0*r0)/sqrt((double)AnzRays);
    P[i]=P[i]+k[i]*r0;
    d=abs(P[i]);
    P[i]=P[i]/d*r0*(1-KORR);
    n=P[i]/r0;
    getKSystem(n,k[i],e0,e1,e2);
   alpha=acos(d/(r0*real(Parms.n0)));
  if (objIndex>-1) e0=-e0;
    trafo(e0,e1,e2,H,R);

  D=rotMatrixA(e2,k[i],alpha);
    x=r0*k0;
//  l=rint(sin(alpha)*x*real(Parms.n0)-0.5);
   l=rint(d*Parms.k0-0.5);
  ts=teneu(l,SENKRECHT,1.0,Parms.n0,x);
  tp=teneu(l,PARALLEL,1.0,Parms.n0,x);

  F(0,0)=tp;
  F(1,1)=tp;
  F(2,2)=ts;


//    D=Dy(alpha);
  k[i]=D*k[i];
  E[i]=H*E[i];
  E[i]=F*E[i];
  E[i]=R*E[i];
  E[i]=D*E[i];
  rp=conj(ri(l,PARALLEL,1.0,Parms.n0,x));
  rs=conj(ri(l,SENKRECHT,1.0,Parms.n0,x));
  if (abs(rp)>abs(rs)) FR=abs(rp); else FR=abs(rs);
  m=0.5*log(Grenze)/log(FR);
  E[i]*=pow(FR,Parms.AnzReflex-m);
   }
  }
}*/


bool tubedRay::schneidePlane (maths::Vector<double> *Erg, const Plane &E)
{
 GOAT::maths::Matrix<double> M1,M2;
 maths::Vector<double> Ps,P1;
 bool invertierbar;
  for (int i=0; i<3; i++)
  {
   M1(i,0)=E.e1[i];
   M1(i,1)=E.e2[i];
   M1(i,2)=-k[4][i];
   }
   M2=invert(M1,invertierbar);
   if (!invertierbar) return false;
   Ps=P[4]-E.P;
   P1=M2*Ps;
   if (P1[2]<0.0) return false;
   Erg[4]=P[4]+P1[2]*k[4];
  return true;
 
}

maths::Vector<double> tubedRay::schneidePlane( const Plane &E, bool &found)
{
    maths::Vector<double> Ps=E.P-P[4];
 double l;
 double kn=k[4]*E.n;
 if (fabs(kn)<1E-10) 
 { 
 /* if (Ps*E.n==0.0) { found=true; return P[4]; }
  else { found=false;  return P[4]; }*/
  found=false;  return P[4];
 }
 l=(Ps*E.n)/kn;
 found=(l>0); 
 return P[4]+l*k[4];    
}

int tubedRay::schneidePlane(const Plane &E, double d, maths::Vector<double> &S1, maths::Vector<double> &S2)
{
 double lp,lm,l;
  if (k[4]*E.n!=0.0)
  {
    lp=-d/2.0-(P[4]-E.P)*E.n;
    lm= d/2.0-(P[4]-E.P)*E.n;
    if ((lp<0) || (lm<0)) return KEIN_SCHNITTPUNKT;
    else
    {
      if (lp<lm) { l=lm; lm=lp; lp=l; }
      S1=P[4]+lm*k[4];
      S2=P[4]+lp*k[4];
      return GESCHNITTEN;
    }
   }
   else
   {
     l=E.n*(P[4] % E.P);
     if (fabs(l)<=d/2.0) return IN_Plane;
   }
   return KEIN_SCHNITTPUNKT;
}
/*
double Strahl::flaeche ()
{
 // Berechne Querschnittsflaeche aus den Randpunkten
 // durch Aufteilung in 2 Teildreiecke 
 Vector<double> cp1,cp2,cp3,P01,P02,P03;
 double F1,F2;

 P01=P[1]-P[0];
 P02=P[2]-P[0];
 P03=P[3]-P[0];
 cp1=P01 % P03;
 cp2=P02 % P03;
 
 if (cp1*cp2>0) {F1=abs(cp1)/2.0; F2=abs(cp2)/2.0; }
 else {cp1=P01 % P02; F1=abs(cp1)/2.0; F2=abs(cp2)/2.0; }
 return F1+F2;
}*/


double tubedRay::flaeche ()
{
    maths::Vector<double> a,b,c,d;
	a=P[0]-P[1];
	b=P[2]-P[1];
	c=P[0]-P[3];
	d=P[2]-P[3];
    
	double F1,F2;
	F1=abs(a % c)/2.0;
	F2=abs(b % d)/2.0;
	return F1+F2;
}


 double tubedRay::normVol (maths::Vector<double> P[5], maths::Vector<double> k, maths::Vector<double> n)
  // Verhaeltnis zwischen Volumen eines Strahls der durch ein Volumenelement
  // geht zum Volumen des Volumenelements
 {
     maths::Vector <double> a10,a12,a20,a23,h;
  double dx,V,V0;
  dx=0;
  a10=P[0]-P[1];
  a12=P[2]-P[1];
  a20=P[2]-P[0];
  a23=P[3]-P[2];
  //dx=r0/Parms.nx;
  h=k/abs(k)*dx;
  V0=1.0;
  //V0=dx*dx*Parms.r0/Parms.ny;
  V=(abs((a10%a12)*h)+abs((a10%a12)*h))/2.0;
  return V/V0;
 }

 /* double Strahl::pjump (Vector<double> P1[5],Vector<double> P2[5])
  // gibt Phasenspruenge auf Grund von Fokussierung zurueck
  // Parameter:
  // P1 : Randpunkte des "Strahls" am Anfang
  // P2 :    "        "      "     am Ende

 {
  double pj=0.0;
  double k;
  k=cross (P1[0],P1[1],P2[0],P2[1]);
  if ((k>0.0) && (k<1.0)) pj=-M_PI/2.0;
  else
  {
   k=cross (P1[2],P1[3],P2[2],P2[3]);
   if ((k>0.0) && (k<1.0)) pj=-M_PI/2.0;
  }

  k=cross (P1[0],P1[3],P2[0],P2[3]);
  if ((k>0.0) && (k<1.0)) pj-=M_PI/2.0;
  else
  {
   k=cross (P1[1],P1[2],P2[1],P2[2]);
   if ((k>0.0) && (k<1.0)) pj-=M_PI/2.0;
  }
  return pj;
 }*/

/* double Strahl::pjump (Vector<double> P1[5],Vector<double> P2[5],const double epsilon)
{
 double l3,l1;
 double d3,d1,pj;
 d3=abs(P2[3]-P1[3]);
 d1=abs(P2[1]-P1[1]);
 Vector<double> h1,h3;
 h3=P1[3]-P1[0];
 pj=0.0;
 if (abs(k[0]*k[3])>epsilon) 
 {
  l3=((h3-(h3*k[0])*k[0])*k[3])/(1-(((k[3]*k[0])*k[0])*k[3]));
  if ( (l3<d3) && (l3>0)) { pj-=M_PI_2; cout << P1[3]+l3*k[3]; }
  else cout << dzero; 
 }
 
 if (abs(k[0]*k[1])>epsilon)
 {
  l1=((h1-(h1*k[0])*k[0])*k[1])/(1-(((k[1]*k[0])*k[0])*k[1]));
  if ( (l1<d1) && (l1>0)) { pj-=M_PI_2; cout << "   " << P1[1]+l1*k[1] << endl; }
  else cout << "   " << dzero << endl;
 }
 return pj;
} */

 maths::Vector<double> Schnittpunkt (maths::Vector<double> P, maths::Vector<double> k, maths::Vector<double> P1, maths::Vector<double> k1)
{
     maths::Vector<double> A,B;
 A=P-P1+((P1-P)*k)*k;
 B=(k1*k)*k-k1;
 double l1=-A*B/(B*B);
 return P1+l1*k1;
}


double tubedRay::pjump (maths::Vector<double> P1[5], maths::Vector<double> P2[5],const double epsilon)
{
 double D1,D2;
 double pj=0;
 double l=1E-10*r0;

 D1=abs(P1[3]-P1[0]+l*(k[3]-k[0]))-abs(P1[3]-P1[0]);
 D2=abs(P2[3]-P2[0]+l*(k[3]-k[0]))-abs(P2[3]-P2[0]);
 if (D1*D2<0) { pj-=M_PI_2; /*cout << Schnittpunkt (P[0],k[0],P[3],k[3]) << "   0 0 0" << endl; */}
 
 D1=abs(P1[1]-P1[0]+l*(k[1]-k[0]))-abs(P1[1]-P1[0]);
 D2=abs(P2[1]-P2[0]+l*(k[1]-k[0]))-abs(P2[1]-P2[0]);
 if (D1*D2<0) { pj-=M_PI_2;/* cout << "0 0 0   " << Schnittpunkt (P[0],k[0],P[1],k[1]) << endl;*/ }
 return pj; 
}



double tubedRay::pjump(void)
{
 double pj=0;
 if (abs(P[0]-P[3]+k[0]-k[3])<abs(P[0]-P[3])) { pj-=M_PI_2; /*cout << Schnittpunkt (P[0],k[0],P[3],k[3]) << "   0 0 0" << endl; */}
 if (abs(P[0]-P[1]+k[0]-k[1])<abs(P[0]-P[1])) { pj-=M_PI_2; /*cout << Schnittpunkt (P[0],k[0],P[1],k[1]) << "   0 0 0" << endl;*/ }
 return pj;
}


double tubedRay::pjump (maths::Vector<double> P1[5], maths::Vector<double> P2[5], maths::Vector<double> *S)
 { 
  double pj=0.0;
  maths::Vector<double> h1,h2;
  double d1,d2,a,b; 
  double l1,l2,l;
  int i;
  S[0]= maths::dzero;
  S[1]= maths::dzero;
  d1=abs(P1[1]-P2[1]);
  d2=abs(P1[2]-P2[2]);

  a=k[0]*k[1];
  b=k[3]*k[2];
  h1=P1[1]-P1[0];
  h2=P1[2]-P1[3];
  l1=(h1*k[1]-(h1*k[0])*a)/(a*a-1.0);
  l2=(h2*k[2]-(h2*k[3])*b)/(b*b-1.0);
  
  if ((l1>0) && (l1<d1)) { l=l1; i=1; }
  else { l=l2; i=2; }
   
  if ( ((l1>0) && (l1<d1)) || ((l2>0) && (l2<d2)) )
   { 
              pj-=M_PI_2; 
	      S[0]=P1[i]+l*k[i];
   }	    	    
  d1=abs(P1[3]-P2[3]);
  d2=abs(P1[1]-P2[1]);	    
  a=k[0]*k[3];
  b=k[1]*k[2];
  h1=P1[3]-P1[0];
  h2=P1[1]-P1[2];
  l1=(h1*k[3]-(h1*k[0])*a)/(a*a-1.0);
  l2=(h2*k[1]-(h2*k[2])*b)/(b*b-1.0);
  if ((l1>0) && (l1<d1)) { l=l1; i=3; }
  else { l=l2; i=1; }
  
  if ( ((l1>0) && (l1<d1)) || ((l2>0) && (l2<d2)) )
            { 
              pj-=M_PI_2; 
	      S[1]=P1[0]+l1*k[0];
            }
   return pj;
  
 }


 
 double tubedRay::cross (const maths::Vector<double> P10, const maths::Vector<double> P11,
                       const maths::Vector<double> P20, const maths::Vector<double> P21)
 {
     maths::Vector<double> a0,a1;
  double k1,k2,k;

  a0=P20-P10;
  a1=P21-P11;
  k1=a0[1]/a0[0]*(P11[0]-P10[0])+P10[1]-P11[1];
  k2=a1[1]-a0[1]/a0[0]*a1[0];
  k=k1/k2;
  return k;
 }

 double tubedRay::crossXAxis (const maths::Vector<double>& P1, const maths::Vector<double>& P2,
   const maths::Vector<double>& k)
 {
  double m=0.0;

  if (k[0]!=0.0) m=-P1[0]/k[0];
  else m=-1000000.0;
  return m;
 }

std::ostream & operator << (std::ostream & os, tubedRay S)
{
 for (int i=0; i<5;i++)
 {
  os << "P[" << i << "]=" << S.P[i] << "   k[" << i << "]=" << S.k[i] << std::endl;
  os << "E[" << i << "]=" << S.E[i] << std::endl;
  os << "Anzahl Einschlsse=" << S.numObj << std::endl;
  os << "imEinschlus=" << S.inObject  << "   objIndex=" << S.objIndex << std::endl;
  os << "OK=" << S.OK << std::endl;
 }
 return os;
}

maths::Vector<double> tubedRay::nextCaustic (double &l)
{
 l=-1;
 if (fabs((k[4]%ka)*P[4])>1E-6) return maths::dzero;
 if (k[4]*ka==1.0) return maths::dzero;

 l=((P[4]*ka)*(ka*k[4])-P[4]*k[4])/(1.0-(ka*k[4])*(ka*k[4]));
 return P[4]+l*k[4];
}

std::ostream& operator << (std::ostream &os, Gauss gs)
{
 os << "Gausstrahl:" << std::endl;
 os << "w0=" << gs.w0 << std::endl;
 os << " F=" << gs.F << std::endl;
 os << " k=" << gs.k << std::endl;
 os << " isGauss=";
 if (gs.isGauss) os << "true" << std::endl;
 else os << "false" << std::endl;
 return os;
}

tubedRay& tubedRay::operator = (const tubedRay& S)
{
  if (this == &S) return *this;
  numObj=S.numObj;
  Obj=S.Obj;
  n0=S.n0;
  for (int i=0; i<5; i++)
  {
   E[i]=S.E[i];
   P[i]=S.P[i];
   k[i]=S.k[i];
   phi[i]=S.phi[i];
  }

  OK=S.OK;
  objIndex=S.objIndex;
  g=S.g;
  getunnelt=S.getunnelt;
  iR=S.iR;
  inObject=S.inObject;
  isValid=S.isValid;
  k0=S.k0;
  ka=S.ka;
  logRay=S.logRay;
  n=S.n;


  pol=S.pol;
  r0=S.r0;
  rc=S.rc;
  return *this;
}

void tubedRay::tunnel(maths::Vector<std::complex<double> > Pol, std::complex<double>  np, std::complex<double>  na, int l)
{
 double s;
 maths::Matrix<double> H,R;
 maths::Matrix<std::complex<double> > T;
 maths::Matrix<double> D;
 maths::Vector <double> n,e0,e1,e2;
 double alpha,gamma;
 std::complex<double>  beta;
 double b;

 // Innerer Reflexionswinkel (von der Tangente an P aus gerechnet)
 gamma = std::acos((l+0.5)/(real(np)*k0*r0));

// cout << "P Anfang:" << P[0] << P[1] << P[2] << P[3] << P[4] << endl;

// cout << "gamma:" << gamma << endl;
 // Rechtsrum oder linksrum tunneln ?

 // noch zu implementieren

 // P liegt auf der Kugeloberflaeche

 for (int i=0;i<5;i++)
 {
 // Transformation in lokales Koordinatensystem
 // P liegt auf der Kugeloberflaeche
 n = P[i]/abs(P[i]);
 getKSystem (n,k[i],e0,e1,e2);
 trafo(e0,e1,e2,H,R);
 // rotMatrix D
 D=rotMatrixA(e2,k[i],gamma);
 // k drehen
// cout << "k:" << k[i] << endl;
 k[i]=D*k[i];
//  cout << "k nach:" << k[i] << endl;

  E[i]=H*E[i];
//  E1=T*E1;
  E[i]=R*E[i];
  E[i]=D*E[i];
 }

//  E2=H*E2;
//  E2=T*E2;
//  E2=R*E2;
//  E2=D*E2;
/*  for (int i=0; i<5; i++)
 {
  E[i]=Pol;
  b=sqrt(P[i][1]*P[i][1]+P[i][2]*P[i][2]);
  n=P[i]/abs(P[i]);
  getKSystem (n,k[i],e0,e1,e2);
  alpha=acos(n*k[i]);
  if (alpha>M_PI/2.0)
  {
   alpha=M_PI-alpha;
   beta=asin(n1/n2*sin(alpha));
   e2=-e2;
  }

  beta=asin(n1/n2*sin(alpha));
  E[i]*=exp(I*n1*k0*P[i][0]);
  P[i][0]=0.0;
  P[i]*=r0/abs(P[i]);
  gamma=-real(acos(b*n1/(r0*n2)));

  s=1.0;
  trafo(e0,e1,e2,H,R);
  D=rotMatrixA(e2,k[4],gamma);
  T=Fresnel_trans(alpha,beta);
  T*=exp(I*k0*cos(beta)*(b-r0)*n1);
  k[i]=D*k[i];
  E[i]=H*E[i];
  E[i]=T*E[i];
  E[i]=R*E[i];
  E[i]=D*E[i];
  P[i]*=(1.0-KORR);
 }*/


}


void tubedRay::tunnel(maths::Vector<std::complex<double> > Pol, std::complex<double>  n1, std::complex<double>  n2)
{
 double s;
 maths::Matrix<double> H,R;
 maths::Matrix<std::complex<double> > T;
 maths::Matrix<double> D;
 maths::Vector <double> n,e0,e1,e2;
 double alpha,gamma;
 double  beta;
 double b;

 for (int i=0; i<5; i++)
 {
  E[i]=Pol;
  b=sqrt(P[i][1]*P[i][1]+P[i][2]*P[i][2]);
  n=P[i]/abs(P[i]);
  getKSystem (n,k[i],e0,e1,e2);
  alpha=std::acos(n*k[i]);
  if (alpha>M_PI/2.0)
  {
   alpha=M_PI-alpha;
   beta=real(asin((std::complex<double>)real(n1)/real(n2)*sin(alpha)));
   e2=-e2;
  }

  beta=real(asin((std::complex<double>)real(n1)/real(n2)*sin(alpha)));
  E[i]*=exp(I*n1*k0*P[i][0]);
  P[i][0]=0.0;
  P[i]*=r0/abs(P[i]);
  gamma=-std::acos(b*real(n1)/(r0*real(n2)));

  s=1.0;
  trafo(e0,e1,e2,H,R);
  D=rotMatrixA(e2,k[4],gamma);

  T=Fresnel_trans(alpha,beta,n1,n2);
  T*=exp(I*k0*cos(beta)*(b-r0)*n1);
  k[i]=D*k[i];
  E[i]=H*E[i];
  E[i]=T*E[i];
  E[i]=R*E[i];
  E[i]=D*E[i];
//  P[i]*=(1.0-KORR);
 }
}

maths::Vector<double> tubedRay::crossPlane (const maths::Vector<double> Pe, const maths::Vector<double> n)
/* Berechnet den Schnittpunkt des zentralen Strahls mit einer Plane, die durch den Aufpunkt Pe und die Normale n beschrieben wird 
 * Die Routine liefert inf-Vektor, falls kein Schnittpunkt */
{
 double l=k[4]*n;
 double h=(Pe-P[4])*n;
 maths::Vector<double> S;

 if ( (l==0) && (h==0) ) return P[4]; // Strahl befindet sich in der Plane und bewegt sich in der Plane
 if ( (l==0) && (h!=0) ) return maths::Vector<double> (-1,-1,-1);
 l=h/l;
 S=P[4]+l*k[4];
 return S;
}

void binWrite (Gauss gs, std::ofstream &os)
{
 os.write((char *) &gs.w0, (char) sizeof(gs.w0));
 gs.F.binWrite(os);
 gs.k.binWrite(os);
 os.write((char *) &gs.isGauss, (char) sizeof (gs.isGauss));
}


void binRead (Gauss &Gs, std::ifstream &is)
{
 Gauss gs; 
 is.read((char *) &gs.w0, (char) sizeof(gs.w0));
 gs.F.binRead(is);
 gs.k.binRead(is);
 is.read((char *) &gs.isGauss, (char) sizeof (gs.isGauss));
 Gs=gs;
}
}
}