#pragma strict_gs_check(on)
#include "surface.h"
#include "misc.h"
#include <stdio.h>
#include <string.h>
#include "box.h"

#define TREE_RECURSIONS 4

//#define TEST_CULL

using namespace std;

/*
surface::surface(int anz)
{
  r0=1.0;
  sf=1.0;
  numTriangles=anz;
  S=new triangle[numTriangles];
  P=dzero;
  currentnorm=Vector<double>(0,0,0);
  currentIndex = -1;
  type=OBJECTSHAPE_SURFACE;
}
*/
surface::surface() 
{
  r0=1.0;
  sf=1.0;
  numTriangles=0;
  S=NULL;
  P = dzero;
  currentnorm=Vector<double>(0,0,0);
  currentIndex = -1;
  type=OBJECTSHAPE_SURFACE;
}
/*
surface::surface (const ObjectShape &F)
{
 r0=1.0;
 sf=1.0;
 type=OBJECTSHAPE_SURFACE;
 S = 0;
}
*/
surface::surface(const surface &Su):ObjectShape(Su)
{ 
    numTriangles=Su.numTriangles;
  sf=Su.sf;
  r0=Su.r0;
  sf=1.0;  
  P = dzero;
  currentnorm=Su.currentnorm;
  currentIndex = Su.currentIndex;
  S=new triangle[numTriangles];
  P=Su.P;
  for (int i=0;i<numTriangles;i++)
  {
   S[i]=Su.S[i];
   // cout << S[i].P[0] << "   " << S[i].P[1] <<  "   " << S[i].P[2] << endl; 
  }

  type=OBJECTSHAPE_SURFACE;  
//  if (FName!=0) delete[] FName;  
    FName="UNBEKANNT";
    
//  cout << "weiter" << endl;
  initQuad(); 
  #ifdef WITH_OCTREE
  Vector<double> pul, por;
  initBounds(pul,por);
  Vector<double> d = por - pul;
  double h = d[0];
  if (d[1] > h) h = d[1];
  if (d[2] > h) h = d[2];
 // cout << "pul=" << pul << "   por=" << por << endl;
/*  cout << "CREATE NEW OCTREE" << endl;
  cout << "d=" << d << endl;
  cout << "numTriangles=" << numTriangles << endl;*/
 /* cube = Octree<triangle>::Section(P, d, 1.0, h);
  Tree.setRootSection(cube);
  Tree.build(S, numTriangles);*/
  Tree.BBox = Box(dzero, d, n);
  Tree.BBox.setOctree(true);
  Tree.createTree(TREE_RECURSIONS);
  for (int i = 0; i < numTriangles; i++)
  {
	  addTriangleToTriangle(Tree, S[i]);
  }
  Tree.trimOctree();
/*  cout << "%TREE Begins ------------------" << endl;
  cout << Tree << endl;
  cout << "%TREE END ---------------------" << endl;*/
#endif 
//     sprintf (FName,"%s",Su.FName);

//   FName="HALLO";  
   
}

surface::surface(Vector<double> Oh) : ObjectShape()
{
  numTriangles=0;
  S=NULL;
  P = Oh;
  r0=1.0;
  type=OBJECTSHAPE_SURFACE;
}


surface::surface(const Vector<double> &P,
                complex<double>  n,
          const Matrix<complex<double> > alpha,
          const Vector<double> &Ex,
          const Vector<double> &Ey,
          const Vector<double> &Ez
         )
: ObjectShape (P, n, alpha, Ex, Ey, Ez, OBJECTSHAPE_SURFACE)
{
 r0=1.0;
  sf=1.0;
  numTriangles=0;
  S=NULL;
  currentnorm=Vector<double>(0,0,0);
  currentIndex = -1;
 trafo (Ex,Ey,Ez,H,R);
 this->P=P;
 type=OBJECTSHAPE_SURFACE;
} 

surface::surface(const Vector<double> &P,
                complex<double>  n,
          int anz, triangle* list,
          const Matrix<complex<double> > alpha,
          const Vector<double> &Ex,
          const Vector<double> &Ey,
          const Vector<double> &Ez)
: ObjectShape (P, n, alpha, Ex, Ey, Ez,OBJECTSHAPE_SURFACE)
{
 trafo (Ex,Ey,Ez,H,R);
 this->P=P;
 this->numTriangles=anz;
 S=new triangle[numTriangles];
 for(int i=0;i<numTriangles;i++)
 this->S[i]=list[i];
 type=OBJECTSHAPE_SURFACE;
 initQuad();
#ifdef WITH_OCTREE
// Vector<double> por, pul;
// initBounds(pul,por); 
 Vector<double> d = por - pul;
 double h = d[0];
 if (d[1] > h) h = d[1];
 if (d[2] > h) h = d[2];
 Tree.createTree(TREE_RECURSIONS);
 Tree.BBox = Box(P, 1.05 * d, n);
 for (int i = 0; i < numTriangles; i++)
	 addTriangleToTriangle(Tree, S[i]);
 // Tree.trimOctree();
 /*cube = Octree<triangle>::Section(P, d, 1.0, h);
 cout << "bounds:" << cube.bounds[0] << "    " << cube.bounds[1] << endl; 
 Tree.setRootSection(cube);
 Tree.build(S, numTriangles);*/
#endif 

}
/*
surface::~surface()
{
 clearS();
}

*/
int surface::createsurface()
{
  double P1x, P1y, P1z, P2x, P2y, P2z, P3x, P3y, P3z;
  Vector<double> P1, P2, P3;
  cout  << "% Anzahl Flaechen:";
  cin >> numTriangles;
  delete[] S;
  S = new triangle[numTriangles];

  for(int i=0;i<numTriangles;i++){
    cout << "% Flaeche" << i+1 << endl;
    cout << "P1x:" << endl;
    cin >> P1x;
    cout << "P1y:" << endl;
    cin >> P1y;
    cout << "P1z:" << endl;
    cin >> P1z;
    cout << "P2x:" << endl;
    cin >> P2x;
    cout << "P2y:" << endl;
    cin >> P2y;
    cout << "P2z:" << endl;
    cin >> P2z;
    cout << "P3x:" << endl;
    cin >> P3x;
    cout << "P3y:" << endl;
    cin >> P3y;
    cout << "P3z:" << endl;
    cin >> P3z;
    P1 = Vector<double>(P1x,P1y,P1z);
    P2 = Vector<double>(P2x,P2y,P2z);
    P3 = Vector<double>(P3x,P3y,P3z);
    S[i]=triangle(P1,P2,P3);
	S[i].setnorm();
  }
  return 0;
}

int surface::createsurface(std::string FName)
{
  Vector<double> P1, P2, P3;

  ifstream is;
//  if (this->FName!=0) delete this->FName;
  int Nc = FName.length() + 1;
 //  this->FName=new char[Nc];
  this->FName = std::string(FName);
  // strcpy (this->FName,FName);
  is.open(FName);
  is >> numTriangles;

  S = new triangle[numTriangles];

  for(int i=0;i<numTriangles;i++)
  {
    is >> P1 >> P2 >> P3;
    S[i] = triangle(P1*r0,P2*r0,P3*r0); 
	S[i].setnorm();
//   cout << P1 << "\t" << P2 << "\t" << P3 << endl; 
   // cout << S[i] << endl;
  //  cout << "r0:" << r0 << ", S[" << i <<"]:" << S[i] << endl;
  }
//  cout << "numTriangles nach:" << numTriangles << endl;

  
 // setCenter2CoM();
  setCenter(P);
initQuad();
#ifdef WITH_OCTREE
        Vector<double> por, pul;
        initBounds(pul,por);

	Vector<double> d = por - pul;
	Vector<double> Ph = (por + pul) / 2.0;
	double h = d[0];
	if (d[1] > h) h = d[1];
	if (d[2] > h) h = d[2];
        Vector<double> hd(h,h,h); 
		Tree.BBox = Box(Ph, 1.05 * d, this->n);
		Tree.createTree();
		
        writeTriangleOctree("octree.log", Tree);
		for (int i = 0; i < numTriangles; i++)
			addTriangleToTriangle(Tree, S[i]);
#endif 
  return 0;
}

int surface::importBinSTL(std::string FName)
{
	ifstream is;
	int dummy;
	int i,j;
	int anz;
	float data;
	char str[255];
	Vector<double> P1,P2,P3,n;
	Vector<double> cm;
	cout << "% ---------------------------- IMPORT STL-FILE --------------------------------" << endl;
        cout << "% Lese:" << FName << endl;
	if (numTriangles!=0) delete[] S;

	is.open (FName,ios::binary);
	if (is.good())
	{
	is.read (str,80);
	anz=readLE_int32(is);
	//if (this->FName!=0) delete[] this->FName;
	//this->FName=new char[strlen(FName)+1];
	// sprintf (this->FName,"%s",FName);
	//strcpy (this->FName,FName);
	this->FName = std::string(FName);
	S=new triangle[anz];
	numTriangles=anz;
	cout << "% Lese " << anz << "  Dreiecke" << endl;

	for (i=0; i<anz && !is.eof(); i++)
	{
		// is.read(str,12);  // Dreiecksnormale --> wird eh neu berechnet 

		for (j=0; j<3; j++)   // Punkt 1
		{  
         n[j]=readLE_float32(is);
		}
		
		for (j=0; j<3; j++)   // Punkt 1
		{  
         P1[j]=readLE_float32(is);
		 cm[j]+=P1[j];
		}

		for (j=0; j<3; j++)   // Punkt 2
		{
		 P2[j]=readLE_float32(is);
		 cm[j]+=P2[j];
		}

		for (j=0; j<3; j++)   // Punkt 3
		{
		 P3[j]=readLE_float32(is);
		 cm[j]+=P3[j];
		}

		cm=cm/anz/3.0;
		is.read(str,2);
		S[i]=triangle(P1,P2,P3);
		S[i].setnorm(n); // Glauben wir mal, dass die Oberfl�chennormale im STL-File richtig ist !
	}
	
	initQuad();
        cout << "% por=" << por << "    pul=" << pul << endl;
	is.close();
/*	setCenter2CoM();
	setCenter(P); */
	P = dzero;
	setCenter(P);
	cout << "% P=" << P << endl;
        initQuad();
#ifdef WITH_OCTREE
        Vector<double> por, pul;
        initBounds(pul,por);
	Vector<double> d = por - pul;
	double h = d[0];
	if (d[1] > h) h = d[1];
	if (d[2] > h) h = d[2];
        Vector<double> hd(h,h,h); 	
		// hd = (pul + por) / 2.0;
//		Tree.BBox = Box(dzero, Vector<double>(20,20,20), this->n);
		Tree.BBox = Box(P, d, this->n);
    	Tree.createTree();
		
		for (int i = 0; i < numTriangles; i++)
			addTriangleToTriangle(Tree, S[i]);
	
#endif 
	std::cout <<  "% ------------------------------- IMPORT ENDE ---------------------------------" << endl;
	return 0;
	}
	else 
        { 
          printf ("loading STL-file failed\n");
          return -1; 
        }
}

bool surface::next(const Vector<double> &r, const Vector<double> &k, Vector<double> &p)
{
// Angenommene Koordinatensysteme:
//
// r     : Globales Koordinatensystem
// S[i].P: Einschlusskoordinatensystem
// rhilf : r im Einschlusskoordinatensystem
// khilf : k im Einschlusskoordinatensystem
// p     : Globales Koordinatensystem
//

#ifndef WITH_OCTREE
  double thilf,t;
  Vector<double> rhilf, philf, khilf;
  bool found=false;
  rhilf = H*(r - P);
  khilf = H*k;
  currentIndex = -1;
  for(int i=0;i<numTriangles;i++)
  {
    if(S[i].calcIntersectionPoint(rhilf,khilf,t,p))
    {
      // found=true;
       if((!found) && (t>0))
       {
        thilf=t; philf=p; currentnorm=S[i].n;
	currentIndex = i;
        found=true;
       }
       else
       {
        if((t<thilf)&&(t>1E-20))
        {
           thilf=t;
           philf=p;
           currentnorm=S[i].n;
	   currentIndex = i;
        }
       }
    }
//    else
//    {
//      cout << i << endl;
//      cout << "nix\n";
//    }
  }
  if (found)
  {
 //  p=R*philf+P;
	  p = r + thilf * k;
 //  cout << "next->currentIndex=" << currentIndex << endl;
   return found;
  }
  else
  { 
//   cout <<"Kein Schnittpunkt!" << endl;
//   p=Vector<double> (-1E20,-1E20,-1E20);
   return found;
  }
#endif
#ifdef WITH_OCTREE
  triangle d;
  Vector<double> rhilf, philf, khilf;
  bool found = false;
  rhilf = H*(r - P);  
  khilf = H*k;

  double erg=rayOctreeIntersection(Tree, rhilf, khilf, d)>0;

  if (erg>0)
  {
	 //  d.calcIntersectionPoint(rhilf, khilf, philf);
	  erg=d.distance(rhilf, khilf);
	  p = r + erg*k;
	  // d.setnorm();
	  currentnorm = R*d.n;
	  // currentnorm = R*d.n;
	  return true;
  }
 /* if (d != nullptr)
  {
	  currentnorm = d->n;
	  // cout << "distance=" << distance << endl;
	  // cout << r << "   " << k << "   " << d->P[0] << "   " << d->P[1] << "   " << d->P[2] << "   " << d->n << endl;
	  // bool  h=d->calcIntersectionPoint(ot.P, ot.k, p);
	   p = r + distance*k;
	  // cout << "p=" << p << endl;
	  return true;
  }*/
  return false;
#endif
}

bool surface::isInside(const Vector<double> &p)
{
  double hilf;
  Vector<double>  hilf2;
  cout << "XXX p:" << p << endl;
  for(int i=0;i<numTriangles;i++)
  {
    hilf2 = S[i].P[2]-p;
    hilf = (S[i].n*hilf2);
    cout << "Dreieck" << i <<": " << S[i] << endl;
    cout << "n:" << S[i].n << endl;
    cout << "hilf2:" << hilf2 << endl;
    cout << "hilf:" << hilf << endl;
    if(hilf<0.0)
      return false;
  }
  return true;
}

void surface::initBounds(Vector<double> &pul, Vector<double> &por)
{
  double xmin,ymin,zmin,xmax,ymax,zmax;
  triangle Sh;

  Sh = S[0];
  xmin=Sh.P[0][0]; ymin = Sh.P[0][1]; zmin=Sh.P[0][2];
  xmax=Sh.P[0][0]; ymax = Sh.P[0][1]; zmax=Sh.P[0][2];

  for(int i=0;i<numTriangles;i++)
  {
   Sh = S[i];
   if(Sh.P[0][0]<xmin)
     xmin = Sh.P[0][0];
   else if(Sh.P[0][0]>xmax)
     xmax = Sh.P[0][0];
   if(Sh.P[0][1]<ymin)
     ymin = Sh.P[0][1];
   else if(Sh.P[0][1]>ymax)
     ymax = Sh.P[0][1];
   if(Sh.P[0][2]<zmin)
     zmin = Sh.P[0][2];
   else if(Sh.P[0][2]>zmax)
     zmax = Sh.P[0][2];

   if(Sh.P[1][0]<xmin)
     xmin = Sh.P[1][0];
   else if(Sh.P[1][0]>xmax)
     xmax = Sh.P[1][0];
   if(Sh.P[1][1]<ymin)
     ymin = Sh.P[1][1];
   else if(Sh.P[1][1]>ymax)
     ymax = Sh.P[1][1];
   if(Sh.P[1][2]<zmin)
     zmin = Sh.P[1][2];
   else if(Sh.P[1][2]>zmax)
     zmax = Sh.P[1][2];

   if(Sh.P[2][0]<xmin)
     xmin = Sh.P[2][0];
   else if(Sh.P[2][0]>xmax)
     xmax = Sh.P[2][0];
   if(Sh.P[2][1]<ymin)
     ymin = Sh.P[2][1];
   else if(Sh.P[2][1]>ymax)
     ymax = Sh.P[2][1];
   if(Sh.P[2][2]<zmin)
     zmin = Sh.P[2][2];
   else if(Sh.P[2][2]>zmax)
     zmax = Sh.P[2][2];
  }

 /* cout << "xmin=" << xmin << "   xmax=" << xmax;
  cout << "  ymin=" << ymin << "   ymax=" << ymax;
  cout << "  zmin=" << zmin << "   zmax=" << zmax << endl; */
  pul=Vector<double>(xmin,ymin,zmin);
  por=Vector<double>(xmax,ymax,zmax);
 // cout << "pul=" << pul << "    por=" << por << endl;
}


void surface::initQuad()
{
  double xmin,ymin,zmin,xmax,ymax,zmax;
  triangle Sh;
  if (numTriangles > 0)
  {
	  Sh = R * S[0] + P;
	  xmin = Sh.P[0][0]; ymin = Sh.P[0][1]; zmin = Sh.P[0][2];
	  xmax = Sh.P[0][0]; ymax = Sh.P[0][1]; zmax = Sh.P[0][2];

	  for (int i = 0; i < numTriangles; i++)
	  {
		  Sh = R * S[i] + P;
		  if (Sh.P[0][0] < xmin)
			  xmin = Sh.P[0][0];
		  else if (Sh.P[0][0] > xmax)
			  xmax = Sh.P[0][0];
		  if (Sh.P[0][1] < ymin)
			  ymin = Sh.P[0][1];
		  else if (Sh.P[0][1] > ymax)
			  ymax = Sh.P[0][1];
		  if (Sh.P[0][2] < zmin)
			  zmin = Sh.P[0][2];
		  else if (Sh.P[0][2] > zmax)
			  zmax = Sh.P[0][2];

		  if (Sh.P[1][0] < xmin)
			  xmin = Sh.P[1][0];
		  else if (Sh.P[1][0] > xmax)
			  xmax = Sh.P[1][0];
		  if (Sh.P[1][1] < ymin)
			  ymin = Sh.P[1][1];
		  else if (Sh.P[1][1] > ymax)
			  ymax = Sh.P[1][1];
		  if (Sh.P[1][2] < zmin)
			  zmin = Sh.P[1][2];
		  else if (Sh.P[1][2] > zmax)
			  zmax = Sh.P[1][2];

		  if (Sh.P[2][0] < xmin)
			  xmin = Sh.P[2][0];
		  else if (Sh.P[2][0] > xmax)
			  xmax = Sh.P[2][0];
		  if (Sh.P[2][1] < ymin)
			  ymin = Sh.P[2][1];
		  else if (Sh.P[2][1] > ymax)
			  ymax = Sh.P[2][1];
		  if (Sh.P[2][2] < zmin)
			  zmin = Sh.P[2][2];
		  else if (Sh.P[2][2] > zmax)
			  zmax = Sh.P[2][2];
	  }

	  /*cout << "xmin=" << xmin << "   xmax=" << xmax;
	  cout << "  ymin=" << ymin << "   ymax=" << ymax;
	  cout << "  zmin=" << zmin << "   zmax=" << zmax << endl; */
	  pul = Vector<double>(xmin, ymin, zmin) + P;
	  por = Vector<double>(xmax, ymax, zmax) + P;
	  // cout << "pul=" << pul << "    por=" << por << endl;
  }
}

Vector<double> surface::norm (const Vector<double> &dummy)
{
 // Dummy Argument ist notwendig aufgrund des Designs der Formklasse
 return currentnorm;
}

void surface::setr0(double rneu)
{
 // cout << "surface::setr0" << endl;
 /* for(int i=0;i<numTriangles;i++)
  { 
    S[i].P[0]=S[i].P[0]*rneu/r0;
    S[i].P[1]=S[i].P[1]*rneu/r0;
    S[i].P[2]=S[i].P[2]*rneu/r0;
  }*/
  r0=rneu;
  initQuad();
}

surface surface::nosurface()
{
  surface h;
  h=*this;
  h.clearS();
  return h; 
}

void surface::clearS()
{
  if (numTriangles!=0)
  {
    delete[] S;
    S = NULL;
    numTriangles=0;
  }
}

void surface::addTriangle(triangle* list,int anz)
{
  if(numTriangles!=0)
   clearS();

  numTriangles=anz;
  
  S = new triangle[numTriangles];
  
  for(int i=0;i<numTriangles;i++)
  {
    S[i]=list[i];
  }

}



ostream& operator << (ostream &os, const surface &su)
{
	os << "Pos=" << su.P << endl;
	os << "Winkel:" << su.Ealpha << "," << su.Ebeta << "," << su.Egamma << endl; 
 os << "anzp:" << su.numTriangles << endl;
 if(su.numTriangles>0)
 {
  os << "[";
  for (int i=0;i<su.numTriangles-1;i++)
    os  << su.S[i] << "\n";
  os << su.S[su.numTriangles-1] << "]";
 }
 else
 {
  os <<"[]";
 } 
 return os;
}


surface& surface::operator = (const surface &s)
{
  this->P=s.P;
  this->H=s.H;
  this->R=s.R;
  this->n=s.n;
  this->alpha=s.alpha;
  this->type=s.type;
  this->pul=s.pul;
  this->por=s.por;
  for(int i=0;i<3;i++)
  this->e[i]=s.e[i];
  this->Ealpha=s.Ealpha;
  this->Ebeta=s.Ebeta;
  this->Egamma=s.Egamma;
  this->r0=s.r0;
  this->numTriangles=s.numTriangles; 
  this->S=new triangle[this->numTriangles];
  for(int i=0;i<this->numTriangles;i++)
  this->S[i]=s.S[i];
  this->sf=s.sf;
  this->currentnorm=s.currentnorm;
  this->currentIndex = s.currentIndex;
  // sprintf (this->FName,"%s",FName);
  //strcpy (this->FName,FName);
  
#ifdef WITH_OCTREE
   Vector<double> d = por - pul;
   double h = d[0];
   if (d[1] > h) h = d[1];
   if (d[2] > h) h = d[2];
    Tree.BBox = Box(dzero, d, n);
   Tree.BBox.setOctree(true);

   Tree.createTree(TREE_RECURSIONS,0);
   
   for (int i = 0; i < numTriangles; i++)
	   addTriangleToTriangle(Tree, S[i]);
  // Tree.trimOctree();
//   cube = Octree<triangle>::Section(P, d, 1.0, h);
//   Tree.setRootSection(cube);
//   Tree.build(S, numTriangles);
#endif 
   initQuad();
  return *this;
}

surface operator + (const surface &s, const Vector<double> &v)
{

 surface h(s);

 for (int i=0;i<s.numTriangles;i++)
 {
   h.S[i]=s.S[i]+v;
 } 
 return h;
}

surface operator - (const surface &s, const Vector<double> &v)
{

 surface h(s);
 
 for (int i=0;i<s.numTriangles;i++)
 {
   h.S[i]=s.S[i]-v;
 } 
 return h;
}

surface operator * (const Matrix<double> &M, const surface &s)
{
 surface h(s);

 for (int i=0;i<s.numTriangles;i++)
 {
   h.S[i]=M*s.S[i];
 }
 return h;
}
/** No descriptions */
void surface::scale (double sf)
{
 for (int i=0; i<numTriangles; i++)
  S[i]=S[i]*sf/this->sf;
 this->sf=sf;
}
/** No descriptions */
/*double surface::isInHost(void)
{
 double D,Pk,l;
 double erg=1;
 double r=0;
 Vector<double> hv,k,ergv;
 if (abs(P)>1.0) return -1;
 for (int i=0; i<numTriangles; i++)
 {
  for (int j=0; j<3; j++)
  {
   hv=R*S[i][j]+P;
   r=abs(hv);
   if ((r>erg) && (r>1.0)) { erg=r; ergv=hv; }
  }
 }
 if (erg>1.0)
 {
  k=hv-P;
  Pk=P*k;
  D=Pk*Pk-abs2(k)*(abs2(P)-1.0);
  l=-Pk+sqrt(D)/(abs2(k));
  return l/abs(ergv);
 }
 return 1.0/erg;
}
*/
void surface::binWrite (ofstream &os)
{
 size_t strl;
 P.binWrite(os);
 H.binWrite(os);
 R.binWrite(os);
 os.write ((char *) &n, (char) sizeof (n)); 
 alpha.binWrite(os);
 pul.binWrite (os);
 por.binWrite (os);
 for (int i=0; i<3; i++)
  e[i].binWrite(os);
 os.write ((char *) &Ealpha,(char) sizeof (Ealpha));
 os.write ((char *) &Ebeta,(char) sizeof (Ebeta));
 os.write ((char *) &Egamma,(char) sizeof (Egamma));
 os.write ((char *) &sf,(char) sizeof (sf));
 os.write ((char *) &r0,(char) sizeof (r0));
 os.write ((char *) &numTriangles, (char) sizeof(numTriangles));
 strl=FName.length();
 os.write ((char *) &strl, sizeof(strl));
 char c;
 for (int i=0; i<=FName.length(); i++)
 {
  c=FName[i];
  cout << "c=" << c << endl;
 os.write ((char *) &c, 1);
 }
 for (int i=0; i<numTriangles; i++)
  S[i].binWrite(os);  
}

void surface::binRead (ifstream &is)
{
 type=OBJECTSHAPE_SURFACE;
 P.binRead(is);
 H.binRead(is);
 R.binRead(is);
 is.read((char *) &n, (char) sizeof(n)); 
 alpha.binRead(is);
 pul.binRead (is);
 por.binRead (is);
 for (int i=0; i<3; i++)
  e[i].binRead(is);
 is.read ((char *) &Ealpha,(char) sizeof (Ealpha));
 is.read ((char *) &Ebeta,(char) sizeof (Ebeta));
 is.read ((char *) &Egamma,(char) sizeof (Egamma));
 is.read ((char *) &sf,(char) sizeof (sf));
 is.read ((char *) &r0,(char) sizeof (r0));
 clearS();
 is.read ((char *) &numTriangles,(char) sizeof (numTriangles));
 size_t strl;
 is.read ((char *) &strl, (char) sizeof(strl));
 char c;
 for (int i=0; i<=strl; i++)
 {
  is.read ((char *) &c, 1);
  FName=FName + c;
 } 
//  cout << "FName=" << FName << endl; 
 S=new triangle[numTriangles];
 for (int i=0; i<numTriangles; i++)
 {
  S[i].binRead(is);
  S[i].setnorm();
 } 
}


void surface::exportSRF (std::string FName)
{
 ofstream os;
 os.open (FName);
 os << numTriangles << endl;
 for (int i=0; i<numTriangles; i++)
  os << S[i] << endl;
 os.close();
}

/*!
    \fn surface::volume()
 */
// double surface::volume() {return volume(); }

double surface::volume()
{
 /* double F,Vg,absn,h;
  Vector<double> n;
  
  if (numTriangles==0) return 0.0; 
  Vg=0.0;
  for (int i=0; i<numTriangles; i++)
  {
   n=(S[i][1]-S[i][0]) % (S[i][2]-S[i][0]);
   absn=abs(n);
   n/=absn;
   F=absn/2.0;
   h=fabs((S[i][0]-P)*n);   
   Vg+=F*h/3.0;
//    cout << "h=" << h << endl;
  }
  return Vg;*/
  double Vg=0;
  double A;
  Vector<double> vsum;
  
  for (int i=0; i<numTriangles; i++)
  {
    A=S[i].area();
	vsum=S[i].P[0]+S[i].P[1]+S[i].P[2];
	Vg+=A/3.0*S[i].getnorm()*vsum;	
  }
  return Vg/3.0;
}
void surface::setP (Vector<double> r)
{
 //Vector<double> dP=P-r;
 P=r;
 /*for (int i=0; i<numTriangles; i++)
   for (int j=0; j<3; j++)
  S[i].P[j]+=dP;*/ 
 initQuad();
}


Vector<double> surface::calcCoM()
{
	Vector<double> n,Fx,Fy,Fz;
	double int_x=0,int_y=0,int_z=0;
	double F,x,y,z,V=volume(),Fges=0;
	for (int i=0; i<numTriangles; i++)
	{
		F=S[i].area();
		Fges+=F;
		x=S[i][0][0];
		y=S[i][0][1];
		z=S[i][0][2];
		n=S[i].getnorm();
		Fx=Vector<double> (0.5*x*x,0,0);
		Fy=Vector<double> (0,0.5*y*y,0);
		Fz=Vector<double> (0,0,0.5*z*z);
		int_x+=F*Fx*n;
		int_y+=F*Fy*n;
		int_z+=F*Fz*n;
//                cout << "i=" << i << "   F=" << F << "   Fx=" << Fx << "   Fy=" << Fy << "   Fz=" << Fz << "   n=" << n << endl;
	}
	cout << "% Gesamtfl�che :" << Fges << endl;
	Vector<double> P=Vector<double>(int_x,int_y,int_z)/V;
	return P;
}

void surface::setCenter(Vector<double> P)
{
	cout << "% P= " << P << "     anzp=" << numTriangles << endl;
	for (int i=0; i<numTriangles; i++)
	{
		S[i].P[0]=S[i].P[0]-P;
		S[i].P[1]=S[i].P[1]-P;
		S[i].P[2]=S[i].P[2]-P;
	}
}

inline void Subexpressions (double w0,double w1, double w2, double &f1, double &f2, double &f3, double &g0, double &g1, double &g2)
{
	double temp0=w0+w1; f1=temp0+w2;
	double temp1=w0*w0;
	double temp2=temp1 + w1*temp0;
	f2=temp2+w2*f1; f3=w0*temp1+w1*temp2+w2*f2;
	g0=f2+w0*(f1+w0); g1=f2+w1*(f1+w1); g2=f2+w2*(f1+w2);
}

Matrix<double> surface::computeInertia()
{
	double f1x,f2x,f3x,g0x,g1x,g2x;
	double f1y,f2y,f3y,g0y,g1y,g2y;
	double f1z,f2z,f3z,g0z,g1z,g2z;
	double x0,y0,z0;
	double x1,y1,z1;
	double x2,y2,z2;
	Vector<double> d;
	Matrix <double> inertia;
		
	const double mult[10]={1.0/6.0,1.0/24.0,1.0/24.0,1.0/24.0,1.0/60.0,1.0/60.0,1.0/60.0,1.0/120.0,1.0/120.0,1.0/120.0};
	double intg[10]={0,0,0,0,0,0,0,0,0,0};
	for (int t=0; t<numTriangles; t++)
	{
		x0=S[t].P[0][0]; y0=S[t].P[0][1]; z0=S[t].P[0][2];
		x1=S[t].P[1][0]; y1=S[t].P[1][1]; z1=S[t].P[1][2];
		x2=S[t].P[2][0]; y2=S[t].P[2][1]; z2=S[t].P[2][2];

		d=(S[t].P[1]-S[t].P[0])%(S[t].P[2]-S[t].P[0]);
//		d=d*(S[t].n*d); // Sicherstellen, dass d auch nach au�en zeigt !
		
		Subexpressions(x0,x1,x2,f1x,f2x,f3x,g0x,g1x,g2x);
		Subexpressions(y0,y1,y2,f1y,f2y,f3y,g0y,g1y,g2y);
		Subexpressions(z0,z1,z2,f1z,f2z,f3z,g0z,g1z,g2z);

		intg[0]+=d[0]*f1x;
		intg[1]+=d[0]*f2x; intg[2]+=d[1]*f2y; intg[3]+=d[2]*f2z;
		intg[4]+=d[0]*f3x; intg[5]+=d[1]*f3y; intg[6]+=d[2]*f3z;
		intg[7]+=d[0]*(y0*g0x+y1*g1x+y2*g2x);
		intg[8]+=d[1]*(z0*g0y+z1*g1y+z2*g2y);
		intg[9]+=d[2]*(x0*g0z+x1*g1z+x2*g2z);
	}

	for (int i=0; i<10; i++) intg[i]*=mult[i];

	double mass=intg[0];

	// center of mass 
	Vector<double> cm=Vector<double> (intg[1],intg[2],intg[3])/mass;
	cout << "% CoM=" << cm << endl;
	cout << "% mass=" << mass << endl;
	inertia(0,0)=intg[5]+intg[6]-mass*(cm[1]*cm[1]+cm[2]*cm[2]);
	inertia(1,1)=intg[4]+intg[6]-mass*(cm[2]*cm[2]+cm[0]*cm[0]);
	inertia(2,2)=intg[4]+intg[5]-mass*(cm[0]*cm[0]+cm[1]*cm[1]);
	inertia(0,1)=-(intg[7]-mass*cm[0]*cm[1]);
	inertia(1,2)=-(intg[8]-mass*cm[0]*cm[2]);
	inertia(0,2)=-(intg[9]-mass*cm[2]*cm[0]);

	return inertia/mass;  
}

void surface::setCenter2CoM()
{
	Vector<double> CoM=calcCoM();
	setCenter(CoM);
}

void surface::getMinMax(double &min, double &max)
{
	double l[3];
	min = -1;
	max = 0;
	if (numTriangles > 0)
	{
		for (int i = 0; i < numTriangles; i++)
		{
			l[0] = abs(S[i].P[0] - S[i].P[1]);
			l[1] = abs(S[i].P[0] - S[i].P[2]);
			l[2] = abs(S[i].P[1] - S[i].P[2]);
			if ((min > l[0]) || (i == 0)) min = l[0];
			if (min > l[1]) min = l[1];
			if (min > l[2]) min = l[2];

			if (max < l[0]) max = l[0];
			if (max < l[1]) max = l[1];
			if (max < l[2]) max = l[2];
		}
	}
}

void getMinMax(int numTriangles, triangle *S, double &min, double &max)
{
	double l[3];
	min = -1;
	max = 0;
	if (numTriangles > 0)
	{
		for (int i = 0; i < numTriangles; i++)
		{
			l[0] = abs(S[i].P[0] - S[i].P[1]);
			l[1] = abs(S[i].P[0] - S[i].P[2]);
			l[2] = abs(S[i].P[1] - S[i].P[2]);
			if ((min > l[0]) || (i == 0)) min = l[0];
			if (min > l[1]) min = l[1];
			if (min > l[2]) min = l[2];

			if (max < l[0]) max = l[0];
			if (max < l[1]) max = l[1];
			if (max < l[2]) max = l[2];
		}
	}
}

surface generatePill (double a, double b, double h, int N, double r0, Matrix<double> M)
{ 
 if (N % 2!=0) N=N+1;
 triangle *D=new triangle [2*N*(N-1)];
 double **x=new double *[N]; for (int l=0;l<N;l++) x[l]=new double  [N]; 
 double **y=new double *[N]; for (int l=0;l<N;l++) y[l]=new double  [N]; 
 double **z=new double *[N]; for (int l=0;l<N;l++) z[l]=new double  [N]; 
 double dr,r;
 double dphi,phi;
 Vector<double> P[3];

 
 dr=2.0*b/((double)N-1.0);
 dphi=2.0*M_PI/((double)N-1.0);
 double sgnr;
 int lr,lphi;
 double zD;
 surface S;



 for (lr=0;lr<N;lr++)
 {
  r=lr*dr-b;
  for (lphi=0;lphi<N;lphi++)
  {
   phi=lphi*dphi;
   x[lr][lphi]=r*cos(phi);
   y[lr][lphi]=r*sin(phi);
  
   sgnr = r<0 ? -1 : 1;

   zD=a*a*b*b-r*r*a*a;
   if (zD<0) zD=0;
   z[lr][lphi]=sgnr*(sqrt(zD)/b+h/2.0);  
  }
 }
 
 int c=0;
 for (lr=0; lr<N/2-1; lr++)
  for (lphi=0; lphi<N-1; lphi++)
 {
  P[0]=M*Vector<double> (x[lr][lphi],y[lr][lphi],z[lr][lphi]);
  P[1]=M*Vector<double> (x[lr+1][lphi],y[lr+1][lphi],z[lr+1][lphi]);
  P[2]=M*Vector<double> (x[lr+1][lphi+1],y[lr+1][lphi+1],z[lr+1][lphi+1]);
  D[c]=triangle (P[0],P[1],P[2]);
  c++;
 
  P[0]=M*Vector<double> (x[lr][lphi],y[lr][lphi],z[lr][lphi]);
  P[1]=M*Vector<double> (x[lr+1][lphi+1],y[lr+1][lphi+1],z[lr+1][lphi+1]);
  P[2]=M*Vector<double> (x[lr][lphi+1],y[lr][lphi+1],z[lr][lphi+1]);
  D[c]=triangle (P[0],P[1],P[2]);
    c++;
 }

 for (lr=0; lr<N/2-1; lr++)
  for (lphi=0; lphi<N-1; lphi++)
 {
  P[0]=M*Vector<double> (x[lr][lphi],y[lr][lphi],-z[lr][lphi]);
  P[1]=M*Vector<double> (x[lr+1][lphi],y[lr+1][lphi],-z[lr+1][lphi]);
  P[2]=M*Vector<double> (x[lr+1][lphi+1],y[lr+1][lphi+1],-z[lr+1][lphi+1]);
  D[c]=triangle (P[0],P[1],P[2]);
  c++;
 
  P[0]=M*Vector<double> (x[lr][lphi],y[lr][lphi],-z[lr][lphi]);
  P[1]=M*Vector<double> (x[lr+1][lphi+1],y[lr+1][lphi+1],-z[lr+1][lphi+1]);
  P[2]=M*Vector<double> (x[lr][lphi+1],y[lr][lphi+1],-z[lr][lphi+1]);
  D[c]=triangle (P[0],P[1],P[2]);
  c++;
 }



 for (lphi=0; lphi<N-1; lphi++)
 {
  P[0]=M*Vector<double> (0,0,-a);
  P[1]=M*Vector<double> (x[N/2-1][lphi],y[N/2-1][lphi],z[N/2-1][lphi+1]);
  P[2]=M*Vector<double> (x[N/2-1][lphi+1],y[N/2-1][lphi+1],z[N/2-1][lphi+1]);
  D[c]=triangle (P[0],P[1],P[2]);
  c++;

  P[0]=M*Vector<double> (0,0,a);
  P[1]=M*Vector<double> (x[N/2-1][lphi],y[N/2-1][lphi],z[N/2][lphi+1]);
  P[2]=M*Vector<double> (x[N/2-1][lphi+1],y[N/2-1][lphi+1],z[N/2][lphi+1]);
  D[c]=triangle (P[0],P[1],P[2]);
  c++;
 }

 lr=N-1;
 for (lphi=0; lphi<N-1; lphi++)
 {
  P[0]=M*Vector<double> (x[lr][lphi],y[lr][lphi],h/2.0);
  P[1]=M*Vector<double> (x[lr][lphi+1],y[lr][lphi+1],h/2.0);
  P[2]=M*Vector<double> (x[lr][lphi],y[lr][lphi],-h/2.0);
  D[c]=triangle (P[0],P[1],P[2]);
  c++;
 
  P[0]=M*Vector<double> (x[lr][lphi+1],y[lr][lphi+1],h/2.0);
  P[1]=M*Vector<double> (x[lr][lphi+1],y[lr][lphi+1],-h/2.0);
  P[2]=M*Vector<double> (x[lr][lphi],y[lr][lphi],-h/2.0);
  D[c]=triangle (P[0],P[1],P[2]);
    c++;
 }
  
       S=surface(dzero, 1.5, c, D);
       S.r0;
  
	return S;

}

surface generateEllipsoid (double a, double b, int N, double r0, Matrix<double> M)
{ 
 if (N % 2!=0) N=N+1;
 triangle *D=new triangle [2*N*(N-1)];
 double **x=new double *[N]; for (int l=0;l<N;l++) x[l]=new double  [N]; 
 double **y=new double *[N]; for (int l=0;l<N;l++) y[l]=new double  [N]; 
 double **z=new double *[N]; for (int l=0;l<N;l++) z[l]=new double  [N]; 
 double dr,r;
 double dphi,phi;
 Vector<double> P[3];
 
 dr=2.0*b/((double)N-1.0);
 dphi=2.0*M_PI/((double)N-1.0);
 double sgnr;
 int lr,lphi;
 double zD;
 surface S;

 for (lr=0;lr<N;lr++)
 {
  r=lr*dr-b;
  for (lphi=0;lphi<N;lphi++)
  {
   phi=lphi*dphi;
   x[lr][lphi]=r*cos(phi);
   y[lr][lphi]=r*sin(phi);
  
   sgnr = r<0 ? -1 : 1;

   zD=a*a*b*b-r*r*a*a;
   if (zD<0) zD=0;
   z[lr][lphi]=sgnr*sqrt(zD)/b;  
  }
 }
 
 int c=0;
 for (lr=0; lr<N/2-1; lr++)
  for (lphi=0; lphi<N-1; lphi++)
 {
  P[0]=M*Vector<double> (x[lr][lphi],y[lr][lphi],z[lr][lphi]);
  P[1]=M*Vector<double> (x[lr+1][lphi],y[lr+1][lphi],z[lr+1][lphi]);
  P[2]=M*Vector<double> (x[lr+1][lphi+1],y[lr+1][lphi+1],z[lr+1][lphi+1]);
  D[c]=triangle (P[0],P[1],P[2]);
  c++;
 
  P[0]=M*Vector<double> (x[lr][lphi],y[lr][lphi],z[lr][lphi]);
  P[1]=M*Vector<double> (x[lr+1][lphi+1],y[lr+1][lphi+1],z[lr+1][lphi+1]);
  P[2]=M*Vector<double> (x[lr][lphi+1],y[lr][lphi+1],z[lr][lphi+1]);
  D[c]=triangle (P[0],P[1],P[2]);
    c++;
 }

 for (lr=0; lr<N/2-1; lr++)
  for (lphi=0; lphi<N-1; lphi++)
 {
  P[0]=M*Vector<double> (x[lr][lphi],y[lr][lphi],-z[lr][lphi]);
  P[1]=M*Vector<double> (x[lr+1][lphi],y[lr+1][lphi],-z[lr+1][lphi]);
  P[2]=M*Vector<double> (x[lr+1][lphi+1],y[lr+1][lphi+1],-z[lr+1][lphi+1]);
  D[c]=triangle (P[0],P[1],P[2]);
  c++;
 
  P[0]=M*Vector<double> (x[lr][lphi],y[lr][lphi],-z[lr][lphi]);
  P[1]=M*Vector<double> (x[lr+1][lphi+1],y[lr+1][lphi+1],-z[lr+1][lphi+1]);
  P[2]=M*Vector<double> (x[lr][lphi+1],y[lr][lphi+1],-z[lr][lphi+1]);
  D[c]=triangle (P[0],P[1],P[2]);
  c++;
 }



 for (lphi=0; lphi<N-1; lphi++)
 {
  P[0]=M*Vector<double> (0,0,-a);
  P[1]=M*Vector<double> (x[N/2-1][lphi],y[N/2-1][lphi],z[N/2-1][lphi+1]);
  P[2]=M*Vector<double> (x[N/2-1][lphi+1],y[N/2-1][lphi+1],z[N/2-1][lphi+1]);
  D[c]=triangle (P[0],P[1],P[2]);
  c++;

  P[0]=M*Vector<double> (0,0,a);
  P[1]=M*Vector<double> (x[N/2-1][lphi],y[N/2-1][lphi],z[N/2][lphi+1]);
  P[2]=M*Vector<double> (x[N/2-1][lphi+1],y[N/2-1][lphi+1],z[N/2][lphi+1]);
  D[c]=triangle (P[0],P[1],P[2]);
  c++;
 }
  
       S=surface(dzero, 1.5, c, D);
       S.r0;
  
	return S;

}


surface generateHexagonCylinder(double a, double h, Matrix<double> M)
{
	int i;
	int N = 24; // 24 Dreiecke
	double s30 = 0.5; // sin 30�
	double c30 = SQRT3 / 2.0;
	triangle *D = new triangle[N];
	Vector<double> P[7],P1,P2,P3;
	Vector<double> hp = 0.5*h*ez;
	Vector<double> hm = -hp;
        surface S;
//		sprintf(S.FName, "hexcyl");

	P[0] = Vector<double>(0, a, 0);
	P[1] = Vector<double>(c30*a, s30*a, 0);
	P[2] = Vector<double>(c30*a, -s30*a, 0);
	P[3] = Vector<double>(0, -a, 0);
	P[4] = Vector<double>(-c30*a, -s30*a, 0);
	P[5] = Vector<double>(-c30*a, s30*a, 0);
	P[6] = P[0];
	for (i = 0; i < 6; i++)
	{
		D[i] = triangle(M*(P[i] + hp), M*(P[i + 1] + hp), M*hp);  // oberer Deckel
                D[i].n = M*ez; 
		D[i + 6] = triangle(M*(P[i] + hm),M*( P[i + 1] + hm), M*hm); // unterer Deckel
                D[i+6].n =-M*ez;
                
                P1=M*(P[i] + hp); 
                P2=M*(P[i] + hm); 
                P3=M*(P[i + 1] + hp);
		D[i + 12] = triangle(P1, P2, P3 );
                D[i + 12].n=(P1-P2)%(P3-P2);
                D[i + 12].n/=abs(D[i + 12].n);

                P1=M*(P[i + 1] + hp);
                P2=M*(P[i + 1] + hm); 
                P3=M*(P[i] + hm);
		D[i + 18] = triangle(P1,P2,P3);
                D[i + 18].n=-(P1-P2)%(P3-P2);
                D[i + 18].n/=abs(D[i + 18].n);
 
	}
        S=surface(dzero, 1.5, N, D);
  
	return S;
}

surface generateHollowHexagonalCylinder (double a, double h, double t, Matrix<double> M)
{
	int i;
	int N = 24; // 24 Dreiecke
	double s30 = 0.5; // sin 30�
	double c30 = SQRT3 / 2.0;
        surface S;
	triangle *D = new triangle[N];
	Vector<double> P[7],P1,P2,P3;
	Vector<double> hp = 0.5*h*ez;
	Vector<double> hm = -hp;

	P[0] = Vector<double>(0, a, 0);
	P[1] = Vector<double>(c30*a, s30*a, 0);
	P[2] = Vector<double>(c30*a, -s30*a, 0);
	P[3] = Vector<double>(0, -a, 0);
	P[4] = Vector<double>(-c30*a, -s30*a, 0);
	P[5] = Vector<double>(-c30*a, s30*a, 0);
	P[6] = P[0];
	for (i = 0; i < 6; i++)
	{
		D[i] = triangle(M*(P[i] + hp), M*(P[i + 1] + hp), M*(hp-t*ez));  // oberer Deckel
                D[i].n = M*ez; 
		D[i + 6] = triangle(M*(P[i] + hm),M*( P[i + 1] + hm), M*(hm+t*ez)); // unterer Deckel
                D[i+6].n =-M*ez;
                
                P1=M*(P[i] + hp); 
                P2=M*(P[i] + hm); 
                P3=M*(P[i + 1] + hp);
		D[i + 12] = triangle(P1, P2, P3 );
                D[i + 12].n=(P1-P2)%(P3-P2);
                D[i + 12].n/=abs(D[i + 12].n);

                P1=M*(P[i + 1] + hp);
                P2=M*(P[i + 1] + hm); 
                P3=M*(P[i] + hm);
		D[i + 18] = triangle(P1,P2,P3);
                D[i + 18].n=-(P1-P2)%(P3-P2);
                D[i + 18].n/=abs(D[i + 18].n);
 
	}
        S=surface(dzero, 1.5, N, D);
        S.FName="UNBEKANNT";

	return S;
}	

surface& generateBullet (double a, double h, double t, Matrix<double> M)
{
	surface obj;
	int i;
	int N = 24; // 24 Dreiecke
	double s30 = 0.5; // sin 30�
	double c30 = SQRT3 / 2.0;
	triangle *D = new triangle[N];
	Vector<double> P[7],P1,P2,P3; 
        double hs=h-t;
	Vector<double> hp = 0.5*hs*ez;
	Vector<double> hm = -hp;
        double dphi=60.0/180.0*M_PI;
        double phi,phip;  
        
        
	for (i = 0; i < 6; i++)
	{
                phi=i*dphi;
                phip=phi+dphi;
                 
                P1=Vector<double> (cos(phi) * a,sin(phi) * a , h/2.0);
                P2=Vector<double> (cos(phip) * a,sin(phip) * a , h/2.0);
                P3=Vector<double> (0,0 , h/2.0+t);
  
                D[i] = triangle(M*P1,M*P2,M*P3);
                D[i].setnorm();

 
                P1=Vector<double> (cos(phi) * a,sin(phi) * a , h/2.0);
                P2=Vector<double> (cos(phip) * a,sin(phip) * a , h/2.0);
                P3=Vector<double> (cos(phip) *a, sin(phip) * a, -h/2.0);
  
                D[i+6] = triangle(M*P1,M*P2,M*P3);
                D[i+6].setnorm();

           
                P1=Vector<double> (cos(phi) * a,sin(phi) * a , h/2.0);
                P2=Vector<double> (cos(phip) * a,sin(phip) * a , -h/2.0);
                P3=Vector<double> (cos(phi) *a, sin(phi) * a, -h/2.0);
  
                D[i+12] = triangle(M*P1,M*P2,M*P3);
                D[i+12].setnorm();

  
                P1=Vector<double> (cos(phi) * a,sin(phi) * a , -h/2.0);
                P2=Vector<double> (cos(phip) * a,sin(phip) * a , -h/2.0);
                P3=Vector<double> (0,0 , -h/2.0);
  
                D[i+18] = triangle(M*P1,M*P2,M*P3);
                D[i+18].setnorm();


  
/*		D[i] = triangle(M*(P[i] + hp), M*(P[i + 1] + hp), M*(hp+t*ez));  // oberer Deckel
                D[i].n = M*ez; 
		D[i + 6] = triangle(M*(P[i] + hm),M*( P[i + 1] + hm), M*hm); // unterer Deckel
                D[i+6].n =-M*ez;
                
                P1=M*(P[i] + hp); 
                P2=M*(P[i] + hm); 
                P3=M*(P[i + 1] + hp);
		D[i + 12] = triangle(P1, P2, P3 );
                D[i + 12].n=(P1-P2)%(P3-P2);
                D[i + 12].n/=abs(D[i + 12].n);

                P1=M*(P[i + 1] + hp);
                P2=M*(P[i + 1] + hm); 
                P3=M*(P[i] + hm);
		D[i + 18] = triangle(P1,P2,P3);
                D[i + 18].n=-(P1-P2)%(P3-P2);
                D[i + 18].n/=abs(D[i + 18].n);
*/
 
	}
        surface *S=new surface(dzero, 1.5, N, D);
        S->FName = "UNBEKANNT"; 
	return *S;
}


surface generatePolyBullet(int n, double a, double h, double t, Matrix<double> M)
{
	//a ist Seitenlänge
	//Koordinatenursprung ist in dem Zentrum des Zylinders vom Körper
	int i;
	int N = 4 * n; // Anzahl Dreiecke für Oberfläche

	double alpha = M_PI / (double)n; //zwei Winkel vom Dreieck im Querschnitt
	double dbeta = 2.0 * M_PI/(double)n; //Winkel vom Zentrum der Grundfläche aus
	double beta, betap;
	double r;
	r = a / (2.0 * sin(alpha));

	triangle* D = new triangle[N];
	Vector<double> P[7], P1, P2, P3;
	double hs = h/2.0 + t; // Abstand Koordinatenursprung bis Spitze
	Vector<double> hp =  hs * ez; // " als Vektor oben
	Vector<double> hm = -hp; // " als Vektor unten
	//double dphi = 60.0 / 180.0 * M_PI;
	
	for (i = 0; i < n; i++)
	{
		beta = i * dbeta;
		betap = beta + dbeta;
		//Mantel des Zylinders mit n-eckiger Grundfläche
		P1 = Vector<double>(cos(beta)*r , sin(beta)*r, h/ 2.0);
		P2 = Vector<double>(cos(beta) * r, sin(beta) * r, -h/ 2.0);
		P3 = Vector<double>(cos(betap)*r, sin(betap)*r, h/2.0); 
		D[i*4] = triangle(M * P1, M * P2, M * P3);
		D[i*4].setnorm();

 		P1 = Vector<double>(cos(beta)*r , sin(beta)*r, -h / 2.0);
		P2 = Vector<double>(cos(betap) * r, sin(betap) * r, -h/ 2.0);
	        P3 = Vector<double>(cos(betap) * r, sin(betap) * r,  h / 2.0);

		D[i*4+1] = triangle(M * P1, M * P2, M * P3);
		D[i*4+1].setnorm();

		//Deckel oben
		P1 = Vector<double>(cos(beta) * r, sin(beta) * r, h / 2.0);
		P2 = Vector<double>(0, 0, hs);
		P3 = Vector<double>(cos(betap) * r, sin(betap) * r, h / 2.0);

		D[i*4+2] = triangle(M * P1, M * P2, M * P3);
		D[i*4+2].setnorm();
		
//Deckel unten... moment, die Oberfläche für beide Deckel ist doch exakt dieselbe, warum berechnen wir die noch mal neu?
		P1 = Vector<double>(cos(beta) * r, sin(beta) * r, -h / 2.0);
		P2 = Vector<double>(0, 0, -h/2.0);
		P3 = Vector<double>(cos(betap) * r, sin(betap) * r, -h / 2.0);


		D[i*4+3] = triangle(M * P1, M * P2, M * P3);
		D[i*4+3].setnorm();
	}
        surface S=surface(dzero, 1.5, N, D);
        S.FName = "UNBEKANNT"; 
	return S;
}
	

/*#ifdef WITH_OCTREE
void surface::initBounds()
{
	Vector<double> Pmin, Pmax;
	int i, j, l;
	for (i = 0; i < numTriangles; i++)  // Dreiecke
	{
		for (j = 0; j < 3; j++)   // Punkte im Dreieck
		{
			for (l = 0; l < 3; l++)    // Koordinaten
			{
				if (S[i].P[j][l] > Pmax[l]) Pmax[l] = S[i].P[j][l];
				if (S[i].P[j][l] < Pmin[l]) Pmin[l] = S[i].P[j][l];
			}
		}
	}
	Tree.BBox.bounds[0] = Pmin;
	Tree.BBox.bounds[1] = Pmax;
}
#endif*/

