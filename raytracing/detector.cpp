#include "detector.h"
#include "matrix.h"




Detector::Detector(void)
{
  n1=0;
  n2=0;
  D=0;
}

Detector::Detector(int n1, int n2)
{
 init(n1,n2);
}

void Detector::init(int n1,int n2)
{
  D=new Vector<std::complex<double> > *[n1];
 for (int i=0; i<n1; i++)
	 D[i]=new Vector<std::complex<double> >[n2];
 this->n1=n1;
 this->n2=n2;
}

Detector::Detector (const Detector& Det)
{
 if (Det.n1>0)
 {
	 init(Det.n1,Det.n2);
     for (int i1=0; i1<n1; i1++)
		 for (int i2=0; i2<n2; i2++)
			 this->D[i1][i2]=Det.D[i1][i2];
 }
 n=Det.n;
 e1=Det.e1;
 e2=Det.e2;
 P=Det.P;
}

Detector& Detector::operator = (const Detector& Det)
{
 if (Det.n1>0)
 {
	 init(Det.n1,Det.n2);
     for (int i1=0; i1<n1; i1++)
		 for (int i2=0; i2<n2; i2++)
			 this->D[i1][i2]=Det.D[i1][i2];
 }
n=Det.n;
 e1=Det.e1;
 e2=Det.e2;
 P=Det.P;
 return *this;
}


Detector::~Detector(void)
{
 clear();
}


void Detector::clean()
{
	if (n1 > 0)
		for (int i1 = 0; i1 < n1; i1++)
			for (int i2 = 0; i2 < n2; i2++)
				D[i1][i2] = czero;
}

void Detector::clear()
{
	 if (n1>0) 
  {
	  for (int i1=0; i1<n1; i1++)
		  delete[] D[i1];
      delete[] D;
	  D=0;
	  n1=0;
	  n2=0;
  }
}

int Detector::N1() { return n1; }
int Detector::N2() { return n2; }

void Detector::save(char* fn)
{
	std::ofstream os;
	os.open(fn);
	for (int i1 = 0; i1 < n1; i1++)
	{
		for (int i2 = 0; i2 < n2; i2++)
		{
	
			os << D[i1][i2] << "   ";
		}
		os << std::endl;
	}
	os.close();
}

void Detector::saveabs (char *fn)
{
  std::ofstream os;
  double h;
  os.open (fn);
   for (int i1=0; i1<n1; i1++)
   {
    for (int i2=0; i2<n2; i2++)
	{
     if (type>=3) h=abs(D[i1][i2]);
	 else
	  h=abs(D[i1][i2]);
	 os << h << "   ";
	}
	os << std::endl;
   }
  os.close();
}


void Detector::savereal (char *fn, int type)
{
  std::ofstream os;
  double h;
  os.open (fn);
   for (int i1=0; i1<n1; i1++)
   {
    for (int i2=0; i2<n2; i2++)
	{
     if (type>=3) h=abs(D[i1][i2]);
	 else
	  h=real(D[i1][i2][type]);
	 os << h << "   ";
	}
	os << std::endl;
   }
  os.close();
}


void Detector::saveimag (char *fn, int type)
{
	std::ofstream os;
  double h;
  os.open (fn);
   for (int i1=0; i1<n1; i1++)
   {
    for (int i2=0; i2<n2; i2++)
	{
     if (type>=3) h=abs(D[i1][i2]);
	 else
	  h=imag(D[i1][i2][type]);
	 os << h << "   ";
	}
	os << std::endl;
   }
  os.close();
}


void Detector::savePhase (char *fn, int type)
{
	std::ofstream os;
  double h;
  os.open (fn);
   for (int i1=0; i1<n1; i1++)
   {
    for (int i2=0; i2<n2; i2++)
	{
	  h=arg(D[i1][i2][type]);
	 os << h << "   ";
	}
	os << std::endl;
   }
  os.close();
}

/*bool Detector::cross(Vector<double> P, Vector<double> k, int& i1, int& i2, double& l)
{ 
  switch (type)
  {
  case DETECTOR_PLANE : return ((DetectorPlane *)this)->cross(P,k,i1,i2,l); break;
  }
  return false;
}
*/

DetectorPlane::DetectorPlane(void)
{
  n1=0;
  n2=0;
  D=0;
  type=DETECTOR_PLANE;
}

DetectorPlane::DetectorPlane(Vector<double> P, Vector<double> n, double d, int N)
{
	this->n = n / abs(n);
	if (abs(this->n % ex)>1E-5)
		e1 = ex - (ex * n) * n;
	else  
		e1 = ey - (ey * n) * n;
	e1 = e1 / abs(e1);
	e2 = n % e1;
	e2 = e2 / abs(e2);
	// d1 = d / (double)N;
	// d2 = d / (double)N;
	d1 = d;
	d2 = d;
	this->P = P;
	n1 = N;
	n2 = N;
	D = new Vector<std::complex<double> > *[N];
	for (int i=0; i<N; i++) 
		D[i]=new Vector<std::complex<double> > [N];
}

DetectorPlane::DetectorPlane (Vector<double> P, Vector<double> e1, Vector<double> e2, int n1, int n2)
{
    D=new Vector<std::complex<double> > *[n1];
 for (int i=0; i<n1; i++)
	 D[i]=new Vector<std::complex<double> >[n2];
 this->n1=n1;
 this->n2=n2;
 this->P=P;
 /*for (int i = 0; i<3; i++)
 {
  if (e1[i]!=0) 
  this->e1[i]=1.0/e1[i]*n1;
  else this->e1[i]=0;
  
  if (e2[i]!=0) 
  this->e2[i]=1.0/e2[i]*n2;
  else this->e2[i]=0;
 }*/
 if (n1 == 1) this->e1 = e1;
 else this->e1 = e1 / (double)(n1 - 1);

 if (n2 == 1) this->e2 = e2;
 else this->e2 = e2 / (double)(n2 - 1);

 this->n=e1%e2;
 this->n=this->n/abs(this->n);
 type=DETECTOR_PLANE;
}

bool DetectorPlane::cross(Vector<double> P, Vector<double> k, int &i1, int &i2, double &l)
{
	/*
	*  Ebene: n*(PE-P)=0 
	*  Strahl: P=k*l+PS 
	*  in Ebenengleichung n*PE-n*(k*l+PS)=0 => n*(PE-PS)=n*k*l => l=n*(PE-PS)/(n*k)
	*/
	Vector<double> dP = this->P - P;
	double kn = k * n;
	if (kn == 0) return false;
	l = (dP * n) / kn;

	/* Index berechnen */
	Vector<double> Ph = P + l * k;
	dP= Ph - this->P;
	i1 = dP * e1 * n1 / d1 + n1 / 2;
	i2 = dP * e2 * n2 / d2 + n2 / 2;

   if ( (i1<0) || (i1>=n1) || (i2<0) || (i2>=n2)) return false;
 return true;
}

std::ostream& operator << (std::ostream &os, Detector& D)
{
 for (int i1=0; i1<D.n1; i1++)
 { 
  for (int i2=0; i2<D.n2; i2++)
     os << D.D[i1][i2] << "  ";
  os << std::endl;
 }
 return os;
}
 
