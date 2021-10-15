#include "goodies.h"
#include <iostream>
#include <math.h>
#include <string>
#include <sstream>

using namespace std;

bool testnan (double x)
{
 stringstream ss;
 string S; 
 ss << x;
 ss >> S;
 return S=="nan"; 
}

double sqr(double x){return x*x;}

double abs (complex<double> x)
{
 return real(sqrt(x*conj(x)));
}



