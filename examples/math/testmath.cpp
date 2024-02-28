#include "goodies.h"
#include <iostream>

double f(double x)
{
	return (x - 2) * (x + 3);
	// x^2 + 3x - 2x - 6=x^2+x-6
}

double df(double x)
{
	return 2 * x + 1;
}

int main (int argc, char **argv)
{
	std::cout << "x^2+x-6:" << std::endl;
	std::cout << GOAT::maths::newton_root(&f, &df) << std::endl;
	std::cout << GOAT::maths::newton_root(&f, &df,-6) << std::endl;
  return 0;
}
