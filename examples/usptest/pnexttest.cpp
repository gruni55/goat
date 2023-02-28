#include "vector.h"
#define SGN(x) (x<0) ? -1 : ((x>0) ? 1 : 0)
 GOAT::maths::Vector<double> pnext(GOAT::maths::Vector<double> p0, GOAT::maths::Vector<double> k0, GOAT::maths::Vector<double> d, GOAT::maths::Vector<int> nges)
{
    double lambdax, lambday, lambdaz, lambda;
    double signx, signy, signz;
    double sx, sy, sz;

    signx = SGN(k0[0]);
    signy = SGN(k0[1]);
    signz = SGN(k0[2]);

    int nx = p0[0] / d[0] + signx;
    int ny = p0[1] / d[1] + signy;
    int nz = p0[2] / d[2] + signz;

    std::cout << "nx=" << nx << "  ny=" << ny << "  nz=" << nz << std::endl;

    lambdax = (nx * d[0] - p0[0]) / k0[0];
    lambday = (ny * d[1] - p0[1]) / k0[1];
    lambdaz = (nz * d[2] - p0[2]) / k0[2];

    // fabs avoids question after k0[i]=0 which can result into -inf !
    if (fabs(lambdax) < fabs(lambday))
        lambda = lambdax;
    else
        lambda = lambday;
    if (fabs(lambdaz) < lambda)
        lambda = lambdaz; 
    return  p0 + lambda * k0;

}

int main(int argc, char** argv)
{
    double r0 = 500E-6;
    GOAT::maths::Vector<double> k(1.0, 1.0, 1.0);
    GOAT::maths::Vector <int> nges(1000, 1000, 1000);
    GOAT::maths::Vector<double> P0(-212.34E-6, 20E-6, -10E-6);
    GOAT::maths::Vector<double> d(2.0 * r0 / (double)nges[0], 2.0 * r0 / (double)nges[1], 2.0 * r0 / (double)nges[2]);
    GOAT::maths::Vector<double> P=P0;
    for (int i = 1; i <= 10; i++)
    {
        P = pnext(P0, k, d, nges);
        std::cout << P << "  " << P-P0 << "  " << d << std::endl;
        P0 = P;
    }
    return 0;
}