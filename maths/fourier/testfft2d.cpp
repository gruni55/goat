#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <stdexcept>

#include "fft2d.h"

using complexd = std::complex<double>;

double maxAbsDiff(const std::vector<complexd>& a, const std::vector<complexd>& b)
{
    if (a.size() != b.size())
        throw std::runtime_error("Size mismatch in maxAbsDiff");

    double maxErr = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i)
    {
        double err = std::abs(a[i] - b[i]);
        if (err > maxErr)
            maxErr = err;
    }
    return maxErr;
}

bool testForwardInverse()
{
    const std::size_t nx = 8;
    const std::size_t ny = 8;
    const std::size_t n = nx * ny;

    GOAT::maths::fourier::fft2D fft(nx, ny);



    std::vector<complexd> in(n);
    std::vector<complexd> spec(n);
    std::vector<complexd> back(n);

    for (int y = 0; y < ny; ++y)
    {
        for (int x = 0; x < nx; ++x)
        {
            std::size_t idx = y * nx + x;
            // in[idx] = complexd(static_cast<double>(x + 2 * y), static_cast<double>(x - y));
			double d = (double)x - (double)y;
            in[idx] = complexd(x + 2*y, x-y);
        }
    }

    fft.forward(in.data(), spec.data());
    fft.inverse(spec.data(), back.data());
    


    double err = maxAbsDiff(in, back);

    std::cout << "[testForwardInverse] max error = " << err << std::endl;
    return err < 1e-10;
}

bool testDeltaSpectrum()
{
    const std::size_t nx = 8;
    const std::size_t ny = 8;
    const std::size_t n = nx * ny;

    GOAT::maths::fourier::fft2D fft(nx, ny);

    std::vector<complexd> in(n, complexd(0.0, 0.0));
    std::vector<complexd> spec(n);

    in[0] = complexd(1.0, 0.0);

    fft.forward(in.data(), spec.data());

    double maxErr = 0.0;
    for (std::size_t i = 0; i < n; ++i)
    {
        double err = std::abs(spec[i] - complexd(1.0, 0.0));
        if (err > maxErr)
            maxErr = err;
    }

    std::cout << "[testDeltaSpectrum] max error = " << maxErr << std::endl;
    return maxErr < 1e-10;
}

bool testConstantField()
{
    const std::size_t nx = 8;
    const std::size_t ny = 8;
    const std::size_t n = nx * ny;

    GOAT::maths::fourier::fft2D fft(nx, ny);

    std::vector<complexd> in(n, complexd(1.0, 0.0));
    std::vector<complexd> spec(n);

    fft.forward(in.data(), spec.data());

    double dcErr = std::abs(spec[0] - complexd(static_cast<double>(n), 0.0));
    double otherMax = 0.0;

    for (std::size_t i = 1; i < n; ++i)
    {
        double val = std::abs(spec[i]);
        if (val > otherMax)
            otherMax = val;
    }

    std::cout << "[testConstantField] dc error = " << dcErr
        << ", max other = " << otherMax << std::endl;

    return dcErr < 1e-10 && otherMax < 1e-10;
}

int main()
{
    bool ok = true;

    ok = testForwardInverse() && ok;
    ok = testDeltaSpectrum() && ok;
    ok = testConstantField() && ok;

    if (ok)
    {
        std::cout << "All fft2D tests passed." << std::endl;
        return 0;
    }

    std::cerr << "At least one fft2D test failed." << std::endl;
    return 1;
}