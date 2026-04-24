#include <fftw3.h>
#include <iostream>

int main()
{
    void* p = fftw_malloc(16);
    fftw_free(p);
    std::cout << "FFTW function link works\n";
    return 0;
}