# GOAT

Geometric optics provides a flexible tool for calculating a variety of optical tasks. If you change from one problem to another, you often have to start again from the beginning. To prevent this, a programming library has been developed that provides a flexible, customizable platform. Geometric optics provides a flexible tool for calculating a variety of optical tasks. If you change from one problem to another, you often have to start again from the beginning. To prevent this, a programming library has been developed that provides a flexible, customizable platform. 

The library consists of an abstract ray tracer that can be adapted to the needs. 
Thereby interfaces are offered, which are called whenever a ray hits an object or emerges from it. To describe the setup in a simple way, a scene consisting of light sources and objects is defined before each calculation. Like the raytracers themselves, base classes are available for this purpose, from which the concrete object and light source types can be derived. As an example, the light source types for Gaussian rays and top hats are available. Ellipsoids are available as object shapes. Furthermore, shapes can be imported from binary STL files. 

To make the mathematical description as simple as possible, a vector as well as a matrix template class for a 3D calculation was developed. 

The library is divided into two parts: 





## Installation instructions

For a proper installation, cmake is highly recommended. You can download it from [](https://cmake.org/)
and of course, a working C++-Compiler (We tested it on Windows with Visual Studio Community Edition 2019 and on Linux with gcc)

open a console window a 

change to your GOAT directory 

type "cmake ." (don't forget the dot!) 

Cmake automatically recognizes the corresponding compiler system. For the Visual Studio, a corresponding project file will be generated and for the gcc a Makefile is created. 





