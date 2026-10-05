# GOAT — Geometrical Optics Application Tool

GOAT is a modular, scene-based C++ library for the simulation of optical fields
and intensity distributions within the framework of geometrical optics with
explicit phase information.

**Author:** Thomas Weigel  
**License:** BSD-3-Clause  
**Initial release:** 2021  

---

## Scientific Reference

GOAT is described in the following peer-reviewed publication:

> Thomas Weigel, Gustav Schweiger, and Andreas Ostendorf  
> *GOAT: a multipurpose optical simulation tool*  
> Journal of the Optical Society of America B **39**, 2061–2065 (2022)  
> https://doi.org/10.1364/JOSAB.457951

If you use GOAT in scientific work, please cite this publication.

---

## Scope and Concept

GOAT is based on a scene description consisting of optical sources,
geometrical objects, and detectors. Ray tracing is performed within this scene,
and complex electric field amplitudes are propagated along the rays.

The library deliberately adopts **geometrical optics with explicit phase
information** as its underlying model. Within this framework, interference
effects arising from phase differences can be described, whereas phenomena
that are inherently wave-optical in nature, such as diffraction, lie outside
the scope of the model.

---

## Key Features

- Scene-based and modular architecture (sources, objects, detectors)
- Ray tracing with propagation of complex electric field amplitudes
- Explicit phase tracking for the analysis of interference effects
- Efficient handling of complex geometries, including STL import and octree acceleration
- Cross-platform C++17 implementation (Linux and Windows)
- Extensive API documentation generated with Doxygen

---

## Installation and Build

GOAT is built using CMake and requires a C++17-compatible compiler.

### Requirements
- CMake
- C++17-compatible compiler (GCC, Clang, or Microsoft Visual Studio)
- TinyXML-2 (included in the repository)

### Linux
```bash
cmake .
make
