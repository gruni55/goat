# GOAT – Geometrical Optics Application Tool

**GOAT** is a modular, scene-based C++ library for flexible simulation of electric fields and intensity distributions using geometrical optics, including explicit phase calculation for the simulation of interference phenomena.

## Why GOAT?

In optical modeling, each new problem often requires starting from scratch. **GOAT** solves this by providing a reusable, extensible library designed for a wide range of optical tasks. Its architecture is scene-based and modular, allowing you to combine and extend components as needed.

## Key Features

- **Abstract, extensible ray tracer:**  
  Easily adapt the ray tracing engine to your problem.
- **Scene architecture:**  
  Define a simulation by specifying light sources and objects in a scene.
- **Flexible interfaces:**  
  Extend base classes for your own object and light source types.
- **Ready-made sources and shapes:**  
  Gaussian beams, top-hat sources, ellipsoids, and import of binary STL shapes.
- **Efficient 3D math classes:**  
  Matrix and Vector templates for mathematical operations (in `/maths`).
- **Parallel processing for large STL files:**  
  The `goat_raytracing_mp.lib` library offers **OpenMP support** to accelerate the loading and processing of large STL files.

## Library Structure

- **goat_maths.lib:**  
  Mathematical basics: 3D vector and matrix classes with operators.
- **goat_raytracing.lib:**  
  Ray tracing engine, light sources, object representations.
- **goat_raytracing_mp.lib:**  
  Same as `goat_raytracing.lib`, but with **OpenMP parallelization** for efficient loading and handling of large STL geometries.

## Documentation

- **Full API documentation (Doxygen):**  
  [https://gruni55.github.io/goat/html/](https://gruni55.github.io/goat/html/)  
  Also included locally in `/docs`.

For more information, see [README.txt](README.txt).

## Installation

**Requirements:**  
- [CMake](https://cmake.org/) (recommended)
- C++17 compatible compiler (tested: Visual Studio 2019+, GCC)
- For `goat_raytracing_mp.lib`: OpenMP support enabled in your compiler.

**Build steps:**  
1. Open a terminal / console.
2. Navigate to your GOAT directory.
3. Run:  cmake .
*(Don’t forget the dot!)*

CMake detects your system and generates project files or a Makefile.
- **Windows (VS 2019+):** Open the GOAT directory directly with Visual Studio, or use CMake as above.
- **Linux:** Run `make` after `cmake`.

Libraries will be created in `/lib`.

## XML Support

Scenes and calculation setups can be saved and loaded as XML files.  
See [XML.md](XML.md) for details.

---

**Citation:**  
If you use GOAT in your work, please cite:  
Thomas Weigel, Gustav Schweiger, and Andreas Ostendorf,  
*"GOAT: a multipurpose optical simulation tool,"*  
J. Opt. Soc. Am. B 39, 2061-2065 (2022).  
[https://doi.org/10.1364/JOSAB.457951](https://opg.optica.org/josab/fulltext.cfm?uri=josab-39-8-2061&id=477834)

