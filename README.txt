GOAT - Geometrical Optics Application Tool
==========================================

Author: Thomas Weigel (C) 2021
License: BSD-3-Clause

GOAT is a modular, scene-based C++ library for the simulation of electric fields and intensity distributions using geometrical optics, including explicit phase calculation for the simulation of interference phenomena.

Homepage / Source code:
    https://github.com/gruni55/goat

Online API documentation (Doxygen):
    https://gruni55.github.io/goat/html/

Scientific reference:
    Thomas Weigel, Gustav Schweiger, and Andreas Ostendorf,
    "GOAT: a multipurpose optical simulation tool,"
    J. Opt. Soc. Am. B 39, 2061-2065 (2022)
    https://doi.org/10.1364/JOSAB.457951

-------------------------------------------------------------------------------

KEY FEATURES

- Scene-based, modular architecture (sources, objects, detectors)
- Explicit phase calculation for interference effects (no diffraction)
- Efficient geometry handling, including octree and STL file import
- Multiple ray types: line rays and tube rays
- Cross-platform, open source (BSD-3-Clause License)
- Extensive Doxygen API documentation

-------------------------------------------------------------------------------

INSTALLATION

REQUIREMENTS
    - CMake (https://cmake.org/)
    - C++17 compatible compiler (GCC, clang, or Visual Studio)
    - TinyXML-2 (included in repository for XML support)

LINUX / macOS
    1. Open a terminal and go to the GOAT directory.
    2. Configure the build:
         cmake .
    3. Compile:
         make

WINDOWS (Visual Studio 2019 or newer)
    1. Open the GOAT directory with Visual Studio (with CMake support)
       OR open a command prompt in the GOAT directory and run:
         cmake .
    Note: Replace '/' with '\' in all folder paths for Windows.

Resulting libraries:
    - goat_maths.lib   (mathematical basics)
    - goat_raytracing.lib   (ray tracing functions and classes)

Libraries can be found in /lib/Release or /lib/Debug depending on build type.
Header files are in the /maths and /raytracing directories.

-------------------------------------------------------------------------------

EXAMPLES

You will find several example programs in the /examples directory:

  /examples/ot         Optical trap examples (single sphere, counterpropagating trap)
  /examples/layers     Simulation of transmission through thin tilted layers
  /examples/paths      Ray tracing from sources, STL object import, saving ray paths
  /examples/axicon     Gaussian beam and axicon example (SRF file required)

Please copy any required SRF files into the folder where the binary is located before running the examples.

Output files will be written to the same directory as the executable.
Visual Studio default output folders are: 
  goat-main/out/build/x64-Debug/bin
  goat-main/out/build/x64-Release/bin

-------------------------------------------------------------------------------

DOCUMENTATION

A complete API documentation generated with Doxygen is available:

    - Online: https://gruni55.github.io/goat/html/
    - **Locally:** in the `/docs` directory of this repository

To generate or update the documentation yourself, run:
    doxygen Doxyfile

-------------------------------------------------------------------------------

CONTACT AND SUPPORT

- For questions, bug reports, or feature requests, please use the GitHub Issues page.
- More information can be found in the online documentation and in the `/docs` folder.

-------------------------------------------------------------------------------

CITATION

If you use GOAT in your work, please cite:

    Thomas Weigel, Gustav Schweiger, and Andreas Ostendorf,
    "GOAT: a multipurpose optical simulation tool,"
    J. Opt. Soc. Am. B 39, 2061-2065 (2022).
    https://doi.org/10.1364/JOSAB.457951

-------------------------------------------------------------------------------
