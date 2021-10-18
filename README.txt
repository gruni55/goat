GOAT - (G)eometrical (O)ptics (A)pplication (T)ool
by Thomas Weigel (C) 2021 

this project requires CMake by Kitware 


Windows (with Visual Studio) 
============================
to create a Visual Studio project just type "cmake ."


Linux with GCC compiler
=======================
type "cmake ." (for Release compilation) or "cmake -DCMAKE_BUILD_TYPE=Debug ." (for Debug compilation)

then start the compilation with "make"


the following libraries will be created 
- goat_maths.lib (for the mathematical basics)
- goat_raytracing.lib (for all function and classes needed for the raytracing)


the libraries can be found in /lib/Release or in /lib/Debug, depending on wether the library is compiled in Release or in Debug mode.

the header files are in the directories /maths and /raytracing

an documentation created by doxygen can be found in the folder /documentation

