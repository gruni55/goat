GOAT - (G)eometrical (O)ptics (A)pplication (T)ool
by Thomas Weigel (C) 2021 

this project requires CMake by Kitware, which can be downloaded from https://cmake.org/ 

Please note that for Windows, the slash "/" must be replaced by the backslash "\" in the specifications for the directories. 


*** INSTALLATION NOTES ***

Windows (with Visual Studio) 
============================
Here, you have two possible ways: 

Just open the GOAT directory with Visual Studio (if your VS version supports CMake, which is true e.g. for VS 2019 Community edition or later) 
or you can use the command line: 
             1. open a command prompt and change to the GOAT directory
             2. create a Visual Studio with typing "cmake ."


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


*** Examples ***

in the folder /examples/ot 
    - Single spherical particle in an optical trap (example_one_sphere.cpp)
    - Particle in a counterpropagating trap, i.e. two opposite traps (counterpropagating_trap.cpp)
    
in /examples/paths
    One example, how to trace the rays from a light source and store it to a file. 
    Here is also shown how to import an object from a STL-file.  
    The generated file contains the start and end points (with the x, y and z coordinates) for each step:
    
    P1Startx P1Starty P1Startz P1Stopx P1Stopy P1Stopz
    P2Startx P2Starty P2Startz P2Stopx P2Stopy P2Stopz
    :        :        :        :       :       :
    
    



