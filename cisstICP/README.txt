
Required Software Dependencies:
-------------------------------
 - CMake (http://www.cmake.org)
 - Visual Studio for Windows: http://www.visualstudio.com/downloads/download-visual-studio-vs
    > application was developed in Visual Studio 2013
 - CISST Libraries (https://trac.lcsr.jhu.edu/cisst)
    > requires: 
    > must compile this from source
 - WildMagic5 Libraries:  http://www.geometrictools.com/Downloads/Downloads.html
    > requires: Wm5Core.lib, Wm5Mathematics.lib
    > must compile these from source


Instructions for Compiling Source Code
--------------------------------------
These instructions are based on the Windows platform with VisualStudio 2013. The output is a static library (cisstICP.lib).

Compile WildMagic5 Libraries
 > compile Wm5Core.lib and Wm5Mathematics.lib
 > see website of WildMagic5 for instructions on this (www.geometrictools.com)
 
Compile CISST Libraries:
 - Get cisst repository from github
    > https://github.com/jhu-cisst/cisst/wiki/Download-cisst-and-SAW
    > https://github.com/jhu-cisst/cisst
    > only jhu-cisst/cisst required (not cisst-saw, etc.)
 - Create Visual Studio solution for CISST Libraries using CMake
    > run CMake on cisst source directory
       -- check boxes for libraries: CISST_cisstCommon, CISST_cisstVector, CISST_cisstNumerical, CISST_cisstOSAbstraction
       -- check box: CISST_HAS_CISSTNETLIB
    > click Configure
       -- check box: CISSTNETLIB_DOWNLOAD_NOW
    > generate project
    > build in Visual Studio
 
Compile cisstICP Library (the IMLP code):
  - run CMake on source directory
     > specify paths to cisstICP source code and build directories
     > specify path to WildMagic5 base directory (other WildMagic5 fields should auto-detect)
     > click Configure
        -- set cisst_Dir to location of CISST Library build (i.e. "<working_dir/cisst_build/cisst")
     > Generate project
     > build in Visual Studio
  - Build cisstICP library in Visual Studio
    
    