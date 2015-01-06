Source Code for Experiments 2-8
 - these experiments are implemented in C++
 - to run the experiments, the required software dependencies must first be installed /
   compiled; the experiments are launched from main.cpp of the ICP_App project


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
These instructions are for the Windows platform running VisualStudio 2013
 
Compile cisstICP Library (the IMLP code):
  - see the "README" file in the "cisstICP" source code folder for instructions
        
Compile ICP_App  (the test code):
  - run CMake on source directory
     > specify paths to ICP_App source code and build directories
     > click Configure
        -- set cisst_Dir: set to main "cisst" folder in the CISST Library build directory
        -- set cisstICP_LIB: set to library file for cisstICP library built above
        -- set cisstICP_LIB_INCLUDE: set to source code directory of cisstICP code above
     > Generate project
  - compile code in Visual Studio
  - run main.cpp
     

Analysing Output Data of ICP_App
 - ICP_App stores data to predefined folders
 - Matlab scripts located in:  <your_path>\PLOSONE\Experiments2-8\ICP_TestData
    have been preconfigured to process this data and output a summary analysis


Other Algorithms Evaluated in the PLOSONE Paper Include:
--------------------------------------------------------

GICP:
 - was downloaded from:  http://www.robots.ox.ac.uk/~avsegal/generalized_icp.html
 - see IMLP paper for minor modifications required for termination condition and to set covariance matrices
 
CPD:
 - was downloaded from: www.bme.ogi.edu/~myron/matlab/cpd
 - now available at: https://sites.google.com/site/myronenko/research/cpd   
 - see IMLP paper for minor modifications required for termination condition    
    
 These algorithms can be run using output from the IMLP algorithm above;
  i.e. point sets, noisy points, covariances, etc. are all saved and can
  be loaded and run for these algorithms.
    
    