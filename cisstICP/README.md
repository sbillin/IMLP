
This is a preliminary README for how to compile and run the cisstICP framework.

Dependencies:
--------------
 - CMake (http://www.cmake.org)
 - Visual Studio for Windows: http://www.visualstudio.com/downloads/download-visual-studio-vs
    > (application was developed in Visual Studio 2013)
 - CISST Libraries (https://trac.lcsr.jhu.edu/cisst)
    > requires: 
    > must compile this from source
 - WildMagic5 Libraries:  http://www.geometrictools.com/Downloads/Downloads.html
    > requires: Wm5Core.lib, Wm5Mathematics.lib
    > must compile these from source
 - Numerical Recipes (needed only for GIMLOP algorithm; other algorithms don't use it and the code can be built without it)
 - PLY_IO Library (for reading/writing PLY files in Matlab)

**NOTES:**

* Visual Studio 12 2013 Win64 was used successfully when creating this README
* Visual Stuido 14 2015 Win64 did not work due to errors such as `cisstNetlib_f2c.lib(endfile.obj) : error LNK2001: unresolved external symbol __imp_sprintf`. This issue has not been solved, thus it is recommended to build the cisstICP application using Visual Studio 2013.

# CISST Libraries

Reference

* [CISST Libraries Git Repo](https://github.com/jhu-cisst/cisst/)

Clone the Git repo for the CISST libraries

    git clone https://github.com/jhu-cisst/cisst.git

Now create a build directory `cisst_build` alongside (not inside) the "cisst" directory that was cloned from git. Then open CMake and point the source and build folders to `cisst` and `cisst_build` respectively. Hit configure and select `Visual Studio 12 2013 Win64` as the compiler and click Finish. Once CMake has finished configuring, apply the following options:

Uncheck: cisstMultitask, cisstMultiTask_*, cisstParameterTypes, cisstRobot_*.
Check: CISST_HAS_CISSTNETLIB

Press "Configure" button again and then check "CISSTNETLIB_DOWNLOAD_NOW" 

Press "Configure" again and then "Generate".

Open the `cisst_build/cisst.sln` file in visual studio. Set the build mode to "Release" and "x64" and then build the solution.

# WildMagic5 Libraries

These libraries are used solely for the closed-form implementation of computing the decomposition of a 3x3 covariance matrix.  The libraries are no longer available online from their author, but they are included in the `dependencies` folder of the cisstICP repo.

Extract the `WildMagic5p13.zip` file located in the dependencies folder. Open the `GeometricTools\WildMagic5\WildMagic5Wgl_VC120.sln` file in Visual Studio.  Set the build configuration to `Release` and `x64`. In the Solution Explorer, right click `LibCore_VC120` and select build. When that build has completed, right click `LibMathematics_VC120` in the Solution Explorer and select build.

# dlib

Extract the `dlib-a8.6.7z` file located in the `dependencies` folder.

# cisstICP Library

Create a folder `cisstICP_build` within the `cisstICP` repo folder (alongside the `cisstICP` folder within the repo). Open CMake, point the source and build directories to these folders, and press Configure.

Set `cisst_DIR` to the path to the `cisst_build` folder and configure again. If the compiled CISST libraries were found correctly by CMake, then `CISST_USE_FILE` should now automatically contain a path to the file `Usecisst.cmake`.

Set `DLIB_INCLUDE` to the path to the extracted dlib folder.
Set `RPLY_DIR` to the path to `cisstICP/dependencies/rply-1.1.4`.
Set 'WM5_BASE_DIR' to the path to 'cisstICP/dependencies/GeometricTools/WildMagic5' 

Press Configure again and then Generate.

Open `cisstICP_build/cisstICP.sln` in Visual Studio, set the build configuration to 'Release' and 'x64' and then build the solution.

# matlabICP Library

This library provides a Matlab interface to the C++ code of the cisstICP Library.

Create a `matlabICP_build` folder alongside the `matlabICP` folder in the cloned repo.  Set CMake's source and build locations to these folders and then press Configure.

Set `cisst_DIR` to the path to the `cisst_build` folder and configure again. If the compiled CISST libraries were found correctly by CMake, then `CISST_USE_FILE` should now automatically contain a path to the file `Usecisst.cmake`.
Set `DLIB_INCLUDE` to the path to the extracted dlib folder.
Set 'WM5_BASE_DIR' to the path to 'cisstICP/dependencies/GeometricTools/WildMagic5' 

Set `MATLAB_ENGINE_INCLUDE_DIR` and `MATLAB_ENGINE_LIB_DIR` to your system paths for Matlab, such as `C:/Program Files/MATLAB/R2015a/extern/include` and `C:/Program Files/MATLAB/R2015a/extern/lib/win64/microsoft` (these paths will change depending on the Matlab version).

Set `cisstICP_LIB` to the path to `/cisstICP/cisstICP_build/Release/cisstICP.lib`
Set `cisstICP_LIB_INCLUDE` to the path to `cisstICP/cisstICP`

Configure again and then Generate.  Now open the `matlabICP_build/matlabICP.sln` file in Visual Studio, set the build options to 'Release' and 'x64', and build the solution.

If you do not have the Numerical Recipes code, then the Matlab interface for GIMLOP will fail to compile (since GIMLOP will not have been compiled in the C++ library), giving an error. This is not a problem, as the Matlab interfaces for the other algorithms won't be affected and can still be used. 

Extract the `/cisstICP/dependencies/MatlabDependencies/mtimesx.zip` and add it to your Matlab path.

Also add `/cisstICP/dependencies/MatlabDependencies/PLY_IO` to your Matlab path.  Note that the [PLY_IO library code](http://people.sc.fsu.edu/~jburkardt/m_src/ply_io/ply_io.html) included here containes a bug fix for `ply_read.m`, changing line 571 from

    if ( ( nargin > 1 & strcmpi(Str,'Tri') ) || nargout > 2 )

to

    if ( ( nargin > 1 && strcmpi(Str,'Tri') ) || nargout > 2 )

# Test Run

To test that everything is working, open `cisstICP/matlabICP/TestApps/App_Test_IMLOP.m` in Matlab and hit run.
