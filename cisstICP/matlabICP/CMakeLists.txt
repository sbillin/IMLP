# See license at http://www.cisst.org/cisst/license.txt

cmake_minimum_required(VERSION 2.6)
project (matlabICP)

IF (CMAKE_CL_64 OR CMAKE_GENERATOR MATCHES Win64)
  set( X86_MODE 0 )  # 64-bit compile
ELSE (CMAKE_CL_64 OR CMAKE_GENERATOR MATCHES Win64)  
  set( X86_MODE 1 )  # 32-bit compile
ENDIF (CMAKE_CL_64 OR CMAKE_GENERATOR MATCHES Win64)

# create a list of libraries needed for this project
set ( REQUIRED_CISST_LIBRARIES cisstCommon cisstVector cisstNumerical cisstOSAbstraction )

# find cisst and make sure the required libraries have been compiled
find_package (cisst REQUIRED ${REQUIRED_CISST_LIBRARIES})

if (cisst_FOUND_AS_REQUIRED)

  # load cisst configuration
  include (${CISST_USE_FILE} )

  #=== dlib ===
  # TODO: figure out how to remove this dependency within vonMisesPrj interface
  find_path( DLIB_INCLUDE dlib
    NAMES "dlib/dlib_include_path_tutorial.txt"
    HINTS "../dependencies/dlib-18.6"
    DOC  "Choose folder that CONTAINS the dlib folder"
    )
  include_directories(${DLIB_INCLUDE})
  
  # #=== GSL ===
  # IF (X86_MODE)  
    # # GSL library
    # find_path( GSL_INCLUDE gsl
      # HINTS "C:/software/gsl-1.14-bin-msvc9/gsl-1.14/include"
      # )
    # find_library( GSL_LIB gsl.lib
      # HINTS "C:/software/gsl-1.14-bin-msvc9/gsl-1.14/lib"
      # )
    # find_library( GSL_CBLAS_LIB gslcblas.lib
      # HINTS "C:/software/gsl-1.14-bin-msvc9/gsl-1.14/lib"
      # )    
    
    # include_directories( 
      # ${GSL_INCLUDE}
      # )    
  # ENDIF (X86_MODE)  
  
  # WildMagic5 Mathematics Library
  find_path (WM5_BASE_DIR
    NAMES wm5path.txt
    PATHS "../dependencies/GeometricTools/WildMagic5"
    )
  if (IS_DIRECTORY ${WM5_BASE_DIR})
    # Release:
    find_path (WM5_LIB_PATH
      NAMES Wm5Core.lib
      #PATHS "C:/software/WildMagic5/SDK/Library/v120/x64/Release"
      HINTS ${WM5_BASE_DIR}"/SDK/Library/v120/x64/Release"
      )    
    find_library( WM5_CORE_LIB Wm5Core.lib PATHS ${WM5_LIB_PATH})
    find_library( WM5_MATH_LIB Wm5Mathematics.lib PATHS ${WM5_LIB_PATH})
    # # Debug:
    # find_path (WM5_LIB_PATH
      # NAMES Wm5CoreD.lib
      # #PATHS "C:/software/WildMagic5/SDK/Library/v120/x64/Release"
      # HINTS ${WM5_BASE_DIR}"/SDK/Library/v120/x64/Debug"
      # )    
    # find_library( WM5_CORE_LIB Wm5CoreD.lib PATHS ${WM5_LIB_PATH})
    # find_library( WM5_MATH_LIB Wm5MathematicsD.lib PATHS ${WM5_LIB_PATH})    
    # # # find_library( WM5_CORE_LIB_DEBUG Wm5CoreD.lib PATHS ${WM5_LIB_PATH_DEBUG})
    # # # find_library( WM5_MATH_LIB_DEBUG Wm5MathematicsD.lib PATHS ${WM5_LIB_PATH_DEBUG})        
    # # # # Debug:
    # # # find_library( WM5MATH_LIB Wm5MathematicsD.lib Wm5CoreD.lib
      # # # HINTS "C:/software/WildMagic5/SDK/Library/v120/x64/Debug"
      # # # )
      
    find_path( WM5_LIB_INCLUDE
      NAMES "Wm5NoniterativeEigen3x3.h"
      PATHS
      HINTS ${WM5_BASE_DIR}"/SDK/Include"
      #HINTS "C:/software/WildMagic5/SDK/Include"
      DOC  ""
      )
    include_directories(${WM5_LIB_INCLUDE})  
  endif (IS_DIRECTORY ${WM5_BASE_DIR})
  
  
  #=== ICP LIBRARY ===
  
  # If compiling in debug set library to cisstICP
  #  otherwise use the cisstICPCombined library which includes all dependencies
  #  compiled in release mode
  find_library(cisstICP_LIB cisstICP
  HINTS "../cisstICP_build/Release/cisstICP.lib")
  find_path(cisstICP_LIB_INCLUDE 
    NAMES cisstICP.h
    HINTS "../cisstICP"
    DOC  "cisstICP source code directory"
    )
  
  #SET(LINK_LIBRARY optimized Foo debug Foo_d)
  #target_link_libraries(MyEXE ${LINK_LIBRARY})
  
  
  #=== Mex Libraries for Matlab Interface (Shared) ===

  # Matlab Engine
  find_path(MATLAB_ENGINE_INCLUDE_DIR
    NAMES engine.h matrix.h
    HINTS "C:/Program Files/MATLAB/R2015a/extern/include"
    DOC  "$MATLABROOT/extern/include"
    )
  set( MEX_BASE_FILES
    matlabClassHandle.h
    matlabParser.h
    matlabExtras.h
    )    
  include_directories(${MATLAB_ENGINE_INCLUDE_DIR})
  find_path(MATLAB_ENGINE_LIB_DIR 
    NAMES libmex.lib
    HINTS "C:/Program Files/MATLAB/R2015a/extern/lib/win64/microsoft"  
    DOC "$MATLABROOT/extern/lib/win64/microsoft"
    )
  set( MATLAB_ENGINE_LIBRARIES
    libmx.lib
    libmat.lib
    libeng.lib
    libmex.lib
    )
  link_directories(${MATLAB_ENGINE_LIB_DIR})    

  #set(CMAKE_SHARED_LINKER_FLAGS "-Wl,--export-all-symbols")  

  
  # ---- mexAlgDirICP_IMLOP_Mesh ---- #
  
  add_library( mexAlgDirICP_IMLOP_Mesh SHARED
    mexInterface_AlgDirICP_IMLOP_Mesh.cpp
    mexInterface_AlgDirICP_IMLOP_Mesh.h
    mexInterface_AlgDirICP_IMLOP_Mesh.def
    ${MEX_BASE_FILES}
    )
  include_directories(
    ${mexAlgDirICP_IMLOP_Mesh_SOURCE_DIR}
    ${cisstICP_LIB_INCLUDE}
    )
  target_link_libraries (mexAlgDirICP_IMLOP_Mesh
    ${MATLAB_ENGINE_LIBRARIES}
    ${cisstICP_LIB}
    #${WM5_MATH_LIB_DEBUG} # only need this if compiling in debug
    ${WM5_MATH_LIB}
    ${WM5_CORE_LIB}
    )
  cisst_target_link_libraries ( mexAlgDirICP_IMLOP_Mesh
    ${REQUIRED_CISST_LIBRARIES}
    )
  IF (X86_MODE)
    MESSAGE("Configuring mex as 32-bit library (mexw32)...\n")
    MESSAGE("32-bit mex file will not run on 64-bit Matlab!\n")
    SET_TARGET_PROPERTIES( mexAlgDirICP_IMLOP_Mesh
      PROPERTIES
      PREFIX ""
      SUFFIX ".mexw32"
      )
  ELSE (X86_MODE)
    MESSAGE("Configuring mex as 64-bit library (mexw64)...\n")
    SET_TARGET_PROPERTIES( mexAlgDirICP_IMLOP_Mesh
      PROPERTIES
      PREFIX ""
      SUFFIX ".mexw64"
      )
  ENDIF (X86_MODE)
  
  
  # ---- mexAlgDirICP_GIMLOP_Mesh ---- #
  
  add_library( mexAlgDirICP_GIMLOP_Mesh SHARED
    mexInterface_AlgDirICP_GIMLOP_Mesh.cpp
    mexInterface_AlgDirICP_GIMLOP_Mesh.h
    mexInterface_AlgDirICP_GIMLOP_Mesh.def
    ${MEX_BASE_FILES}
    )
  include_directories(
    ${mexAlgDirICP_GIMLOP_Mesh_SOURCE_DIR}
    ${cisstICP_LIB_INCLUDE}
    )
  target_link_libraries (mexAlgDirICP_GIMLOP_Mesh
    ${MATLAB_ENGINE_LIBRARIES}
    ${cisstICP_LIB}
    #${WM5_MATH_LIB_DEBUG} # only need this if compiling in debug
    ${WM5_MATH_LIB}
    ${WM5_CORE_LIB}
    )
  cisst_target_link_libraries ( mexAlgDirICP_GIMLOP_Mesh
    ${REQUIRED_CISST_LIBRARIES}
    )
  IF (X86_MODE)
    MESSAGE("Configuring mex as 32-bit library (mexw32)...\n")
    MESSAGE("32-bit mex file will not run on 64-bit Matlab!\n")
    SET_TARGET_PROPERTIES( mexAlgDirICP_GIMLOP_Mesh
      PROPERTIES
      PREFIX ""
      SUFFIX ".mexw32"
      )
  ELSE (X86_MODE)
    MESSAGE("Configuring mex as 64-bit library (mexw64)...\n")
    SET_TARGET_PROPERTIES( mexAlgDirICP_GIMLOP_Mesh
      PROPERTIES
      PREFIX ""
      SUFFIX ".mexw64"
      )
  ENDIF (X86_MODE)
  
  
  # ---- mexAlgDirICP_PIMLOP_Mesh ---- #
  
  add_library( mexAlgDirICP_PIMLOP_Mesh SHARED
    mexInterface_AlgDirICP_PIMLOP_Mesh.cpp
    mexInterface_AlgDirICP_PIMLOP_Mesh.h
    mexInterface_AlgDirICP_PIMLOP_Mesh.def
    ${MEX_BASE_FILES}
    )
  include_directories(
    ${mexAlgDirICP_PIMLOP_Mesh_SOURCE_DIR}
    ${cisstICP_LIB_INCLUDE}
    )
  target_link_libraries (mexAlgDirICP_PIMLOP_Mesh
    ${MATLAB_ENGINE_LIBRARIES}
    ${cisstICP_LIB}
    ${WM5_MATH_LIB}
    ${WM5_CORE_LIB}
    )
  cisst_target_link_libraries ( mexAlgDirICP_PIMLOP_Mesh
    ${REQUIRED_CISST_LIBRARIES}
    )
  IF (X86_MODE)
    MESSAGE("Configuring mex as 32-bit library (mexw32)...\n")
    MESSAGE("32-bit mex file will not run on 64-bit Matlab!\n")
    SET_TARGET_PROPERTIES( mexAlgDirICP_PIMLOP_Mesh
      PROPERTIES
      PREFIX ""
      SUFFIX ".mexw32"
      )
  ELSE (X86_MODE)
    MESSAGE("Configuring mex as 64-bit library (mexw64)...\n")
    SET_TARGET_PROPERTIES( mexAlgDirICP_PIMLOP_Mesh
      PROPERTIES
      PREFIX ""
      SUFFIX ".mexw64"
      )
  ENDIF (X86_MODE)
    
    
  # ---- mexAlgICP_StdICP_Mesh ---- #
  
  add_library( mexAlgICP_StdICP_Mesh SHARED
    mexInterface_AlgICP_StdICP_Mesh.cpp
    mexInterface_AlgICP_StdICP_Mesh.h
    mexInterface_AlgICP_StdICP_Mesh.def 
    ${MEX_BASE_FILES}
    )  
  include_directories(
    ${mexAlgICP_StdICP_Mesh_SOURCE_DIR}
    ${cisstICP_LIB_INCLUDE}
    )  
  target_link_libraries (mexAlgICP_StdICP_Mesh
    ${MATLAB_ENGINE_LIBRARIES}
    ${cisstICP_LIB}
    ${WM5_MATH_LIB}
    ${WM5_CORE_LIB}
    )
  # IF (X86_MODE)
    # target_link_libraries( mexAlgDirICP_StdICP_Mesh 
      # ${GSL_LIB}
      # ${GSL_CBLAS_LIB}
      # )
  # ENDIF (X86_MODE)
  cisst_target_link_libraries ( mexAlgICP_StdICP_Mesh
    ${REQUIRED_CISST_LIBRARIES}
    )
  IF (X86_MODE)
    MESSAGE("Configuring mex as 32-bit library (mexw32)...\n")
    MESSAGE("32-bit mex file will not run on 64-bit Matlab!\n")
    SET_TARGET_PROPERTIES( mexAlgICP_StdICP_Mesh
      PROPERTIES
      PREFIX ""
      SUFFIX ".mexw32"
      )   
  ELSE (X86_MODE)
    MESSAGE("Configuring mex as 64-bit library (mexw64)...\n")
    SET_TARGET_PROPERTIES( mexAlgICP_StdICP_Mesh
      PROPERTIES
      PREFIX ""
      SUFFIX ".mexw64"
      )  
  ENDIF (X86_MODE)
	    
    
  # ---- mexAlgDirICP_StdICP_Mesh ---- #
  
  add_library( mexAlgDirICP_StdICP_Mesh SHARED
    #stdint.h   # not needed for Visual Studio 2010 and later
    mexInterface_AlgDirICP_StdICP_Mesh.cpp
    mexInterface_AlgDirICP_StdICP_Mesh.h
    mexInterface_AlgDirICP_StdICP_Mesh.def 
    ${MEX_BASE_FILES}
    )  
  include_directories(
    ${mexAlgDirICP_StdICP_Mesh_SOURCE_DIR}
    ${cisstICP_LIB_INCLUDE}
    )  
  target_link_libraries (mexAlgDirICP_StdICP_Mesh
    ${MATLAB_ENGINE_LIBRARIES}
    ${cisstICP_LIB}
    ${WM5_MATH_LIB}
    ${WM5_CORE_LIB}   
    )    
  # IF (X86_MODE)
    # target_link_libraries( mexAlgDirICP_StdICP_Mesh 
      # ${GSL_LIB}
      # ${GSL_CBLAS_LIB}
      # )
  # ENDIF (X86_MODE)
  cisst_target_link_libraries ( mexAlgDirICP_StdICP_Mesh
    ${REQUIRED_CISST_LIBRARIES}
    )
  IF (X86_MODE)
    MESSAGE("Configuring mex as 32-bit library (mexw32)...\n")
    MESSAGE("32-bit mex file will not run on 64-bit Matlab!\n")
    SET_TARGET_PROPERTIES( mexAlgDirICP_StdICP_Mesh
      PROPERTIES
      PREFIX ""
      SUFFIX ".mexw32"
      )   
  ELSE (X86_MODE)
    MESSAGE("Configuring mex as 64-bit library (mexw64)...\n")
    SET_TARGET_PROPERTIES( mexAlgDirICP_StdICP_Mesh
      PROPERTIES
      PREFIX ""
      SUFFIX ".mexw64"
      )  
  ENDIF (X86_MODE)
	

  # ---- mexAlgPDTree_MLP_Mesh ---- #
  
  add_library( mexAlgPDTree_MLP_Mesh SHARED
    #stdint.h   # not needed for Visual Studio 2010 and later
    mexInterface_AlgPDTree_MLP_Mesh.cpp
    mexInterface_AlgPDTree_MLP_Mesh.h
    mexInterface_AlgPDTree_MLP_Mesh.def    
    ${MEX_BASE_FILES}
    )  
  include_directories(
    ${mexAlg_PDTree_MLP_Mesh_SOURCE_DIR}
    ${cisstICP_LIB_INCLUDE}
    )  
  target_link_libraries (mexAlgPDTree_MLP_Mesh
    ${MATLAB_ENGINE_LIBRARIES}
    ${cisstICP_LIB}
    ${WM5_MATH_LIB}
    ${WM5_CORE_LIB}
    )
  # IF (X86_MODE)
    # target_link_libraries( mexAlgPDTree_MLP_Mesh 
      # ${GSL_LIB}
      # ${GSL_CBLAS_LIB}
      # )
  # ENDIF (X86_MODE)
  cisst_target_link_libraries ( mexAlgPDTree_MLP_Mesh
    ${REQUIRED_CISST_LIBRARIES}
    )
  IF (X86_MODE)
    MESSAGE("Configuring mex as 32-bit library (mexw32)...\n")
    MESSAGE("32-bit mex file will not run on 64-bit Matlab!\n")
    SET_TARGET_PROPERTIES( mexAlgPDTree_MLP_Mesh
      PROPERTIES
      PREFIX ""
      SUFFIX ".mexw32"
      )   
  ELSE (X86_MODE)
    MESSAGE("Configuring mex as 64-bit library (mexw64)...\n")
    SET_TARGET_PROPERTIES( mexAlgPDTree_MLP_Mesh
      PROPERTIES
      PREFIX ""
      SUFFIX ".mexw64"
      )  
  ENDIF (X86_MODE)      
  
  
  # ---- mexAlgPDTree_CP_Mesh ---- #
  
  add_library( mexAlgPDTree_CP_Mesh SHARED
    #stdint.h   # not needed for Visual Studio 2010 and later
    mexInterface_AlgPDTree_CP_Mesh.cpp
    mexInterface_AlgPDTree_CP_Mesh.h
    mexInterface_AlgPDTree_CP_Mesh.def    
    ${MEX_BASE_FILES}
    )  
  include_directories(
    ${mexAlg_PDTree_CP_Mesh_SOURCE_DIR}
    ${cisstICP_LIB_INCLUDE}
    )  
  target_link_libraries (mexAlgPDTree_CP_Mesh
    ${MATLAB_ENGINE_LIBRARIES}
    ${cisstICP_LIB}
    ${WM5_MATH_LIB}
    ${WM5_CORE_LIB}
    )
  # IF (X86_MODE)
    # target_link_libraries( mexAlgPDTree_CP_Mesh 
      # ${GSL_LIB}
      # ${GSL_CBLAS_LIB}
      # )
  # ENDIF (X86_MODE)
  cisst_target_link_libraries ( mexAlgPDTree_CP_Mesh
    ${REQUIRED_CISST_LIBRARIES}
    )
  IF (X86_MODE)
    MESSAGE("Configuring mex as 32-bit library (mexw32)...\n")
    MESSAGE("32-bit mex file will not run on 64-bit Matlab!\n")
    SET_TARGET_PROPERTIES( mexAlgPDTree_CP_Mesh
      PROPERTIES
      PREFIX ""
      SUFFIX ".mexw32"
      )   
  ELSE (X86_MODE)
    MESSAGE("Configuring mex as 64-bit library (mexw64)...\n")
    SET_TARGET_PROPERTIES( mexAlgPDTree_CP_Mesh
      PROPERTIES
      PREFIX ""
      SUFFIX ".mexw64"
      )  
  ENDIF (X86_MODE)     
  
  
  # ---- mexAlg2D_DirPDTree_CP_Edges ---- #
  
  add_library( mexAlg2D_DirPDTree_CP_Edges SHARED
    #stdint.h   # not needed for Visual Studio 2010 and later
    mexInterface_Alg2D_DirPDTree_CP_Edges.cpp
    mexInterface_Alg2D_DirPDTree_CP_Edges.h
    mexInterface_Alg2D_DirPDTree_CP_Edges.def    
    ${MEX_BASE_FILES}
    )  
  include_directories(
    ${mexAlg2D_DirPDTree_CP_Edges_SOURCE_DIR}
    ${cisstICP_LIB_INCLUDE}
    )  
  target_link_libraries (mexAlg2D_DirPDTree_CP_Edges
    ${MATLAB_ENGINE_LIBRARIES}
    ${cisstICP_LIB}
    ${WM5_MATH_LIB}
    ${WM5_CORE_LIB}  
    )
  # IF (X86_MODE)
    # target_link_libraries( mexAlg2D_DirPDTree_CP_Edges 
      # ${GSL_LIB}
      # ${GSL_CBLAS_LIB}
      # )
  # ENDIF (X86_MODE)
  cisst_target_link_libraries ( mexAlg2D_DirPDTree_CP_Edges
    ${REQUIRED_CISST_LIBRARIES}
    )
  IF (X86_MODE)
    MESSAGE("Configuring mex as 32-bit library (mexw32)...\n")
    MESSAGE("32-bit mex file will not run on 64-bit Matlab!\n")
    SET_TARGET_PROPERTIES( mexAlg2D_DirPDTree_CP_Edges
      PROPERTIES
      PREFIX ""
      SUFFIX ".mexw32"
      )   
  ELSE (X86_MODE)
    MESSAGE("Configuring mex as 64-bit library (mexw64)...\n")
    SET_TARGET_PROPERTIES( mexAlg2D_DirPDTree_CP_Edges
      PROPERTIES
      PREFIX ""
      SUFFIX ".mexw64"
      )  
  ENDIF (X86_MODE)    
  
  
  # ---- mexAlg2D_DirPDTree_CP_PointCloud ---- #
  
  add_library( mexAlg2D_DirPDTree_CP_PointCloud SHARED
    #stdint.h   # not needed for Visual Studio 2010 and later
    mexInterface_Alg2D_DirPDTree_CP_PointCloud.cpp
    mexInterface_Alg2D_DirPDTree_CP_PointCloud.h
    mexInterface_Alg2D_DirPDTree_CP_PointCloud.def    
    ${MEX_BASE_FILES}
    )  
  include_directories(
    ${mexAlg2D_DirPDTree_CP_PointCloud_SOURCE_DIR}
    ${cisstICP_LIB_INCLUDE}
    )  
  target_link_libraries (mexAlg2D_DirPDTree_CP_PointCloud
    ${MATLAB_ENGINE_LIBRARIES}
    ${cisstICP_LIB}
    ${WM5_MATH_LIB}
    ${WM5_CORE_LIB}  
    )
  # IF (X86_MODE)
    # target_link_libraries( mexAlg2D_DirPDTree_CP_PointCloud
      # ${GSL_LIB}
      # ${GSL_CBLAS_LIB}
      # )
  # ENDIF (X86_MODE)
  cisst_target_link_libraries ( mexAlg2D_DirPDTree_CP_PointCloud
    ${REQUIRED_CISST_LIBRARIES}
    )
  IF (X86_MODE)
    MESSAGE("Configuring mex as 32-bit library (mexw32)...\n")
    MESSAGE("32-bit mex file will not run on 64-bit Matlab!\n")
    SET_TARGET_PROPERTIES( mexAlg2D_DirPDTree_CP_PointCloud
      PROPERTIES
      PREFIX ""
      SUFFIX ".mexw32"
      )   
  ELSE (X86_MODE)
    MESSAGE("Configuring mex as 64-bit library (mexw64)...\n")
    SET_TARGET_PROPERTIES( mexAlg2D_DirPDTree_CP_PointCloud
      PROPERTIES
      PREFIX ""
      SUFFIX ".mexw64"
      )  
  ENDIF (X86_MODE)      
  
  
  # ---- mexAlg2D_DirPDTree_vonMises_Edges ---- #
  
  add_library( mexAlg2D_DirPDTree_vonMises_Edges SHARED
    #stdint.h   # not needed for Visual Studio 2010 and later
    mexInterface_Alg2D_DirPDTree_vonMises_Edges.cpp
    mexInterface_Alg2D_DirPDTree_vonMises_Edges.h
    mexInterface_Alg2D_DirPDTree_vonMises_Edges.def    
    ${MEX_BASE_FILES}
    )  
  include_directories(
    ${mexAlg2D_DirPDTree_vonMises_Edges_SOURCE_DIR}
    ${cisstICP_LIB_INCLUDE}
    )  
  target_link_libraries (mexAlg2D_DirPDTree_vonMises_Edges
    ${MATLAB_ENGINE_LIBRARIES}
    ${cisstICP_LIB}
    ${WM5_MATH_LIB}
    ${WM5_CORE_LIB}  
    )
  # IF (X86_MODE)
    # target_link_libraries( mexAlg2D_DirPDTree_vonMises_Edges 
      # ${GSL_LIB}
      # ${GSL_CBLAS_LIB}
      # )
  # ENDIF (X86_MODE)
  cisst_target_link_libraries ( mexAlg2D_DirPDTree_vonMises_Edges
    ${REQUIRED_CISST_LIBRARIES}
    )
  IF (X86_MODE)
    MESSAGE("Configuring mex as 32-bit library (mexw32)...\n")
    MESSAGE("32-bit mex file will not run on 64-bit Matlab!\n")
    SET_TARGET_PROPERTIES( mexAlg2D_DirPDTree_vonMises_Edges
      PROPERTIES
      PREFIX ""
      SUFFIX ".mexw32"
      )   
  ELSE (X86_MODE)
    MESSAGE("Configuring mex as 64-bit library (mexw64)...\n")
    SET_TARGET_PROPERTIES( mexAlg2D_DirPDTree_vonMises_Edges
      PROPERTIES
      PREFIX ""
      SUFFIX ".mexw64"
      )  
  ENDIF (X86_MODE)      
  
  # ---- mexAlg2D_DirPDTree_vonMises_PointCloud ---- #
  
  add_library( mexAlg2D_DirPDTree_vonMises_PointCloud SHARED
    #stdint.h   # not needed for Visual Studio 2010 and later
    mexInterface_Alg2D_DirPDTree_vonMises_PointCloud.cpp
    mexInterface_Alg2D_DirPDTree_vonMises_PointCloud.h
    mexInterface_Alg2D_DirPDTree_vonMises_PointCloud.def    
    ${MEX_BASE_FILES}
    )  
  include_directories(
    ${mexAlg2D_DirPDTree_vonMises_PointCloud_SOURCE_DIR}
    ${cisstICP_LIB_INCLUDE}
    )  
  target_link_libraries (mexAlg2D_DirPDTree_vonMises_PointCloud
    ${MATLAB_ENGINE_LIBRARIES}
    ${cisstICP_LIB}
    ${WM5_MATH_LIB}
    ${WM5_CORE_LIB}  
    )
  # IF (X86_MODE)
    # target_link_libraries( mexAlg2D_DirPDTree_vonMises_Edges 
      # ${GSL_LIB}
      # ${GSL_CBLAS_LIB}
      # )
  # ENDIF (X86_MODE)
  cisst_target_link_libraries ( mexAlg2D_DirPDTree_vonMises_PointCloud
    ${REQUIRED_CISST_LIBRARIES}
    )
  IF (X86_MODE)
    MESSAGE("Configuring mex as 32-bit library (mexw32)...\n")
    MESSAGE("32-bit mex file will not run on 64-bit Matlab!\n")
    SET_TARGET_PROPERTIES( mexAlg2D_DirPDTree_vonMises_PointCloud
      PROPERTIES
      PREFIX ""
      SUFFIX ".mexw32"
      )   
  ELSE (X86_MODE)
    MESSAGE("Configuring mex as 64-bit library (mexw64)...\n")
    SET_TARGET_PROPERTIES( mexAlg2D_DirPDTree_vonMises_PointCloud
      PROPERTIES
      PREFIX ""
      SUFFIX ".mexw64"
      )  
  ENDIF (X86_MODE)       
  
  
  # # ---- mexAlg2D_DirICP_StdICP_Edges ---- #
  
  # add_library( mexAlg2D_DirICP_StdICP_Edges SHARED
    # #stdint.h   # not needed for Visual Studio 2010 and later
    # alg2D_DirICP_StdICP_Edges_MexInterface.cpp
    # alg2D_DirICP_StdICP_Edges_MexInterface.h
    # alg2D_DirICP_StdICP_Edges_MexInterface.def    
    # ${MEX_BASE_FILES}
    # )  
  # include_directories(
    # ${mexAlg2D_DirICP_StdICP_Edge_SOURCE_DIR}
    # ${cisstICP_LIB_INCLUDE}
    # )  
  # target_link_libraries (mexAlg2D_DirICP_StdICP_Edges
    # ${MATLAB_ENGINE_LIBRARIES}
    # ${cisstICP_LIB}
    # # #${WM5_MATH_LIB_DEBUG} # only need this if compiling in debug
    # # #${WM5_MATH_LIB_DEBUG}    
    # # ${WM5_MATH_LIB} # only need this if compiling in debug
    # # ${WM5_CORE_LIB}  
    # )
  # # IF (X86_MODE)
    # # target_link_libraries( mexAlg2D_DirICP_StdICP_Edges 
      # # ${GSL_LIB}
      # # ${GSL_CBLAS_LIB}
      # # )
  # # ENDIF (X86_MODE)
  # cisst_target_link_libraries ( mexAlg2D_DirICP_StdICP_Edges
    # ${REQUIRED_CISST_LIBRARIES}
    # )
  # IF (X86_MODE)
    # MESSAGE("Configuring mex as 32-bit library (mexw32)...\n")
    # MESSAGE("32-bit mex file will not run on 64-bit Matlab!\n")
    # SET_TARGET_PROPERTIES( mexAlg2D_DirICP_StdICP_Edges
      # PROPERTIES
      # PREFIX ""
      # SUFFIX ".mexw32"
      # )   
  # ELSE (X86_MODE)
    # MESSAGE("Configuring mex as 64-bit library (mexw64)...\n")
    # SET_TARGET_PROPERTIES( mexAlg2D_DirICP_StdICP_Edges
      # PROPERTIES
      # PREFIX ""
      # SUFFIX ".mexw64"
      # )  
  # ENDIF (X86_MODE)  
  
  
else (cisst_FOUND_AS_REQUIRED)
  message ("Information: code in ${CMAKE_CURRENT_SOURCE_DIR} will not be compiled, it requires ${REQUIRED_CISST_LIBRARIES}")
endif (cisst_FOUND_AS_REQUIRED)			
