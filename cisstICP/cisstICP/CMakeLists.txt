# See license at http://www.cisst.org/cisst/license.txt

cmake_minimum_required(VERSION 2.8)
project (cisstICP)

# NOTES:
#  To get the last component of a path, use
#   get_filename_component(... NAME)
#  To get the parent directory, use
#   get_filename_component(... PATH)

IF (CMAKE_CL_64 OR CMAKE_GENERATOR MATCHES Win64)
  set( X86_MODE 0 )  # 64-bit compile
ELSE (CMAKE_CL_64 OR CMAKE_GENERATOR MATCHES Win64)
  set( X86_MODE 1 )  # 32-bit compile
ENDIF (CMAKE_CL_64 OR CMAKE_GENERATOR MATCHES Win64)


# create a list of CISST libraries needed for this project
set ( REQUIRED_CISST_LIBRARIES cisstCommon cisstVector cisstNumerical cisstOSAbstraction )

# find cisst and make sure the required libraries have been compiled
find_package (cisst REQUIRED ${REQUIRED_CISST_LIBRARIES})

IF (cisst_FOUND_AS_REQUIRED)

  # load cisst configuration
  include (${CISST_USE_FILE})

  # rply library
  find_path( RPLY_DIR
    NAMES rply.c
    HINTS ../dependencies/rply-1.1.4
    DOC  "path to the rply library directory"
    )
  IF (RPLY_DIR)
    include_directories(${RPLY_DIR})
  ENDIF (RPLY_DIR)  
  
  # dlib library
  #  (this library consists only of header files)
  find_path( DLIB_INCLUDE
    NAMES dlib/dlib_include_path_tutorial.txt
    HINTS ../dependencies/dlib-18.6
    DOC  "path to the directory that CONTAINS the 'dlib' folder"
    )
  IF (DLIB_INCLUDE)
    include_directories(${DLIB_INCLUDE})
  ENDIF (DLIB_INCLUDE)

  # Numerical Recipes
  #  (this code consists only of header files)
  find_path( NUMERICAL_RECIPES_INCLUDE
    NAMES nr3.h mins.h 
    HINTS ../dependencies/NumericalRecipes
    DOC  "path to the directory that contains the Numerical Recipes header files"
    )
  IF (NUMERICAL_RECIPES_INCLUDE)
    include_directories(${NUMERICAL_RECIPES_INCLUDE})  
  ENDIF (NUMERICAL_RECIPES_INCLUDE)

  # WildMagic5 Mathematics Library
  # Release
  find_path (WM5_BASE_DIR
    NAMES wm5path.txt
    HINTS C:/software/WildMagic5
    )
  find_path (WM5_LIB_PATH
    NAMES Wm5Core.lib libWm5Core.a
    HINTS ${WM5_BASE_DIR}/SDK/Library/v120/x64/Release
    )
  find_library( WM5_CORE_LIB Wm5Core PATHS ${WM5_LIB_PATH})
  find_library( WM5_MATH_LIB Wm5Mathematics PATHS ${WM5_LIB_PATH})
  # # Debug:
  # find_library( WM5MATH_LIB Wm5MathematicsD.lib Wm5CoreD.lib
    # HINTS C:/software/WildMagic5/SDK/Library/v120/x64/Debug
    # )
  find_path( WM5_LIB_INCLUDE
    NAMES Wm5NoniterativeEigen3x3.h
    HINTS ${WM5_BASE_DIR}/SDK/Include
    )
  IF (WM5_LIB_INCLUDE)
    include_directories(${WM5_LIB_INCLUDE})  
  ENDIF (WM5_LIB_INCLUDE)
  
  # # GSL: has a non-linear optimization library that may be used as
  # #      an alternative to dlib; the quasi-newton method of
  # #      dlib has better performance, however
  # IF (X86_MODE)
    # # GSL library
    # # (only have pre-compiled 32-bit versions of this library)
    # find_path( GSL_INCLUDE gsl
     # HINTS C:/software/gsl-1.14-bin-msvc9/gsl-1.14/include
    # )
    # find_library( GSL_LIB gsl.lib
     # HINTS C:/software/gsl-1.14-bin-msvc9/gsl-1.14/lib
    # )
    # find_library( GSL_CBLAS_LIB gslcblas.lib
     # HINTS C:/software/gsl-1.14-bin-msvc9/gsl-1.14/lib
    # )
    # include_directories(${GSL_INCLUDE})
  # ENDIF (X86_MODE)
  

  #=== ICP Library ===

  set( cisstICP_FILES
    cisstICP.cpp
    cisstICP.h
    cisstMesh.cpp
    cisstMesh.h
    cisstPointCloud.h
    cisstPointCloud.cpp
    TriangleClosestPointSolver.cpp
    TriangleClosestPointSolver.h    
    BoundingBox.cpp
    BoundingBox.h
    utilities.h
    utilities.cpp
    PointProjectionRoutines.cpp
    PointProjectionRoutines.h 
    RegisterP2P.cpp
    RegisterP2P.h
    Ellipsoid_OBB_Intersection_Solver.cpp
    Ellipsoid_OBB_Intersection_Solver.h
    cisstException.h
    ply_io.cpp
    ply_io.h
    # 3D Registration
    # non-oriented point routines
    PDTreeBase.cpp
    PDTreeBase.h
    PDTreeNode.cpp
    PDTreeNode.h
    PDTree_Mesh.cpp
    PDTree_Mesh.h
    PDTree_PointCloud.cpp
    PDTree_PointCloud.h
    algPDTree.cpp
    algPDTree.h
    algPDTree_CP.cpp
    algPDTree_CP.h
    algPDTree_CP_Mesh.cpp
    algPDTree_CP_Mesh.h
    algPDTree_CP_PointCloud.cpp
    algPDTree_CP_PointCloud.h
    algPDTree_MLP.cpp
    algPDTree_MLP.h
    algPDTree_MLP_Mesh.cpp
    algPDTree_MLP_Mesh.h
    algPDTree_MLP_PointCloud.cpp
    algPDTree_MLP_PointCloud.h    
    algICP.cpp
    algICP.h
    algICP_StdICP.cpp
    algICP_StdICP.h
    algICP_StdICP_Mesh.h
    algICP_StdICP_PointCloud.h
    algICP_IMLP.cpp
    algICP_IMLP.h
    algICP_IMLP_Mesh.cpp
    algICP_IMLP_Mesh.h
    algICP_IMLP_PointCloud.cpp
    algICP_IMLP_PointCloud.h
    # oriented point routines
    DirPDTreeBase.cpp
    DirPDTreeBase.h
    DirPDTreeNode.cpp
    DirPDTreeNode.h
    DirPDTree_Mesh.cpp
    DirPDTree_Mesh.h
    DirPDTree_PointCloud.cpp
    DirPDTree_PointCloud.h
    algDirPDTree.cpp
    algDirPDTree.h
    algDirPDTree_BoundedAngle.cpp
    algDirPDTree_BoundedAngle.h
    algDirPDTree_BoundedAngle_Mesh.cpp    
    algDirPDTree_BoundedAngle_Mesh.h  
    algDirPDTree_BoundedAngle_PointCloud.cpp
    algDirPDTree_BoundedAngle_PointCloud.h
    algDirPDTree_vonMisesPrj.cpp
    algDirPDTree_vonMisesPrj.h
    algDirPDTree_vonMisesPrj_Mesh.cpp
    algDirPDTree_vonMisesPrj_Mesh.h    
    algDirICP.cpp
    algDirICP.h
    algDirICP_StdICP.cpp
    algDirICP_StdICP.h
    algDirICP_StdICP_Mesh.cpp
    algDirICP_StdICP_Mesh.h
    algDirICP_StdICP_PointCloud.cpp
    algDirICP_StdICP_PointCloud.h
    algDirICP_IMLOP.cpp
    algDirICP_IMLOP.h
    algDirICP_IMLOP_Mesh.cpp
    algDirICP_IMLOP_Mesh.h
    # 2D Registration
    DirPDTree2DBase.cpp
    DirPDTree2DBase.h
    DirPDTree2DNode.cpp
    DirPDTree2DNode.h
    DirPDTree2D_Edges.cpp
    DirPDTree2D_Edges.h    
    DirPDTree2D_PointCloud.cpp
    DirPDTree2D_PointCloud.h
    alg2D_DirPDTree.cpp
    alg2D_DirPDTree.h
    alg2D_DirPDTree_CP.cpp
    alg2D_DirPDTree_CP.h
    alg2D_DirPDTree_CP_Edges.cpp
    alg2D_DirPDTree_CP_Edges.h    
    alg2D_DirPDTree_CP_PointCloud.cpp
    alg2D_DirPDTree_CP_PointCloud.h    
    alg2D_DirPDTree_vonMises.cpp
    alg2D_DirPDTree_vonMises.h
    alg2D_DirPDTree_vonMises_Edges.cpp
    alg2D_DirPDTree_vonMises_Edges.h
    alg2D_DirPDTree_vonMises_PointCloud.cpp
    alg2D_DirPDTree_vonMises_PointCloud.h    
    alg2D_DirICP.cpp
    alg2D_DirICP.h
    alg2D_DirICP_StdICP.cpp
    alg2D_DirICP_StdICP.h
    alg2D_DirICP_StdICP_Edges.cpp
    alg2D_DirICP_StdICP_Edges.h
    BoundingBox2D.cpp
    BoundingBox2D.h
    cisstEdge2D.h
    cisstEdgeList2D.cpp
    cisstEdgeList2D.h
    utilities2D.h
    utilities2D.cpp
    cisstPointCloud2D.cpp
    cisstPointCloud2D.h
    )
    
  # Extra Algorithms
  #  these algorithms implement or nearly implement other 
  #  prior works and were used for performance comparison
  set( USE_EXTRA_ALGORITHMS false CACHE BOOL "Enable this option in order to include extra algorithms that were written to emulate prior works" )
  IF (USE_EXTRA_ALGORITHMS)
    set(EXTRA_ALG_DIR ${cisstICP_SOURCE_DIR}/ExtraAlgorithms)
    include_directories(${EXTRA_ALG_DIR})
    
    message("Including the following extra algorithms: RobustICP, IMLP_MahalDist & IMLP_ClosestPoint")
    
    set( cisstICP_FILES
      ${cisstICP_FILES}
      ExtraAlgorithms/algICP_IMLP_MahalDist.cpp
      ExtraAlgorithms/algICP_IMLP_MahalDist.h 
      ExtraAlgorithms/algICP_IMLP_MahalDist_PointCloud.cpp
      ExtraAlgorithms/algICP_IMLP_MahalDist_PointCloud.h
      ExtraAlgorithms/algICP_IMLP_MahalDist_Mesh.cpp
      ExtraAlgorithms/algICP_IMLP_MahalDist_Mesh.h
      ExtraAlgorithms/algICP_IMLP_ClosestPoint.cpp
      ExtraAlgorithms/algICP_IMLP_ClosestPoint.h 
      ExtraAlgorithms/algICP_IMLP_ClosestPoint_PointCloud.cpp
      ExtraAlgorithms/algICP_IMLP_ClosestPoint_PointCloud.h
      ExtraAlgorithms/algICP_IMLP_ClosestPoint_Mesh.cpp
      ExtraAlgorithms/algICP_IMLP_ClosestPoint_Mesh.h
      ExtraAlgorithms/algICP_RobustICP.cpp
      ExtraAlgorithms/algICP_RobustICP.h
      ExtraAlgorithms/algICP_RobustICP_Mesh.h
      ExtraAlgorithms/algICP_RobustICP_PointCloud.h
      )
  ENDIF (USE_EXTRA_ALGORITHMS)
    
  IF (RPLY_DIR)
    set( cisstICP_FILES
      ${cisstICP_FILES}
      ${RPLY_DIR}/rply.c
      #rply.h
      #rplyfile.h
      )
  ENDIF (RPLY_DIR)
    
  # Algorithms requiring dlib (for quasi-Newton optimizer)
  IF (DLIB_INCLUDE)
    set( cisstICP_FILES
      ${cisstICP_FILES}
      wrapper_dlib.cpp
      wrapper_dlib.h
      algDirICP_PIMLOP.cpp
      algDirICP_PIMLOP.h
      algDirICP_PIMLOP_Mesh.cpp
      algDirICP_PIMLOP_Mesh.h
      algDirICP_PIMLOP_dlibWrapper.cpp
      algDirICP_PIMLOP_dlibWrapper.h
      )
  
    # GIMLOP also requires Numerical Recipes (for 1D Brent Search)
    IF (NUMERICAL_RECIPES_INCLUDE)
      set( cisstICP_FILES
        ${cisstICP_FILES}
        algDirICP_GIMLOP.cpp
        algDirICP_GIMLOP.h
        algDirICP_GIMLOP_Mesh.cpp
        algDirICP_GIMLOP_Mesh.h      
        )
    ELSE (NUMERICAL_RECIPES_INCLUDE)
      message(WARNING "The G-IMLOP algorithm is being excluded due to Numerical Recipes unfound")
    ENDIF (NUMERICAL_RECIPES_INCLUDE)  
    
  ELSE (DLIB_INCLUDE)
    message(WARNING "The G-IMLOP and P-IMLOP algorithms are being excluded due to dlib unfound")
  ENDIF (DLIB_INCLUDE)
    
  # IF (X86_MODE)
    # set (cisstICP_FILES
      # ${cisstICP_FILES}
      # wrapper_gsl.cpp
      # wrapper_gsl.h
      # )
  # ENDIF (X86_MODE)

  # cisstICP Library
  add_library (cisstICP
    ${cisstICP_FILES}    
    )
  include_directories( ${cisstICP_SOURCE_DIR} )

  IF (WM5_LIB_PATH)
    target_link_libraries( cisstICP
      ${WM5_MATH_LIB}
      ${WM5_CORE_LIB}
      )
  ELSE (WM5_LIB_PATH)
    message(WARNING "WildMagic5 is not found; must uncomment the line for define REDIRECT_COV_EIGENDECOMP_NONITER_TO_SVD in utilities.cpp prior to compiling in order to redirect the calls for non-iterative covariance eigen-decomposition to the iterative SVD-based alternative.")
  ENDIF (WM5_LIB_PATH)

  # IF (GSL_LIB)
    # IF (X86_MODE)
      # target_link_libraries( cisstICP
        # ${GSL_LIB}
        # ${GSL_CBLAS_LIB}
        # )
    # ENDIF (X86_MODE)
  # ENDIF (GSL_LIB)
  
  cisst_target_link_libraries ( cisstICP ${REQUIRED_CISST_LIBRARIES})
  
ELSE (cisst_FOUND_AS_REQUIRED)
  message ("Information: code in ${CMAKE_CURRENT_SOURCE_DIR} will not be compiled, it requires ${REQUIRED_CISST_LIBRARIES}")
ENDIF (cisst_FOUND_AS_REQUIRED)
