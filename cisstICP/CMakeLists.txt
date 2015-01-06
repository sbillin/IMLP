# See license at http://www.cisst.org/cisst/license.txt

cmake_minimum_required(VERSION 2.8)
project (cisstICP)

# NOTES:
#
# To get the last component of a path, use
# get_filename_component(... NAME)
# To get the parent directory, use
# get_filename_component(... PATH)

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
  include (${CISST_USE_FILE})

  # WildMagic5 Mathematics Library
  # Release
  find_path (WM5_BASE_DIR
    NAMES wm5path.txt
    PATHS "C:/software/WildMagic5"
    )
  find_path (WM5_LIB_PATH
    NAMES Wm5Core.lib
    #PATHS "C:/software/WildMagic5/SDK/Library/v120/x64/Release"
    HINTS ${WM5_BASE_DIR}"/SDK/Library/v120/x64/Release"
    )
  find_library( WM5_CORE_LIB Wm5Core.lib PATHS ${WM5_LIB_PATH})
  find_library( WM5_MATH_LIB Wm5Mathematics.lib PATHS ${WM5_LIB_PATH})
  # # Debug:
  # find_library( WM5MATH_LIB Wm5MathematicsD.lib Wm5CoreD.lib
    # HINTS "C:/software/WildMagic5/SDK/Library/v120/x64/Debug"
    # )
  find_path( WM5_LIB_INCLUDE
    NAMES "Wm5NoniterativeEigen3x3.h"
    PATHS
    HINTS ${WM5_BASE_DIR}"/SDK/Include"
    #HINTS "C:/software/WildMagic5/SDK/Include"
    DOC  ""
    )
  include_directories(${WM5_LIB_INCLUDE})  

  
  #=== ICP Library ===

  set (cisstICP_FILES
    cisstICP.cpp
    cisstICP.h
    cisstCovTreeBase.cpp
    cisstCovTreeBase.h
    cisstCovTreeNode.cpp
    cisstCovTreeNode.h
    cisstCovTree_Mesh.cpp
    cisstCovTree_Mesh.h
    cisstCovTree_PointCloud.cpp
    cisstCovTree_PointCloud.h
    cisstAlgorithmCovTree.cpp
    cisstAlgorithmCovTree.h
    cisstAlgorithmCovTree_CP.cpp
    cisstAlgorithmCovTree_CP.h
    cisstAlgorithmCovTree_CP_Mesh.cpp
    cisstAlgorithmCovTree_CP_Mesh.h
    cisstAlgorithmCovTree_CP_PointCloud.cpp
    cisstAlgorithmCovTree_CP_PointCloud.h
    cisstAlgorithmICP.cpp
    cisstAlgorithmICP.h
    cisstAlgorithmICP_StdICP.cpp
    cisstAlgorithmICP_StdICP.h
    cisstAlgorithmICP_StdICP_Mesh.h
    cisstAlgorithmICP_StdICP_PointCloud.h
    cisstAlgorithmICP_RobustICP.cpp
    cisstAlgorithmICP_RobustICP.h
    cisstAlgorithmICP_RobustICP_Mesh.h
    cisstAlgorithmICP_RobustICP_PointCloud.h
    cisstAlgorithmICP_IMLP.cpp
    cisstAlgorithmICP_IMLP.h
    cisstAlgorithmICP_IMLP_Mesh.cpp
    cisstAlgorithmICP_IMLP_Mesh.h
    cisstAlgorithmICP_IMLP_PointCloud.cpp
    cisstAlgorithmICP_IMLP_PointCloud.h
    cisstAlgorithmICP_IMLP_MahalDist.cpp
    cisstAlgorithmICP_IMLP_MahalDist.h 
    cisstAlgorithmICP_IMLP_MahalDist_PointCloud.cpp
    cisstAlgorithmICP_IMLP_MahalDist_PointCloud.h
    cisstAlgorithmICP_IMLP_MahalDist_Mesh.cpp
    cisstAlgorithmICP_IMLP_MahalDist_Mesh.h
    cisstAlgorithmICP_IMLP_ClosestPoint.cpp
    cisstAlgorithmICP_IMLP_ClosestPoint.h 
    cisstAlgorithmICP_IMLP_ClosestPoint_PointCloud.cpp
    cisstAlgorithmICP_IMLP_ClosestPoint_PointCloud.h
    cisstAlgorithmICP_IMLP_ClosestPoint_Mesh.cpp
    cisstAlgorithmICP_IMLP_ClosestPoint_Mesh.h
    cisstMesh.cpp
    cisstMesh.h
    cisstTriangle.cpp
    cisstTriangle.h
    cisstTriangleClosestPointSolver.cpp
    cisstTriangleClosestPointSolver.h
    cisstBoundingBox.cpp
    cisstBoundingBox.h
    cisstPointCloud.h
    cisstPointCloud.cpp
    cisstPointCloudDir.h
    cisstPointCloudDir.cpp
    utilities.h
    utilities.cpp
    RegisterP2P.cpp
    RegisterP2P.h
    Ellipsoid_OBB_Intersection_Solver.cpp
    Ellipsoid_OBB_Intersection_Solver.h
    )

  # cisstICP Library
  add_library (cisstICP
    ${cisstICP_FILES}    
    )
  include_directories(
    ${cisstICP_SOURCE_DIR}
    ${cisstICP_SOURCE_DIR}/NumericalRecipes
    )
  target_link_libraries( cisstICP
    ${WM5_MATH_LIB}
    ${WM5_CORE_LIB}
    )
  cisst_target_link_libraries ( cisstICP ${REQUIRED_CISST_LIBRARIES})

else (cisst_FOUND_AS_REQUIRED)
  message ("Information: code in ${CMAKE_CURRENT_SOURCE_DIR} will not be compiled, it requires      ${REQUIRED_CISST_LIBRARIES}")
endif (cisst_FOUND_AS_REQUIRED)
