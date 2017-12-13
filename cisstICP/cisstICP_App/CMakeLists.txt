# See license at http://www.cisst.org/cisst/license.txt

cmake_minimum_required(VERSION 2.6)
project (ICP_App)

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

  # # Boost
  # find_package(Boost COMPONENTS system filesystem REQUIRED)
  # include_directories(${Boost_INCLUDE_DIRS})
  # link_directories(${Boost_LIBRARY_DIRS})  

 #=== dlib ===
 find_path( dlib_include dlib
   names "dlib/dlib_include_path_tutorial.txt"
   hints "c:/software/dlib-18.6"
   doc  "choose folder that contains the dlib folder"
 )
 include_directories(${dlib_include})
 
#  #=== GSL ===
#  IF (X86_MODE)  
#    # GSL library
#    find_path( GSL_INCLUDE gsl
#     HINTS "C:/software/gsl-1.14-bin-msvc9/gsl-1.14/include"
#    )
#    find_library( GSL_LIB gsl.lib
#     HINTS "C:/software/gsl-1.14-bin-msvc9/gsl-1.14/lib"
#    )
#    find_library( GSL_CBLAS_LIB gslcblas.lib
#     HINTS "C:/software/gsl-1.14-bin-msvc9/gsl-1.14/lib"
#    )    
#    
#    include_directories( 
#      ${GSL_INCLUDE}
#      )    
#  ENDIF (X86_MODE)  
 
  # WildMagic5 Mathematics Library
  find_path (WM5_BASE_DIR
    NAMES wm5path.txt
    PATHS "C:/software/WildMagic5"
    )
  if (IS_DIRECTORY ${WM5_BASE_DIR})
    find_path (WM5_LIB_PATH
      NAMES Wm5Core.lib
      #PATHS "C:/software/WildMagic5/SDK/Library/v120/x64/Release"
      HINTS ${WM5_BASE_DIR}"/SDK/Library/v120/x64/Release"
      )
    # Release:
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
  endif (IS_DIRECTORY ${WM5_BASE_DIR})
  
  #=== ICP LIBRARY ===
  
  find_library(cisstICP_LIB cisstICP)
  find_path(cisstICP_LIB_INCLUDE 
    NAMES cisstICP.h
    )

  #SET(LINK_LIBRARY optimized Foo debug Foo_d)
  #target_link_libraries(MyEXE ${LINK_LIBRARY})
  
  #=== EXECUTABLE PROGRAM ===
  
  # Add a new program
  add_executable(ICP_App
    main.cpp
    utility.h
    utility.cpp
    testICP.h
    testICPNormals.h
  )

  include_directories( 
    ${ICP_App_SOURCE_DIR}
    ${cisstICP_LIB_INCLUDE}
    )
    
  target_link_libraries( ICP_App 
    ${cisstICP_LIB}
    ${Boost_LIBRARIES}
    ${WM5_MATH_LIB}
    ${WM5_CORE_LIB}
    )
      
#  IF (X86_MODE)
#    target_link_libraries( ICP_App 
#      ${GSL_LIB}
#      ${GSL_CBLAS_LIB}
#      )
#  ENDIF (X86_MODE)  

  # CISST_REQUIRES will check that the libraries are compiled
  # and set the correct link options
  cisst_target_link_libraries ( ICP_App 
    ${REQUIRED_CISST_LIBRARIES}
    )  
	
else (cisst_FOUND_AS_REQUIRED)
  message ("Information: code in ${CMAKE_CURRENT_SOURCE_DIR} will not be compiled, it requires ${REQUIRED_CISST_LIBRARIES}")
endif (cisst_FOUND_AS_REQUIRED)			
