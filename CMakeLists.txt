#
# mkdir build && mkdir build
# cmake -DCMAKE_INSTALL_PREFIX=<your-favorite-installation-directory> -Wno-dev ..
# make install
#

PROJECT(BMF)

# Install under /usr/local per default;
IF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  SET(CMAKE_INSTALL_PREFIX "/usr/local" CACHE PATH "Installation directory" FORCE)
ENDIF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

# Yes, help my buggy assert() statements; FIXME: do it better later;
IF( NOT CMAKE_BUILD_TYPE )
   SET( CMAKE_BUILD_TYPE Debug)
ENDIF()

# Should suffice?;
cmake_minimum_required(VERSION 2.8.10)

if(NOT DEFINED CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 11)
endif()

SET( CMAKE_CXX_FLAGS "-g -Wall" )

include_directories( ${PROJECT_SOURCE_DIR}/include )

SET(BMF_LIB bmf)
SET(BMF_EXE bmf-main)

SET( BMF_SRC  source/BeastMagneticField.cc )
SET( MAIN_SRC source/main.cc )

# Well, ROOT is optional;
IF( DEFINED ENV{ROOTSYS})
  FIND_PACKAGE(ROOT 6.14.0)

  #---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
  include(${ROOT_USE_FILE})

  #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_BMF_ROOT_")
  SET(BMF_ROOTMAP ${CMAKE_CURRENT_BINARY_DIR}/libbmf_rdict.pcm ${CMAKE_CURRENT_BINARY_DIR}/libbmf.rootmap)

  SET(BMF_SRC ${BMF_SRC} G__bmf.cxx) 

  ROOT_GENERATE_DICTIONARY(G__bmf BeastMagneticField.h LINKDEF include/bmfLinkDef.h)
endif()

ADD_LIBRARY( ${BMF_LIB} SHARED ${BMF_SRC} )
# ${ROOT_LIBRARIES} will be a blank string if building with no ROOT; fine;
TARGET_LINK_LIBRARIES(${BMF_LIB} ${ROOT_LIBRARIES} )

ADD_EXECUTABLE( ${BMF_EXE} ${MAIN_SRC} )
TARGET_LINK_LIBRARIES( ${BMF_EXE} ${BMF_LIB} )

INSTALL(FILES   ${PROJECT_SOURCE_DIR}/include/BeastMagneticField.h DESTINATION include)
IF( DEFINED ENV{ROOTSYS})
  INSTALL(FILES ${BMF_ROOTMAP}                                     DESTINATION lib)
endif()
INSTALL(TARGETS ${BMF_LIB}                                         DESTINATION lib)
INSTALL(TARGETS ${BMF_EXE}                                         DESTINATION bin)
