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

SET( CMAKE_CXX_FLAGS "-g -Wall -std=c++11" )

include_directories( ${PROJECT_SOURCE_DIR}/include )

SET(BMF_LIB bmf)

SET( BMF_SRC  source/BeastMagneticField.cc )
SET( MAIN_SRC source/main.cc )

# No dependencies whatsoever for the library itself, just ADD_LIBRARY;
ADD_LIBRARY( ${BMF_LIB} SHARED ${BMF_SRC} )

ADD_EXECUTABLE( bmf-main ${MAIN_SRC} )
TARGET_LINK_LIBRARIES( bmf-main ${BMF_LIB} )

INSTALL(FILES   ${PROJECT_SOURCE_DIR}/include/BeastMagneticField.h DESTINATION include)
INSTALL(TARGETS ${BMF_LIB}                                         DESTINATION lib)
