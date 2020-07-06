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

# No dependencies whatsoever for the library itself, just ADD_LIBRARY;
ADD_LIBRARY( ${BMF_LIB} SHARED ${BMF_SRC} )

ADD_EXECUTABLE( ${BMF_EXE} ${MAIN_SRC} )
TARGET_LINK_LIBRARIES( ${BMF_EXE} ${BMF_LIB} )

INSTALL(FILES   ${PROJECT_SOURCE_DIR}/include/BeastMagneticField.h DESTINATION include)
INSTALL(TARGETS ${BMF_LIB}                                         DESTINATION ${CMAKE_INSTALL_LIBDIR})
INSTALL(TARGETS ${BMF_EXE}                                         DESTINATION ${CMAKE_INSTALL_BINDIR})
