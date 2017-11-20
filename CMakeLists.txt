CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)

FIND_PACKAGE(deal.II 9.0.0 QUIET
  HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
  )
IF(NOT ${deal.II_FOUND})
  MESSAGE(FATAL_ERROR "\n"
    "*** Could not locate a (sufficiently recent) version of deal.II. ***\n\n"
    "You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake\n"
    "or set an environment variable \"DEAL_II_DIR\" that contains this path."
    )
ENDIF()

#
# Elemental. Set with -DElemental_DIR=/path/to/prefix/CMake/Elemental
#
# use FIND_PACKAGE in "Config" mode to utilize <package>Config.cmake
FIND_PACKAGE(Elemental CONFIG)
IF (Elemental_FOUND)
  MESSAGE("-- Use Elemental from ${Elemental_DIR}")
  INCLUDE_DIRECTORIES(${Elemental_INCLUDE_DIRS})
  TARGET_LINK_LIBRARIES(${TARGETLIB}
  	"${Elemental_LIBRARIES}"
  )
  add_definitions(-DWITH_ELEMENTAL)
ENDIF (Elemental_FOUND)

DEAL_II_INITIALIZE_CACHED_VARIABLES()

SET(PROJECT_NAME "dense_benchmark")
PROJECT(${PROJECT_NAME})

function(ADD_BENCHMARK _name)
  ADD_EXECUTABLE(${_name}
      ${_name}.cc
  )
  DEAL_II_SETUP_TARGET(${_name})
  IF (Elemental_FOUND)
    TARGET_LINK_LIBRARIES(${_name}
      ${Elemental_LIBRARIES}
    )
  ENDIF()
endfunction()

#
# add benchmarks:
#
ADD_BENCHMARK("cholesky")
