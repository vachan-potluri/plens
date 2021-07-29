CMAKE_MINIMUM_REQUIRED(VERSION 2.8.12)
FIND_PACKAGE(deal.II 9.3.0 QUIET
    HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
)
IF(NOT ${deal.II_FOUND})
    MESSAGE(FATAL_ERROR "\n"
    "*** Could not locate a (sufficiently recent) version of deal.II. ***\n\n"
    "You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake\n"
    "or set an environment variable \"DEAL_II_DIR\" that contains this path."
)
ENDIF()

DEAL_II_INITIALIZE_CACHED_VARIABLES()
PROJECT(plens)

# Introduce a preprocessor macro SPACETIME to tweak some settings for compilation on spacetime
# To use this, invoke cmake with '-DBUILD_SPACETIME=0N' and modify the code where required
# See WJ-22-Jan-2021
OPTION(BUILD_SPACETIME "Build the code on spacetime" OFF)
IF(BUILD_SPACETIME)
    ADD_DEFINITIONS(-DSPACETIME) # for older cmake versions
ENDIF(BUILD_SPACETIME)

# SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14") # required for 9.3.0

ADD_SUBDIRECTORY(utilities)
ADD_SUBDIRECTORY(manifolds)
ADD_SUBDIRECTORY(modelling)
ADD_SUBDIRECTORY(IBCs)
ADD_SUBDIRECTORY(dgsem)

ADD_EXECUTABLE(plens.out
    main.cc
)

DEAL_II_SETUP_TARGET(plens.out)

TARGET_LINK_LIBRARIES(plens.out utilities manifolds modelling IBCs dgsem)

# for make debug
ADD_CUSTOM_TARGET(debug
    COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Debug ${CMAKE_SOURCE_DIR}
    COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
    COMMENT "Switch CMAKE_BUILD_TYPE to Debug"
)

# for make release
ADD_CUSTOM_TARGET(release
    COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Release ${CMAKE_SOURCE_DIR}
    COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
    COMMENT "Switch CMAKE_BUILD_TYPE to Release"
)

SET(CMAKE_EXPORT_COMPILE_COMMANDS ON)
