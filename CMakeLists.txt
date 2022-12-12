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

# To profile the release version of code, a flag '-g' has to be added. This option will take care
# of that. To use, invoke cmake with '-DPROFILE=ON'. This option works irrespective of the build
# type, although it is redundant for a debug build (since dealii adds related flags automatically).
OPTION(PROFILE "Add '-g' flag to enable profiling" OFF)
IF(PROFILE)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
ENDIF(PROFILE)

# Switch to RANS SA turbulent solver
OPTION(RANS_SA "Use RANS Spalart-Allmaras turbulence modelling" OFF)
IF(RANS_SA)
    ADD_DEFINITIONS(-DRANS_SA)
ENDIF(RANS_SA)

# Use conservative variable gradients for stress and heat flux calculation
OPTION(
    CVAR_GRAD
    "Use conservative variable gradients for shear stress and heat flux calculation"
    OFF
)
IF(CVAR_GRAD)
    ADD_DEFINITIONS(-DCVAR_GRAD)
ENDIF(CVAR_GRAD)

# For 2d simulations
OPTION(
    PLENS_2D
    "Manually set RHS of z-direction momentum to 0. Useful for 2d order of convergence tests"
    OFF
)
IF(PLENS_2D)
    ADD_DEFINITIONS(-DPLENS_2D)
ENDIF(PLENS_2D)

# Constant viscosity
OPTION(
    CONST_VISCOSITY
    "Use constant viscosity for simulations"
)
IF(CONST_VISCOSITY)
    ADD_DEFINITIONS(-DCONST_VISCOSITY)
ENDIF(CONST_VISCOSITY)

# For using Persson's indicator, without Hennemenn's modification
OPTION(
    PERSSON_INDICATOR
    "Use purely Persson's version of trouble energy calculation, without Hennemenn's modification"
    OFF
)
IF(PERSSON_INDICATOR)
    ADD_DEFINITIONS(-DPERSSON_INDICATOR)
ENDIF(PERSSON_INDICATOR)

# SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14") # required for 9.3.0

ADD_SUBDIRECTORY(utilities)
ADD_SUBDIRECTORY(manifolds)
ADD_SUBDIRECTORY(modelling)
ADD_SUBDIRECTORY(IBCs)
ADD_SUBDIRECTORY(dgsem)

ADD_EXECUTABLE(plens.out
    main.cc
)

# Gets git information
# This is only executed when the ".h.in" file is updated
# A hack to do that is to call "touch utilities/plens_git_revision.h.in"
# This will change the timestamp and force cmake to rerun
# See https://stackoverflow.com/questions/50102893/rebuild-configuration-file-when-there-are-no-cmake-changes
DEAL_II_QUERY_GIT_INFORMATION(PLENS)
IF(GIT_FOUND AND EXISTS ${CMAKE_SOURCE_DIR}/.git/HEAD)
    # git exists on the system and this directory is a git repo
    MESSAGE("-- plens-branch:       " ${PLENS_GIT_BRANCH})
    MESSAGE("-- plens-hash (full):  " ${PLENS_GIT_REVISION})
    MESSAGE("-- plens-hash (short): " ${PLENS_GIT_SHORTREV})
    cmake_policy(SET CMP0053 OLD)
    SET(PLENS_GIT_BRANCH "@PLENS_GIT_BRANCH@")
    SET(PLENS_GIT_REVISION "@PLENS_GIT_REVISION@")
    SET(PLENS_GIT_SHORTREV "@PLENS_GIT_SHORTREV@")
    CONFIGURE_FILE(
        ${CMAKE_SOURCE_DIR}/utilities/plens_git_revision.h.in
        ${CMAKE_SOURCE_DIR}/utilities/plens_git_revision.h
    )
ENDIF()
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
