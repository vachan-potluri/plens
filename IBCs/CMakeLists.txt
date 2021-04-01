ADD_LIBRARY(IBCs SHARED
    BC.cc
    free.cc
    outflow.cc
    uniform_inflow.cc
    uniform_temp_wall.cc
    symmetry.cc
    periodic.cc
    IC.cc
)

TARGET_INCLUDE_DIRECTORIES(IBCs PUBLIC "..")

DEAL_II_SETUP_TARGET(IBCs)

TARGET_LINK_LIBRARIES(IBCs utilities modelling)

SET(CMAKE_EXPORT_COMPILE_COMMANDS ON)

