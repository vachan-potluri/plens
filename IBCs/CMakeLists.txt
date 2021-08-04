ADD_LIBRARY(IBCs SHARED
    BC.cc
    free.cc
    outflow.cc
    uniform_inflow.cc
    uniform_temp_wall.cc
    symmetry.cc
    empty.cc
    periodic.cc
    IC.cc
    piecewise_function.cc
    from_archive.cc
    from_archive_restart.cc
)

TARGET_INCLUDE_DIRECTORIES(IBCs PUBLIC "..")

DEAL_II_SETUP_TARGET(IBCs)

TARGET_LINK_LIBRARIES(IBCs utilities manifolds modelling)

SET(CMAKE_EXPORT_COMPILE_COMMANDS ON)

