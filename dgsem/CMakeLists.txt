ADD_LIBRARY(dgsem SHARED
    metric_terms.cc
    blender_calculator.cc
    rk4_stage5_register3.cc
    subcell_interpolator.cc
    plens.cc
    # plens_test.cc
)

TARGET_INCLUDE_DIRECTORIES(dgsem PUBLIC "..")

DEAL_II_SETUP_TARGET(dgsem)

TARGET_LINK_LIBRARIES(dgsem utilities manifolds modelling IBCs)

SET(CMAKE_EXPORT_COMPILE_COMMANDS ON)

