ADD_LIBRARY(manifolds SHARED
    manifold_description.cc
    cylinder.cc
    nose_cylinder.cc
)

TARGET_INCLUDE_DIRECTORIES(manifolds PUBLIC "..")

DEAL_II_SETUP_TARGET(manifolds)

TARGET_LINK_LIBRARIES(manifolds utilities)

SET(CMAKE_EXPORT_COMPILE_COMMANDS ON)

