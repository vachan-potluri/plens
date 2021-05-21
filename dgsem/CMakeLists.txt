ADD_LIBRARY(dgsem SHARED
    plens.cc
    metric_terms.cc
    plens_test.cc
)

TARGET_INCLUDE_DIRECTORIES(dgsem PUBLIC "..")

DEAL_II_SETUP_TARGET(dgsem)

TARGET_LINK_LIBRARIES(dgsem utilities modelling IBCs)

SET(CMAKE_EXPORT_COMPILE_COMMANDS ON)

