/**
 * @file rk4_stage5_register3.cc
 * @brief A class for 5 stage, 3 register RK4 time integration.
 */

#include "rk4_stage5_register3.h"

/**
 * Constructor. Sets the coefficients according to RK4(3)5[3R+]N scheme from Kennedy et al (2000).
 * See WJ-09-Jun-2021.
 */
RK4Stage5Register3::RK4Stage5Register3()
{
    b = {{
        314199625218.0/7198350928319,
        6410344372641.0/17000082738695,
        292278564125.0/5593752632744,
        5010207514426.0/21876007855139,
        5597675544274.0/18784428342765
    }};

    a_outer = {{
        4745337637855.0/22386579876409,
        6808157035527.0/13197844641179,
        4367509502613.0/10454198590847,
        1236962429870.0/3429868089329,
        b[4]
    }};

    a_inner = {{
        546509042554.0/9152262712923,
        625707605167.0/5316659119056,
        582400652113.0/7078426004906,
        b[3]
    }};
}
