/**
 * @file plens_test.h
 * @brief A class to test PLENS class.
 */

#ifndef PLENS_TEST_H
#define PLENS_TEST_H

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include <deal.II/base/exceptions.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/fe/mapping_q_generic.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/base/utilities.h>

#include "plens.h"
#include <utilities/testing.h>
#include <modelling/navier_stokes.h>
#include <modelling/var_enums.h>

using namespace dealii;
/**
 * @class plens_test
 * @brief A class to test PLENS class.
 */
class plens_test
{
    private:
    utilities::Testing t;

    public:
    plens_test();
    ~plens_test();
    void read_mesh_test() const;
    void set_NS_test() const;
    void set_IC_test() const;
};

#endif
