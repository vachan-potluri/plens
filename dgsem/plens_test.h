/**
 * @file plens_test.h
 * @brief A class to test PLENS class.
 */

#ifndef PLENS_TEST_H
#define PLENS_TEST_H

#include "plens.h"
#include <utilities/testing.h>

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
};

#endif
