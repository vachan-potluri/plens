/**
 * @file testing.h
 * @brief A class that automates some output generation while testing
 */

#ifndef TESTING_H
#define TESTING_H

#include <string>
#include <iostream>

namespace utilities{

/**
 * @class Testing
 * @brief Automates certain output generation
 *
 * The functionalities this provides are
 * 1. Print "Testing xyz class/function" when Testing::begin() is invoked
 * 2. Add some visible separating content when Testing::new_block() is invoked
 *
 * All output is directed to std::cout
 */
class Testing
{
    private:
    std::string name_; // name of the function/class being tested
    std::string type_; // the type (e.g. function/class/namespace)
    
    public:
    Testing(const std::string &name, const std::string &type);
    ~Testing();
    void begin() const;
    void end() const;
    void new_block() const;
    void new_block(const std::string &msg) const;
};

}// namespace utilities

#endif

