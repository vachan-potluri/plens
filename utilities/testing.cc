/**
 * @file testing.cc
 * @brief A class that automates some output generation while testing
 */

#include "testing.h"

using namespace utilities;

/**
 * @brief Constuctor. Calls Testing::begin()
 */
Testing::Testing(const std::string &name, const std::string &type)
: name_(name), type_(type)
{
    begin();
}



/**
 * @brief Destructor. Calls Testing::end()
 */
Testing::~Testing()
{
    end();
}



/**
 * @brief Print some statements to say that a new class/function is being tested.
 *
 * This is generally not invoked separately, but through the constructor.
 */
void Testing::begin() const
{
    std::cout << "\n\n\n\n"
    << "# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #\n\n"
    << "Testing " << name_ << " " << type_ << "\n\n"
    << "# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #\n\n";
}



/**
 * @brief Print some statements to say that testing of the current class/function is finished.
 *
 * This is generally not invoked separately, but through the destructor.
 */
void Testing::end() const
{
    std::cout
    << "\n\n+ + + + + + + + + + + + END TESTING BLOCK + + + + + + + + + + + + + + + +\n\n";
}



/**
 * @brief Add some visible separation content to indicate a new set of functionality of the class/
 * function in question is now being tested.
 */
void Testing::new_block() const
{
    std::cout << "\n\n----------------------------------------------------------------\n";
}



/**
 * @brief Calls Testing::new_block() and prints `msg` in the end.
 */
void Testing::new_block(const std::string &msg) const
{
    new_block();
    std::cout << msg << "\n";
}

