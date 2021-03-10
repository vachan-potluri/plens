/**
 * @file testing.h
 * @brief A class that automates some output generation while testing
 */

#include "testing.h"

using namespace utilities;

Testing::Testing(const std::string &name, const std::string &type)
: name_(name), type_(type)
{
    begin();
}

Testing::~Testing()
{
    end();
}

void Testing::begin() const
{
    std::cout << "\n\n\n\n"
    << "# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #\n\n"
    << "Testing " << type_ << " '" << name_ << "'\n\n"
    << "# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #\n\n";
}

void Testing::end() const
{
    std::cout
    << "\n\n+ + + + + + + + + + + + END TESTING BLOCK + + + + + + + + + + + + + + + +\n\n";
}

void Testing::new_block() const
{
    std::cout << "\n\n----------------------------------------------------------------\n";
}

void Testing::new_block(const std::string &msg) const
{
    new_block();
    std::cout << msg << "\n";
}
