/**
 * @file split_string.h
 * @brief A c++ parallel of python's split function
 */

#include <string>
#include <vector>

#ifdef DEBUG
#include <iostream>
#include <testing.h>
#endif

#ifndef SPLIT_STRING_H
#define SPLIT_STRING_H

namespace utilities{

/**
 * @brief Splits a given string based on a given delimiter
 * 
 * See https://stackoverflow.com/questions/14265581/parse-split-a-string-in-c-using-string-delimiter-standard-c
 * @param[in] s The string to be split
 * @param[in] delim The delimiter to be used for splitting
 * @param[out] splits A vector of splitted sub-strings. This variable will be overwritten here and
 * hence must preferably be empty.
 * 
 * @warning No checks of any kind are done. This is purely meant to be a supporting function.
 */
inline void split_string(
    const std::string &s,
    const std::string &delim,
    std::vector<std::string> &splits
)
{
    // clear the vector
    splits.clear();
    
    std::size_t pos;
    std::string token;
    std::string s_(s);
    while( (pos = s_.find(delim)) != std::string::npos ){
        token = s_.substr(0, pos);
        if(token.size() > 0) splits.emplace_back(token); // add non-empty strings
        s_.erase(0, pos + delim.length());
    }
    if(s_.size() > 0) splits.emplace_back(s_); // add final remaining part if non-empty
}

#ifdef DEBUG
inline void split_string_test()
{
    std::string test_string = "My name  is vac ha n  potluri? yo@#$.><{}[]()\/-+*";
    std::vector<std::string> splits;
    split_string(test_string, " ", splits);
    Testing t("split_string", "function");
    std::cout << "Input string: '" << test_string << "'\n";
    std::cout << "Obtained splits:\n";
    for(auto str: splits) std::cout << "'" << str << "'\n";
}
#endif

} // namespace utilities

#endif
