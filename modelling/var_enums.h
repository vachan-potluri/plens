/**
 * @file var_enums.h
 * @brief Defines certain enums for variables in NS equations
 * 
 * These enums may not be actually used in the code. But they serve as reference for ordering 
 * convention of quantities like in a conservative state. The names are useful during data output.
 * 
 * The ordering of enum avar is as follows, considering the stresses represented by a matrix.
@verbatim
    0  1  2
           
    *  3  4
           
    *  *  5
@endverbatim
 */

#include <array>
#include <string>

#ifndef VAR_ENUMS_H
#define VAR_ENUMS_H

enum cvar {rho, rhou, rhov, rhow, rhoE};
// list of cvars for looping
const std::array<cvar, 5> cvar_list = {cvar::rho, cvar::rhou, cvar::rhov, cvar::rhow, cvar::rhoE};
// cvar names, for data output
const std::array<std::string, 5> cvar_names = {"rho", "rhou", "rhov", "rhow", "rhoE"};

enum pvar {p, u, v, w, T};
const std::array<pvar, 5> pvar_list = {pvar::p, pvar::u, pvar::v, pvar::w, pvar::T};
const std::array<std::string, 5> pvar_names = {"p", "u", "v", "w", "T"};

enum avar {txx, txy, txz, tyy, tyz, tzz, qx, qy, qz};
const std::array<avar, 9> avar_list = {
    avar::txx, avar::txy, avar::txz, avar::tyy, avar::tyz, avar::tzz, avar::qx, avar::qy, avar::qz
};
const std::array<std::string, 9> avar_names = {
    "txx", "txy", "txz", "tyy", "tyz", "tzz", "qx", "qy", "qz"
};

#endif

