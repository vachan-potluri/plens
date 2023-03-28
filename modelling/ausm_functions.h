/**
 * @file ausm_functions.h
 * @brief Contains AUSM family related functions. Most functions will not be documented as they are
 * straight forward to understand.
 */

#ifndef AUSM_FUNCTIONS_H
#define AUSM_FUNCTIONS_H

namespace ausm
{

constexpr double Minfty = 0.3;

constexpr double beta = 0.125;

constexpr double Kp = 0.25;

constexpr double Ku = 0.75;

constexpr double sigma = 1.0;

// eq 18 of Liou (2006)
inline double mach_split_1_pos(const double M) {return 0.5*(M+fabs(M));}

inline double mach_split_1_neg(const double M) {return 0.5*(M-fabs(M));}

// eq 19 of Liou (2006)
inline double mach_split_2_pos(const double M) {return 0.25*(M+1)*(M+1);}

inline double mach_split_2_neg(const double M) {return -0.25*(M-1)*(M-1);}

// eq (20) of Liou (2006)
inline double mach_split_4_pos(const double M)
{
    if(fabs(M) >= 1){
        return mach_split_1_pos(M);
    }
    else{
        // expanded form
        // return 0.25*(M+1)*(M+1) + beta*(M*M-1)*(M*M-1);
        // original form
        return mach_split_2_pos(M)*(1-16*beta*mach_split_2_neg(M));
    }
}

inline double mach_split_4_neg(const double M)
{
    if(fabs(M) >= 1){
        return mach_split_1_neg(M);
    }
    else{
        // expanded form
        // return -(0.25*(M-1)*(M-1) + beta*(M*M-1)*(M*M-1));
        // original form
        return mach_split_2_neg(M)*(1+16*beta*mach_split_2_pos(M));
    }
}

// eq (24) of Liou (2006)
inline double pressure_split_5_pos(const double M, const double alpha)
{
    if(fabs(M) >= 1){
        return mach_split_1_pos(M)/M;
    }
    else{
        // expanded form
        // return 0.25*(M+1)*(M+1)*(2-M) + alpha*M*(M*M-1)*(M*M-1);
        // original form
        return mach_split_2_pos(M)*( (2-M) - 16*alpha*M*mach_split_2_neg(M) );
    }
}

inline double pressure_split_5_neg(const double M, const double alpha)
{
    if(fabs(M) >= 1){
        return mach_split_1_neg(M)/M;
    }
    else{
        // expanded form
        // return 0.25*(M-1)*(M-1)*(2+M) - alpha*M*(M*M-1)*(M*M-1);
        // original form
        return mach_split_2_neg(M)*( (-2-M) + 16*alpha*M*mach_split_2_pos(M) );
    }
}

} // namespace ausm

#endif
