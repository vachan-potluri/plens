/**
 * @file piecewise_function.h
 * @brief Sets initial condition based on piecewise defined function.
 */

#include "piecewise_function.h"

using namespace ICs;

/**
 * Constructor. Calls the base constructor and parses all the functions in the file @p filename.
 */
PiecewiseFunction::PiecewiseFunction(
    const DoFHandler<dim> &dh,
    std::array<LA::MPI::Vector, 5> &gcv,
    const std::string &filename
)
: IC(dh, gcv)
{
    std::ifstream fn_file(filename);
    AssertThrow(
        fn_file.good(),
        StandardExceptions::ExcMessage(
            "Unable to open file provided for reading IC functions."
        )
    );

    // Note the choice: primitive or conservative functions
    std::string line;
    try{
        std::getline(fn_file, line);
    } catch(...){
        AssertThrow(
            false,
            StandardExceptions::ExcMessage(
                "Unable to read first line of the IC file."
            )
        );
    }
    if(line == "c"){
        prim_fns_ = false;
    } else if(line == "p"){
        prim_fns_ = true;
    } else{
        AssertThrow(
            false,
            StandardExceptions::ExcMessage(
                "The first line in the file containing functions must be 'p' or 'c' to indicate "
                "whether the functions are for primitive or conservative variables."
            )
        );
    }

    // get the number of pieces
    for(int dir=0; dir<dim; dir++){
        try{
            std::getline(fn_file, line);
        } catch(...){
            AssertThrow(
                false,
                StandardExceptions::ExcMessage(
                    "Unable to read number of pieces in each direction from IC file."
                )
            );
        }
        np_[dir] = stod(line); // assuming line contains a number
        AssertThrow(
            np_[dir] >= 1,
            StandardExceptions::ExcMessage(
                "The number of pieces provided in each direction must be >= 1."
            )
        );
    }

    // get the interface locations
    std::vector<std::string> splits;
    for(int dir=0; dir<dim; dir++){
        try{
            std::getline(fn_file, line);
        } catch(...){
            AssertThrow(
                false,
                StandardExceptions::ExcMessage(
                    "Unable to read the interface locations in each direction from IC file."
                )
            );
        }
        utilities::split_string(line, " ", splits); // split at spaces
        AssertThrow(
            splits.size() == np_[dir]-1,
            StandardExceptions::ExcMessage(
                "The number of interface locations must be exactly one less than the number of "
                "pieces provided earlier for each direction."
            )
        );
        for(const auto &s: splits) ilocs_[dir].emplace_back(stod(s)); // assuming s is a number
    }

    // resize fpps_. This automatically initialises the pointed objects using default ctor of
    // FunctionParser which is a scalar function ctor
    fpps_.resize(np_[0]*np_[1]*np_[2]);

    // parse the functions
    std::string variables("x,y,z");
    std::map<std::string, double> constants; // empty

    for(usi pid=0; pid<fpps_.size(); pid++){
        for(cvar var: cvar_list){
            try{
                std::getline(fn_file, line);
            } catch(...){
                AssertThrow(
                    false,
                    StandardExceptions::ExcMessage(
                        "Unable to read functions from IC file. Probably insufficient number of "
                        "functions were provided in the file."
                    )
                );
            }
            fpps_[pid][var].reset(new FunctionParser<dim>());
            fpps_[pid][var]->initialize(variables, line, constants);
        }
    }
}



#ifdef DEBUG
void PiecewiseFunction::test()
{
    Testing t("PiecewiseFunction", "class");
    utilities::ICTestData ictd(5,2); // divisions, degree

    // run this test in a folder where there is an IC file
    {
        t.new_block("testing construction");
    }
}
#endif
