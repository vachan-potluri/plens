/**
 * @file cavars.h
 * @brief A data structure for combined storage of @p state and @p avars. 'c' in this name is for
 * conservative and 'a' is for auxiliary.
 */

#ifndef CAVARS_H
#define CAVARS_H

#include <deal.II/base/exceptions.h>

#include <array>
#include <memory>

#include "state.h"
#include "avars.h"
#include "utilities/testing.h"

#ifdef DEBUG
#include <iostream>
#endif

/**
 * @class cavars
 * @brief A container like class for storing conservative variables and auxiliary variables
 * together.
 *
 * The main usage of this class is in viscous flux (surface and volume) calculation. This will also
 * be useful for BC functions. It stores shared pointers to a @p state variable and a @p avars
 * variable and provides read and write access to their data. Shared pointers are used because in
 * general, a @p state object and a @p avars object are constructed independently and then combined
 * through this class thus reducing the number of arguments required for passing in functions that
 * require both these quantities.
 *
 * The read write access is not provided directly, but indirectly through cavars::get_state and
 * cavars::get_avars. This is the safest way. And really, I don't think any more fancy functionality
 * is required for this, because all the other classes internally use @p state and @p avars
 * directly.
 */
class cavars
{
    private:
    std::shared_ptr<state> sp_; // pointer to state
    std::shared_ptr<avars> ap_; // pointer to avars
    bool only_state; // true: only pointer to a state is set, false: both pointers set
    
    public:
    /**
     * @brief Constructor. Assigns shared pointers cavars::sp_ and cavars::ap_ to @p sp and @p ap.
     * cavars::only_state is set to false
     */
    cavars(state *sp, avars *ap)
    {
        set_state(sp);
        set_avars(ap);
        only_state = false;
    }
    
    
    
    /**
     * @brief Constructor using state alone.
     *
     * cavars::ap_ is left unset. cavars::only_state is set to true
     */
    cavars(state *sp)
    {
        set_state(sp);
        only_state = true;
    }
    
    
    
    /**
     * @brief Desctructor. Releases the ownership of pointers
     */
    ~cavars()
    {
        sp_.reset();
        ap_.reset();
    }
    
    
    
    /**
     * @brief (Re)Setter for state variable pointer
     */
    void set_state(state *sp)
    {
        sp_.reset(sp); // reset requires address
    }
    
    
    
    /**
     * @brief (Re)Setter for avars variable pointer
     */
    void set_avars(avars *ap)
    {
        ap_.reset(ap); // reset requires address
        only_state = false;
    }
    
    
    
    /**
     * @brief Returns reference to state object held by cavars::sp_
     */
    state& get_state()
    {
        return *sp_;
    }
    
    
    
    /**
     * @brief Returns reference to state object held by cavars::sp_ (const version)
     */
    const state& get_state() const
    {
        return *sp_;
    }
    
    
    
    /**
     * @brief Returns reference to avars object held by cavars::ap_
     */
    avars& get_avars()
    {
        AssertThrow(
            !only_state,
            dealii::StandardExceptions::ExcMessage(
                "Calling get_avars() without setting avars invalid!"
            )
        );
        return *ap_;
    }
    
    
    
    /**
     * @brief Returns reference to avars object held by cavars::ap_ (const version)
     */
    const avars& get_avars() const
    {
        AssertThrow(
            !only_state,
            dealii::StandardExceptions::ExcMessage(
                "Calling get_avars() without setting avars invalid!"
            )
        );
        return *ap_;
    }
    
    
    
    #ifdef DEBUG
    static void test()
    {
        utilities::Testing t("cavars", "class");
        {
            t.new_block("Testing construction");
            state s={1,2,3,4};
            avars a={1,2,3,4,5,6,7,8,9};
            cavars ca(&s, &a);
            std::cout << "Hi\n";
        }
    }
    #endif
};

#endif

