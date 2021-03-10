/**
 * @file cavars.h
 * @brief A data structure for combined storage of @p state and @p avars. 'c' in this name is for
 * conservative and 'a' is for auxiliary.
 */

#ifndef CAVARS_H
#define CAVARS_H

#include <deal.II/base/exceptions.h>

#include <array>

#include "state.h"
#include "avars.h"

#ifdef DEBUG
#include <iostream>
#include "utilities/testing.h"
#include "utilities/printing.h"
#endif

/**
 * @class cavars
 * @brief A container like class for storing conservative variables and auxiliary variables
 * together.
 *
 * The main usage of this class is in viscous flux (surface and volume) calculation. This will also
 * be useful for BC functions. It stores pointers to a @p state variable and a @p avars variable
 * and provides read and write access to their data. Since pointers are stored, the behaviour
 * of this class becomes undefined if the original objects go out of scope. However, in the specific
 * setting of plens, this is not an issue as the usage of this class is only like a wrapper.
 *
 * The read write access is not provided directly, but indirectly through cavars::get_state() and
 * cavars::get_avars(). This is the safest way. And really, I don't think any more fancy
 * functionality is required for this, because all the other classes internally use @p state and
 * @p avars directly.
 */
class cavars
{
    private:
    state* sp_; // pointer to state
    avars* ap_; // pointer to avars
    bool only_state; // true: only pointer to a state is set, false: both pointers set
    
    
    
    /**
     * @brief (Re)Setter for state variable pointer
     */
    void set_state(state *sp)
    {
        sp_ = sp;
    }
    
    
    
    /**
     * @brief (Re)Setter for avars variable pointer
     */
    void set_avars(avars *ap)
    {
        ap_ = ap;
        only_state = false;
    }
    
    public:
    /**
     * @brief Constructor. Assigns pointers cavars::sp_ and cavars::ap_ to @p sp and @p ap.
     * cavars::only_state is set to false (internally, in cavars::set_avars())
     */
    cavars(state *sp, avars *ap)
    {
        set_state(sp);
        set_avars(ap);
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
                "Calling get_avars() without setting avars (in constructor) invalid!"
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
                "Calling get_avars() without setting avars (in constructor) invalid!"
            )
        );
        return *ap_;
    }
    
    
    
    #ifdef DEBUG
    static void test()
    {
        utilities::Testing t("cavars", "class");
        {
            t.new_block("Testing construction (both)");
            state s={1,2,3,4};
            avars a={1,2,3,4,5,6,7,8,9};
            cavars ca(&s, &a);
            std::cout << "OK\n";
        }
        
        {
            t.new_block("Testing construction (only state)");
            state s={1,2,3,4};
            cavars ca(&s);
            std::cout << "OK\n";
        }
        
        {
            t.new_block("Value modification testing");
            state s={1,2,3,4,5};
            avars a={1,2,3,4,5,6,7,8,9};
            cavars ca(&s, &a);
            
            std::cout << "\nStage 1";
            utilities::print_state(ca.get_state());
            utilities::print_avars(ca.get_avars());
            
            s[2] = 300;
            a[6] = 700;
            std::cout << "\nStage 2";
            utilities::print_state(ca.get_state());
            utilities::print_avars(ca.get_avars());
            
            state &s2 = ca.get_state();
            avars &a2 = ca.get_avars();
            s2[1] = -100;
            a2[3] = -400;
            std::cout << "\nStage 3";
            utilities::print_state(ca.get_state());
            utilities::print_avars(ca.get_avars());
        }
        
        {
            t.new_block("Out of scope behaviour");
            state s={1,2,3,4,5};
            cavars ca(&s);
            {
                avars a={1,2,3,4,5,6,7,8,9};
                ca.set_avars(&a);
                std::cout << "avars in scope";
                utilities::print_avars(ca.get_avars());
            }
            std::cout << "avars out of scope";
            utilities::print_avars(ca.get_avars());
        }
        
        {
            t.new_block("Unset avars behaviour");
            state s={1,2,3,4,5};
            cavars ca(&s);
            // utilities::print_avars(ca.get_avars()); // throw exception
        }
        
        {
            t.new_block("Const behaviour test");
            state s={1,2,3,4,5};
            avars a={1,2,3,4,5,6,7,8,9};
            cavars ca(&s, &a);
            // const state &s2 = ca.get_state();
            // s2[0] = 0; // error
            
            const cavars ca2(&s, &a);
            const state &s2 = ca.get_state(); // OK
        }
    }
    #endif
};

#endif

