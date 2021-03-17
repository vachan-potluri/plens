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
 * @class CAvars
 * @brief A container like class for storing conservative variables and auxiliary variables
 * together.
 *
 * The main usage of this class is in viscous flux (surface and volume) calculation. This will also
 * be useful for BC functions. It stores pointers to a @p state variable and a @p avars variable
 * and provides read and write access to their data. Since pointers are stored, the behaviour
 * of this class becomes undefined if the original objects go out of scope. However, in the specific
 * setting of plens, this is not an issue as the usage of this class is only like a wrapper.
 *
 * The read write access is not provided directly, but indirectly through CAvars::get_state() and
 * CAvars::get_avars(). This is the safest way. And really, I don't think any more fancy
 * functionality is required for this, because all the other classes internally use @p state and
 * @p avars directly.
 */
class CAvars
{
    private:
    State* sp_; // pointer to state
    Avars* ap_; // pointer to avars
    bool only_state; // true: only pointer to a state is set, false: both pointers set
    
    
    
    /**
     * @brief (Re)Setter for state variable pointer
     */
    void set_state(State *sp)
    {
        sp_ = sp;
    }
    
    
    
    /**
     * @brief (Re)Setter for avars variable pointer
     */
    void set_avars(Avars *ap)
    {
        ap_ = ap;
        only_state = false;
    }
    
    public:
    /**
     * @brief Constructor. Assigns pointers CAvars::sp_ and CAvars::ap_ to @p sp and @p ap.
     * CAvars::only_state is set to false (internally, in CAvars::set_avars())
     */
    CAvars(State *sp, Avars *ap)
    {
        set_state(sp);
        set_avars(ap);
    }
    
    
    
    /**
     * @brief Constructor using state alone.
     *
     * CAvars::ap_ is left unset. CAvars::only_state is set to true
     */
    CAvars(State *sp)
    {
        set_state(sp);
        only_state = true;
    }
    
    
    
    /**
     * @brief Assignment operator overload
     *
     * @pre `cav.only_state` must be `false` because this function internally calls
     * CAvars::get_avars() of @p cav
     */
    void operator=(const CAvars &cav)
    {
        *sp_ = cav.get_state();
        *ap_ = cav.get_avars();
    }
    
    
    
    /**
     * @brief Returns reference to state object held by CAvars::sp_
     */
    State& get_state()
    {
        return *sp_;
    }
    
    
    
    /**
     * @brief Returns reference to state object held by CAvars::sp_ (const version)
     */
    const State& get_state() const
    {
        return *sp_;
    }
    
    
    
    /**
     * @brief Returns reference to Avars object held by CAvars::ap_
     */
    Avars& get_avars()
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
     * @brief Returns reference to avars object held by CAvars::ap_ (const version)
     */
    const Avars& get_avars() const
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
        utilities::Testing t("CAvars", "class");
        {
            t.new_block("Testing construction (both)");
            State s={1,2,3,4};
            Avars a={1,2,3,4,5,6,7,8,9};
            CAvars ca(&s, &a);
            std::cout << "OK\n";
        }
        
        {
            t.new_block("Testing construction (only state)");
            State s={1,2,3,4};
            CAvars ca(&s);
            std::cout << "OK\n";
        }
        
        {
            t.new_block("Value modification testing");
            State s={1,2,3,4,5};
            Avars a={1,2,3,4,5,6,7,8,9};
            CAvars ca(&s, &a);
            
            std::cout << "\nStage 1";
            utilities::print_state(ca.get_state());
            utilities::print_avars(ca.get_avars());
            
            s[2] = 300;
            a[6] = 700;
            std::cout << "\nStage 2";
            utilities::print_state(ca.get_state());
            utilities::print_avars(ca.get_avars());
            
            State &s2 = ca.get_state();
            Avars &a2 = ca.get_avars();
            s2[1] = -100;
            a2[3] = -400;
            std::cout << "\nStage 3";
            utilities::print_state(ca.get_state());
            utilities::print_avars(ca.get_avars());
        }
        
        {
            t.new_block("Out of scope behaviour");
            State s={1,2,3,4,5};
            CAvars ca(&s);
            {
                Avars a={1,2,3,4,5,6,7,8,9};
                ca.set_avars(&a);
                std::cout << "avars in scope";
                utilities::print_avars(ca.get_avars());
            }
            std::cout << "avars out of scope";
            utilities::print_avars(ca.get_avars());
        }
        
        {
            t.new_block("Unset avars behaviour");
            State s={1,2,3,4,5};
            CAvars ca(&s);
            // utilities::print_avars(ca.get_avars()); // throw exception
        }
        
        {
            t.new_block("Const behaviour test");
            State s={1,2,3,4,5};
            Avars a={1,2,3,4,5,6,7,8,9};
            CAvars ca(&s, &a);
            // const state &s2 = ca.get_state();
            // s2[0] = 0; // error
            
            const CAvars ca2(&s, &a);
            const State &s2 = ca.get_state(); // OK
        }
        
        {
            t.new_block("operator= test");
            State s1={1,2,3,4,5};
            Avars a1={1,2,3,4,5,6,7,8,9};
            CAvars ca1(&s1, &a1);
            
            State s2={10,20,30,40,50};
            Avars a2={10,20,30,40,50,60,70,80,90};
            CAvars ca2(&s2, &a2);
            
            std::cout << "Before";
            utilities::print_state(ca1.get_state());
            utilities::print_avars(ca1.get_avars());
            ca1 = ca2;
            std::cout << "\nAfter";
            utilities::print_state(ca1.get_state());
            utilities::print_avars(ca1.get_avars());            
        }
    }
    #endif
};

#endif

