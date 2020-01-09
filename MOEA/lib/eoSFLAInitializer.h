#ifndef EOSFLAINITIALIZER_H
#define EOSFLAINITIALIZER_H


#include <eoPop.h>


/// VER FICHEIRO eoInitializer !!!



///**
//        Base (name) class for Initialization of algorithm SFLA

//        @see eoInitializerBase eoUF apply
//*/
//template <class POT> class eoSFLAInitializer : public eoInitializerBase<POT>
//{
//public:

//    //!	Constructor
//    //! @param _proc Evaluation function

//    //! @param _pop Population
//    eoSFLAInitializer(
//        eoUF<POT&, void>& _proc,
//        eoPop<POT>& _pop
//    ) : proc(_proc), pop(_pop)
//    {}

//    //! Give the name of the class
//    //! @return The name of the class
//    virtual std::string className (void) const
//    {
//        return "eoSFLAInitializer";
//    }

//    virtual void operator() ()
//    {
//        // Global exploration
//        //
//        // Step 0 Initialize. Select m and n, where m is the number of memeplexes and n is the number of frogs in each
//        // memeplex. Therefore, the total sample size, F, in the swamp is given by F = mn.
//        //
//        // Step 1 Generate a virtual population. Sample F virtual frogs, U(1),...,U(F).
//        // Compute the performance value f(i) for each frog U.


//        // Step 2 Rank frogs. Sort the F frogs in order of decreasing performance value. Store them in an array X={U(i),
//        // f(i), i=1,...,F} so that i=1 represents the frog with the best performance value. Record the best frog’s
//        // position, PX, in the entire population (F frogs; where PX =U(1)).


//    }

//private :


//};

#endif // EOSFLAINITIALIZER_H
