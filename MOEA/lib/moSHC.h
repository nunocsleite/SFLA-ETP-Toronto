#ifndef MOSHC_H
#define MOSHC_H

#include <algo/moLocalSearch.h>
#include "moSHCexplorer.h"
#include <continuator/moTrueContinuator.h>
#include <eval/moEval.h>
#include <eoEvalFunc.h>

/**
 * moSHC - Stochastic hill climber (minimization of f)
 */
template<class Neighbor>
class moSHC: public moLocalSearch<Neighbor>
{
public:

    typedef typename Neighbor::EOT EOT;
    typedef moNeighborhood<Neighbor> Neighborhood ;

    /**
     * Basic constructor for Stochastic Hill Climber
     * @param _neighborhood the neighborhood
     * @param _fullEval the full evaluation function
     * @param _eval neighbor's evaluation function
     * @param _tmax Number of iterations
     * @param _tempT Temperature
     * @param _sense Problem sense (minimization or maximization)
     */
    moSHC(Neighborhood& _neighborhood, eoEvalFunc<EOT>& _fullEval, moEval<Neighbor>& _eval,
          double _tmax, double _tempT, ProblSense _sense) :
            moLocalSearch<Neighbor>(explorer, trueCont, _fullEval),
            explorer(_neighborhood, _eval, defaultSolNeighborComp, _tmax, _tempT, _sense)
    { }


/// TODO: other constructors

//    /**
//     * Simple constructor for a simulated annealing
//     * @param _neighborhood the neighborhood
//     * @param _fullEval the full evaluation function
//     * @param _eval neighbor's evaluation function
//     * @param _cool a cooling schedule
//     */
//    moSA(Neighborhood& _neighborhood, eoEvalFunc<EOT>& _fullEval, moEval<Neighbor>& _eval, moCoolingSchedule<EOT>& _cool):
//            moLocalSearch<Neighbor>(explorer, trueCont, _fullEval),
//            defaultCool(0, 0, 0, 0),
//            explorer(_neighborhood, _eval, defaultSolNeighborComp, _cool)
//    {}

//    /**
//     * General constructor for a simulated annealing
//     * @param _neighborhood the neighborhood
//     * @param _fullEval the full evaluation function
//     * @param _eval neighbor's evaluation function
//     * @param _cool a cooling schedule
//     * @param _cont an external continuator
//     */
//    moSA(Neighborhood& _neighborhood, eoEvalFunc<EOT>& _fullEval, moEval<Neighbor>& _eval, moCoolingSchedule<EOT>& _cool, moContinuator<Neighbor>& _cont):
//            moLocalSearch<Neighbor>(explorer, _cont, _fullEval),
//            defaultCool(0, 0, 0, 0),
//            explorer(_neighborhood, _eval, defaultSolNeighborComp, _cool)
//    {}

//    /**
//     * General constructor for a simulated annealing
//     * @param _neighborhood the neighborhood
//     * @param _fullEval the full evaluation function
//     * @param _eval neighbor's evaluation function
//     * @param _cool a cooling schedule
//     * @param _comp a solution vs neighbor comparator
//     * @param _cont an external continuator
//     */
//    moSA(Neighborhood& _neighborhood, eoEvalFunc<EOT>& _fullEval, moEval<Neighbor>& _eval, moCoolingSchedule<EOT>& _cool, moSolNeighborComparator<Neighbor>& _comp, moContinuator<Neighbor>& _cont):
//            moLocalSearch<Neighbor>(explorer, _cont, _fullEval),
//            defaultCool(0, 0, 0, 0),
//            explorer(_neighborhood, _eval, _comp, _cool)
//    {}

private:
    moTrueContinuator<Neighbor> trueCont;
    moSolNeighborComparator<Neighbor> defaultSolNeighborComp;
    moSHCexplorer<Neighbor> explorer;
};


#endif // MOSHC_H
