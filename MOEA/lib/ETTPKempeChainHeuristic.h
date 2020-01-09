#ifndef ETTPKEMPECHAINHEURISTIC_H
#define ETTPKEMPECHAINHEURISTIC_H

#include "eoFunctor.h"
#include "Chromosome.h"

using namespace std;
using namespace boost;


//////////////////////////////////////////////////////////
//
// (rf. Demeester paper)
// Concerning the uncapacitated Toronto examination time-
// tabling problem, we apply the Kempe chain based heuris-
// tics. These low-level heuristics perturb feasible solutions
// to the uncapacitated examination timetabling problem,
// without making them infeasible. Suppose a partial solu-
// tion that corresponds to the left hand side of Fig. 1. If
// we want to move exam e1 to time slot tj, the solution
// becomes infeasible, since exam e1 has students in com-
// mon with exams e6, e7, and e8 that are assigned to time
// slot tj. To overcome this, exams e6, e7, and e8 should be
// moved to time slot ti. This process is repeated until all the
// exams that have students in common are assigned to dif-
// ferent time slots. The result is depicted at the right hand
// side of Fig. 1. The Kempe Chain heuristic Uncap1 can be
// seen as moving an exam to another time slot whilst main-
// taining feasibility by repairing any infeasibilities that may
// have been introduced.
//
class ETTPKempeChainHeuristic : public eoUF<eoChromosome&, void> {

public:
    virtual void operator()(eoChromosome& _chrom);

private:
    void kempeChainOperator(eoChromosome& _chrom, unordered_set<int>& ti_period,
                            unordered_set<int>& tj_period, int exami);

};


#endif // ETTPKEMPECHAINHEURISTIC_H




























