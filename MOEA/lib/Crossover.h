#ifndef CROSSOVER_H
#define CROSSOVER_H

#include "Chromosome.h"


using namespace std;


class Crossover : public eoQuadOp<moeoChromosome>
{

public:

    /**
    * the class name (used to display statistics)
    */
    string className() const;

    /**
    * eoQuad crossover - _chromosome1 and _chromosome2 are the (future) offspring, i.e. _copies_ of the parents
    * @param _chromosome1 the first parent
    * @param _chromosome2 the second parent
    */
    bool operator()(moeoChromosome& _chromosome1, moeoChromosome& _chromosome2);

private:

    /**
    * generation of an offspring
    * @param _parent1 the first parent
    * @param _parent2 the second parent
    */
    void generateOffspring(moeoChromosome &_parent1, moeoChromosome &_parent2);

    /**
    * best day computation
    * @param chromosome
    */
    int bestDay(moeoChromosome const& chromosome);

    /**
    * insert day into chromosome
    * @param _parent1
    * @param _parent2
    * @param day2
    */
    void insertDay(moeoChromosome& _parent1, moeoChromosome& _parent2, int day2);
};




#endif // CROSSOVER_H
