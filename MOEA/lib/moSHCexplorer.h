#ifndef MOSHCEXPLORER_H
#define MOSHCEXPLORER_H

#include <cstdlib>

#include <explorer/moNeighborhoodExplorer.h>
#include <comparator/moSolNeighborComparator.h>
#include <neighborhood/moNeighborhood.h>

#include <utils/eoRNG.h>


enum ProblSense { minimize, maximize };

/**
 * Explorer for the Stochastic hill climber
 * Fitness must be > 0
 *
 */
template <class Neighbor>
class moSHCexplorer : public moNeighborhoodExplorer<Neighbor>
{
public:
    typedef typename Neighbor::EOT EOT ;
    typedef moNeighborhood<Neighbor> Neighborhood ;

    using moNeighborhoodExplorer<Neighbor>::neighborhood;
    using moNeighborhoodExplorer<Neighbor>::eval;
    using moNeighborhoodExplorer<Neighbor>::selectedNeighbor;

    /**
     * Constructor
     * @param _neighborhood the neighborhood
     * @param _eval the evaluation function
     * @param _solNeighborComparator a solution vs neighbor comparator
     * @param _tmax Number of iterations
     * @param _tempT Temperature
     * @param _sense Problem sense (minimization or maximization)
     */
  moSHCexplorer(Neighborhood& _neighborhood, moEval<Neighbor>& _eval,
                moSolNeighborComparator<Neighbor>& _solNeighborComparator,
                double _tmax, double _tempT, ProblSense _sense)
        : moNeighborhoodExplorer<Neighbor>(_neighborhood, _eval), solNeighborComparator(_solNeighborComparator),
          tmax(_tmax), tempT(_tempT), sense(_sense)
    {
        isAccept = false;
        if (!neighborhood.isRandom()) {
            std::cout << "moSHCexplorer::Warning -> the neighborhood used is not random" << std::endl;
        }
    }

    /**
     * Destructor
     */
    ~moSHCexplorer() {}

    /**
    * @return true if the neighborhood is random (default false)
    */
    virtual bool isRandom() {
        return true;
    }

    /**
     * initialization of the initial temperature
     * @param _solution the solution
     */
    virtual void initParam(EOT & _solution) {
        t = 0;
        isAccept = false;
    };

    /**
     * decrease the temperature if necessary
     * @param _solution unused solution
     */
    virtual void updateParam(EOT & _solution) {
        ++t;
    };

    /**
     * terminate: NOTHING TO DO
     * @param _solution unused solution
     */
    virtual void terminate(EOT & _solution) {};

    /**
     * Explore one random solution in the neighborhood
     * @param _solution the solution
     */
    virtual void operator()(EOT & _solution) {

//        cout << endl << "moSHCexplorer.operator()" << endl << endl;
       // cin.get();

        // Test if _solution has a Neighbor
        if (neighborhood.hasNeighbor(_solution)) {
            // Init on the first neighbor: supposed to be random solution in the neighborhood
            neighborhood.init(_solution, selectedNeighbor);
            // Eval the _solution moved with the neighbor and stock the result in the neighbor
            eval(_solution, selectedNeighbor);
        }
        else {
            // If _solution hasn't neighbor,
            isAccept = false;
        }
    };

    /**
     * continue if not reached tmax iterations
     * @param _solution the solution
     * @return true if not reached tmax iterations
     */
    virtual bool isContinue(EOT & _solution) {
        return t < tmax;
    };

//    /**
//     * acceptance criterion according to the boltzmann criterion
//     * @param _solution the solution
//     * @return true if better neighbor or rnd < exp(delta f / T)
//     */
//    virtual bool accept(EOT & _solution) {
//        if (neighborhood.hasNeighbor(_solution)) {
//            if (solNeighborComparator(_solution, selectedNeighbor)) // accept if the current neighbor is better than the solution
//                isAccept = true;
//            else {
//            double alpha=0.0;
//            double fit1, fit2;
//            fit1=(double)selectedNeighbor.fitness();
//                fit2=(double)_solution.fitness();
//                if (fit1 < fit2) // this is a maximization
//                    alpha = exp((fit1 - fit2) / temperature );
//                else // this is a minimization
//                    alpha = exp((fit2 - fit1) / temperature );
//                isAccept = (rng.uniform() < alpha) ;
//            }
//        }
//        return isAccept;
//    };

    /**
    * Stochastic Hill Climber acceptance criterion
    * @param _solution the solution
    * @return true if rnd < 1 / (1 + exp(factor*(fv-fu) / (tempT*fu)))
    */
    virtual bool accept(EOT & _solution) {
        // SHC - Stochastic hill climber (minimization of f)
        // Choose a solution u (at random) and compute fu = f(u)
        // Repeat step 1 to 3 (until t = tmax)
        //  Step 1  Find (at random) a neighbor of u, say v.
        //  Step 2  Evaluate v, i.e., compute fv = f(v)
        //  Step 3  Make u = v with probability p(fu,fv), see Eq.(1);
        //  Make t = t+1.
        if (neighborhood.hasNeighbor(_solution)) {
            int factor;
            if (sense == maximize) // maximization problem
                factor = -1;
            else if (sense == minimize) // minimization problem
                factor = 1;

            double alpha = 0.0;
            double fu, fv;
            fu = (double)_solution.fitness();                     /// ATTENTION: COMPARE SAME KIND OF FITNESS!!!
            fv = (double)selectedNeighbor.fitness();

//            cout << (sense == minimize) << endl;
//            cout << "fu = " << fu << ", fv = " << fv << endl;


            alpha = 1 / (1 + exp(factor*(fv-fu) / (tempT*fu)));

            if (fv > fu)
                isAccept = false;
            else
                isAccept = (rng.uniform() < alpha);

//            double p = rng.uniform();
//            isAccept = (rng.uniform() < alpha);
//            isAccept = (p < alpha);
//            cout << "p = " << p << ", alpha = " << alpha << " isAccept = " << isAccept << endl;

        }
        return isAccept;
    };


private:
    // comparator betwenn solution and neighbor
    moSolNeighborComparator<Neighbor>& solNeighborComparator;
    // Optimization problem sense (maximization or minimization)
    ProblSense sense;
    // Iteration counter
    int t;
    // Max number of iterations
    int tmax;
    // SHC Temperature
    double tempT;
    // true if the move is accepted
    bool isAccept;
};

#endif // MOSHCEXPLORER_H


//function Res = SHC(tmax, T, data, initialSolution, getRandomNeigh, ...
//    evalFunc)
//    % SHC - Stochastic hill climber (minimization of f)
//    %
//    % Choose a solution u (at random) and compute fu = f(u)
//    % Repeat step 1 to 3 (until t = tmax)
//    %   Step 1  Find (at random) a neighbor of u, say v.
//    %   Step 2  Evaluate v, i.e., compute fv = f(v)
//    %   Step 3  Make u = v with probability p(fu,fv), see Eq.(1);
//    %   Make t = t+1.
//    %
//    sense = 'minimize';

//	% Number of evaluations
//    numEvaluations = 0;
//    % Choose a solution u (at random) and compute fu = f(u)
//    u = initialSolution;
//    fu = evalFunc(u, data);
//    % Increment number of evaluations
//    numEvaluations = numEvaluations + 1;
//    z = 1;
//    F(z) = fu;
//    z = z+1;
//    orig = fu;

//    % Repeat step 1 to 3 (until t = tmax)
//    t = 1;
//    while (t <= tmax)
//        % Step 1  Find (at random) a neighbor of u, say v.
//        v = getRandomNeigh(data, u);
//        % Step 2 Evaluate v, i.e., compute fv = f(v)
//	    fv = evalFunc(v, data);
//	    % Increment number of evaluations
//        numEvaluations = numEvaluations + 1;

//        % Step 4 Make u = v with probability p(fu,fv)
//        prob = p(fu, fv, T, sense);

//        x = myRand();
//        if (x < prob)
//            % Accept this solution
//            u = v;
//            fu = fv;
//        end
//        F(z) = fu;
//        z = z+1;
//        % Make t = t+1.
//        t = t + 1;
//    end

//%     orig
//%     fu

//    Res = struct('T', T, 'tmax', tmax, 'NumEvaluations', numEvaluations, ...
//        'Cost', fu, 's', u, 'u', u);
//end


//%/////////////////////////////////////////////////////

//% Random number generator
//function r = myRand()
//    r = rand;
//end

//% Probability function
//function x = p(fu, fv, T, sense)
//    % Minimization problem
//    %x = 1 / (1 + exp((fv-fu) / (T*fu)));
//    % Maximization problem
//    %x = 1 / (1 + exp((fu-fv) / (T*fu)));
//    if strcmp(sense, 'maximize')
//        factor = -1;
//    elseif strcmp(sense, 'minimize')
//        factor = 1;
//    end
//    x = 1 / (1 + exp(factor*(fv-fu) / (T*fu)));
//end
