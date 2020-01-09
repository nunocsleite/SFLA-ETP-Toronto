#ifndef MOLASAEXPLORER_H
#define MOLASAEXPLORER_H


/*
  <moLASAexplorer.h>

  Based on <moSAexplorer.h> written by
     DOLPHIN Project-Team, INRIA Lille - Nord Europe, 2006-2010

     Sébastien Verel, Arnaud Liefooghe, Jérémie Humeau

     ParadisEO WebSite : http://paradiseo.gforge.inria.fr
     Contact: paradiseo-help@lists.gforge.inria.fr
*/

#include <cstdlib>

#include <explorer/moNeighborhoodExplorer.h>
#include <comparator/moSolNeighborComparator.h>
#include <coolingSchedule/moCoolingSchedule.h>
#include <neighborhood/moNeighborhood.h>
#include <utils/eoRNG.h>

#include <boost/circular_buffer.hpp> // Optimized Boost's Queue
#include <queue>

/**
 * Explorer for the Simulated Annealing
 * Only the symetric case is considered when Q(x,y) = Q(y,x)
 * Fitness must be > 0
 *
 */
template< class Neighbor >
class moLASAexplorer : public moNeighborhoodExplorer<Neighbor>
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
     * @param _coolingSchedule the cooling schedule
     */
  moLASAexplorer(Neighborhood& _neighborhood, moEval<Neighbor>& _eval, moSolNeighborComparator<Neighbor>& _solNeighborComparator, moCoolingSchedule<EOT>& _coolingSchedule)
      :
//        queue_size(1),
//                queue_size(25),
        queue_size(500),
//        queue_size(14000),

        moNeighborhoodExplorer<Neighbor>(_neighborhood, _eval), solNeighborComparator(_solNeighborComparator),
        coolingSchedule(_coolingSchedule)/*, solutionsQueue(queue_size)*/ {


        isAccept = false;

        if (!neighborhood.isRandom()) {
            std::cout << "moLASAexplorer::Warning -> the neighborhood used is not random" << std::endl;
        }
    }

    /**
     * Destructor
     */
    ~moLASAexplorer() {
    }

    /**
     * initialization of the initial temperature
     * @param _solution the solution
     */
    virtual void initParam(EOT & _solution) {
      temperature = coolingSchedule.init(_solution);
      isAccept = false;

//      cout << "Solutions queue size:" << endl;
//      cout << solutionsQueue.size() << endl;

      // Initialize solutions queue: at the beginning, all the elements
      // of the queue are the same and are equal to the initial cost function
      for (size_t i = 0; i < queue_size; ++i)
//          solutionsQueue.push_back(_solution.fitness());
          solutionsQueue.push(_solution.fitness());


//      cin.get();
    };

    /**
     * decrease the temperature if necessary
     * @param _solution unused solution
     */
    virtual void updateParam(EOT & _solution) {
//        cout << "Before Temperature = " << temperature << endl;
        coolingSchedule.update(temperature, this->moveApplied());
//        cout << "After Temperature = " << temperature << endl;
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
        //Test if _solution has a Neighbor
        if (neighborhood.hasNeighbor(_solution)) {
            //init on the first neighbor: supposed to be random solution in the neighborhood
            neighborhood.init(_solution, selectedNeighbor);

            //eval the _solution moved with the neighbor and stock the result in the neighbor
            eval(_solution, selectedNeighbor);
        }
        else {
            //if _solution hasn't neighbor,
            isAccept = false;
        }
    };

    /**
     * continue if the temperature is not too low
     * @param _solution the solution
     * @return true if the criteria from the cooling schedule is true
     */
    virtual bool isContinue(EOT & _solution) {
        return coolingSchedule(temperature);
    };

    /**
     * acceptance criterion according to the boltzmann criterion
     * @param _solution the solution
     * @return true if better neighbor or rnd < exp(delta f / T)
     */
    //  % Simulated Annealing (algorithm 2 for minimization of f)
    //  %
    //  % Step 1  Make T = Tmax  and Choose a solution u (at random)
    //  %
    //  % Step 2  Select a neighbor of u, say v
    //  %         If f(v) < f(u) make u = v;
    //  %         Else make u = v with probability
    //  %              p = exp((fu-fv)/(fu T)))
    //  %
    //  %         Repeat Step 2 k times
    //  %
    //  % Step 3  Make t = t+1; Set T = T(t) see Eq.(4)
    //  %
    //  %         If  T < Tmin  go to Step 2;
    //  %         Else go to Step 1.
    //  %
    virtual bool accept(EOT & _solution) {
//        cout << "accept method" << endl;

        if (neighborhood.hasNeighbor(_solution)) {
//            if (solNeighborComparator(_solution, selectedNeighbor)) { // accept if the current neighbor is better than the solution

            // DEBUG
//            cout << "Solutions queue contents:" << endl;

//            for (boost::circular_buffer<double>::const_iterator it = solutionsQueue.begin(); it != solutionsQueue.end(); ++it) {
//                cout << *it << " " << endl;
//            }


            double fitNeigh, fitSol;
//            fit1 = (double)selectedNeighbor.fitness();
            fitNeigh = (double)selectedNeighbor.fitness();
//            fitSol = (double)_solution.fitness();
            fitSol = solutionsQueue.back();
//            solutionsQueue.pop(); // hum
//            solutionsQueue.pop_back();

//            cout << "fit1 (neighbour) = " << fitNeigh << endl;
//            cout << "fit2 (current solution) = " << fitSol << endl;


//            if (selectedNeighbor.fitness() < _solution.fitness()) { // Minimization problem
            if (fitNeigh < fitSol) { // Minimization problem
                isAccept = true;
//                cout << "accept because the current neighbor is better than the solution" << endl;
//                solutionsQueue.push_front(selectedNeighbor.fitness());
                solutionsQueue.push(selectedNeighbor.fitness());

//                  solutionsQueue.push_front(_solution.fitness());



            }
            else {
                double alpha = 0.0;
                alpha = exp((fitSol - fitNeigh) / (fitSol*temperature) );

                double r = rng.uniform();

//                cout << "temperature = " << temperature << endl;
//                cout << "alpha = " << alpha << endl;
//                cout << "rand = " << r << endl;

                isAccept = (r < alpha);

//                cout << "isAccept = " << isAccept << endl;

                if (isAccept) {
                  solutionsQueue.push(selectedNeighbor.fitness());

//                    solutionsQueue.push_front(selectedNeighbor.fitness());
//                    solutionsQueue.push_front(_solution.fitness());
                }

            }

            // Insert new neighbor in the queue
//        solutionsQueue.push_front(selectedNeighbor.fitness());

//            solutionsQueue.push_front(_solution.fitness());

/*
            double fit1, fit2;
            fit1 = (double)selectedNeighbor.fitness();
            fit2 = (double)_solution.fitness();

//            cout << "fit1 (neighbour) = " << fit1 << endl;
//            cout << "fit2 (current solution) = " << fit2 << endl;


            if (selectedNeighbor.fitness() < _solution.fitness()) { // Minimization problem
                isAccept = true;
//                cout << "accept because the current neighbor is better than the solution" << endl;
            }
            else {
                double alpha = 0.0;
//                double fit1, fit2;
//                fit1 = (double)selectedNeighbor.fitness();
//                fit2 = (double)_solution.fitness();

//                cout << "fit1 = " << fit1 << endl;
//                cout << "fit2 = " << fit2 << endl;

//                if (fit1 < fit2) // this is a maximization
//                    alpha = exp((fit1 - fit2) / temperature );
//                else // this is a minimization
//                    alpha = exp((fit2 - fit1) / temperature );
//                isAccept = (rng.uniform() < alpha);

//                alpha = exp((fit2 - fit1) / (fit2*0.01) ); // 40.56
//                alpha = exp((fit2 - fit1) / (fit2*2) ); //
                alpha = exp((fit2 - fit1) / (fit2*temperature) );
//                alpha = exp((fit2 - fit1) / temperature );

                double r = rng.uniform();

//                cout << "temperature = " << temperature << endl;
//                cout << "alpha = " << alpha << endl;
//                cout << "rand = " << r << endl;


                isAccept = (r < alpha);
//                isAccept = (rng.uniform() < alpha);

//                cout << "isAccept = " << isAccept << endl;

            }
*/
        }
        return isAccept;
    };

    /**
     * Getter
     * @return the temperature
     */
    double getTemperature() {
        return temperature;
    }

private:
    // comparator betwenn solution and neighbor
    moSolNeighborComparator<Neighbor>& solNeighborComparator;

    moCoolingSchedule<EOT>& coolingSchedule;

    // temperature of the process
    double temperature;

    // true if the move is accepted
    bool isAccept;

    // Optimized queue implemented as circular buffer
    queue<double> solutionsQueue;
//    boost::circular_buffer<double> solutionsQueue;

    const size_t queue_size;
};


#endif // MOLASAEXPLORER_H


