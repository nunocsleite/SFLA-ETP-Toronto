#ifndef ETTPNEIGHBORHOOD_H
#define ETTPNEIGHBORHOOD_H

#include <neighborhood/moNeighborhood.h>
#include <neighborhood/moNeighbor.h>
#include "ETTPneighbor.h"
#include <utils/eoRNG.h>

using namespace boost;

class ETTPneighborhood : public moNeighborhood<ETTPneighbor> {
    //
    // Auxiliary private methods
    //

    // Builds the clash list, consisting of all exams that are involved
    // in at least one clash
    void buildClashList(eoChromosome& timetable, unordered_map<int, Exam> & clashList) const {

        /// TODO: Exclude saturday clashes

//        vector<vector<int> >& periods = timetable.getPeriods();
        vector<unordered_set<int> >& periods = timetable.getPeriods();
        // Get exam's list pointer
//        vector<int>* examListP1 = &periods[0];
//        vector<int>* examListP2 = 0;
        unordered_set<int>* examListP1 = &periods[0];
        unordered_set<int>* examListP2 = 0;
        for (int p = 1; p < timetable.getNumPeriods(); ++p) {
            examListP2 = &periods[p];
//            for (int i = 0; i < examListP1->size(); ++i) {
//                for (int j = 0; j < examListP2->size(); ++j) {

            for (unordered_set<int>::const_iterator it_i = examListP1->begin(); it_i != examListP1->end(); ++it_i) {
                for (unordered_set<int>::const_iterator it_j = examListP2->begin(); it_j != examListP2->end(); ++it_j) {
                    int ei = *it_i;
                    int ej = *it_j;
                    // int ei = (*examListP1)[i], ej = (*examListP2)[j];

                    // Exam Periods
                    int pi = p-1, pj = p;
                    // Verify if there's a clash between exams
                    int numStudents = timetable.getConflictMatrix()->getVal(ei, ej);
                    if (numStudents > 0) {
                        // If not already inserted, insert exams ei and ej into the clash list
                        if (clashList.find(ei) == clashList.end()) {
                            clashList.insert(pair<int, Exam>(ei, Exam(ei, pi)));
                        }
                        if (clashList.find(ej) == clashList.end()) {
                            clashList.insert(pair<int, Exam>(ej, Exam(ej, pj)));
                        }
                    }
                }
            }
            examListP1 = examListP2;
        }
    }

//    struct NumberClashesCmp {
//        int operator()(ETTPneighbor const& n1, ETTPneighbor const& n2) {
//            return n1.getNumClashes() < n2.getNumClashes();
//        }
//    };

    struct ProximityCostCmp {
        int operator()(ETTPneighbor const& n1, ETTPneighbor const& n2) {
            return n1.getProximityCost() < n2.getProximityCost();
        }
    };


public:
    void randomNeighbor(eoChromosome& timetable, ETTPneighbor& neighbor) {

/*

        // In order to identify the most promising moves, a clash list,
        // like the one used in the period expansion operator, is maintained.
        unordered_map<int, Exam> clashList;
        buildClashList(timetable, clashList);

//        cout << endl << endl << "randomNeighbor" << endl;

//        cout << "timetable.fitness() = " << timetable.fitness() << endl;
//        cin.get();

        // Reserve a vector for keys in the clash list map,
        // for selecting entries randomly
        std::vector<int> keys;
        size_t numClashExams = clashList.size();
        keys.reserve(numClashExams);
        for (auto kv : clashList) {
            keys.push_back(kv.first);
        }

        feasibleNeigh = false;
        for (int i = 0; i < numClashExams && !feasibleNeigh; ++i) {
            // Randomly select a clashed exam to move
            int clashedExam = keys[eo::rng.random(numClashExams)];
            Exam ex = clashList.at(clashedExam);
            int clashPeriod = ex.getPeriod();
            // Determine feasible periods for clashed exam
            vector<int> feasiblePeriods = timetable.getFeasiblePeriods(clashedExam, clashPeriod);
            int numFeasiblePeriods = feasiblePeriods.size();

//            cout << "numFeasiblePeriods = " << numFeasiblePeriods << endl;

            // If there are any feasible periods, then evaluate all moves in order to determine the best move
            if (numFeasiblePeriods > 0) {
                // The move which leads to the greatest decrease in the number of clashes is selected
                // Create vector of neighbors
                vector<ETTPneighbor> neighbors(numFeasiblePeriods);
                for (int p = 0; p < numFeasiblePeriods; ++p) {
                    // Compute move cost
                    int periodToMove = feasiblePeriods[p];
                    neighbors[p].computeMoveCost(timetable, clashPeriod, clashedExam, periodToMove);

//                    movesCell{p} = computeMoveClashes(Data, timetable, clashPeriod, ...
//                                            clashedExam, periodToMove);

                }
                // Determine the move which produces the lowest number of clashes
//                vector<ETTPneighbor>::iterator minIt = min_element(neighbors.begin(), neighbors.end(), NumberClashesCmp());
                vector<ETTPneighbor>::iterator minIt =
                        min_element(neighbors.begin(), neighbors.end(), ProximityCostCmp());

                // Set neighbor
                neighbor = *minIt;
                feasibleNeigh = true;


        // Set neighbor
        neighbor = *minIt;
        feasibleNeigh = true;






//                cout << "neighbor.fitness() = " << neighbor.fitness() << endl;

                /// TODO: the evaluation and acceptance of this neighbor is done in SHC...

                // If the min move reduces the number of clashes of the actual solution, then accept this move.
//                if ((*minIt).getNumClashes() < timetable.getProximityCost())
//                    feasibleNeigh = true;

//                if (idxMin == -1 || movesCell{p}.numClashes < movesCell{idxMin}.numClashes)
//                    idxMin = p;
//                end

                // If the min move reduces the number of clashes (of the actual
                // solution), then accept this move.
//                randomNeigh = movesCell{idxMin}.timetable;
//                t = computeNumClashes(Data, randomNeigh);
//                if (randomNeigh.NumClashes ~= t.NumClashes)
//                    disp('Error')
//                    t.NumClashes
//                    randomNeigh.NumClashes
//                    pause
//                end

//                if (~verifyTimetable(Data, randomNeigh))
//                    fprintf('randomNeigh is not valid\n');
//                    pause
//                else
//                    return;
//                end
//            end
            }
        }

        */

        neighbor.setSolution(timetable);
        feasibleNeigh = true;
    }


    //        // In order to identify the most promising moves, a clash list,
    //        // like the one used in the period expansion operator, is maintained.
    //        [clashList, periodList] = buildClashList(_sol);

    //        numClashExams = length(clashList);
    //        randExamsClashList = randi([1, numClashExams], 1, numClashExams);
    //        feasibleNeigh = false;
    //        i = 1;
    //        while (i <= numClashExams && ~feasibleNeigh)
    //            // Clashed exam and respective period
    //            clashedExam = clashList(randExamsClashList(i));
    //            clashPeriod = periodList(randExamsClashList(i));
    //            // Determine feasible periods for clashed exam i
    //            feasiblePeriods = getFeasiblePeriods(Data, timetable, clashedExam, clashPeriod);
    //            numFeasiblePeriods = length(feasiblePeriods);
    //            if (numFeasiblePeriods ~= 0)
    //                // Evaluate moves
    //                // The move which leads to the greatest decrease in the number of clashes
    //                // is selected and the exam is removed from the clash list.
    //                movesCell = cell(1, numFeasiblePeriods);
    //                idxMin = -1;
    //                for p = 1 : numFeasiblePeriods
    //                    // compute move cost
    //                    periodToMove = feasiblePeriods(p);
    //                    movesCell{p} = computeMoveClashes(Data, timetable, clashPeriod, ...
    //                                            clashedExam, periodToMove);
    //                    // Determine the move which produces the lowest number of clashes
    //                    if (idxMin == -1 || movesCell{p}.numClashes < movesCell{idxMin}.numClashes)
    //                        idxMin = p;
    //                    end
    //                end
    //                // If the min move reduces the number of clashes (of the actual
    //                // solution), then accept this move.
    //                randomNeigh = movesCell{idxMin}.timetable;
    //                t = computeNumClashes(Data, randomNeigh);
    ////                if (randomNeigh.NumClashes ~= t.NumClashes)
    ////                    disp('Error')
    ////                    t.NumClashes
    ////                    randomNeigh.NumClashes
    ////                    pause
    ////                end

    //                if (~verifyTimetable(Data, randomNeigh))
    //                    fprintf('randomNeigh is not valid\n');
    //                    pause
    //                else
    //                    return;
    //                end
    //            end
    //            i = i+1;
    //        end
    ////        disp('Can''t find feasible random neighbor')
    ////    %     pause
    //    end


        //    function feasiblePeriods = getFeasiblePeriods(Data, timetable, ...
        //        clashedExam, clashPeriod)
        //        feasiblePeriods = [];
        //        % Determine possible periods for clashed exam
        //        for period = 1 : timetable.NumPeriods
        //            if (period ~= clashPeriod && isFeasibleExamPeriod(Data, timetable, period, clashedExam))
        //                feasiblePeriods = [feasiblePeriods period];
        //            end
        //        end
        //    end

        //    function moveInfo = computeMoveClashes(Data, timetable, clashPeriod, ...
        //                                        clashedExam, period)
        //        moveInfo.timetable = [];
        //        moveInfo.numClashes = -1;
        //        % Verify if clashedExam can move to period (if it's feasible)
        //        % If it can move, then compute the resulting timetable and number of
        //        % clashes.
        //        if (isFeasibleExamPeriod(Data, timetable, period, clashedExam))
        //            % Add clashed exam to period
        //            timetable.Periods{period} = [timetable.Periods{period} clashedExam];
        //            % Remove clashed exam from original period
        //            timetable.Periods{clashPeriod} = remove(clashedExam, timetable.Periods{clashPeriod});
        //            % Compute number of clashes of new timetable
        //            moveInfo.timetable = computeNumClashes(Data, timetable);
        //            moveInfo.numClashes = moveInfo.timetable.NumClashes;
        //        end
        //    end

public:

    ETTPneighborhood() : feasibleNeigh(true) {}

    /**
     * @return true if the neighborhood is random (default false)
     */
    virtual bool isRandom() {
        return true;
    }


    /**
     * Test if a solution has (again) a Neighbor
     * @param _solution the related solution
     * @return true if _solution has a Neighbor
     */
    virtual bool hasNeighbor(EOT & _solution) {
        // In the beginning, this variable is true in order to call randomNeighbor method
        return feasibleNeigh;
    }

    /**
     * Initialization of the neighborhood
     * @param _solution the solution to explore
     * @param _current the first neighbor
     */
    virtual void init(EOT & _solution, ETTPneighbor & _current) {
        randomNeighbor(_solution, _current);
    }

    /**
     * Give the next neighbor
     * @param _solution the solution to explore
     * @param _current the next neighbor
     */
    virtual void next(EOT & _solution, ETTPneighbor & _current) {
        randomNeighbor(_solution, _current);
    }

    /**
     * Test if there is again a neighbor
     * @param _solution the solution to explore
     * @return true if there is again a neighbor not explored
     */
    virtual bool cont(EOT & _solution) {
        return true;                                /// TODO: eager or lazy?
    }

    /**
     * Return the class Name
     * @return the class name as a std::string
     */
    virtual std::string className() const {
        return "ETTPNeighborhood";
    }

private:
    bool feasibleNeigh;
};

#endif // ETTPNEIGHBORHOOD_H

//function randomNeigh = getRandomNeigh(Data, timetable)
//    randomNeigh = [];
//    % In order to identify the most promising moves,
//    % a clash list, like the one used in the period expansion operator,
//    % is maintained.
//    [clashList, periodList] = buildClashList(Data, timetable);

//    numClashExams = length(clashList);
//    randExamsClashList = randi([1, numClashExams], 1, numClashExams);
//    feasibleNeigh = false;
//    i = 1;
//    while (i <= numClashExams && ~feasibleNeigh)
//        % Clashed exam and respective period
//        clashedExam = clashList(randExamsClashList(i));
//        clashPeriod = periodList(randExamsClashList(i));
//        % Determine feasible periods for clashed exam i
//        feasiblePeriods = getFeasiblePeriods(Data, timetable, clashedExam, clashPeriod);
//        numFeasiblePeriods = length(feasiblePeriods);
//        if (numFeasiblePeriods ~= 0)
//            % Evaluate moves
//            % The move which leads to the greatest decrease in the number
//            % of clashes is selected and the exam is removed from the
//            % clash list.
//            movesCell = cell(1, numFeasiblePeriods);
//            idxMin = -1;
//            for p = 1 : numFeasiblePeriods
//                % compute move cost
//                periodToMove = feasiblePeriods(p);
//                movesCell{p} = computeMoveClashes(Data, timetable, clashPeriod, ...
//                                        clashedExam, periodToMove);
//                % Determine the move which produces the lowest number of clashes
//                if (idxMin == -1 || movesCell{p}.numClashes < movesCell{idxMin}.numClashes )
//                    idxMin = p;
//                end
//            end
//            % If the min move reduces the number of clashes (of the actual
//            % solution), then accept this move.
//            randomNeigh = movesCell{idxMin}.timetable;
//            t = computeNumClashes(Data, randomNeigh);
//            if (randomNeigh.NumClashes ~= t.NumClashes)
//                disp('Error')
//                t.NumClashes
//                randomNeigh.NumClashes
//                pause
//            end

//            if (~verifyTimetable(Data, randomNeigh))
//                fprintf('randomNeigh is not valid\n');
//                pause
//            else
//                return;
//            end
//        end
//        i = i+1;
//    end
//    disp('Can''t find feasible random neighbor')
//%     pause
//end

//function feasiblePeriods = getFeasiblePeriods(Data, timetable, ...
//    clashedExam, clashPeriod)
//    feasiblePeriods = [];
//    % Determine possible periods for clashed exam
//    for period = 1 : timetable.NumPeriods
//        if (period ~= clashPeriod && isFeasibleExamPeriod(Data, timetable, period, clashedExam))
//            feasiblePeriods = [feasiblePeriods period];
//        end
//    end
//end

//function [clashList, periodList] = buildClashList(Data, timetable)
//    clashList = [];
//    periodList = [];
//    examListP1 = timetable.Periods{1};
//    for p = 2 : timetable.NumPeriods
//        examListP2 = timetable.Periods{p};
//        if (~(p == 7 || p == 13 || p == 19))
//            numExamsP1 = length(examListP1);
//            numExamsP2 = length(examListP2);
//            marked1 = zeros(1, length(examListP1));
//            marked2 = zeros(1, length(examListP2));
//            for i = 1 : numExamsP1
//                for j = 1 : numExamsP2
//                    numStudents = Data.ConflictMatrix(examListP1(i), examListP2(j));
//                    if (numStudents > 0) % There's a clash between exams
//                        marked1(i) = 1;
//                        marked2(j) = 1;
//                    end
//                end
//            end
//            % If marked exams aren't yet in the clash list, then copy them.
//            [clashList, periodList] = copyExamsClashList(p-1, ...
//                periodList, clashList, examListP1(find(marked1 == 1)));
//            [clashList, periodList] = copyExamsClashList(p, ...
//                periodList, clashList, examListP2(find(marked2 == 1)));
//        end
//        examListP1 = examListP2;
//    end
//end

//function [clashList periodList] = copyExamsClashList(period, periodList, ...
//                                                        clashList, examList)
//    if (~isempty(examList))
//        for i = 1 : length(examList)
//            if (isempty(find(clashList == examList(i))))
//                clashList = [clashList examList(i)];
//                periodList = [periodList period];
//            end
//        end
//    end
//end

//function moveInfo = computeMoveClashes(Data, timetable, clashPeriod, ...
//                                    clashedExam, period)
//    moveInfo.timetable = [];
//    moveInfo.numClashes = -1;
//    % Verify if clashedExam can move to period (if it's feasible)
//    % If it can move, then compute the resulting timetable and number of
//    % clashes.
//    if (isFeasibleExamPeriod(Data, timetable, period, clashedExam))
//        % Add clashed exam to period
//        timetable.Periods{period} = [timetable.Periods{period} clashedExam];
//        % Remove clashed exam from original period
//        timetable.Periods{clashPeriod} = remove(clashedExam, timetable.Periods{clashPeriod});
//        % Compute number of clashes of new timetable
//        moveInfo.timetable = computeNumClashes(Data, timetable);
//        moveInfo.numClashes = moveInfo.timetable.NumClashes;
//    end
//end
//%////////////////////////////////////////////////////////////////////////
