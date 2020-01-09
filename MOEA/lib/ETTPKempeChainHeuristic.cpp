
#include "ETTPKempeChainHeuristic.h"
#include <boost/unordered_set.hpp>


void ETTPKempeChainHeuristic::operator()(eoChromosome& _chrom) {

//    cout << "kempe chain" << endl;


    int exami;
    int ti, tj;

    unordered_set<int>* ti_period;
    unordered_set<int>* tj_period;


    if (true) {
//    if (false) {

    /// Original implementation
    // *
    // *
    // Select randomly two timeslots, ti and tj.
    int ti, tj;

    do {
        ti = rng.random(_chrom.getNumPeriods());
        do {
            tj = rng.random(_chrom.getNumPeriods());
        }
        while (ti == tj);
    }
    while (_chrom.getPeriods()[ti].size() == 0);


//    do {
//        ti = rng.random(_chrom.getNumPeriods());

////        int delta = rng.uniform(1, _chrom.getNumPeriods()); // 37.95
//        int delta = rng.uniform(1, 6); //


//        if (ti == _chrom.getNumPeriods()-1)
//            tj = ti-delta;
//        else if (ti == 0)
//            tj = ti+delta;
//        else {
//            if (rng.flip() < 0.5)
//                tj = ti-delta < 0 ? 0 : ti-delta;
//            else
//                tj = ti+delta >= _chrom.getNumPeriods()-1 ? _chrom.getNumPeriods()-1 : ti+delta;
//        }
//    }
//    while (_chrom.getPeriods()[ti].size() == 0);


//    cout << "\t" << abs(ti-tj) << endl;

    // Select randomly an exam ei from timeslot ti and move it
    // to timeslot tj. If the solution becomes infeasible,
    // because exam ei has students in common with exams ej, ek, ...,
    // that are assigned to time slot tj, one have to move exams
    // ej, ek, ..., to time slot ti. This process is repeated until all the
    // exams that have students in common are assigned to dif-
    // ferent time slots.


//    cout << "ti = " << ti << ", tj = " << tj << ", _chrom.getNumPeriods() = " << _chrom.getNumPeriods() << endl;


    vector<unordered_set<int> >& periods = _chrom.getPeriods();
//    unordered_set<int>& ti_period = periods[ti];
//    unordered_set<int>& tj_period = periods[tj];

    ti_period = &periods[ti];
    tj_period = &periods[tj];


    // Generate random index
//    int idx = rng.random(ti_period.size());
    int idx = rng.random(ti_period->size());



//    cout << "idx = " << idx << ", ti_period.size() = " << ti_period.size() << endl;

    /// TODO - NOT EFFICIENT
    // Get exam id
    unordered_set<int>::const_iterator it_i;
//    for (it_i = ti_period.begin(); idx-- > 0; ++it_i)
//        ;
//    it_i = ti_period.begin();
    it_i = ti_period->begin();

    for (int i = 0; i < idx; ++i, ++it_i)
        ;

    exami = *it_i;

//    ****************************************************************

//}

//else {



/*
    // New developments

    ////////////////////////////////////////////////////////////////////////
    // Kempe Chain variant: move the exam with highest penalty
    //

    // Timeslots ti and tj
    int ti, tj;
    int examHighCost = -1, period = -1;
    double highProxCost = 0.0;
    // Determine exam with highest penalty along with respective timeslot ti
    _chrom.getInfoExamHighestProximityCost(examHighCost, period, highProxCost); // ref parameters
    ti = period;

    // Select randomly a timeslot tj
    do {
        tj = rng.random(_chrom.getNumPeriods());
    }
    while (ti == tj);

    // Move highest penalty exam from timeslot ti
    // to timeslot tj. If the solution becomes infeasible,
    // because exam ei has students in common with exams ej, ek, ...,
    // that are assigned to time slot tj, one have to move exams
    // ej, ek, ..., to time slot ti. This process is repeated until all the
    // exams that have students in common are assigned to dif-
    // ferent time slots.


//    cout << "ti = " << ti << ", tj = " << tj << ", _chrom.getNumPeriods() = " << _chrom.getNumPeriods() << endl;

    vector<unordered_set<int> >& periods = _chrom.getPeriods();
//    unordered_set<int>& ti_period = periods[ti];
//    unordered_set<int>& tj_period = periods[tj];

    ti_period = &periods[ti];
    tj_period = &periods[tj];

//    // Generate random index
//    int idx = rng.random(ti_period.size());

////    cout << "idx = " << idx << ", ti_period.size() = " << ti_period.size() << endl;

//    /// TODO - NOT EFFICIENT
//    // Get exam id
//    unordered_set<int>::const_iterator it_i;
////    for (it_i = ti_period.begin(); idx-- > 0; ++it_i)
////        ;
//    it_i = ti_period.begin();
//    for (int i = 0; i < idx; ++i, ++it_i)
//        ;

//    int exami = *it_i;

    //int exami = examHighCost;
    exami = examHighCost;

    ////////////////////////////////////////////////////////////////////////
*/
//}


//    cout << "Before kempeChainOperator(_chrom, ti_period, tj_period, exami);" << endl;

//    cout << "ti = " << ti << ", tj = " << tj << endl;

//    cout << "#exams in slot ti = " << ti_period.size() << endl;
//    cout << "#exams in slot tj = " << tj_period.size() << endl;

//    cout << "Exams's list in slot ti = " << endl;
//    copy(ti_period.begin(), ti_period.end(), ostream_iterator<int>(cout, "\n"));

//    cout << "Exams's list in slot tj = " << endl;
//    copy(tj_period.begin(), tj_period.end(), ostream_iterator<int>(cout, "\n"));

    // Kempe chain operator
//    kempeChainOperator(_chrom, ti_period, tj_period, exami);
    kempeChainOperator(_chrom, *ti_period, *tj_period, exami);

}


//    cout << "After kempeChainOperator(_chrom, ti_period, tj_period, exami);" << endl;


    _chrom.computeProximityCosts();

//    _chrom.validate();

//    cout << "Validate Ok" << endl;

//    cout << "#exams in slot ti = " << ti_period.size() << endl;
//    cout << "#exams in slot tj = " << tj_period.size() << endl;

//    cout << "Exams's list in slot ti = " << endl;
//    copy(ti_period.begin(), ti_period.end(), ostream_iterator<int>(cout, "\n"));

//    cout << "Exams's list in slot tj = " << endl;
//    copy(tj_period.begin(), tj_period.end(), ostream_iterator<int>(cout, "\n"));


    // BUG:

//    kempe chain
//    ti = 2, tj = 7, _chrom.getNumPeriods() = 24
//    idx = 0, ti_period.size() = 0


}


void ETTPKempeChainHeuristic::kempeChainOperator(eoChromosome& _chrom, unordered_set<int>& ti_period,
                                                 unordered_set<int>& tj_period, int exami) {

//    cout << "kempe chain" << endl;
//    cout << "ei = " << exami << endl;


    /// TODO > SEE OTHER WAYS: WORKING WITH THE GRAPH, FIND & MOVE,...

    // Verify exam ei feasibility with exams in slot tj and move exams
    // if necessary
    int n;
    vector<int> conflictingExams_tj;
    for (unordered_set<int>::const_iterator it_j = tj_period.begin(); it_j != tj_period.end(); ++it_j) {
        int examj = *it_j;
        n = _chrom.getConflictMatrix()->getVal(exami, examj);
        if (n > 0) {
            // Exam examj have conflicts with exam ei
            conflictingExams_tj.push_back(examj);
        }
    }

//    cout << "#exams in slot ti = " << ti_period.size() << endl;
//    cout << "#exams in slot tj = " << tj_period.size() << endl;

//    cout << "Exams's list in slot ti = " << endl;
//    copy(ti_period.begin(), ti_period.end(), ostream_iterator<int>(cout, "\n"));

//    cout << "Exams's list in slot tj = " << endl;
//    copy(tj_period.begin(), tj_period.end(), ostream_iterator<int>(cout, "\n"));

//    cout << "Conflicting exams in slot tj = " << endl;
//    copy(conflictingExams_tj.begin(), conflictingExams_tj.end(), ostream_iterator<int>(cout, "\n"));

    vector<int> conflictingExams_ti;
    for (vector<int>::const_iterator it_conflict = conflictingExams_tj.begin(); it_conflict != conflictingExams_tj.end(); ++it_conflict) {
        for (unordered_set<int>::const_iterator it_i = ti_period.begin(); it_i != ti_period.end(); ++it_i) {
            int exam_conflict = *it_conflict;
            if (*it_i != exami) {
                n = _chrom.getConflictMatrix()->getVal(exam_conflict, *it_i);
                if (n > 0) {
                    // Verify if it's inserted already
                    if (find_if(conflictingExams_ti.begin(), conflictingExams_ti.end(),
                                bind2nd(equal_to<int>(), *it_i)) == conflictingExams_ti.end()) {
                        // Exam exam_conflict have conflicts with exam in slot ti
                        conflictingExams_ti.push_back(*it_i);
                    }
                }
            }
        }
    }
    // Consider also exam exami in slot ti
    conflictingExams_ti.push_back(exami);

//    cout << "Conflicting exams in slot ti = " << endl;
//    copy(conflictingExams_ti.begin(), conflictingExams_ti.end(), ostream_iterator<int>(cout, "\n"));


    // Move conflicting exams from slot ti
//    for (vector<int>::const_iterator it_conflict = conflictingExams_ti.begin(); it_conflict != conflictingExams_ti.end(); ++it_conflict) {
        // Remove from slot ti
//        cout << "*it_conflict = " << *it_conflict << " found: " << (ti_period.find(*it_conflict) != ti_period.end()) <<  endl;
//        ti_period.erase(ti_period.find(*it_conflict));

        // Insert in slot tj
//        tj_period.insert(*it_conflict);
//    }

    // Remove from slot ti
    ti_period.erase(exami);
    // Insert just exami in slot tj and then repeat the process for the remainder conflicting exams
    tj_period.insert(exami);
    // Remove exami from conflicting exams
    conflictingExams_ti.erase(find_if(conflictingExams_ti.begin(), conflictingExams_ti.end(), bind2nd(equal_to<int>(), exami)));

    // Move conflicting exams from slot tj
    for (vector<int>::const_iterator it_conflict = conflictingExams_tj.begin(); it_conflict != conflictingExams_tj.end(); ++it_conflict) {
        // Remove from slot tj
        tj_period.erase(tj_period.find(*it_conflict));
        // Insert in slot ti
        ti_period.insert(*it_conflict);
    }

//    cout << "After Kempe chain movement" << endl;
//    cout << "#exams in slot ti = " << ti_period.size() << endl;
//    cout << "#exams in slot tj = " << tj_period.size() << endl;

//    cout << "Exams's list in slot ti = " << endl;
//    copy(ti_period.begin(), ti_period.end(), ostream_iterator<int>(cout, "\n"));

//    cout << "Exams's list in slot tj = " << endl;
//    copy(tj_period.begin(), tj_period.end(), ostream_iterator<int>(cout, "\n"));

    for (vector<int>::const_iterator it_conflict = conflictingExams_ti.begin(); it_conflict != conflictingExams_ti.end(); ++it_conflict) {
        kempeChainOperator(_chrom, ti_period, tj_period, *it_conflict);
    }
}







