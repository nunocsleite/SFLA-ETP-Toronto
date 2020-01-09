
#include "Chromosome.h"
#include "moeoETTPEval.h"


// TODO
//
//
//% % Compute the number of clashes of solution x
//% % Sum_{i=1}^{|E|-1} Sum_{j=i+1}^{|E|} Sum_{p=1}^{|P|-1} aip aj(p+1) cij
//% function x = computeNumConflicts(Data, x)
//%     x.NumConflicts = 0;
//%     examListP1 = x.Periods{1};
//%     for p = 2 : x.NumPeriods
//%         examListP2 = x.Periods{p};
//% %         if (~(p == 7 || p == 13 || p == 19))
//%             numExamsP1 = length(examListP1);
//%             numExamsP2 = length(examListP2);
//%             for i = 1 : numExamsP1
//%                 for j = 1 : numExamsP2
//%                     numStudents = Data.ConflictMatrix(examListP1(i), examListP2(j));
//%                     if (numStudents > 0)
//%                         x.NumConflicts = x.NumConflicts + numStudents;
//%                     end
//%                 end
//%             end
//% %         end
//%         examListP1 = examListP2;
//%     end
//% end


void Chromosome::init(Data const& data) {
    // Set range
    range = data.getRange();
    // Set original number of Periods
    originalNumPeriods = data.getOriginalNumPeriods();
    // Set conflict matrix
    conflictMatrix = data.getProblemData().getConflictMatrix();
    // Set adjacency list
    graph = data.getProblemData().getGraph();
    // Set number of students
    numStudents = data.getProblemData().getNumStudents();
    // Set number of exams
    numExams = data.getProblemData().getNumExams();
    // Set number of enrolments
    numEnrolments = data.getProblemData().getNumEnrolments();
    // Set student counts
    courseStudentCounts = data.getProblemData().getCourseStudentCounts();
    proximityCost = 0;


    // Generate a random number of periods in the range interval
    int numPeriods = rand() % (data.getRange()[1] - data.getRange()[0] + 1) + data.getRange()[0];
//    int numPeriods = 19;
/// SEE

    //boost::random::mt19937 rng;  // random generator
    //boost::random::uniform_int_distribution<> distribution(data.getRange()[0], data.getRange()[1]); // distribution that maps to getRange()[0]..getRange()[1]
    //numPeriods = distribution(rng);

//    periods = vector<vector<int> >(numPeriods);
    periods = vector<unordered_set<int> >(numPeriods);

}


void Chromosome::baseComputeProximityCosts() {
	// The problem has one hard constraint where conflicting exams cannot be assigned to the same time slot. 
	// In addition, a soft constraint is present where conflicting exams should be spread throughout the timetable 
	// as far as possible from each other. The objective here is to minimise the sum of proximity costs given as:
	//    sum_{i=0}^{4} (wi x n) / S
	// where
	// - wi = 2^{4-i} is the cost of assigning two exams with i time slots apart. Only exams with
	//   common students and are four or less time slots apart are considered as violations
	// - n is the number of students involved in the conflict
	// - S is the total number of students in the problem
//    int pj;
//	int wi, n;
//    // For each period do
//    for (int pi = 0; pi < periods.size(); ++pi) {
//        // For each exam do
//        for (int ei = 0; ei < periods[pi].size(); ++ei) {
//            // Consider exams four periods apart
//			for (int i = 0; i <= 4; ++i) {
//                pj = pi+i+1;
//				if (pj < periods.size()) {
//                    for (int ej = 0; ej < periods[pj].size(); ++ej) {
////                        int exami = periods[pi][ei];
////                        int examj = periods[pj][ej];
//                        int exami = *periods[pi].find(ei);
//                        int examj = *periods[pj].find(ej);
//                        n = conflictMatrix->getVal(exami, examj);
//						if (n > 0) {
//							wi = (int)pow(2.0, 4.0-i);
//                            proximityCost += wi*n;
//						}
//					}
//				}
//			}
//		}
//	}

    proximityCost = 0;

    int pj;
    int wi, n;
    // For each period do
    for (int pi = 0; pi < periods.size(); ++pi) {
        // For each exam do
//        for (int ei = 0; ei < periods[pi].size(); ++ei) {
        for (unordered_set<int>::const_iterator it_i = periods[pi].begin(); it_i != periods[pi].end(); ++it_i) {
            // Consider exams four periods apart
            for (int i = 0; i <= 4; ++i) {
                pj = pi+i+1;
                if (pj < periods.size()) {
//                    for (int ej = 0; ej < periods[pj].size(); ++ej) {
                    for (unordered_set<int>::const_iterator it_j = periods[pj].begin(); it_j != periods[pj].end(); ++it_j) {
//                        int exami = periods[pi][ei];
//                        int examj = periods[pj][ej];
                        int exami = *it_i;
                        int examj = *it_j;
                        n = conflictMatrix->getVal(exami, examj);
                        if (n > 0) {
                            wi = (int)pow(2.0, 4.0-i);
                            proximityCost += wi*n;
                        }
                    }
                }
            }
        }
    }

    proximityCost /= numStudents;

//    cout << "numStudents = " << numStudents << endl;
//    cout << "proximityCost = " << proximityCost << endl;

//    cin.get();
}



//// Get number of clashes of exam in the given period.
//// A clash is the conflict of one student having an exam in two consecutive periods.
//int Chromosome::getNumClashesExamPeriod(int exami, int pi) const {
//    // If it's the last period, there's no clashes
//    if (pi == getNumPeriods())
//        return 0;

//    int pj, numClashes = 0, n;
//    pj = pi+1;
//    if (pj < periods.size()) {
//        // Compute conflicts with exams in period pi+1
////        for (int ej = 0; ej < periods[pj].size(); ++ej) {
//        for (unordered_set<int>::const_iterator it_j = periods[pj].begin(); it_j != periods[pj].end(); ++it_j) {
////            int examj = periods[pj][ej];
//            int examj = *it_j;

//            n = conflictMatrix->getVal(exami, examj);
//            if (n > 0)
//                numClashes += n;
//        }
//    }
//    return numClashes;
//}



// Get number of clashes of exam in the given period.
// A clash is the conflict of one student having an exam in two consecutive periods.
double Chromosome::getProximityCostExamPeriod(int exami, int pi) const {
    // If it's the last period, there's no clashes
    if (pi == getNumPeriods())
        return 0;

    int pj, wi, n, numStudentsPeriod = 0;

    /// TODO: SHOULD INCLUDE PROXIMITY COST OF PERIOD PI?
    double proxCost = 0.0;

    // Consider exams four periods apart
    for (int i = 0; i <= 4; ++i) {
        pj = pi+i+1;
        if (pj < periods.size()) {
            // Compute conflicts with exams in period pi+1
    //        for (int ej = 0; ej < periods[pj].size(); ++ej) {
            for (unordered_set<int>::const_iterator it_j = periods[pj].begin(); it_j != periods[pj].end(); ++it_j) {
    //            int examj = periods[pj][ej];
                int examj = *it_j;

                n = conflictMatrix->getVal(exami, examj);
                if (n > 0) {
                    wi = (int)pow(2.0, 4.0-i);
                    proxCost += wi*n;
//                    numStudentsPeriod += getNumStudentsPeriod(pj); // ACHO QUE NAO ESTA BEM
                }
            }
            // ADDED
            numStudentsPeriod += getNumStudentsPeriod(pj); // AQUI SIM
        }
    }
    if (numStudentsPeriod > 0)
        proxCost /= numStudentsPeriod;
    return proxCost;
}

/*

// Determine the exam with the highest proximity cost
void Chromosome::getInfoExamHighestProximityCost(int& examHighCost, int& period, double& highProxCost) {
    highProxCost = 0.0;

    period = -1;
    examHighCost = -1;
    // For each period do
    for (int pi = 0; pi < periods.size(); ++pi) {
        int exam = 1;
        double proxCost = 0.0;

        getExamHighestProximityCost(pi, exam, proxCost); // ref examHighCost, ref highProxCost
        // Update highest proximity cost exam
        if (proxCost >= highProxCost) {
            highProxCost = proxCost;
            examHighCost = exam;
            period = pi;
        }
    }
}


// Determine the exam with the highest proximity cost within given period
void Chromosome::getExamHighestProximityCost(int pi, int& exam, double& highProxCost)  {
    int pj, wi, n, numStudentsPeriod = 0;
    double proxCost;
    for (unordered_set<int>::const_iterator it_i = periods[pi].begin(); it_i != periods[pi].end(); ++it_i) {
        int exami = *it_i;
        proxCost = 0.0;
        //
        // Compute proximity cost for exami
        //
        // Consider exams four periods apart
        for (int i = 0; i <= 4; ++i) {
            pj = pi+i+1;
            if (pj < periods.size()) {
                // Compute conflicts between exami and exams in period pj
                for (unordered_set<int>::const_iterator it_j = periods[pj].begin(); it_j != periods[pj].end(); ++it_j) {
                    int examj = *it_j;
                    n = conflictMatrix->getVal(exami, examj);
                    if (n > 0) {
                        wi = (int)pow(2.0, 4.0-i);
                        proxCost += wi*n;
                    }
                }
                numStudentsPeriod += getNumStudentsPeriod(pj);
            }
        }
        if (numStudentsPeriod > 0)
            proxCost /= numStudentsPeriod; // Proximity cost for exami

        // Update highest proximity cost exam
        if (proxCost >= highProxCost) {
            highProxCost = proxCost;
            exam = exami;
        }
    }
}

*/


//// Number of clashes computation
//void Chromosome::computeNumClashes() {
//    int numberClashes = 0;
//    // For each period do
//    for (int pi = 0; pi < periods.size(); ++pi) {
//        numberClashes += getNumClashesPeriod(pi);
//    }
//    updateNumClashes(numberClashes);
//}


// Get number of clashes of the current "day" (or period).
// A clash is the conflict of one student having an exam in two consecutive periods.
int Chromosome::getNumClashesPeriod(int pi) const {
    // If it's the last period, there's no clashes
    if (pi == getNumPeriods())
        return 0;

    int pj, numClashes = 0, n;
    // For each exam in period pi
//    for (int ei = 0; ei < periods[pi].size(); ++ei) {
    for (unordered_set<int>::const_iterator it_i = periods[pi].begin(); it_i != periods[pi].end(); ++it_i) {
        pj = pi+1;
        if (pj < periods.size()) {
            // Compute conflicts with exams in period pi+1
//            for (int ej = 0; ej < periods[pj].size(); ++ej) {
//                int exami = periods[pi][ei];
//                int examj = periods[pj][ej];
            for (unordered_set<int>::const_iterator it_j = periods[pj].begin(); it_j != periods[pj].end(); ++it_j) {
                int exami = *it_i;
                int examj = *it_j;

                n = conflictMatrix->getVal(exami, examj);
                if (n > 0)
                    numClashes += n;
            }
        }
    }
    return numClashes;
}


// Get number of students of the current "day" (or period)
int Chromosome::getNumStudentsPeriod(int pi) const {
    int numberStudents = 0;
    int exami;
    // For each exam in period pi
//    for (int ei = 0; ei < periods[pi].size(); ++ei) {
    for (unordered_set<int>::const_iterator it_i = periods[pi].begin(); it_i != periods[pi].end(); ++it_i) {
//        exami = periods[pi][ei];
        int exami = *it_i;
        numberStudents += (*courseStudentCounts).at(exami);
    }
    return numberStudents;
}



/////////////////////////////////////////////////////////////////////
// Timetable packing auxiliary methods
/////////////////////////////////////////////////////////////////////

vector<int> Chromosome::getFeasiblePeriods(int exam, int period) const {
    vector<int> feasiblePeriods;
    // Determine possible feasible periods where exam can be placed
    for (int p = 0; p < getNumPeriods(); ++p) {
        if (p != period && isFeasiblePeriodExam(p, exam))
            feasiblePeriods.push_back(p);
    }
    return feasiblePeriods;
}


bool Chromosome::isFeasiblePeriodExam(int period, int exam) const {
    // Verify exam period feasibility.
    // Constraint that no student is to be scheduled
    // to take two exams at any one time:
    //   Sum_{i=1}^{|E|-1} Sum_{j=i+1}^{|E|} Sum_{p=1}^{|P|} aip ajp cij = 0
//    vector<int> const& examList = periods[period];
    unordered_set<int> const& examList = periods[period];

    bool feasible = true;
//    for (int i = 0; i < examList.size(); ++i) {
    for (unordered_set<int>::const_iterator it_i = examList.begin(); it_i != examList.end(); ++it_i) {
//        int numStudents = (*conflictMatrix).getVal(examList[i], exam);
        int numStudents = (*conflictMatrix).getVal(*it_i, exam);
        if (numStudents > 0) {
            feasible = false;
            break;
        }
    }
    return feasible;
}


// Get the period with the smallest number of students
int Chromosome::getMinPeriod() {
    int minPeriod = 0, count;
    int minNumSudents = getNumStudentsPeriod(minPeriod);
    int numPeriods = getNumPeriods();
    for (int p = 1; p < numPeriods; ++p) {
        count = getNumStudentsPeriod(p);
        if (count < minNumSudents) {
            minNumSudents = count;
            minPeriod = p;
        }
    }
    return minPeriod;
}

//multiset<Period> Chromosome::getExamsAvailablePeriods(vector<int> const& sourceExams, int sourcePeriod) {
multiset<Period> Chromosome::getExamsAvailablePeriods(unordered_set<int> const& sourceExams, int sourcePeriod) {

    // Initialize vector in order to determine period capacity
    vector<Period> availablePeriodsVec(getNumPeriods());
    int i = 0;
    // Initialize period's Id
    for (vector<Period>::iterator it = availablePeriodsVec.begin(); it != availablePeriodsVec.end(); ++it, ++i) {
        (*it).setId(i);
    }

//    cout << "Initial periods" << endl;
//    cout << "size = " << availablePeriodsVec.size() << endl;
//    for (vector<Period>::iterator it = availablePeriodsVec.begin(); it != availablePeriodsVec.end(); ++it)
//        cout << *it << endl;

    // Determine available period capacity by filling 'availablePeriodsVec' vector
    computeAvailablePeriods(sourceExams, sourcePeriod, availablePeriodsVec);
    // Create sorted set
    multiset<Period> availablePeriodsSet;
    // Insert periods into the sorted set
    for (vector<Period>::iterator it = availablePeriodsVec.begin(); it != availablePeriodsVec.end(); ++it) {
        if ((*it).getCapacity() > 0)
            availablePeriodsSet.insert((*it));
    }
    return availablePeriodsSet;
}


// Compute available periods.
// A given p period's capacity is defined as the number of exams from
// the set 'sourceExams', the 'sourcePeriod' exams, which can scheduled into p
// without causing any clashes while maintaining feasibility.
//void Chromosome::computeAvailablePeriods(vector<int> const& sourceExams, int sourcePeriod,
//                                         vector<Period>& availablePeriodsVec) {

void Chromosome::computeAvailablePeriods(unordered_set<int> const& sourceExams, int sourcePeriod,
                                         vector<Period>& availablePeriodsVec) {

    int numClashes, exam;
    int numPeriods = getNumPeriods();
    // Compute each period availability
    for (int p = 0; p < numPeriods; ++p) {
        if (p != sourcePeriod) {
//            for (int i = 0; i < sourceExams.size(); ++i) {
            for (unordered_set<int>::const_iterator it_i = sourceExams.begin(); it_i != sourceExams.end(); ++it_i) {
//                exam = sourceExams[i];
                exam = *it_i;
                // Verify if period 'p' is a feasible period for the current exam
                if (isFeasiblePeriodExam(p, exam)) {
                    // Compute number of clashes between this exam and other exams
                    // in periods p-1 and p+1, excluding sourcePeriod
                    numClashes = computeNumClashesExamPeriod(exam, p, sourcePeriod);
                    if (numClashes == 0) {
                        // Add current exam to period vector which increments by one the period capacity
                        boost::shared_ptr<Exam> pExam(new Exam(exam));
                        availablePeriodsVec[p].addExam(pExam);
                    }
                }
            }
        }
    }
}


// Compute number of clashes between 'exam'
// and other exams in period p
int Chromosome::computeNumClashesPeriod(int exam, int p) const {
    int numClashes = 0;
    // Get exam's list
//    vector<int> const& examList = periods[p];
    unordered_set<int> const& examList = periods[p];
//    for (int i = 0; i < examList.size(); ++i) {
    for (unordered_set<int>::const_iterator it_i = examList.begin(); it_i != examList.end(); ++it_i) {
//        int numStudents = conflictMatrix->getVal(examList[i], exam);
        int numStudents = conflictMatrix->getVal(*it_i, exam);
        if (numStudents > 0)
            numClashes += numStudents;
    }
    return numClashes;
}


// Compute number of clashes between 'exam' and other exams
// in periods p-1 and p+1, excluding sourcePeriod
int Chromosome::computeNumClashesExamPeriod(int exam, int p, int sourcePeriod) const {

    /// TODO -> In Saturdays there's no clashes with the following period

    int numClashes0 = 0, numClashes1 = 0;
    if (p-1 >= 0 && p-1 != sourcePeriod) {
        numClashes0 = computeNumClashesPeriod(exam, p-1);
    }
    if (p+1 < getNumPeriods() && p+1 != sourcePeriod) {
        numClashes1 = computeNumClashesPeriod(exam, p+1);
    }
    return numClashes0 + numClashes1;
}


struct size_pred {
//    bool operator()(vector<int> const& periodExams) { return periodExams.size() == 0; }
    bool operator()(unordered_set<int> const& periodExams) { return periodExams.size() == 0; }
};

void Chromosome::removeEmptyPeriods() {

    // TODO -> OPTIMIZAR

//    bool found;
//    do {
//        found = false;
//        for (vector<vector<int> >::iterator it = periods.begin(); it != periods.end(); ++it) {
//            if ((*it).size() == 0) {
//                periods.erase(it);
//                found = true;
//                break;
//            }
//        }
//    }
//    while (found);
//    cout << "Remove empty periods" << endl;

    periods.erase(std::remove_if(periods.begin(), periods.end(), size_pred()), periods.end());


}



// Timetable period packing
void Chromosome::pack() {
    // Period packing: Starting from the period with the
    // smallest number of students, the operation searches in order
    // of available period capacity, starting from the smallest,
    // for a period which can accommodate exams from the former
    // without causing any clashes while maintaining feasibility.
    // The operation stops when it goes one cycle through all periods
    // without rescheduling any exam or when the timetable
    // length is reduced to a random number within the desired
    // range.
    //////

//    cout << "Packing timetable" << endl;


    // Remove minPeriod
//    removeEmptyPeriods();  /// TODO: VER EFEITO DESTA LINHA

    bool packed = false;
    while (!packed) {
        bool scheduleAnyExam = false;
        // Random timetable length
        int timetableLength = eo::random(getRange()[0], getRange()[1]);
        // Get the period with the smallest number of students
        int minPeriod = getMinPeriod();
//        vector<int> minExams = periods[minPeriod];
        unordered_set<int> minExams = periods[minPeriod];
        // Determine exams possible moves
//        multiset<Period>& availablePeriodsSet = *(getExamsAvailablePeriods(minExams, minPeriod).get());
        multiset<Period> availablePeriodsSet = getExamsAvailablePeriods(minExams, minPeriod);
        // The operation searches in order of available period capacity,
        // starting from the smallest, for a period which can accommodate exams
        // from the former without causing any clashes while maintaining
        // feasibility.
        if (availablePeriodsSet.size() > 0) {
//            cout << endl << "Reschedule exams" << endl << endl;



            // Reschedule exams
            //
            // Randomly select a period with available period capacity, starting from the smallest
            // 1. Get range of equal capacity periods
//            pair<multiset<Period>::iterator, multiset<Period>::iterator> rangeInfo;
//            rangeInfo = availablePeriodsSet.equal_range(*availablePeriodsSet.begin());
//            // 2. Select one period randomly, which returns an index in the interval [0, upperLim-lowerLim)
//            int dist = distance(rangeInfo.first, rangeInfo.second);
//            int randPeriodIdx = eo::random(dist);
//            // Reschedule exams and remove them from minPeriod
//            int i = 0;
            Period const* targetPeriod = 0;
//            for (multiset<Period>::iterator it = availablePeriodsSet.begin(); it != availablePeriodsSet.end(); ++it, ++i) {
//                if (i == randPeriodIdx) {
//                   targetPeriod = &*it;
//                   break;
//                }
//            }
//            assert(targetPeriod == 0);
//            if (targetPeriod == 0) {
//                cout << "targetPeriod == 0"  << endl;
//            }
            targetPeriod = &*availablePeriodsSet.begin();

//            cout << *targetPeriod << endl;


            vector<boost::shared_ptr<Exam> > const& examsToSchedule = targetPeriod->getExams();
            for (int i = 0; i < examsToSchedule.size(); ++i) {
                Exam* exam = examsToSchedule[i].get();

//                cout << exam->getId() << endl;

                if (!exam->isAlreadyScheduled()) {

//                    cout << "Schedule exam" << endl;
                    // Remove exam from minPeriod

                    // TODO -> OPTIMIZAR
//                    vector<int>::iterator removeItr =
//                            find_if(periods[minPeriod].begin(), periods[minPeriod].end(),
//                                    bind2nd(equal_to<int>(), exam->getId()));
                    unordered_set<int>::iterator removeItr =
                            find_if(periods[minPeriod].begin(), periods[minPeriod].end(),
                                    bind2nd(equal_to<int>(), exam->getId()));

                    if (removeItr == periods[minPeriod].end()) {
                        cout << "could not find exam" << endl;
                        cout << exam->getId() << endl;
                        copy(periods[minPeriod].begin(), periods[minPeriod].end(), ostream_iterator<int>(cout, " "));
                        cout << endl;
                        cin.get();
                    }
                    periods[minPeriod].erase(removeItr);

                    // Reschedule exam into target period
//                    periods[targetPeriod->getId()].push_back(exam->getId());
                    periods[targetPeriod->getId()].insert(exam->getId());
                    scheduleAnyExam = true;
                }
            }
            if (periods[minPeriod].size() == 0) {
                cout << "Remove empty periods" << endl;
//                cin.get();
                // Remove minPeriod
                removeEmptyPeriods();
            }
        }

//        scheduleAnyExam??

        // The operation stops when it goes one cycle through all periods
        // without rescheduling any exam or when the timetable
        // length is reduced to a random number within the desired range.
//        if (!scheduleAnyExam || getNumPeriods() == timetableLength)
        if (!scheduleAnyExam)
            packed = true;

    }
    // Recompute number of clashes
    baseComputeProximityCosts();
}


void Chromosome::expand() {
    cout << "Expanding timetable" << endl;
}

void Chromosome::validate() const {
    int numExams = 0, numEnrolments = 0;
    unordered_map<int, int> exams;
    bool unique = true;
    for (ConstPeriodIterator it = periods.begin(); it != periods.end(); ++it) {
//        for (vector<int>::const_iterator examIt = (*it).begin(); examIt != (*it).end(); ++examIt) {
        for (unordered_set<int>::const_iterator examIt = (*it).begin(); examIt != (*it).end(); ++examIt) {
            ++numExams;
            numEnrolments += (*courseStudentCounts).at(*examIt);
            // Verify unicity of exams
            pair<unordered_map<int, int>::iterator, bool> p =
                    exams.insert(pair<int, int>(*examIt, (*courseStudentCounts).at(*examIt)));
            if (p.second == false) {
                unique = false;
                break;
            }
        }
    }

//    cout << "numExams = " << numExams << endl;
//    cout << "numEnrolements = " << numEnrolments << endl;

//    cout << "Original values" << endl;
//    cout << "numExams = " << this->numExams << endl;
//    cout << "numEnrolements = " << this->numEnrolments << endl;

    if (!unique) {
        cout << "There are duplicated exams" << endl;
        cin.get();
    }
    else if (numExams != this->numExams) {
        cout << "numExams is different" << endl;
        cout << "counted: " << numExams << "; original: " <<  this->numExams << endl;

        cin.get();
    }
    else if (numEnrolments != this->numEnrolments) {
        cout << "numEnrolments is different" << endl;
        int count = 0, enrol = 0;
        // Exams indexed from [1..numExams]
        vector<int> sorted_exams(this->numExams+1);
        for (unordered_map<int, int>::iterator it = exams.begin(); it != exams.end(); ++it) {
//            cout << (*it).first << " " << (*it).second << endl;
            ++count;
            enrol += (*it).second;
            sorted_exams[(*it).first] = (*it).second;
        }
        cout << endl << "count = " << count << ", enrol = " << enrol << endl;
        copy(sorted_exams.begin(), sorted_exams.end(), ostream_iterator<int>(cout, "\n"));

        cin.get();
    }

    // Verify feasibility
    int i = 0;
    for (ConstPeriodIterator it = periods.begin(); it != periods.end(); ++it, ++i) {
        for (unordered_set<int>::const_iterator examIt = (*it).begin(); examIt != (*it).end(); ++examIt) {
            if (!isFeasiblePeriodExam(i, *examIt)) {
                cout << endl << "Not FeasiblePeriodExam" << endl;
                cin.get();
            }
        }
    }

}


ostream& operator<<(ostream& os, const Chromosome& timetable) {
    os << endl << "Timetable" << endl;
//    int i = 1;
//    for (vector<unordered_set<int> >::const_iterator it = timetable.getConstPeriods().begin(); it != timetable.periods.end(); ++it) {
//        os << endl << "Period " << i << endl;
//        for (unordered_set<int>::const_iterator examit = (*it).begin(); examit != (*it).end(); ++examit) {
//            //os << timetable.deetc_ucs[*examit-1] << " ";
//            os << *examit << " ";
//        }
//        os << endl;
//        i++;
//    }


//    os << "NumClashes = " << timetable.getNumClashes() << " - Periods = " << timetable.getNumPeriods() << endl;
//    os << "Proximity cost = " << timetable.getProximityCost() << " - Periods = " << timetable.getNumPeriods();
//    if (timetable.getNumPeriods() == timetable.getOriginalNumPeriods()) {
//        os << " <-- " << endl;
//    }
//    else
//        os << endl;

    timetable.validate();

    os << endl;

//    if (timetable.getNumPeriods() == timetable.getOriginalNumPeriods())
//        os << " & " << timetable.getProximityCost() << " & " << timetable.getNumPeriods() << " \\\\ " << endl;

    os << " cost = " << timetable.getProximityCost() << " - periods = " << timetable.getNumPeriods() << endl;



    os << endl << ".Sol Data" << endl;
    int period = 0;
    for (vector<unordered_set<int> >::const_iterator it = timetable.getConstPeriods().begin(); it != timetable.periods.end(); ++it) {
        for (unordered_set<int>::const_iterator examit = (*it).begin(); examit != (*it).end(); ++examit) {
            os << *examit << " " << period << endl;
        }
        ++period;
    }

    return os;
}

ostream& operator<<(ostream& os, const eoChromosome& timetable) {
    os << timetable.getChromosome();
    return os;
}


ostream& operator<<(ostream& os, const moeoChromosome& timetable) {
    os << timetable.getChromosome();
    return os;
}

















