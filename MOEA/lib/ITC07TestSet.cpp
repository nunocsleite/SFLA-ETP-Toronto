

#include "ITC07TestSet.h"
#include <boost/regex.hpp>
#include <boost/unordered_map.hpp>

/// TODO

void ITC07TestSet::load(string _testSetName, string _rootDir) {
    /*
    Input Format

    The problem instance files have the following format;

    Number of Exams:
    Number of exams to be timetabled e.g. [Exams:2176]. As with all other entities the numbering of Exams starts at 0.
    Sequence of lines detailing information regarding each exam:
    The first number is the duration of the examination (specified in minutes). This is followed by the student numbers of
    those students taking the exam. Students numbers are integers starting with 0.
    The line ends with a return character and line feed and is comma separated.

    Number of Periods:
    Number of Periods specified within the Timetabling Session e.g. [Periods:42]
    Sequence of lines detailing Period Dates, Times, Durations associated Penalty:
    E.g.  31:05:2005, 09:00:00, 180, 0.   The Date is in a standard date format. The time is in 24 hour format, the duration
    is in minutes and the penalty has a value which is a positive integer. It should be noted that 0 represents no penalty.
    In relation to other periods a relative value may be present here which indicates that a penalty is added when this period
    is used by placing an exam. The higher the relative value the least the period should be used.

    Number of Rooms:
    Number of Rooms to be used e.g.. [Rooms:10]
    Sequence of lines detailing room capacity and associated penalty:
    Room penalty has a value which is a positive integer. It should be noted that 0 represents no penalty. In relation to other
    rooms a relative value may be present here which indicates that a penalty is added when this room is used by placing an exam.
    The higher the relative value the least the room should be used.
    E.g. 500, 7
    Rooms are numbered in the order from the problem file, starting at zero

    Period Related Hard Constraints
    This section begins with the tag [PeriodHardConstraints] and provides data on conditions which are necessary for a feasible solution.
    These are EXAM_COINCIDENCE, EXCLUSION, and AFTER.  It should be noted that not all these may be present for any one data set.
    E.g.
    0, EXAM_COINCIDENCE, 1   Exam ‘0' and Exam ‘1' should be timetabled in the same period. If two exams are set associated in
                             this manner yet 'clash' with each other due to student enrolment, this hard constraint is ignored.
    0, EXCLUSION, 2          Exam ‘0' and Exam ‘2' should be not be timetabled in the same period
    0, AFTER, 3              Exam ‘0' should be timetabled after Exam ‘3'

    Room Related Hard Constraints
    This section begins with the line [RoomHardConstraints] and provides data on conditions which are necessary for a feasible solution.
    This is
    ROOM_EXCLUSIVE.  An exam must be timetabled in a room by itself e.g.
    2, ROOM_EXCLUSIVE      Exam ‘2' must be timetabled in a room by itself.

    Institutional Model Weightings
    This section begins with the line [InstitutionalWeightings] and provides information on values given to 'global' soft constraints.
    TWOINAROW, 7
    TWOINADAY, 5
    PERIODSPREAD, 3
    NONMIXEDDURATIONS, 10
    FRONTLOAD, 100, 30, 5

    These are all fully explained and illustrated in the evaluation section.

    It should be noted that the format of the input file should not be altered in any way.
    In addition, it is recomended that competitors should ignore unknown lines in the provided format.
    */
    string filename = _rootDir + "/" + _testSetName;
    string wholeFile = read_from_file(filename.c_str());
    // Get line tokenizer
    tokenizer<escaped_list_separator<char> > tok(wholeFile, escaped_list_separator<char>('\\', '\r')); // Lines terminate with \r\n
    // Get line iterator
    tokenizer<escaped_list_separator<char> >::iterator it = tok.begin();
    // Read exams and students
    readExams(it, tok);
    // Read periods
    readPeriods(it, tok);
    // Read rooms
    readRooms(it, tok);
    // Read constraints and weightings
    readConstraints(it, tok);

}


/////////////////////////////////////
// Read exams and students
/////////////////////////////////////
void ITC07TestSet::readExams(tokenizer<escaped_list_separator<char> >::iterator& it,
                             tokenizer<escaped_list_separator<char> > const& tok) {
    //
    // The problem instance files have the following format;
    // Number of Exams:
    // Number of exams to be timetabled e.g. [Exams:2176]. As with all other entities the numbering of Exams starts at 0.
    // Sequence of lines detailing information regarding each exam:
    // The first number is the duration of the examination (specified in minutes). This is followed by the student numbers of
    // those students taking the exam. Students numbers are integers starting with 0.
    // The line ends with a return character and line feed and is comma separated.
    ////

    // Number of exams
    int numExams = 0;
    if (it != tok.end()) {
        // Read number of exams
        std::string s(*it);
        boost::cmatch matches;
        boost::regex re("(\\[Exams:)(\\d+)(\\])");

        if (boost::regex_match(s.c_str(), matches, re)) {
//            cout << "\tmatches" << endl;
//            // matches[0] contains the original string.  matches[n]
//            // contains a sub_match object for each matching
//            // subexpression
//            for (int i = 1; i < matches.size(); i++)
//            {
//               // sub_match::first and sub_match::second are iterators that
//               // refer to the first and one past the last chars of the
//               // matching subexpression
//               string match(matches[i].first, matches[i].second);
//               cout << "\tmatches[" << i << "] = " << match << endl;
//            }

            // Extract digits from string. The digits are the third subexpression, index = 2.
            string match(matches[2].first, matches[2].second);
            numExams = atoi(match.c_str());
            cout << "\tnumExams: " << numExams << endl;

            // Create shared ptr to manage Conflict Matrix
            boost::shared_ptr<Matrix> ptrMatrix(new Matrix(numExams, numExams));
//            ptrConflictMatrix = ptrMatrix;
            // Instantiate graph with ncols vertices
//            boost::shared_ptr<AdjacencyList> ptrGraphAux(new AdjacencyList(numExams+1));
//            ptrGraph = ptrGraphAux

            // Conflict matrix
            boost::shared_ptr<Matrix> ptrConflictMatrix = ptrMatrix;
            // Graph representing exam relations
//            boost::shared_ptr<AdjacencyList> ptrGraph = ptrGraphAux;
            // Obtain internal vector
            vector<int>& v = ptrConflictMatrix->getVec();
            /////////////////////////////////////////////////////////////////////
            // Define student map containing the list of exams for each student
            /////////////////////////////////////////////////////////////////////
            unordered_map<int, vector<int> > studentMap;
            int examDuration;
            // Go to the first line of exams
            ++it;
            for (int i = 0; i < numExams && it != tok.end(); ++i, ++it) {
                string line(*it);
                cout << i << ": " << line << endl; // Lines start with \n because it wasn't removed
                // Get number tokenizer
                tokenizer<escaped_list_separator<char> > numberTok(line, escaped_list_separator<char>('\\', ','));
                // Get number iterator
                tokenizer<escaped_list_separator<char> >::iterator beg = numberTok.begin();
                examDuration = atoi((*beg).c_str());
                cout << "examDuration: " << examDuration << " - Students: ";

                ++beg; // Go to the first student
                for (; beg != numberTok.end(); ++beg) {
                    cout << *beg << " ";
                    //////////////////////////////////////////
                    // Build student map
                    //////////////////////////////////////////
                    // Current student
                    int student = atoi((*beg).c_str());
                    // Exam number start at 1
                    int exam = i+1;
                    // Insert current exam, given by the index i, in the student map
                    unordered_map<int, vector<int> >::iterator mapIt = studentMap.find(student);
                    if (mapIt == studentMap.end()) {
                        vector<int> exams;
                        exams.push_back(exam);
                        studentMap.insert(pair<int, vector<int> >(student, exams));
                    }
                    else {
                        mapIt->second.push_back(exam);
                    }
                }
                cout << endl << endl;
            }

            // Print student info
            cout << "Print student info" << endl << endl;
            int numStudents = 0, smallest = INT_MAX, greatest = INT_MIN;
            for (auto s : studentMap) {
//                cout << "Student: " << s.first << " - total exams: " << s.second.size() << endl;
                ++numStudents;
                smallest = min(smallest, s.first);
                greatest = max(greatest, s.first);

            }
            cout << "numStudents = " << numStudents << endl;
            cout << "smallest = " << smallest << endl;
            cout << "greatest = " << greatest << endl;

            //////////////////////////////////////////
            // Build Conflict matrix and graph
            //////////////////////////////////////////
            // Reserve space for matrix
            ptrConflictMatrix->reserve(ptrConflictMatrix->nLines(), ptrConflictMatrix->nCols());

            for (auto entry : studentMap) {
                vector<int>& exams = entry.second;
                int examListSize = exams.size();
                for (int i = 0; i < examListSize; ++i) {
                    for (int j = 0; j < examListSize; ++j) {
                        if (j != i) {
                            ++v[exams[j]*ptrConflictMatrix->nCols() + exams[i]];
                        }
                    }
                }
            }
            // Count the number of non-zero elements
            int nonZeroElements = 0;
            for (auto elem : v) {
                if (elem != 0)
                    ++nonZeroElements;
            }

//            cout << "nlines = " << ptrConflictMatrix->nLines() << endl;
//            cout << "ncols = " << ptrConflictMatrix->nCols() << endl;

//            // Verify if it's symmetric
//            for (int i = 0; i < ptrConflictMatrix->nLines(); ++i) {
//                for (int j = 0; j < ptrConflictMatrix->nCols(); ++j) {
//                    if (ptrConflictMatrix->getVal(i, j) != ptrConflictMatrix->getVal(j, i))
//                        throw runtime_error("Not symmetric");
//                }
//            }

            // The ‘conflict’ density is the ratio of the number of non-zero elements
            // in the conflict matrix to the total number of conflict matrix elements.
            double numMatrixElements = ptrConflictMatrix->nCols() * ptrConflictMatrix->nCols();
            double conflictDensity = nonZeroElements / (numMatrixElements - ptrConflictMatrix->nCols());

            cout << "conflictDensity = " << conflictDensity << endl;
            cout << "conflictDensity [%] = " << setprecision(3) << (conflictDensity * 100) << endl;

            // Build exam graph representing exam relations
            boost::shared_ptr<AdjacencyList> ptrGraph = buildExamGraph(ptrConflictMatrix);

        }
        else
        {
            cout << "The regexp \"" << re << "\" does not match \"" << s << "\"" << endl;
            cin.get();
        }
    }
    else
        throw runtime_error("Error while reading ITC07 exams");
}



/////////////////////////////////////
// Read periods
/////////////////////////////////////
void ITC07TestSet::readPeriods(tokenizer<escaped_list_separator<char> >::iterator& it,
                               tokenizer<escaped_list_separator<char> > const& tok) {
    //
    // The problem instance files have the following format;
    // Number of Periods:
    // Number of Periods specified within the Timetabling Session e.g. [Periods:42]
    // Sequence of lines detailing Period Dates, Times, Durations and associated Penalty:
    // E.g.  31:05:2005, 09:00:00, 180, 0.   The Date is in a standard date format. The time is in 24 hour format, the duration
    // is in minutes and the penalty has a value which is a positive integer. It should be noted that 0 represents no penalty.
    // In relation to other periods a relative value may be present here which indicates that a penalty is added when this period
    // is used by placing an exam. The higher the relative value the least the period should be used.
    ////

    if (it != tok.end()) {
        // Match period number
        int numPeriods = matchPeriods(it);
        ++it; // Read next line
        // Read periods info
        for (int i = 0; i < numPeriods; ++i) {
            cout << "line " << (i+1) << ":" << endl;
            // Match a sequence line
            SeqLineInfo seqLineInfo;
            matchPeriodSequenceLine(it, seqLineInfo);
            // Do something with seqLineInfo....

            ++it; // Read next line
        }
    }
    else
        throw runtime_error("Error while reading ITC07 periods");
}



int ITC07TestSet::matchPeriods(tokenizer<escaped_list_separator<char> >::iterator& it) {
    // TODO - FIX

    /////////////////////////////////////////////////
    // Lines start with \n because it wasn't removed
    //
    // Remove first \n in the regexp
    /////////////////////////////////////////////////

    // Number of periods
    int numPeriods = 0;
    // Read number of periods specified within the Timetabling Session e.g. [Periods:42]
    std::string s(*it);
    boost::cmatch matches;
//        boost::regex re("(\\[Periods:)(\\d+)(\\])");
    boost::regex re("(\\\n)(\\[Periods:)(\\d+)(\\])");

    if (boost::regex_match(s.c_str(), matches, re)) {
        // Extract digits from string. The digits are the third subexpression, index = 3.
        string match(matches[3].first, matches[3].second);

        numPeriods = atoi(match.c_str());
        cout << "\tnumPeriods: " << numPeriods << endl;
    }
    else {
        cout << "The regexp \"" << re << "\" does not match \"" << s << "\"" << endl;
        cin.get();
    }

    return numPeriods;
}



int ITC07TestSet::matchPeriodSequenceLine(tokenizer<escaped_list_separator<char> >::iterator& it, SeqLineInfo& seqLineInfo) {
    // Read sequence of lines detailing Period Dates, Times, Durations and associated Penalty:
    //  E.g.  31:05:2005, 09:00:00, 180, 0.
    std::string s(*it);
//    std::string s("\n15:04:2005, 09:30:00");

    boost::cmatch matches;
//    boost::regex re("(\\\n)(\\d{2}):(\\d{2}):(\\d{4}), (\\d{2}):(\\d{2}):(\\d{2})");

    boost::regex re("(\\\n)(\\d{2}):(\\d{2}):(\\d{4}), (\\d{2}):(\\d{2}):(\\d{2}), (\\d+), (\\d+)");

    if (boost::regex_match(s.c_str(), matches, re)) {
        // Extract values from string. E.g., the day part is the second subexpression, index = 2.
        string matchDay(matches[2].first, matches[2].second);
        int day = atoi(matchDay.c_str());
        string matchMonth(matches[3].first, matches[3].second);
        int month = atoi(matchMonth.c_str());
        string matchYear(matches[4].first, matches[4].second);
        int year = atoi(matchYear.c_str());

        string matchHour(matches[5].first, matches[5].second);
        int hour = atoi(matchHour.c_str());
        string matchMinute(matches[6].first, matches[6].second);
        int minute = atoi(matchMinute.c_str());
        string matchSecond(matches[7].first, matches[7].second);
        int second = atoi(matchSecond.c_str());

        string matchDuration(matches[8].first, matches[8].second);
        int duration = atoi(matchDuration.c_str());

        string matchPenalty(matches[9].first, matches[9].second);
        int penalty = atoi(matchPenalty.c_str());

        cout << "\tday: " << day << endl;
        cout << "\tmonth: " << month << endl;
        cout << "\tyear: " << year << endl;
        cout << "\thour: " << hour << endl;
        cout << "\tminute: " << minute << endl;
        cout << "\tsecond: " << second << endl;
        cout << "\tduration: " << duration << endl;
        cout << "\tpenalty: " << penalty << endl;
        cout << endl;

        // Set struct fields
        seqLineInfo.day = day;
        seqLineInfo.month = month;
        seqLineInfo.year = year;

        seqLineInfo.hour = hour;
        seqLineInfo.minute = minute;
        seqLineInfo.second = second;

        seqLineInfo.duration = duration;
        seqLineInfo.penalty = penalty;

    }
    else {
        cout << "The regexp \"" << re << "\" does not match \"" << s << "\"" << endl;
        cin.get();
    }
}



/////////////////////////////////////
// Read rooms
/////////////////////////////////////
void ITC07TestSet::readRooms(tokenizer<escaped_list_separator<char> >::iterator& it,
                             tokenizer<escaped_list_separator<char> > const& tok) {
    //
    // The problem instance files have the following format;
    // Number of Rooms:
    // Number of Rooms to be used e.g.. [Rooms:10]
    // Sequence of lines detailing room capacity and associated penalty:
    // Room penalty has a value which is a positive integer. It should be noted that 0 represents no penalty. In relation to other
    // rooms a relative value may be present here which indicates that a penalty is added when this room is used by placing an exam.
    // The higher the relative value the least the room should be used.
    // E.g. 500, 7
    // Rooms are numbered in the order from the problem file, starting at zero
    ////

    if (it != tok.end()) {
        // Match rooms number
        int numRooms = matchRooms(it);
        ++it; // Read next line
        // Read rooms info
        for (int i = 0; i < numRooms; ++i) {
            cout << "line " << (i+1) << ":" << endl;
            // Match a sequence line
            SeqLineInfo seqLineInfo;
            matchRoomSequenceLine(it, seqLineInfo);
            // Do something with seqLineInfo....

            ++it; // Read next line
        }
    }
    else
        throw runtime_error("Error while reading ITC07 rooms");
}



int ITC07TestSet::matchRooms(tokenizer<escaped_list_separator<char> >::iterator& it) {
    // TODO - FIX

    /////////////////////////////////////////////////
    // Lines start with \n because it wasn't removed
    //
    // Remove first \n in the regexp
    /////////////////////////////////////////////////

    // Number of rooms
    int numRooms = 0;
    // Read number of rooms used e.g. [Rooms:10]
    std::string s(*it);
    boost::cmatch matches;
//        boost::regex re("(\\[Rooms:)(\\d+)(\\])");
    boost::regex re("(\\\n)(\\[Rooms:)(\\d+)(\\])");

    if (boost::regex_match(s.c_str(), matches, re)) {
        // Extract digits from string. The digits are the third subexpression, index = 3.
        string match(matches[3].first, matches[3].second);

        numRooms = atoi(match.c_str());
        cout << "\tnumRooms: " << numRooms << endl;
    }
    else {
        cout << "The regexp \"" << re << "\" does not match \"" << s << "\"" << endl;
        cin.get();
    }

    return numRooms;
}



int ITC07TestSet::matchRoomSequenceLine(tokenizer<escaped_list_separator<char> >::iterator& it, SeqLineInfo& seqLineInfo) {
    // Read sequence of lines detailing room capacity and associated penalty:
    //  E.g.  260, 0
    std::string s(*it);

    boost::cmatch matches;

    boost::regex re("(\\\n)(\\d+), (\\d+)");

    if (boost::regex_match(s.c_str(), matches, re)) {
        // Extract values from string.
        string matchCapacity(matches[2].first, matches[2].second);
        int roomCapacity = atoi(matchCapacity.c_str());

        string matchPenalty(matches[3].first, matches[3].second);
        int penalty = atoi(matchPenalty.c_str());

        cout << "\tcapacity: " << roomCapacity << endl;
        cout << "\tpenalty: " << penalty << endl;
        cout << endl;

        // Set struct fields
        seqLineInfo.roomCapacity = roomCapacity;
        seqLineInfo.penalty = penalty;
    }
    else {
        cout << "The regexp \"" << re << "\" does not match \"" << s << "\"" << endl;
        cin.get();
    }
}



/////////////////////////////////////
// Read constraints and weightings
/////////////////////////////////////
void ITC07TestSet::readConstraints(tokenizer<escaped_list_separator<char> >::iterator& it,
                                   tokenizer<escaped_list_separator<char> > const& tok) {

    //
    // The problem instance files have the following format:
    // Period Related Hard Constraints
    // This section begins with the tag [PeriodHardConstraints] and provides data on conditions which are necessary for a feasible solution.
    // These are EXAM_COINCIDENCE, EXCLUSION, and AFTER.  It should be noted that not all these may be present for any one data set.
    // E.g.
    // 0, EXAM_COINCIDENCE, 1   Exam ‘0' and Exam ‘1' should be timetabled in the same period. If two exams are set associated in
    //                          this manner yet 'clash' with each other due to student enrolment, this hard constraint is ignored.
    // 0, EXCLUSION, 2          Exam ‘0' and Exam ‘2' should be not be timetabled in the same period
    // 0, AFTER, 3              Exam ‘0' should be timetabled after Exam ‘3'

    // Room Related Hard Constraints
    // This section begins with the line [RoomHardConstraints] and provides data on conditions which are necessary for a feasible solution.
    // This is
    // ROOM_EXCLUSIVE.  An exam must be timetabled in a room by itself e.g.
    // 2, ROOM_EXCLUSIVE      Exam ‘2' must be timetabled in a room by itself.

    // Institutional Model Weightings
    // This section begins with the line [InstitutionalWeightings] and provides information on values given to 'global' soft constraints.
    // TWOINAROW, 7
    // TWOINADAY, 5
    // PERIODSPREAD, 3
    // NONMIXEDDURATIONS, 10
    // FRONTLOAD, 100, 30, 5
    ////

    vector<BinaryHardConstraint> periodConstraints;

    if (it != tok.end()) {
        readPeriodHardConstraints(it, periodConstraints);


        // TODO

        // Room Related Hard Constraints

        // Institutional Model Weightings
    }
}


void ITC07TestSet::readPeriodHardConstraints(tokenizer<escaped_list_separator<char> >::iterator& it,
                                             vector<BinaryHardConstraint>& periodConstraints) {
    std::string s(*it);
    boost::cmatch matches;
    // Section begins with the tag [PeriodHardConstraints]
    boost::regex re("(\\\n)(\\[PeriodHardConstraints\\])");
    if (!boost::regex_match(s.c_str(), matches, re))
    {
        cout << "The regexp \"" << re << "\" does not match \"" << s << "\"" << endl;
        cin.get();

        // throw exception
        throw runtime_error("No PeriodHardConstraints specified");
    }

    cout << "Read period hard constraints header" << endl;
    ++it; // Proceed to read the next line

    // Read constraints if any
    while (readPeriodConstraint(it, periodConstraints))
        ++it;
}

bool ITC07TestSet::readPeriodConstraint(tokenizer<escaped_list_separator<char> >::iterator& it,
                                             vector<BinaryHardConstraint>& periodConstraints) {

    std::string s(*it);

    boost::cmatch matches;

    boost::regex re("(\\\n)(\\d+), (AFTER|EXAM_COINCIDENCE|EXCLUSION), (\\d+)");

    if (boost::regex_match(s.c_str(), matches, re)) {
        // Extract values from string.
        string matchExam1(matches[2].first, matches[2].second);
        int exam1 = atoi(matchExam1.c_str());

        string constraintType(matches[3].first, matches[3].second);

        string matchExam2(matches[4].first, matches[4].second);
        int exam2 = atoi(matchExam2.c_str());

        cout << "\texam1: " << exam1 << endl;
        cout << "\ttype of constraint: " << constraintType << endl;
        cout << "\texam2: " << exam2 << endl;
        cout << endl;

        // TODO
        // Create constraint object and insert in the vector...

    }
    else {
        cout << "The regexp \"" << re << "\" does not match \"" << s << "\"" << endl;
        cin.get();
    }

}


boost::shared_ptr<AdjacencyList> ITC07TestSet::buildExamGraph(boost::shared_ptr<Matrix> const& ptrConflictMatrix) {
    // Instantiate graph with ncols vertices
    boost::shared_ptr<AdjacencyList> ptrGraphAux(new AdjacencyList(ptrConflictMatrix->nCols()+1));
    int cost;
    for (int v1 = 1; v1 <= ptrConflictMatrix->nLines(); ++v1) {
        for (int v2 = 1; v2 <= ptrConflictMatrix->nCols(); ++v2) {
            cost = ptrConflictMatrix->getVal(v1, v2);
            if (cost != 0) {
                if (v2 > v1) {
                    add_edge(v1, v2, *ptrGraphAux.get());
                }
            }
        }
    }
    return ptrGraphAux;
}





















