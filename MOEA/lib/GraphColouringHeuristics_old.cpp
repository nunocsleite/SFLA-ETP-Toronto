
#include "GraphColouringHeuristics.h"


typedef int exam_t;


//#define _DEBUG



void GCHeuristics::SD(Matrix const& conflictMatrix, AdjacencyList const& graph, Chromosome& timetable,
                      vector<int>& fixedExams, unordered_map<int, int> &examsPeriods, vector<int>& remainderExams) {
    //===
    // SD (Saturation Degree) graph colouring heuristic
    //
    // 1. Initially, create a Priority queue with all exams inserted randomly.
    //    Set available periods of all unscheduled exams equal to T, where T is the max. number of periods.
    // 2. Pop exam 'ei' from the top of exam's priority queue and schedule it into a random available period 'tj'.
    //    2.1. Mark period 'tj' of all the exams connected to 'ei' as unvailable, because they have common registered students.
    //         The number of available periods is the priority used to sort the priority queue of unscheduled exams.
    // 3. Goto step 2.
    //===

    // random generator
    boost::random::mt19937 rng;
    // Timetable's periods
    vector<unordered_set<exam_t> >& periods = timetable.getPeriods();

    // 1. Initially, create a Priority queue with all exams inserted randomly.
    //    Set available periods of all unscheduled exams equal to T, where T is the max. number of periods.
    //
    // N is the number of vertices
    const int N = conflictMatrix.nCols();
    VertexPriorityQueue pq(N);

    // Max number of periods
    const int numPeriods = timetable.getNumPeriods();
#ifdef _DEBUG
    cout << "N = " << N << endl;
    cout << "numPeriods = " << numPeriods << endl;
#endif

    // Exam's vector
    vector<exam_t> exams;
    for (exam_t ei = 1; ei <= N; ++ei)
        exams.push_back(ei);
    // Removed exams
    vector<bool> removed_exams(N+1); // Initially all false

    // available period's list
    vector<vector<int> > examsAvailablePeriodsList(N+1);
    vector<int> allPeriods(numPeriods);
    for (int p = 0; p < numPeriods; ++p)
        allPeriods[p] = p+1;

    //#ifdef _DEBUG
    //    copy(allPeriods.begin(), allPeriods.end(), ostream_iterator<int>(cout, " "));
    //#endif

    for (int exam_i = 1; exam_i <= N; ++exam_i)
        examsAvailablePeriodsList[exam_i] = allPeriods;


    //
    // Insert fixed exams
    //
    GCHeuristics::insertFixedExams(periods, exams, removed_exams, numPeriods, pq, examsAvailablePeriodsList,
                                   allPeriods, fixedExams, examsPeriods, remainderExams, N, graph);

    // shuffle exams
//    random_shuffle(exams.begin(), exams.end());
    //copy(exams.begin(), exams.end(), ostream_iterator<int>(cout, " "));
    // Insert them in the priority queue
//    for (int i = 0; i < N; ++i) {
//        pq.push(exams[i], numPeriods);
//    }
//    // available period's list
//    vector<vector<int> > examsAvailablePeriodsList(N+1);
//    vector<int> allPeriods(numPeriods);
//    for (int p = 0; p < numPeriods; ++p)
//        allPeriods[p] = p+1;

//#ifdef _DEBUG
//    copy(allPeriods.begin(), allPeriods.end(), ostream_iterator<int>(cout, " "));
//#endif

//    for (int exam_i = 1; exam_i <= N; ++exam_i)
//        examsAvailablePeriodsList[exam_i] = allPeriods;

    while (!pq.empty()) {
        // 2. Pop exam 'ei' from the top of exam's priority queue and schedule it into a random period 'tj'.
        VertexPriorityQueue::heap_data data = pq.top();
        exam_t ei = data.vertex;
        int availablePeriods = data.priority;
        pq.pop();

        cout << "ei = " << ei << ", availablePeriods = " << availablePeriods << endl;

#ifdef _DEBUG
        cout << "ei = " << ei << ", availablePeriods = " << availablePeriods << endl;
#endif
        int tj;
        if (availablePeriods > 0) {
            // schedule it into a random available period 'tj'.
            boost::random::uniform_int_distribution<> periodsDist(1, availablePeriods); // period distribution
            int index = periodsDist(rng);
            // get available period's list
            vector<int> availablePeriodsList = examsAvailablePeriodsList[ei];
            tj = availablePeriodsList[index-1];


            cout << "tj = " << tj << endl;


#ifdef _DEBUG
            cout << "index-1 = " << index-1 << ", tj = " << tj << endl;
#endif
            periods[tj-1].insert(ei);
        }
        else {
            cout << "Add a new period to the end of timetable" << endl;
            cout << "SHOULD NOT HAPPEN..." << endl;

            cin.get();

            // Add a new period to the end of timetable
            unordered_set<int> newPeriod;
            newPeriod.insert(ei);
            periods.push_back(newPeriod);

            tj = numPeriods;

#ifdef _DEBUG
            cout << "add new period" << endl;
#endif
        }
        //    2.1. Mark period 'tj' of all the exams connected to 'ei' as unvailable, because they have common registered students.
        //         The number of available periods is the priority used to sort the priority queue of unscheduled exams.

        // get adjacent vertices
        property_map<AdjacencyList, vertex_index_t>::type index_map = get(vertex_index, graph);

        graph_traits<AdjacencyList>::adjacency_iterator ai, a_end;
        boost::tie(ai, a_end) = adjacent_vertices(ei, graph);
        for (; ai != a_end; ++ai) {
            // get adjacent exam
            exam_t ej = get(index_map, *ai);

            cout << "ej = " << ej << endl;

#ifdef _DEBUG
            //cout << "ej = " << ej << ", tj = " << tj << endl;
            cout << "ej = " << ej << endl;

#endif
            // If exam ej was not scheduled yet, mark adjacent exam's tj period as unvailable
            if (!removed_exams[ej])  {
                bool removed = false;

                // Remove period tj from ej's exams periods list
                // OPTIMIZATION: USE a bool array
                vector<int>::iterator it;
                for (it = examsAvailablePeriodsList[ej].begin(); it != examsAvailablePeriodsList[ej].end(); ++it) {
                    if (*it == tj) {
#ifdef _DEBUG
                        cout << "size before removing = " << examsAvailablePeriodsList[ej].size() << endl;
#endif

                        examsAvailablePeriodsList[ej].erase(it);

#ifdef _DEBUG
                        cout << "size after removing = " << examsAvailablePeriodsList[ej].size() << endl;
#endif

                        cout << "Before update" << endl;

                        pq.update(ej, examsAvailablePeriodsList[ej].size());
                        removed = true;
                        break;
                    }
                }
            }
            removed_exams[ei] = true;
        }
        // 3. Goto step 2.
    }

    // Compute timetable conflicts
    timetable.computeProximityCosts();

    //cin.get();
}



void GCHeuristics::insertFixedExams(vector<unordered_set<exam_t> >& periods, vector<exam_t>& exams,
                                    vector<bool> removed_exams, int numPeriods,
                                    VertexPriorityQueue& pq, vector<vector<int> >& examsAvailablePeriodsList,
                                    vector<int>& allPeriods, vector<int>& fixedExams, unordered_map<int, int>& examsPeriods,
                                    vector<int>& remainderExams, int N, AdjacencyList const& graph) {


    // random generator
    boost::random::mt19937 rng;

    // shuffle exams
    //random_shuffle(exams.begin(), exams.end());
    //copy(exams.begin(), exams.end(), ostream_iterator<int>(cout, " "));

    // Insert them in the priority queue
//    for (int i = 0; i < N; ++i) {
//        pq.push(exams[i], numPeriods);
//    }
    // Adaptation:
    // Push only fixed exams first and then push the remainder
    for (int i = 0; i < fixedExams.size(); ++i) {
        pq.push(fixedExams[i], -1); // Set high priority
    }
    // Push remainder
    for (int i = 0; i < remainderExams.size(); ++i) {
        pq.push(remainderExams[i], numPeriods); // Set low priority
    }


//    // available period's list
//    vector<vector<int> > examsAvailablePeriodsList(N+1);
//    vector<int> allPeriods(numPeriods);
//    for (int p = 0; p < numPeriods; ++p)
//        allPeriods[p] = p+1;

//#ifdef _DEBUG
//    copy(allPeriods.begin(), allPeriods.end(), ostream_iterator<int>(cout, " "));
//#endif

//    for (int exam_i = 1; exam_i <= N; ++exam_i)
//        examsAvailablePeriodsList[exam_i] = allPeriods;



//    while (!pq.empty()) {

    // Schedule only the fixed exams
    int z = 0;
    while (z < fixedExams.size()) {
        // 2. Pop exam 'ei' from the top of exam's priority queue and schedule it into a random period 'tj'.
        VertexPriorityQueue::heap_data data = pq.top();
        exam_t ei = data.vertex;
        int availablePeriods = data.priority;
        pq.pop();

        cout << "availablePeriods = " << availablePeriods << endl;

        cout << "ei = " << ei;

        // Insert exam into fixed position
        if (examsPeriods.find(ei) == examsPeriods.end()) {
            throw runtime_error("Could not find key in examsPeriods");

        }

        int tj = examsPeriods.at(ei); // tj must be >= 1

        // Insert exam into the timetable
        //periods[tj-1].insert(ei);
        //
        //   The exam is already inserted!!!
        //

        cout << ", tj = " << tj << endl;


//#ifdef _DEBUG
//        cout << "ei = " << ei << ", availablePeriods = " << availablePeriods << endl;
//#endif
//        int tj;
//        if (availablePeriods > 0) {
//            // schedule it into a random available period 'tj'.
//            boost::random::uniform_int_distribution<> periodsDist(1, availablePeriods); // period distribution
//            int index = periodsDist(rng);
//            // get available period's list
//            vector<int> availablePeriodsList = examsAvailablePeriodsList[ei];
//            tj = availablePeriodsList[index-1];

//#ifdef _DEBUG
//            cout << "index-1 = " << index-1 << ", tj = " << tj << endl;
//#endif
//            periods[tj-1].insert(ei);
//        }
//        else {
//            // Add a new period to the end of timetable
//            unordered_set<int> newPeriod;
//            newPeriod.insert(ei);
//            periods.push_back(newPeriod);

//            tj = numPeriods;

//#ifdef _DEBUG
//            cout << "add new period" << endl;
//#endif
//        }


        //    2.1. Mark period 'tj' of all the exams connected to 'ei' as unvailable, because they have common registered students.
        //         The number of available periods is the priority used to sort the priority queue of unscheduled exams.

        // get adjacent vertices
        property_map<AdjacencyList, vertex_index_t>::type index_map = get(vertex_index, graph);

        graph_traits<AdjacencyList>::adjacency_iterator ai, a_end;
        boost::tie(ai, a_end) = adjacent_vertices(ei, graph);
        for (; ai != a_end; ++ai) {
            // get adjacent exam
            exam_t ej = get(index_map, *ai);


#ifdef _DEBUG
            //cout << "ej = " << ej << ", tj = " << tj << endl;
            cout << "ej = " << ej << endl;

#endif
            // If exam ej was not scheduled yet, mark adjacent exam's tj period as unvailable
            if (!removed_exams[ej])  {
                bool removed = false;

                // Remove period tj from ej's exams periods list
                // OPTIMIZATION: USE a bool array
                vector<int>::iterator it;
                for (it = examsAvailablePeriodsList[ej].begin(); it != examsAvailablePeriodsList[ej].end(); ++it) {
                    if (*it == tj) {
#ifdef _DEBUG
                        cout << "size before removing = " << examsAvailablePeriodsList[ej].size() << endl;
#endif

                        examsAvailablePeriodsList[ej].erase(it);

#ifdef _DEBUG
                        cout << "size after removing = " << examsAvailablePeriodsList[ej].size() << endl;
#endif
                        // If not a fixed exam, update it
                        if (examsPeriods.find(ej) == examsPeriods.end())
                            pq.update(ej, examsAvailablePeriodsList[ej].size());

                        removed = true;
                        break;
                    }
                }
            }
            removed_exams[ei] = true;
        }

        // 3. Goto step 2.

        ++z;
    }

    cout << "Finished inserting fixed exams " << endl;

//    cin.get();

}


/////////////////////////////////////////////////////////////////////////////////////////

void GCHeuristics::SD(Matrix const& conflictMatrix, AdjacencyList const& graph, Chromosome& timetable) {
	//===
	// SD (Saturation Degree) graph colouring heuristic
	//
	// 1. Initially, create a Priority queue with all exams inserted randomly.
	//    Set available periods of all unscheduled exams equal to T, where T is the max. number of periods.
	// 2. Pop exam 'ei' from the top of exam's priority queue and schedule it into a random available period 'tj'.
	//    2.1. Mark period 'tj' of all the exams connected to 'ei' as unvailable, because they have common registered students.
	//         The number of available periods is the priority used to sort the priority queue of unscheduled exams.
	// 3. Goto step 2.
	//===

	// random generator
	boost::random::mt19937 rng;  
	// Timetable's periods
//	vector<vector<exam_t> >& periods = timetable.getPeriods();
    vector<unordered_set<exam_t> >& periods = timetable.getPeriods();

	// 1. Initially, create a Priority queue with all exams inserted randomly.
	//    Set available periods of all unscheduled exams equal to T, where T is the max. number of periods.
	//
	// N is the number of vertices
	const int N = conflictMatrix.nCols(); 
	VertexPriorityQueue pq(N);

	// Max number of periods
	const int numPeriods = timetable.getNumPeriods();
#ifdef _DEBUG
	cout << "N = " << N << endl;
	cout << "numPeriods = " << numPeriods << endl;
#endif

	// Exam's vector
	vector<exam_t> exams;
	for (exam_t ei = 1; ei <= N; ++ei) 
		exams.push_back(ei); 
	// Removed exams
	vector<bool> removed_exams(N+1); // Initially all false


	// shuffle exams
	random_shuffle(exams.begin(), exams.end());
	//copy(exams.begin(), exams.end(), ostream_iterator<int>(cout, " "));
	// Insert them in the priority queue
	for (int i = 0; i < N; ++i) {
        pq.push(exams[i], numPeriods);
	}
	// available period's list
	vector<vector<int> > examsAvailablePeriodsList(N+1);
	vector<int> allPeriods(numPeriods);
	for (int p = 0; p < numPeriods; ++p)
		allPeriods[p] = p+1;

#ifdef _DEBUG
	copy(allPeriods.begin(), allPeriods.end(), ostream_iterator<int>(cout, " "));
#endif


	for (int exam_i = 1; exam_i <= N; ++exam_i)
		examsAvailablePeriodsList[exam_i] = allPeriods;

	while (!pq.empty()) {
		// 2. Pop exam 'ei' from the top of exam's priority queue and schedule it into a random period 'tj'.
		VertexPriorityQueue::heap_data data = pq.top();
		exam_t ei = data.vertex;
		int availablePeriods = data.priority;
		pq.pop();
		
#ifdef _DEBUG
		cout << "ei = " << ei << ", availablePeriods = " << availablePeriods << endl;
#endif
		int tj;
		if (availablePeriods > 0) {
			// schedule it into a random available period 'tj'.
			boost::random::uniform_int_distribution<> periodsDist(1, availablePeriods); // period distribution 
			int index = periodsDist(rng);
			// get available period's list
			vector<int> availablePeriodsList = examsAvailablePeriodsList[ei];
			tj = availablePeriodsList[index-1];

#ifdef _DEBUG
			cout << "index-1 = " << index-1 << ", tj = " << tj << endl;
#endif
//			periods[tj-1].push_back(ei);
            periods[tj-1].insert(ei);
		}
		else {
			// Add a new period to the end of timetable
//			vector<exam_t> newPeriod;
//			newPeriod.push_back(ei);
            unordered_set<int> newPeriod;
            newPeriod.insert(ei);
            periods.push_back(newPeriod);

			tj = numPeriods;

//			timetable.setNumPeriods(numPeriods+1);

#ifdef _DEBUG
			cout << "add new period" << endl;
#endif
		}
		//    2.1. Mark period 'tj' of all the exams connected to 'ei' as unvailable, because they have common registered students.
		//         The number of available periods is the priority used to sort the priority queue of unscheduled exams.
		
		// get adjacent vertices
		property_map<AdjacencyList, vertex_index_t>::type index_map = get(vertex_index, graph);
//
//	for (boost::tie(i, end) = vertices(graph); i != end; ++i) {
//		std::cout << get(index_map, *i);
		graph_traits<AdjacencyList>::adjacency_iterator ai, a_end;
		boost::tie(ai, a_end) = adjacent_vertices(ei, graph);
		for (; ai != a_end; ++ai) {
			// get adjacent exam
			exam_t ej = get(index_map, *ai);
			
		
#ifdef _DEBUG
			//cout << "ej = " << ej << ", tj = " << tj << endl;
			cout << "ej = " << ej << endl;

#endif
			// If exam ej was not scheduled yet, mark adjacent exam's tj period as unvailable
			if (!removed_exams[ej])  {
				bool removed = false;
			
				// Remove period tj from ej's exams periods list
				// OPTIMIZATION: USE a bool array
				vector<int>::iterator it;
				for (it = examsAvailablePeriodsList[ej].begin(); it != examsAvailablePeriodsList[ej].end(); ++it) {
					if (*it == tj) {
#ifdef _DEBUG
						cout << "size before removing = " << examsAvailablePeriodsList[ej].size() << endl;
#endif

						examsAvailablePeriodsList[ej].erase(it);

#ifdef _DEBUG
						cout << "size after removing = " << examsAvailablePeriodsList[ej].size() << endl;
#endif

						pq.update(ej, examsAvailablePeriodsList[ej].size());
						removed = true;
						break;
					}
				}

				//if (!removed) {
				//	cout << "not removed" << endl;
				//	//cin.get();
				//}

				//copy(examsAvailablePeriodsList[ej].begin(), examsAvailablePeriodsList[ej].end(), ostream_iterator<int>(cout, " "));
			

				//examsAvailablePeriodsList[ej].erase(examsAvailablePeriodsList[ej].begin()+(tj-1));

				// Decrement available periods (the priority) and update priority queue

				//cout << "examsAvailablePeriodsList[ej].size() = " << examsAvailablePeriodsList[ej].size() << endl;

				//cout << "updating..." << endl;

				/*if (removed) {
					pq.update(ej, examsAvailablePeriodsList[ej].size());
				}*/
			}

			removed_exams[ei] = true;

			//cin.get();

		}
		//cin.get();
		/*if (ai == a_end)
			std::cout << " empty";
		else
			std::cout << " links: ";
		for (; ai != a_end; ++ai) {
			std::cout << get(index_map, *ai);
			if (boost::next(ai) != a_end)
			std::cout << ", ";
		}

		cin.get();*/

		// 3. Goto step 2.
	}

//    cout << "Good" << endl;

    // Compute timetable conflicts
    timetable.computeProximityCosts();

	//cin.get();
}


//	graph_traits<AdjacencyList>::vertex_iterator i, end;
//	graph_traits<AdjacencyList>::adjacency_iterator ai, a_end;
//	property_map<AdjacencyList, vertex_index_t>::type index_map = get(vertex_index, graph);
//
//	for (boost::tie(i, end) = vertices(graph); i != end; ++i) {
//		std::cout << get(index_map, *i);
//		boost::tie(ai, a_end) = adjacent_vertices(*i, graph);
//		if (ai == a_end)
//			std::cout << " empty";
//		else
//			std::cout << " links: ";
//		for (; ai != a_end; ++ai) {
//			std::cout << get(index_map, *ai);
//			if (boost::next(ai) != a_end)
//			std::cout << ", ";
//		}
