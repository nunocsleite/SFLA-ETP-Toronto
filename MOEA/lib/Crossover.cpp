
#include "Crossover.h"



string Crossover::className() const
{
    return "ETTP Crossover";
}


bool Crossover::operator()(moeoChromosome& _chromosome1, moeoChromosome& _chromosome2)
{
//    cout << "Crossover" << endl;

    bool oneAtLeastIsModified = true;
    //
    // Day-exchange crossover
    //
    // In day-exchange crossover, only the best days (excluding
    // Saturdays, since exams scheduled on Saturdays are always
    // clash-free) of chromosomes, selected based on the crossover
    // rate, are eligible for exchange. The best day consists of three
    // periods and is the day with the lowest number of clashes per student.

    // Computation of the offspring
    generateOffspring(_chromosome1, _chromosome2);
    generateOffspring(_chromosome2, _chromosome1);


//    // does at least one genotype has been modified ?
//    if ((_chromosome1 != offspring1) || (_chromosome2 != offspring2))
//    {
//        // update
//        _chromosome1.value(offspring1);
//        _chromosome2.value(offspring2);
//        // at least one genotype has been modified
//        oneAtLeastIsModified = true;
//    }
//    else
//    {
//        // no genotype has been modified
//        oneAtLeastIsModified = false;
//    }
    // return 'true' if at least one genotype has been modified
    return oneAtLeastIsModified;
}


void Crossover::generateOffspring(moeoChromosome &_parent1, moeoChromosome &_parent2)
{
    // We change directly _parent1 to be the new offspring because it's already a copy
    // of the original parent.
    moeoChromosome& offspring = _parent1;
    // Compute the "Best day" of the second parent.
    int day2 = bestDay(_parent2);
    // Insert best day in the first parent and remove duplicated exams
    insertDay(offspring, _parent2, day2);
    // Actualize proximity costs
    offspring.computeProximityCosts();

//    cout << "offspring" << endl;

//    cout << offspring << endl;


//    cin.get();

//  std::vector<bool> taken_values(result.size(), false);
//  if (_point1 > _point2)
//    std::swap(_point1, _point2);
//  /* first parent */
//  for (unsigned int i=0 ; i<=_point1 ; i++)
//    {
//      // result[i] == _parent1[i]
//      taken_values[_parent1[i]] = true;
//    }
//  for (unsigned int i=_point2 ; i<result.size() ; i++)
//    {
//      // result[i] == _parent1[i]
//      taken_values[_parent1[i]] = true;
//    }
//  /* second parent */
//  unsigned int i = _point1+1;
//  unsigned int j = 0;
//  while (i<_point2 && j<_parent2.size())
//    {
//      if (! taken_values[_parent2[j]])
//        {
//          result[i] = _parent2[j];
//          i++;
//        }
//      j++;
//    }
//  return result;
}

// Compute the "Best day" to exchange between chromosomes.
// The best day consist of three periods (we exclude saturdays because
// saturdays are always clash-free) and is the day with the lowest number
// of clashes per student.
int Crossover::bestDay(const moeoChromosome &chromosome) {
    int numPeriods = chromosome.getNumPeriods();
    int indexBestDay = 0, numClashesDay = 0, numStudentsDay = 0;
    double numClashesPerStudent = 0.0, bestNumClashes = 0.0;

//    for (int i = 0; i < numPeriods; ++i) {

    // Except last day which has no conflicts
    for (int i = 0; i < numPeriods-1; ++i) {


        // Get number of clashes of the current "day" (or period)
        numClashesDay = chromosome.getNumClashesPeriod(i);
        // Get number of students of the current "day" (or period)
        numStudentsDay = chromosome.getNumStudentsPeriod(i);
        // Compute number of clashes per student
        numClashesPerStudent = (double)numClashesDay / numStudentsDay;
//        numClashesPerStudent = numStudentsDay;

//        cout << "numClashesDay = " << numClashesDay << endl;
//        cout << "numStudentsDay = " << numStudentsDay << endl;
//        cout << "numClashesPerStudent = " << numClashesPerStudent << endl;
//        cin.get();

        // Actualize best day
        if (i == 0 || numClashesPerStudent < bestNumClashes) {
            bestNumClashes = numClashesPerStudent;
            indexBestDay = i;
        }
    }
    return indexBestDay;
}

class IsInExamsToRemove {
    vector<int> const& examsToRemove;
public:
    IsInExamsToRemove(vector<int> const& examsToRemove) : examsToRemove(examsToRemove) { }
    bool operator()(int exam) {
        return find_if(examsToRemove.begin(), examsToRemove.end(), bind2nd(equal_to<int>(), exam)) != examsToRemove.end();
    }
};


void Crossover::insertDay(moeoChromosome &_parent1, moeoChromosome &_parent2, int day2) {

//    cout << endl << endl << "Before inserting day" << endl;
    _parent1.validate();

    // The newly inserted day takes the place of a randomly chosen day
    // which is pushed to the end of the timetable.
//    int numPeriods = _parent1.getNumPeriods();
    // Get periods
//    vector<vector<int> >& periods1 = _parent1.getPeriods();
//    vector<vector<int> >& periods2 = _parent2.getPeriods();
    vector<unordered_set<int> >& periods1 = _parent1.getPeriods();
    vector<unordered_set<int> >& periods2 = _parent2.getPeriods();

    // Generate random period
//    int randPeriod = rng.random(numPeriods);
    int randPeriod = rng.random(periods1.size());


//    cout << "randPeriod = " << randPeriod << endl;
//    cout << "periods1 exams before inserting day from periods2" << endl;
//    int i = 1;
//    for (vector<vector<int> >::iterator it = periods1.begin(); it != periods1.end(); ++it) {
//        cout << endl << "Period " << i << endl;
//        for (vector<int>::const_iterator examit = (*it).begin(); examit != (*it).end(); ++examit) {
//            cout << *examit << " ";
//        }
//        cout << endl;
//        ++i;
//    }

    // Push existing day to the end of the timetable
    periods1.push_back(periods1[randPeriod]);
    // Insert new day from second parent
    periods1[randPeriod] = periods2[day2];
    // Increment number of periods
//    ++numPeriods;
//    _parent1.setNumPeriods(numPeriods);

//    cout << "periods1 exams after inserting day from periods2" << endl;
//    i = 1;
//    for (vector<vector<int> >::iterator it = periods1.begin(); it != periods1.end(); ++it) {
//        cout << endl << "Period " << i << endl;
//        for (vector<int>::const_iterator examit = (*it).begin(); examit != (*it).end(); ++examit) {
//            cout << *examit << " ";
//        }
//        cout << endl;
//        ++i;
//    }



    // To ensure the feasibility of chromosomes after the
    // crossover, duplicated exams are deleted. These exams are
    // removed from the original periods, while the newly inserted
    // periods are left intact.
//    vector<int> const& examsToRemove = periods2[day2];
    unordered_set<int> const& examsToRemove = periods2[day2];

//    cout << "periods2 exams to remove:" << endl;
//    for (vector<int>::const_iterator examit = examsToRemove.begin(); examit != examsToRemove.end(); ++examit) {
//        cout << *examit << " ";
//    }

//    IsInExamsToRemove pred(examsToRemove);

    int numRemovedExams = 0;

    /// TODO: USAR unordored_map<INT> PARA GUARDAR EXAMES

    // Remove each exam in examsToRemove vector
//    for (vector<int>::const_iterator itRem = examsToRemove.begin(); itRem != examsToRemove.end(); ++itRem) {
    for (unordered_set<int>::const_iterator itRem = examsToRemove.begin(); itRem != examsToRemove.end(); ++itRem) {

        int i;
//        for (i = 0; i < numPeriods; ++i) {
        for (i = 0; i < periods1.size(); ++i) {
            if (i != randPeriod) {
//                vector<int>& exams = periods1[i];
                unordered_set<int>& exams = periods1[i];
    //            vector<int>::iterator newEnd = remove_if(exams.begin(), exams.end(), pred);
    //            numRemovedExams += (exams.end()-newEnd);
    //            exams.erase(newEnd, exams.end());


                if (exams.find(*itRem) != exams.end()) {
                    exams.erase(*itRem);
                    break;
                }



//                vector<int>::iterator newEnd = remove_if(exams.begin(), exams.end(), bind2nd(equal_to<int>(), *itRem));

                // ESTA LINHA NAO FUNCIONA
//                unordered_set<int>::iterator newEnd = remove_if(exams.begin(), exams.end(),
//                                                                     bind2nd(equal_to<int>(), *itRem));

//                if (newEnd != exams.end()) {
//                    // Remove exam
//                    exams.erase(newEnd, exams.end());
//                    ++numRemovedExams;
//                    break;
//                }
//                int numExamsToRemove = (exams.end()-newEnd);
//                cout << endl << "numExamsToRemove = " << numExamsToRemove << endl;
//                cout << endl << "exams before remove" << endl;
//                copy(exams.begin(), exams.end(), ostream_iterator<int>(cout, " "));
//                cout << endl;
//                numRemovedExams += numExamsToRemove;
//                exams.erase(newEnd, exams.end());
//                cout << endl << "exams after remove" << endl;
//                copy(exams.begin(), exams.end(), ostream_iterator<int>(cout, " "));
//                cout << endl;
            }
        }
        if (i == periods1.size()) {
            cout << endl << "numPeriods = " << periods1.size() << endl;
            cout << endl << "Couldn't find exam " << *itRem << endl;
            cin.get();
        }
    }

//    cout << endl << endl << "numRemovedExams = " << numRemovedExams << endl;
//    cout << endl << endl << "examsToRemove.size() = " << examsToRemove.size() << endl;

//    cout << endl << endl << "After inserting day and remove duplicated exams" << endl;

//    cout << "periods1 exams" << endl;
//    i = 1;
//    for (vector<vector<int> >::iterator it = periods1.begin(); it != periods1.end(); ++it) {
//        cout << endl << "Period " << i << endl;
//        for (vector<int>::const_iterator examit = (*it).begin(); examit != (*it).end(); ++examit) {
//            cout << *examit << " ";
//        }
//        cout << endl;
//        ++i;
//    }

    _parent1.validate();
}

/////////////////////////////////////////////////////////////////////////////////////////////////





/////////////////////////////////////////////////////////////////////////////////////////////////

//% Day-exchange crossover
//function [timetable1, timetable2] = dayExchangeCrossover(Data, timetable1, timetable2)
//    % In day-exchange crossover, only the best days (excluding
//    % Saturdays, since exams scheduled on Saturdays are always
//    % clash-free) of chromosomes, selected based on the crossover
//    % rate, are eligible for exchange. The best day consists of three
//    % periods and is the day with the lowest number of clashes per student.
//    %
//    % Compute the "Best day" to exchange between chromosomes.
//    day1 = bestDay(Data, timetable1);
//    day2 = bestDay(Data, timetable2);

//    if (~verifyTimetable(Data, timetable1))
//        fprintf('[Before exchangeDays] Timetable 1 is not valid\n');
//        printTimetable(Data, timetable1);
//        pause
//    end
//    if (~verifyTimetable(Data, timetable2))
//        fprintf('[Before exchangeDays] Timetable 2 is not valid\n');
//        printTimetable(data, timetable2);
//        pause
//    end

//    % Exchange best days between the two chromossomes.
//    [timetable1, timetable2] = exchangeDays(Data, timetable1, day1, timetable2, day2);

//    if (~verifyTimetable(Data, timetable1))
//        fprintf('[After exchangeDays] Timetable 1 is not valid\n');
//        printTimetable(data, timetable1 );
//        pause
//    end
//    if (~verifyTimetable(Data, timetable2))
//        fprintf('[After exchangeDays] Timetable 2 is not valid\n');
//        printTimetable(data, timetable2);
//        pause
//    end

//    % Actualize number of clashes
//    timetable1 = computeNumClashes(Data, timetable1);
//    timetable2 = computeNumClashes(Data, timetable2);
//end


//% Compute the "Best day" to exchange between chromosomes.
//% The best day consist of three periods (we exclude saturdays because
//% saturdays are always clash-free) and is the day with the lowest number
//% of clashes per student.
//function bday = bestDay(Data, timetable)
//    NumPeriods = timetable.NumPeriods;
//    clashesPerStudent = zeros(1, NumPeriods);
//    % For each period do:
//    %
//    %   - Get number of clashes in each period divided by
//    %     the number of students in that period.
//    numClashes0 = 0;

//    for p = 1 : NumPeriods
//        numStudents = getNumberStudentsPeriod(Data, timetable, p);
//        if (p > 1)
//            numClashes0 = getNumberClashesPeriod(Data, timetable, p-1);
//        end
//        if (p == NumPeriods)
//            numClashes = numClashes0;
//        else
//            numClashes = numClashes0 + getNumberClashesPeriod(Data, timetable, p);
//        end
//%         examList = timetable.Periods{p};
//%         Data.Classes(examList)
//%         clashesPerStudent(p) = numClashes/numStudents;
//        clashesPerStudent(p) = numClashes;

//        % VER MELHOR MEDIDA

//    end

//    [V, I] = min(clashesPerStudent); % Considering just one period per day
//    bday = I(1);
//end


//% Exchange best days between the two chromossomes.
//function [timetable1, timetable2] = exchangeDays(Data, timetable1, day1, timetable2, day2)
//    % Exchange days between timetables.
//    % The best day and the day after are exchanged -> SEE THIS
//    %
//    % The newly inserted day takes the place of a randomly chosen day
//    % which is pushed to the end of the timetable.
//    periodsToInsertTimetable1Copy = cell(1, 1);
//    periodsToInsertTimetable1Copy{1} = timetable1.Periods{day1};
//    % Insert day2 in timetable1
//    idx1 = randi([1, int32(timetable1.NumPeriods)]);
//    fromPeriodIdx = idx1;
//    periodsToInsert = cell(1, 1);
//    periodsToInsert{1} = timetable2.Periods{day2};
//    timetable1 = insertDay(Data, timetable1, fromPeriodIdx, periodsToInsert);
//    % To ensure the feasibility of chromosomes after the
//    % crossover, duplicated exams are deleted. These exams are
//    % removed from the original periods, while the newly inserted
//    % periods are left intact.
//    timetable1 = removeDuplicatedExams(Data, timetable1, fromPeriodIdx, periodsToInsert);

//    numExams = 0;
//    for j = 1 : timetable1.NumPeriods
//        examList = timetable1.Periods{j};
//        numExams = numExams + length(examList);
//    end
//    if (numExams ~= 80)
//        pause
//    end
//    % Insert day1 in timetable2
//    idx2 = randi([1, int32(timetable2.NumPeriods)]);
//    fromPeriodIdx = idx2;
//    periodsToInsert{1} = periodsToInsertTimetable1Copy{1};
//    timetable2 = insertDay(Data, timetable2, fromPeriodIdx, periodsToInsert);
//    % To ensure the feasibility of chromosomes after the
//    % crossover, duplicated exams are deleted. These exams are
//    % removed from the original periods, while the newly inserted
//    % periods are left intact.
//    timetable2 = removeDuplicatedExams(Data, timetable2, fromPeriodIdx, periodsToInsert);
//    numExams = 0;
//    for j = 1 : timetable2.NumPeriods
//        examList = timetable2.Periods{j};
//        numExams = numExams + length(examList);
//    end

//    if (numExams ~= 80)
//        pause
//    end
//end

//function showClassNames(Data, periodExamList)
//    Data.Classes(periodExamList)
//end

//function timetable = insertDay(Data, timetable, fromPeriodIdx, periodsToInsert)
//    periodsToRemove = cell(1, 1);
//    periodsToRemove{1} = timetable.Periods{fromPeriodIdx};
//    timetable.Periods{fromPeriodIdx} = periodsToInsert{1};
//    % Add two periods to the timetable
//    timetable = addPeriodToTimetable(timetable);
//    timetable.Periods{timetable.NumPeriods} = periodsToRemove{1};
//end

//% To ensure the feasibility of chromosomes after the
//% crossover, duplicated exams are deleted. These exams are
//% removed from the original periods, while the newly inserted
//% periods are left intact.
//function timetable = removeDuplicatedExams(Data, timetable, fromPeriodIdx, periodsToInsert)
//    % Form an exam list with the inserted periods
//    insertedExams = [periodsToInsert{1}];
//    numRepeatedExams = 0;

//    for i = 1 : timetable.NumPeriods
//        % Except the periods just inserted
//         if (i ~= fromPeriodIdx)
//            % Remove duplicated exams
//            examList = timetable.Periods{i};
//            numExams = length(examList);
//            % Find exams from this list which are repeated and remove them
//            exams = [];
//            for j = 1 : numExams
//                I = find(insertedExams == examList(j));
//                if (~isempty(I))
//                    exams = [exams examList(j)];
//                    numRepeatedExams = numRepeatedExams+1;
//                end
//            end
//            % Remove exams
//            if (~isempty(exams))
//                for j = 1 : length(exams)
//                    examList = remove(exams(j), examList);
//                end
//                timetable.Periods{i} = examList;
//            end
//        end
//    end
//    if (numRepeatedExams ~= length(insertedExams))
//        disp('Bug...')
//        pause
//    end
//end

/////////////////////////////////////////////////////////////////////////////////////////////////

//std::string FlowShopOpCrossoverQuad::className() const
//  {
//    return "FlowShopOpCrossoverQuad";
//  }


//bool FlowShopOpCrossoverQuad::operator()(FlowShop & _flowshop1, FlowShop & _flowshop2)
//{
//  bool oneAtLeastIsModified;
//  // computation of the 2 random points
//  unsigned int point1, point2;
//  do
//    {
//      point1 =  rng.random(std::min(_flowshop1.size(), _flowshop2.size()));
//      point2 =  rng.random(std::min(_flowshop1.size(), _flowshop2.size()));
//    }
//  while (fabs((double) point1-point2) <= 2);
//  // computation of the offspring
//  FlowShop offspring1 = generateOffspring(_flowshop1, _flowshop2, point1, point2);
//  FlowShop offspring2 = generateOffspring(_flowshop2, _flowshop1, point1, point2);
//  // does at least one genotype has been modified ?
//  if ((_flowshop1 != offspring1) || (_flowshop2 != offspring2))
//    {
//      // update
//      _flowshop1.value(offspring1);
//      _flowshop2.value(offspring2);
//      // at least one genotype has been modified
//      oneAtLeastIsModified = true;
//    }
//  else
//    {
//      // no genotype has been modified
//      oneAtLeastIsModified = false;
//    }
//  // return 'true' if at least one genotype has been modified
//  return oneAtLeastIsModified;
//}


//FlowShop FlowShopOpCrossoverQuad::generateOffspring(const FlowShop & _parent1, const FlowShop & _parent2, unsigned int _point1, unsigned int _point2)
//{
//  FlowShop result = _parent1;
//  std::vector<bool> taken_values(result.size(), false);
//  if (_point1 > _point2)
//    std::swap(_point1, _point2);
//  /* first parent */
//  for (unsigned int i=0 ; i<=_point1 ; i++)
//    {
//      // result[i] == _parent1[i]
//      taken_values[_parent1[i]] = true;
//    }
//  for (unsigned int i=_point2 ; i<result.size() ; i++)
//    {
//      // result[i] == _parent1[i]
//      taken_values[_parent1[i]] = true;
//    }
//  /* second parent */
//  unsigned int i = _point1+1;
//  unsigned int j = 0;
//  while (i<_point2 && j<_parent2.size())
//    {
//      if (! taken_values[_parent2[j]])
//        {
//          result[i] = _parent2[j];
//          i++;
//        }
//      j++;
//    }
//  return result;
//}

