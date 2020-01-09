#ifndef REPAIRING_H
#define REPAIRING_H

#include <utils/eoUpdater.h>
#include <eoPop.h>
#include "Chromosome.h"


using namespace std;


/**
   an eoUpdater that repairs the timetables
*/
template <typename EOT> // Added
class eoRepair : public eoUpdater
{
public :
    /** Default Ctor - requires a reference to the chromosome to repair */
//    eoRepair(eoPop<Chromosome>& _pop) : pop(_pop) { }
    eoRepair(eoPop<EOT>& _pop) : pop(_pop) { }


    virtual string className(void) const { return "eoRepair"; }

    /** Repair chromosome */
    virtual void operator()() {
//        cout << endl << endl << "Packing" << endl;
        int infeasibleSolsAbove = 0, infeasibleSolsBelow = 0;

//        for (eoPop<Chromosome>::iterator it = pop.begin(); it != pop.end(); ++it) {
        for (typename eoPop<EOT>::iterator it = pop.begin(); it != pop.end(); ++it) {
            // Verify if numPeriods exceeds Max range
            if ((*it).getNumPeriods() > (*it).getRange()[1]) {
                // Timetable period packing
                (*it).pack();
                infeasibleSolsAbove++;
            }
            // Verify if numPeriods is below Min range
            else if ((*it).getNumPeriods() < (*it).getRange()[0]) {
                // Timetable period expansion
                (*it).expand();
                infeasibleSolsBelow++;
            }
        }
//        cout << "Found " << infeasibleSolsAbove << " infeasible Sols Above, and " << infeasibleSolsBelow << " infeasible Sols Below" << endl;
    }

private:
//    eoPop<Chromosome>& pop;
    eoPop<EOT>& pop;
};

#endif // REPAIRING_H



/*

%////////////////////////////////////////////////////////////////////////
%
% Timetable packing
%
%////////////////////////////////////////////////////////////////////////
% Pack timetables
function Pop = packTimetables(Data, Pop)
%     fprintf('\nPacking timetables...\n');
    for i = 1 : length(Pop)
        x = Pop{i};
        if (x.NumPeriods > Data.MaxPeriods)
            % Timetable period packing
%             fprintf('Packing timetable %d\n', i);
            x = periodPacking(Data, x);
        elseif (x.NumPeriods < Data.MinPeriods)
            % Timetable period expansion
%             fprintf('Expanding timetable %d\n', i);
            x = periodExpansion(Data, x);
        end
        Pop{i} = x;
    end
end


% Timetable period expansion
function timetable = periodExpansion(Data, timetable)
    % Period expansion: The operation first adds empty periods
    % to the end of the timetable such that the timetable
    % length is equal to a random number within the desired range.
    % A clash list, consisting of all exams that are involved in at
    % least one clash, is also maintained. An exam is randomly
    % selected from the clash list and the operation searches in a
    % random order for a period which the selected exam can be
    % rescheduled without causing any clashes while maintaining
    % feasibility. The exam remains intact if no such period exists.
    % The operation ends after one cycle through all exams in the
    % clash list.

    % Feasible number of periods
    numPeriods = randi([Data.MinPeriods Data.MaxPeriods]);
    numExtraPeriods = numPeriods - timetable.NumPeriods;
    % add empty periods to the end of the timetable such that the timetable
    % length is equal to a random number within the desired range.
    for i = 1 : numExtraPeriods
        timetable = addPeriodToTimetable(timetable);
    end
    % A clash list, consisting of all exams that are involved in at
    % least one clash, is also maintained.
    [clashList, periodList] = buildClashList(Data, timetable);

    while (~isempty(clashList))
        % An exam is randomly
        % selected from the clash list and the operation searches in a
        % random order for a period which the selected exam can be
        % rescheduled without causing any clashes while maintaining
        % feasibility. The exam remains intact if no such period exists.
        % The operation ends after one cycle through all exams in the
        % clash list.

        % Select an exam randomly from the clash list
        examToScheduleIdx = randi([1 length(clashList)]);
        examToSchedule = clashList(examToScheduleIdx);
        period = periodList(examToScheduleIdx);

        % Remove exam from clash list
        clashList = remove(examToSchedule, clashList);
        periodList = remove(period, periodList);

        % Get available periods (set of periods where the exam can be
        % rescheduled without causing any clashes while maintaining
        % feasibility)
        numPeriods = timetable.NumPeriods;
        AllPeriodsCapacity = zeros(1, numPeriods);
        % Get available periods capacity
        for p = 1 : numPeriods
            if (p ~= period)
                AllPeriodsCapacity(p) = getPeriodCapacity(Data, ...
                    timetable, p, examToSchedule)
            end
        end
        ValidPeriodsCapacityIdxs = find(AllPeriodsCapacity);
        ValidPeriodsCapacity = AllPeriodsCapacity(ValidPeriodsCapacityIdxs);
        if (~isempty(ValidPeriodsCapacity))
            % the operation searches in a
            % random order for a period which the selected exam can be
            % rescheduled without causing any clashes while maintaining
            % feasibility. The exam remains intact if no such period exists.
            % The operation ends after one cycle through all exams in the
            % clash list.
            idx = randi([1 length(ValidPeriodsCapacity)]);
            % Period
            p = ValidPeriodsCapacityIdxs(idx);
            % Reschedule exam
            timetable.Periods{p} = [timetable.Periods{p} examToSchedule];
            % Remove exam from original period
            timetable.Periods{period} = remove(examToSchedule, ...
                                            timetable.Periods{period});

            if (~verifyTimetable(Data, timetable))
                printTimetable(Data, timetable);
                fprintf('[After expansion] Timetable is not valid\n');
                pause
            end
            if (~isFeasibleTimetable(Data, timetable))
                printTimetable(Data, timetable);
                fprintf('[After expansion] Timetable is not feasible\n');
                pause
            end
        else
%             disp('no valid periods')
%             pause
        end
    end
    % Recompute number of clashes
    timetable = computeNumClashes(Data, timetable)
end


% Get period's available capacity
function capacity = getPeriodCapacity(Data, timetable, period, examToSchedule)
    capacity = 0;
    % Verify if period 'p' is a feasible period for current exam
    if (isFeasibleExamPeriod(Data, timetable, period, examToSchedule))
        % Compute number of clashes between 'exam' and other exams
        % in periods p-1 and p+1
        numClashes = computeNumClashesExamPeriod(Data, timetable, period, examToSchedule);
        if (numClashes == 0)
            capacity = 1;
        end
    end
 end


function feasible = isFeasibleTimetable(Data, timetable)
    feasible = true;
    for p = 1 : timetable.NumPeriods
        examList = timetable.Periods{p};
        for i = 1 : length(examList)
            feasible = isFeasibleExamPeriod(Data, timetable, p, examList(i));
            if (~feasible)
                disp('Not feasible')
                Data.Classes{examList(i)}
                p
                return;
            end
        end
    end
end


% Timetable period packing
function timetable = periodPacking(Data, timetable)
    % Period packing: Starting from the period with the
    % smallest number of students, the operation searches in order
    % of available period capacity, starting from the smallest,
    % for a period which can accommodate exams from the former
    % without causing any clashes while maintaining feasibility.
    % The operation stops when it goes one cycle through all periods
    % without rescheduling any exam or when the timetable
    % length is reduced to a random number within the desired
    % range.

    timetable = removeEmptyPeriods(Data, timetable);
    numPeriods = timetable.NumPeriods;
    % Get period with the smallest number of students.
    minPeriod = 1;
    min = getNumberStudentsPeriod(Data, timetable, 1);
    for p = 2 : numPeriods
        count = getNumberStudentsPeriod(Data, timetable, p);
        if (count < min)
            min = count;
            minPeriod = p;
        end
    end
    examList = timetable.Periods{minPeriod};
    % The operation searches in order of available period capacity,
    % starting from the smallest, for a period which can accommodate exams
    % from the former without causing any clashes while maintaining
    % feasibility.
    % The operation stops when it goes one cycle through all periods
    % without rescheduling any exam or when the timetable
    % length is reduced to a random number within the desired range.
    AllPeriodsCapacity = zeros(1, numPeriods);
    % Get available periods capacity
    for p = 1 : numPeriods
        if (p ~= minPeriod)
            AllPeriodsCapacity(p) = getPeriodAvailableCapacity(Data, timetable, p, minPeriod);
        end
    end
    ValidPeriodsCapacityIdxs = find(AllPeriodsCapacity);
    ValidPeriodsCapacity = AllPeriodsCapacity(ValidPeriodsCapacityIdxs);

    if (~isempty(ValidPeriodsCapacity))
       % Sort periods by descending order of available period capacity
       [V, I] = sort(ValidPeriodsCapacity);
       % Reschedule an exam
       LeastAvailablePeriodsIdxs = find(ValidPeriodsCapacity == ValidPeriodsCapacity(1));
       idx = randi([1, length(LeastAvailablePeriodsIdxs)]);

       p = ValidPeriodsCapacityIdxs(LeastAvailablePeriodsIdxs(idx));

       [numExams, ExamList] = getPeriodAvailableCapacity(Data, timetable, ...
           p, minPeriod);
       % Reschedule exams in the returned exam list from period 'minPeriod'
       % and insert them in period 'p'
       examListMinPeriod = timetable.Periods{minPeriod};
       for i = 1 : length(ExamList)
           exam = ExamList(i);
           % Remove exam from minPeriod
           examListMinPeriod = remove(exam, examListMinPeriod);
           % Reschedule exam
           timetable.Periods{p} = [timetable.Periods{p} exam];
       end
       timetable.Periods{minPeriod} = examListMinPeriod;

       if (isempty(examListMinPeriod))
           % Remove minPeriod
            timetable = removeEmptyPeriods(Data, timetable);
        end
    end
    % Recompute number of clashes
    timetable = computeNumClashes(Data, timetable);
end

function timetable = removeEmptyPeriods(Data, timetable)
    timetable.Periods = timetable.Periods(~cellfun('isempty', timetable.Periods));
    timetable.NumPeriods = length(timetable.Periods);
end


% Get period's available capacity
function [capacity, exams] = getPeriodAvailableCapacity(Data, timetable, p, minPeriod)
    % p period's capacity is defined as the number of exams from
    % period 'minPeriod' which can accommodated without causing
    % any clashes while maintaining feasibility.
    capacity = 0;
    exams = [];
    examList = timetable.Periods{minPeriod};
    for i = 1 : length(examList)
        exam = examList(i);
        % Verify if period 'p' is a feasible period for current exam
        if (isFeasibleExamPeriod(Data, timetable, p, exam))
            % Compute number of clashes between 'exam' and other exams
            % in periods p-1 and p+1
            numClashes = computeNumClashesExamPeriod(Data, timetable, p, exam);
            if (numClashes == 0)
                capacity = capacity+1;
                exams = [exams exam];
            end
        end
    end
 end

% Compute number of clashes between 'exam' and other exams
% in periods p-1 and p+1
function numClashes = computeNumClashesExamPeriod(Data, timetable, p, exam)
    if (p == 6 || p == 12 || p == 18)
        numClashes1 = 0;
    else
        numClashes1 = 0;
        if (p+1 <= timetable.NumPeriods)
            examListP1 = timetable.Periods{p+1};
        else
            examListP1 = [];
        end
        numExamsP1 = length(examListP1);
        for i = 1 : numExamsP1
            numStudents = Data.ConflictMatrix(examListP1(i), exam);
            if (numStudents > 0)
                numClashes1 = numClashes1 + numStudents;
            end
        end
    end

    if (p == 7 || p == 13 || p == 19)
        numClashes0 = 0;
    else
        numClashes0 = 0;
        if (p-1 >= 1)
            examListP0 = timetable.Periods{p-1};
        else
            examListP0 = [];
        end
        numExamsP0 = length(examListP0);
        for i = 1 : numExamsP0
            numStudents = Data.ConflictMatrix(examListP0(i), exam);
            if (numStudents > 0)
                numClashes0 = numClashes0 + numStudents;
            end
        end
    end
    numClashes = numClashes0 + numClashes1;
end


% Get number of students in a period
function studentCount = getNumberStudentsPeriod(Data, timetable, period)
    examList = timetable.Periods{period};
    studentCount = 0;
    for i = 1 : length(examList)
        studentCount = studentCount + Data.ExamCounts(examList(i));
    end
end
*/

