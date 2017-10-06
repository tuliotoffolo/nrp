/**
 * This class creates and keeps the sub-problems derived from the
 * decomposition approaches, ensuring each one is solved at most once
 * unless an improvement may be obtained.
 *
 * @author Tulio Toffolo
 */

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <string>
#include "decomposition.h"

// decompositions
#define DAYS   0
#define NURSES 2
#define SHIFTS 1


/********************************
 * Decomposition Manager
 ********************************/

DecompositionsManager::DecompositionsManager(Instance &instance, int nDays, int daysStep, int maxDays, int nNurses, int nursesStep, int maxNurses)
{
    improved = false;

    // creating decomposition subproblems
    decompositions.push_back(new DaysDecomposition(instance, nDays, daysStep, maxDays));
    decompositions.push_back(new NursesDecomposition(instance, nNurses, nursesStep, maxNurses));
    decompositions.push_back(new ShiftsDecomposition(instance));

    sortedDecompositions.insert(sortedDecompositions.begin(), decompositions.begin(), decompositions.end());
}

Subproblem *DecompositionsManager::getSubproblem(ostream &output)
{
    // checking number of available subproblems
    int nSubproblems = 0;
    for (auto decomposition : sortedDecompositions)
        nSubproblems += decomposition->getNSubproblems();
    if (nSubproblems == 0) {
        for (auto decomposition : sortedDecompositions) {
            decomposition->reset(improved);
            nSubproblems += decomposition->getNSubproblems();
        }
        improved = false;
    }
    if (nSubproblems == 0) return lastSubproblem = nullptr;

    // selecting a random subproblem (in fact, a random decomposition)
    int randomNumber = rand() % nSubproblems;

    // sortedDecompositions.insert(sortedDecompositions.begin(), decompositions.begin(), decompositions.end());
    // std::random_shuffle(sortedDecompositions.begin(), sortedDecompositions.end());

    for (auto decomposition : decompositions) {
        if (randomNumber >= decomposition->getNSubproblems()) {
            randomNumber -= decomposition->getNSubproblems();
            continue;
        }

        Subproblem *subproblem = decomposition->getSubproblem();
        if (subproblem == nullptr) continue;

        if (decomposition->index == DAYS) {
            output << "#################################" << endl;
            output << "# Fix DAYS (" << subproblem->name << ")" << endl;
            output << "# index: " << subproblem->index << ", fixed vars: " << subproblem->varsToFix.size() << endl;
            output << "#################################" << endl << endl;
        }
        else if (decomposition->index == NURSES) {
            output << "#################################" << endl;
            output << "# Fix NURSES (" << subproblem->index << ")" << endl;
            output << "# index: " << subproblem->index << ", fixed vars: " << subproblem->varsToFix.size() << endl;
            output << "#################################" << endl << endl;
        }
        else if (decomposition->index == SHIFTS) {
            output << "#################################" << endl;
            output << "# Fix SHIFTS (" << subproblem->index << ")" << endl;
            output << "# index: " << subproblem->index << ", fixed vars: " << subproblem->varsToFix.size() << endl;
            output << "#################################" << endl << endl;
        }

        lastDecomposition = decomposition->index;
        return lastSubproblem = subproblem;
    }
    return lastSubproblem = nullptr;
}

void DecompositionsManager::reset()
{
    improved = true;
    //for (auto decomposition:sortedDecompositions) {
    //    if (lastSubproblem != nullptr && decomposition->index == lastDecomposition)
    //        decomposition->reset(lastSubproblem);
    //    else
    //        decomposition->reset();
    //}
}


/********************************
 * Common to all decompositions
 ********************************/

int Decomposition::getVarIndex(int nurse, int shift, int day)
{
    int index = nurse * instance.nShifts * instance.nDays
                + shift * instance.nDays
                + day;
#ifdef DEBUG
    if (index > x.size()) {
    cerr << "Index error concerning vector X within class Solution:" << endl;
    cerr << "  size = " << x.size() << " and requested index = " << index << endl;
    exit(0);
}
#endif
    return index;
}


/********************************
 * DAYS Decomposition
 ********************************/

DaysDecomposition::DaysDecomposition(Instance &instance, int nDays, int daysStep, int maxDays)
        : Decomposition(DAYS, instance), nDays(nDays), daysStep(daysStep), maxDays(maxDays)
{
    makeSubproblems();
}

Subproblem *DaysDecomposition::getSubproblem()
{
    if (nSubproblems == 0 && nDays < maxDays) {
        nDays++;
        daysStep = nDays / 2;
        makeSubproblems();
    }

    if (nSubproblems == 0)
        return nullptr;

    int p = rand() % nSubproblems;
    Subproblem *subproblem = subproblems[p];
    swap(subproblems[p], subproblems[--nSubproblems]);
    subproblems[p]->sortedIndex = p;
    subproblems[nSubproblems]->sortedIndex = nSubproblems;
    return subproblem;
}

void DaysDecomposition::makeSubproblems()
{
    for (Subproblem *subproblem : subproblems)
        delete subproblem;
    subproblems.clear();
    sortedSubproblems.clear();

    for (int day = 0; day < instance.nDays; day += daysStep) {
        // creating sub-problem {i, ..., i + nDays}
        Subproblem *subproblem = new Subproblem();

        for (int d = 0; d < instance.nDays; d++) {
            if (d >= day && d < day + nDays) continue;
            if (day + nDays > instance.nDays && (d >= day || d < (day + nDays) % instance.nDays)) continue;
            for (int nurse = 0; nurse < instance.nNurses; nurse++)
                for (int shift = 0; shift < instance.nShifts; shift++)
                    subproblem->varsToFix.push_back(getVarIndex(nurse, shift, d));
        }

        subproblem->index = (int) subproblems.size();
        sprintf(subproblem->name, "Days %d-%d", day, day + nDays - 1);
        subproblems.push_back(subproblem);
    }
    nSubproblems = (int) subproblems.size();

    // creating 'sorted' subproblems vector
    sortedSubproblems.insert(sortedSubproblems.begin(), subproblems.begin(), subproblems.end());
}

void DaysDecomposition::reset(bool improved)
{
    if (nDays == maxDays) {
        nSubproblems = (int) subproblems.size();
    }
    else if (nSubproblems == 0 && nDays < maxDays) {
        nDays += 2;
        daysStep++;
        makeSubproblems();
    }
}

void DaysDecomposition::reset(bool improved, Subproblem *ignore)
{
    if (nDays == maxDays) {
        nSubproblems = (int) subproblems.size();
        int p = ignore->sortedIndex;
        swap(subproblems[p], subproblems[--nSubproblems]);
    }
    else {
        reset(improved);
    }
}


/********************************
 * NURSES Decomposition
 ********************************/

NursesDecomposition::NursesDecomposition(Instance &instance, int nNurses, int nursesStep, int maxNurses)
        : Decomposition(NURSES, instance), nNurses(nNurses), nursesStep(nursesStep), maxNurses(maxNurses)
{
    makeSubproblems();
}

Subproblem *NursesDecomposition::getSubproblem()
{
    if (nSubproblems == 0 && nNurses < maxNurses) {
        nNurses++;
        nursesStep = nNurses / 2;
        makeSubproblems();
    }

    if (nSubproblems == 0)
        return nullptr;

    int p = rand() % nSubproblems;
    Subproblem *subproblem = subproblems[p];
    swap(subproblems[p], subproblems[--nSubproblems]);
    subproblems[p]->sortedIndex = p;
    subproblems[nSubproblems]->sortedIndex = nSubproblems;
    return subproblem;
}

void NursesDecomposition::makeSubproblems()
{
    if (nNurses == 0) return;

    // creating random nurse sorting
    vector<int> sortedNurses((size_t) instance.nNurses);
    for (int i = 0; i < sortedNurses.size(); i++)
        sortedNurses[i] = i;
    std::random_shuffle(sortedNurses.begin(), sortedNurses.end());

    for (Subproblem *subproblem : subproblems)
        delete subproblem;
    subproblems.clear();
    sortedSubproblems.clear();

    for (int nurse = 0; nurse < instance.nDays; nurse += nursesStep) {
        // creating sub-problem {i, ..., i + nDays}
        Subproblem *subproblem = new Subproblem();

        for (int n = 0; n < instance.nNurses; n++) {
            if (n >= nurse && n < nurse + nNurses) continue;
            if (nurse + nNurses > instance.nNurses && (n >= nurse || n < (nurse + nNurses) % instance.nNurses)) continue;
            for (int day = 0; day < instance.nDays; day++)
                for (int shift = 0; shift < instance.nShifts; shift++)
                    subproblem->varsToFix.push_back(getVarIndex(sortedNurses[n], shift, day));
        }

        subproblem->index = (int) subproblems.size();
        sprintf(subproblem->name, "Nurses %d-%d", nurse, nurse + nNurses - 1);
        subproblems.push_back(subproblem);
    }
    nSubproblems = (int) subproblems.size();

    // creating 'sorted' subproblems vector
    sortedSubproblems.insert(sortedSubproblems.begin(), subproblems.begin(), subproblems.end());
}

void NursesDecomposition::reset(bool improved)
{
    if (nSubproblems == 0 && nNurses < maxNurses) {
        nNurses++;
        nursesStep++;
    }
    makeSubproblems();
}

void NursesDecomposition::reset(bool improved, Subproblem *ignore)
{
    reset(improved);
}

//Subproblem *NursesDecomposition::getSubproblem()
//{
//    if (nSubproblems == 0)
//        return nullptr;
//
//    int p = rand() % nSubproblems;
//    swap(combinations[p], combinations[--nSubproblems]);
//
//    int *tuple = combinations[nSubproblems];
//
//    refSubproblem->varsToFix.clear();
//
//    for (int nurse = 0; nurse < instance.nNurses; nurse++) {
//        if (nurse == tuple[0] || nurse == tuple[1] || nurse == tuple[2])
//            continue;
//
//        for (int day = 0; day < instance.nDays; day++)
//            for (int shift = 0; shift < instance.nShifts; shift++)
//                refSubproblem->varsToFix.push_back(getVarIndex(nurse, shift, day));
//    }
//
//    refSubproblem->index = p;
//    refSubproblem->sortedIndex = p;
//    sprintf(refSubproblem->name, "Nurses %d,%d,%d", tuple[0], tuple[1], tuple[2]);
//    return refSubproblem;
//}
//
//void NursesDecomposition::makeCombinations()
//{
//    while (!combinations.empty()) {
//        delete combinations.back();
//        combinations.pop_back();
//    }
//
//    for (int n1 = 0; n1 < instance.nNurses; n1++)
//        for (int n2 = n1 + 1; n2 < instance.nNurses; n2++)
//            for (int n3 = n2 + 1; n3 < instance.nNurses; n3++)
//                combinations.push_back(new int[3]{n1, n2, n3});
//
//    nSubproblems = combinations.size();
//}
//
//void NursesDecomposition::reset()
//{
//    nSubproblems = combinations.size();
//}
//
//void NursesDecomposition::reset(Subproblem *ignore)
//{
//    reset();
//    int p = ignore->sortedIndex;
//    swap(combinations[p], combinations[--nSubproblems]);
//}


/********************************
 * SHIFTS Decomposition
 ********************************/

ShiftsDecomposition::ShiftsDecomposition(Instance &instance)
        : Decomposition(SHIFTS, instance)
{
    for (int shift = 0; shift < instance.nShifts; shift++) {
        Subproblem *subproblem = new Subproblem();

        for (int s = 0; s < instance.nShifts; s++) {
            if (s == shift) continue;

            for (int day = 0; day < instance.nDays; day++) {
                for (int nurse = 0; nurse < instance.nNurses; nurse++)
                    subproblem->varsToFix.push_back(getVarIndex(nurse, s, day));
            }
        }

        subproblem->index = (int) subproblems.size();
        subproblems.push_back(subproblem);
    }
    nSubproblems = (int) subproblems.size();

    // creating 'sorted' subproblems vector
    sortedSubproblems.insert(sortedSubproblems.begin(), subproblems.begin(), subproblems.end());
    std::reverse(sortedSubproblems.begin(), sortedSubproblems.end());
}

Subproblem *ShiftsDecomposition::getSubproblem()
{
    if (nSubproblems == 0)
        return nullptr;

    int p = rand() % nSubproblems;
    swap(subproblems[p], subproblems[--nSubproblems]);
    subproblems[p]->sortedIndex = p;
    subproblems[nSubproblems]->sortedIndex = nSubproblems;
    return sortedSubproblems[nSubproblems];
}

void ShiftsDecomposition::reset(bool improved)
{
    nSubproblems = (int) subproblems.size();
}

void ShiftsDecomposition::reset(bool improved, Subproblem *ignore)
{
    reset(improved);
    int p = ignore->sortedIndex;
    swap(subproblems[p], subproblems[--nSubproblems]);
}
