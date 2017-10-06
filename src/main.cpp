/**
 * This file contains the decomposition-based algorithms developed for the NRP.
 *
 * Note that this version uses CPLEX and GUROBI directly instead of employing the
 * Open Solver Interface (OSI) as it reduced the computational overhead and
 * resulted in less memory usage.
 *
 * @author Tulio Toffolo
 */

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <vector>
#include <cmath>
#include <algorithm>
#include <ctime>
#include <climits>
#include <string>

#include "decomposition.h"
#include "instance.h"
#include "solution.h"


using namespace std;

// --------------------------------------------------------------------------

#define ORIG_PROBLEM  0
#define DECOMP        1
#define FIX_CONTRACTS 2
#define FIX_WEEKS     3
#define FIX_SHIFTS    4
#define FIX_DAYS      5
#define FIX_WEEKENDS  6
#define FIX_NURSES    7
#define VND           10

// --------------------------------------------------------------------------

char *LPFILE = nullptr;         // input LP file
char *DATFILE = nullptr;        // input DAT file
char *INISOLFILE = nullptr;        // input solution(s) file
char *OUTFILE = nullptr;        // prefix for output files: '.sol', '.log' and '.timelog'

int RSEED = 9;                  // random seed
int METHOD = DECOMP;            // selected method
int K = 0;                      // parameter k
int RANDOM_CHOICE = 1;
int NURSE_CHOICE = 0;
int MIP_EMPHASIS = 2;           // parameter mip_emphasis (CPLEX)
int MIP_FOCUS = 1;              // parameter mip_focus (GUROBI)
int THREADS = 0;                // number of threads
int SOL_COUNT = 1;              // number of solutions in pool
int LOCAL_TIME_LIMIT = 0;       // local time limit (per iteration)
int TIME_LIMIT = 0;             // global time limit
int ITER_LIMIT = 0;             // iteration limit
int MAX_NODES = 0;              // maximum number of nodes
int OBJ_CUT = 0;                // add objective function constraint
int CPLEX_CUTS = 1;             // use CPLEX cuts
int CALLBACKS = 0;              // enable/disable callbacks

int N_DAYS = 4;
int DAYS_STEP = 3;
int MAX_NDAYS = 30;

int N_NURSES = 2;
int NURSES_STEP = 4;
int MAX_NNURSES = 5;

int LOWER_BOUND = 0;

int global_counter = 0;
int best_objval = INT_MAX;
vector<Solution *> solutions;

// --------------------------------------------------------------------------

bool solve(Instance &instance);
static void usage(const char *progname);

// --------------------------------------------------------------------------

/********************************
 * Implementation for CPLEX
 ********************************/


#ifdef CPLEX

#include "ilcplex/cplex.h"
#include "ilcplex/ilocplex.h"

#undef GUROBI

IloEnv env;

IloExtractableArray fixVars_All(Solution *sol, IloCplex &cplex, IloModel &model, IloNumVarArray &vars)
{
    IloExtractableArray constraints(env);
    for (int i = 0; i < sol->x.size(); i++)
        if (sol->x[i] == 1)
            constraints.add(vars[sol->idx[i]] == 1.0);

    return constraints;
}

IloExtractableArray fixVars_List(Solution *sol, IloCplex &cplex, IloModel &model, IloNumVarArray &vars, vector<int> &varsToFix)
{
    IloExtractableArray constraints(env);

    for (int var : varsToFix) {
        if (round(sol->getX(var)) == 1) {
            constraints.add(vars[sol->idx[var]] == 1.0);
        }
    }

    return constraints;
}

ILOMIPINFOCALLBACK7(infoCallback, IloCplex, cplex, IloNumVarArray &, vars, Solution*, sol, double &, objval,
                    IloNum &, timeStart, IloNum &, timeImproved, ostream &, output)
{
    if (hasIncumbent() && (round(objval) > round(getIncumbentObjValue()) || timeImproved == -1)) {
        objval = getIncumbentObjValue();

        double exe_time;
        if (timeImproved == -1) {
            exe_time = INISOLFILE ? 0 : cplex.getCplexTime() - timeStart + sol->totaltime;
            cout << setprecision(0) << "--> First (initial) integer solution [ time = " << exe_time << "s ; obj = " << objval << " ] <--" << endl;
        }
        else {
            exe_time = cplex.getCplexTime() - timeStart + sol->totaltime;
            cout << setprecision(0) << "--> Improvement on the best solution [ time = " << exe_time << "s ; obj = " << objval << " ] <--" << endl;
        }
        output << setprecision(0) << exe_time << ";" << objval << ";" << setprecision(2) << sol->getGap(LOWER_BOUND) << endl;
        output.flush();

        timeImproved = cplex.getCplexTime();

        IloNumArray x(env);
        getIncumbentValues(x, vars);
        sol->read(x, round(getIncumbentObjValue()));
        sol->save((string(OUTFILE) + ".sol").c_str(), exe_time);
        //sol->saveXML((string(OUTFILE) + ".sol").c_str());
    }

    IloNum timeUsed = cplex.getCplexTime() - timeImproved;
    IloNum totalTime = cplex.getCplexTime() - timeStart - sol->getTime();
    if (LOCAL_TIME_LIMIT > 0 && timeImproved > 0 && timeUsed > LOCAL_TIME_LIMIT) {
        cout << "------ Local time limit exceeded ------" << endl;
        timeImproved = 0;
        abort();
    }
    if (TIME_LIMIT > 0 && totalTime > TIME_LIMIT) {
        cout << "------ Total time limit exceeded ------" << endl;
        timeImproved = 0;
        abort();
    }
}

bool validateSolution(IloCplex &cplex, IloModel &model, IloNumVarArray &vars, Solution *sol, IloNum &objval,
                      IloNum &timeStart, IloNum &timeImproved, ostream &output)
{
    // printing a report of the parameters in the output file
    output << "##############################" << endl;
    output << "# Run Information" << endl;
    output << "##################" << endl;
    output << endl;
    output << "*** VALIDATING INTEGER SOLUTION ***" << endl;
    output << endl;
    output << "LPFILE  = " << LPFILE << endl;
    output << "DATFILE = " << DATFILE << endl;
    output << "OUTFILE = " << OUTFILE << endl;
    if (INISOLFILE) output << "SOLFILE = " << INISOLFILE << endl;
    output << endl;

    // executing CPLEX
    timeStart = cplex.getCplexTime();
    IloExtractableArray fixVars = fixVars_All(sol, cplex, model, vars);
    model.add(fixVars);
    cplex.solve();
    sol->totaltime += cplex.getCplexTime() - timeStart;

    // saving solution
    string xml(DATFILE);
    xml[xml.length() - 3] = 'x';
    xml[xml.length() - 2] = 'm';
    xml[xml.length() - 1] = 'l';
    cout << xml << endl;
    IloNumArray x(env);
    cplex.getValues(x, vars);
    cplex.writeSolution((string(OUTFILE) + ".xml").c_str());
    sol->read(x, cplex.getObjValue());
    sol->save((string(OUTFILE) + ".sol").c_str());
    sol->saveXML((string(OUTFILE) + ".xml").c_str(), Instance::readDateFromTXT(xml.c_str()));

    // printing a report of the final result in the output file
    output << "##############################" << endl;
    output << "# Solution" << endl;
    output << "################" << endl;
    output << endl;
    output << "Obj. = " << sol->objval << endl;
    output << "Time = " << setprecision(0) << sol->getTime() << endl;
    output << endl;
    int nurse, day;
    char shift[10];
    for (IloInt j = 0; j < vars.getSize(); j++)
        if (sscanf(vars[j].getName(), "x(N%d,%[a-zA-Z],%d)", &nurse, shift, &day) == 3 && x[j] == 1)
            output << left << setw(12) << vars[j].getName() << " = " << setprecision(0) << x[j] << endl;
    output << endl;

    output.flush();

    string outfile = OUTFILE;
    ofstream gurobi((outfile + ".gurobi").c_str(), ios::out);
    for (IloInt j = 0; j < vars.getSize(); j++)
        if (round(x[j]) == 1)
            gurobi << vars[j].getName() << " " << round(x[j]) << endl;
    gurobi.close();

    return 1;
}

bool runOrigProblem(IloCplex &cplex, IloNumVarArray &vars, Solution *sol, IloNum &objval, IloNum &timeStart,
                    IloNum &timeImproved, ostream &output)
{
    // printing a report of the parameters in the output file
    output << "##############################" << endl;
    output << "# Run Information" << endl;
    output << "##################" << endl;
    output << endl;
    output << "LPFILE  = " << LPFILE << endl;
    output << "DATFILE = " << DATFILE << endl;
    output << "OUTFILE = " << OUTFILE << endl;
    if (INISOLFILE) output << "SOLFILE = " << INISOLFILE << endl;
    output << endl;
    output << "METHOD       = ORIGINAL LP" << endl;
    output << "MIP_EMPHASIS = " << MIP_EMPHASIS << endl;
    output << "THREADS      = " << THREADS << endl;
    if (TIME_LIMIT) output << "TIME_LIMIT   = " << TIME_LIMIT << endl;
    if (LOCAL_TIME_LIMIT) output << "LTIME_LIMIT  = " << LOCAL_TIME_LIMIT << endl;
    if (MAX_NODES) output << "MAX_NODES    = " << MAX_NODES << endl;
    if (ITER_LIMIT) output << "ITER_LIMIT   = " << ITER_LIMIT << endl;
    output << endl;

    // executing CPLEX
    timeStart = cplex.getCplexTime();
    cplex.solve();
    sol->totaltime += cplex.getCplexTime() - timeStart;

    // printing a report of the final result in the output file
    output << "##############################" << endl;
    output << "# Result report " << endl;
    output << "################" << endl;
    output << endl;
    output << setiosflags(ios::fixed);
    output << "Solution status  = " << cplex.getStatus() << endl;
    output << "Lower bound      = " << setprecision(2) << cplex.getBestObjValue() << endl;
    output << "Upper bound      = " << setprecision(0) << cplex.getObjValue() << endl;
    output << "Initial solution = " << sol->objval << endl;
    output << "Solution gap (%) = " << setprecision(2) << 100 * (cplex.getObjValue() - cplex.getBestObjValue()) / cplex.getObjValue() << endl;
    output << "Iterations       = " << cplex.getNiterations() << endl;
    output << "Total time (s)   = " << setprecision(2) << sol->getTime() << endl;
    output << setfill(' ') << endl;
    output.flush();

    // saving solution
    IloNumArray x(env);
    cplex.getValues(x, vars);
    sol->read(x, cplex.getObjValue());

    return 1;
}

bool runDecompositionHeuristic(Instance &instance, IloCplex &cplex, IloModel &model, IloNumVarArray &vars, Solution *sol,
                               IloNum &objval, IloNum &timeStart, IloNum &timeImproved, ostream &output)
{
    // auxiliary variables
    timeImproved = cplex.getCplexTime();
    IloNum lastObjval = sol->objval;
    IloNumArray x(env);
    IloExtractableArray constraints;
    int counter = 0;

    DecompositionsManager manager(instance, N_DAYS, DAYS_STEP, MAX_NDAYS, N_NURSES, NURSES_STEP, MAX_NNURSES);

    while (1) {
        Subproblem *subproblem = manager.getSubproblem(output);
        if (subproblem == nullptr) break;

        if (TIME_LIMIT) cplex.setParam(IloCplex::TiLim, TIME_LIMIT - sol->getTime());

        if (counter == 0) {
            IloExtractableArray fixVars = fixVars_All(sol, cplex, model, vars);
            model.add(fixVars);
            cplex.solve();
            lastObjval = cplex.getObjValue();
            constraints = fixVars_List(sol, cplex, model, vars, subproblem->varsToFix);
            model.remove(fixVars);
            model.add(constraints);
        }
        else {
            IloExtractableArray oldConstraints = constraints;
            constraints = fixVars_List(sol, cplex, model, vars, subproblem->varsToFix);
            model.remove(oldConstraints);
            model.add(constraints);
        }

        // executando o cplex
        timeStart = cplex.getCplexTime();
        cplex.solve();
        sol->totaltime += cplex.getCplexTime() - timeStart;
        counter++;
        global_counter++;

        // checking for improvement
        if (cplex.getObjValue() < lastObjval)
            manager.reset();

        // printing a report concerning the current iteration
        output << setiosflags(ios::fixed);
        output << "r" << setfill('0') << setw(3) << global_counter << ": Subproblem        = " << subproblem->name << endl;
        output << "r" << setfill('0') << setw(3) << global_counter << ": Improvement (%)   = " << 100 * (lastObjval - cplex.getObjValue()) / lastObjval << endl;
        output << "r" << setfill('0') << setw(3) << global_counter << ": Initial solution  = " << setprecision(2) << lastObjval << endl;
        output << "r" << setfill('0') << setw(3) << global_counter << ": Solution obtained = " << setprecision(2) << cplex.getObjValue() << endl;
        output << "r" << setfill('0') << setw(3) << global_counter << ": Time (s)          = " << setprecision(2) << cplex.getCplexTime() - timeStart << endl;
        output << "r" << setfill('0') << setw(3) << global_counter << ": Total time (s)    = " << setprecision(2) << sol->getTime() << endl;
        output << setfill(' ') << endl;
        output.flush();

        // updating the solution
        cplex.getValues(x, vars);
        sol->read(x, cplex.getObjValue());

        if (LOWER_BOUND && sol->objval == LOWER_BOUND) {
            output << "***\nGlobal optimum solution found." << endl;
            cout << endl << "***\nGlobal optimum solution found." << endl;
            break;
        }

        if (TIME_LIMIT > 0 && sol->getTime() > TIME_LIMIT) break;

        lastObjval = sol->objval;
    }
    return 0;
}

bool solve(Instance &instance)
{
    srand(RSEED);

    try {
        Solution *sol = solutions[0];
        string outfile = OUTFILE;
        ofstream output((outfile + ".log").c_str(), ios::app);
        ofstream timelog((outfile + ".timelog").c_str(), ios::app);
        cout << setiosflags(ios::fixed);
        output << setiosflags(ios::fixed);
        timelog << setiosflags(ios::fixed) << "Time(s);Obj;Gap(%)" << endl;
        timelog.precision(2);

        IloModel model(env);
        IloCplex cplex(env);
        IloObjective obj;
        IloNumVarArray vars(env);
        IloRangeArray rng(env);
        IloNum timeStart = cplex.getCplexTime();
        IloNum timeImproved = -1;
        IloNum objval = INISOLFILE ? sol->objval : DBL_MAX;

        cplex.setOut(env.getNullStream());

        // importing folrmulation into cplex
        cplex.importModel(model, LPFILE, obj, vars, rng);
        cplex.extract(model);

        // alterando callbacks
        cplex.use(infoCallback(env, cplex, vars, sol, objval, timeStart, timeImproved, timelog));

        // alterando o nro de threads na execucao do cplex
        cplex.setParam(IloCplex::MIPEmphasis, MIP_EMPHASIS);
        if (THREADS > 0) cplex.setParam(IloCplex::Threads, IloInt(THREADS));
        if (TIME_LIMIT > 0) cplex.setParam(IloCplex::TiLim, IloNum(TIME_LIMIT));
        if (ITER_LIMIT > 0) cplex.setParam(IloCplex::ItLim, IloNum(ITER_LIMIT));

        // desabilitando cortes do cplex
        if (!CPLEX_CUTS) {
            cplex.setParam(IloCplex::Cliques, -1);
            cplex.setParam(IloCplex::Covers, -1);
            cplex.setParam(IloCplex::FlowCovers, -1);
            cplex.setParam(IloCplex::GUBCovers, -1);
            cplex.setParam(IloCplex::FracCuts, -1);
            cplex.setParam(IloCplex::MIRCuts, -1);
            cplex.setParam(IloCplex::FlowPaths, -1);
            cplex.setParam(IloCplex::ImplBd, -1);
            cplex.setParam(IloCplex::DisjCuts, -1);
            cplex.setParam(IloCplex::ZeroHalfCuts, -1);
        }

        // adicionando os indices de conversao
        int nurse, day;
        char shift[10];
        for (int j = 0; j < vars.getSize(); j++) {
            if (sscanf(vars[j].getName(), "x(N%d,%[a-zA-Z],%d)", &nurse, shift, &day) == 3) {
                for (int s = 0; s < solutions.size(); s++) {
                    solutions[s]->setIdx(nurse, string(shift), day - 1, j);
                }
                sol->setIdx(nurse, string(shift), day - 1, j);
            }
        }

        // adicionando a solucao inicial ao cplex
        if (INISOLFILE) {
            IloNumArray newValues(env);
            IloNumVarArray newVars(env);
            for (int i = 0; i < sol->x.size(); i++) {
                //if (sol->idx[i] == -1) continue;
                newVars.add(vars[sol->idx[i]]);
                newValues.add(sol->x[i]);
            }
            cplex.addMIPStart(newVars, newValues);
        }

        // imprimindo no arquivo de saida um pequeno relatorio do resultado final
        {
            output << "##########################################################################################" << endl;
            output << endl;
            output << "########################" << endl;
            output << "# EXECUTION PARAMETERS #" << endl;
            output << "########################" << endl;
            output << endl;
            output << "LPFILE  = " << LPFILE << endl;
            output << "DATFILE = " << DATFILE << endl;
            output << "OUTFILE = " << OUTFILE << ".*" << endl;
            if (INISOLFILE) output << "INISOLFILE = " << INISOLFILE << endl;
            output << endl;
            output << "METHOD        = " << METHOD << endl;
            output << "K             = " << K << endl;
            output << "RANDOM_CHOICE = " << RANDOM_CHOICE << endl;
            output << "NURSE_CHOICE  = " << NURSE_CHOICE << endl;
            output << "N_DAYS         = " << N_DAYS << endl;
            output << "DAYS_STEP     = " << DAYS_STEP << endl;
            output << "MAX_NDAYS     = " << MAX_NDAYS << endl;
            output << "MIP_EMPHASIS  = " << MIP_EMPHASIS << endl;
            output << "THREADS       = " << THREADS << endl;
            if (TIME_LIMIT) output << "TIME_LIMIT   = " << TIME_LIMIT << endl;
            if (LOCAL_TIME_LIMIT) output << "LTIME_LIMIT  = " << LOCAL_TIME_LIMIT << endl;
            if (MAX_NODES) output << "MAX_NODES    = " << MAX_NODES << endl;
            if (ITER_LIMIT) output << "ITER_LIMIT   = " << ITER_LIMIT << endl;
            if (LOWER_BOUND) output << "LOWER_BOUND  = " << LOWER_BOUND << endl;
            output << endl;
        }

        int ini_objval = sol->objval;

        timelog << setprecision(0) << sol->getTime() << ";" << sol->objval << ";" << setprecision(2) << sol->getGap(LOWER_BOUND) << endl;

        if (METHOD == -1)
            validateSolution(cplex, model, vars, sol, objval, timeStart, timeImproved, output);
        else if (METHOD == ORIG_PROBLEM)
            runOrigProblem(cplex, vars, sol, objval, timeStart, timeImproved, output);
        //else if (METHOD == FIX_SHIFTS)
        //    runFixShifts(cplex, model, vars, sol, objval, timeStart, timeImproved, output, K, false);
        //else if (METHOD == FIX_DAYS)
        //    runFixDays(cplex, model, vars, sol, objval, timeStart, timeImproved, output, K, false);
        //else if (METHOD == VND)
        //    runVND(cplex, model, vars, sol, objval, timeStart, timeImproved, output);
        else if (METHOD == DECOMP)
            runDecompositionHeuristic(instance, cplex, model, vars, sol, objval, timeStart, timeImproved, output);

        timelog << setprecision(0) << sol->getTime() << ";" << sol->objval << ";" << setprecision(2) << sol->getGap(LOWER_BOUND) << endl;

        output.close();
        timelog.close();

        sol->save((string(OUTFILE) + ".sol").c_str());
    }
    catch (IloException &e) {
        cerr << "Concert exception caught: " << e << endl;
        return 0;
    }
    catch (...) {
        cerr << "Unknown exception caught" << endl;
        return 0;
    }

    return 1;
}

#endif

/********************************
 * Implementation for GUROBI
 ********************************/

#ifdef GUROBI

#include "gurobi_c++.h"

void fixVars_All(Solution *sol, GRBModel &gurobi, GRBVar *vars)
{
    for (int i = 0; i < sol->x.size(); i++)
        if (sol->x[i] == 1)
            vars[sol->idx[i]].set(GRB_DoubleAttr_LB, sol->x[i]);
}

void unfixVars_All(Solution *sol, GRBModel &gurobi, GRBVar *vars)
{
    for (int i = 0; i < sol->x.size(); i++)
        vars[sol->idx[i]].set(GRB_DoubleAttr_LB, 0.0);
}

void fixVars_List(Solution *sol, GRBModel &gurobi, GRBVar *vars, vector<int> &varsToFix)
{
    for (int var : varsToFix) {
        if (round(sol->getX(var)) == 1) {
            vars[sol->idx[var]].set(GRB_DoubleAttr_LB, 1.0);
        }
    }
}

bool runDecompositionHeuristic(Instance &instance, GRBModel &gurobi, GRBVar *vars, int nVars, Solution *sol,
                               double &objval, ostream &output, ostream &timelog)
{
    // auxiliary variables
    double lastObjval = sol->objval;
    double *x;
    int counter = 0;

    DecompositionsManager manager(instance, N_DAYS, DAYS_STEP, MAX_NDAYS, N_NURSES, NURSES_STEP, MAX_NNURSES);

    while (1) {
        Subproblem *subproblem = manager.getSubproblem(output);
        if (subproblem == nullptr) break;

        if (TIME_LIMIT) gurobi.set(GRB_DoubleParam_TimeLimit, TIME_LIMIT - sol->getTime());

        if (counter == 0) {
            fixVars_All(sol, gurobi, vars);
            gurobi.optimize();
            lastObjval = gurobi.get(GRB_DoubleAttr_ObjVal);
            unfixVars_All(sol, gurobi, vars);
            fixVars_List(sol, gurobi, vars, subproblem->varsToFix);
        }
        else {
            fixVars_List(sol, gurobi, vars, subproblem->varsToFix);
        }

        // executando o cplex
        gurobi.update();
        clock_t start = clock();
        gurobi.optimize();
        double runtime = (double) (clock() - start) / (double) CLOCKS_PER_SEC;
        sol->totaltime += runtime;
        counter++;
        global_counter++;

        // checking for improvement
        if (gurobi.get(GRB_DoubleAttr_ObjVal) < lastObjval) {
            manager.reset();
            timelog << setprecision(0) << sol->getTime() << ";" << sol->objval << ";" << setprecision(2) << sol->getGap(LOWER_BOUND) << endl;
        }

        // printing a report concerning the current iteration
        output << setiosflags(ios::fixed);
        output << "r" << setfill('0') << setw(3) << global_counter << ": Subproblem        = " << subproblem->name << endl;
        output << "r" << setfill('0') << setw(3) << global_counter << ": Improvement (%)   = "
               << 100 * (lastObjval - gurobi.get(GRB_DoubleAttr_ObjVal)) / lastObjval << endl;
        output << "r" << setfill('0') << setw(3) << global_counter << ": Initial solution  = " << setprecision(2) << lastObjval << endl;
        output << "r" << setfill('0') << setw(3) << global_counter << ": Solution obtained = " << setprecision(2) << gurobi.get(GRB_DoubleAttr_ObjVal) << endl;
        output << "r" << setfill('0') << setw(3) << global_counter << ": Time (s)          = " << setprecision(2) << runtime << endl;
        output << "r" << setfill('0') << setw(3) << global_counter << ": Total time (s)    = " << setprecision(2) << sol->getTime() << endl;
        output << setfill(' ') << endl;
        output.flush();

        // updating the solution
        x = gurobi.get(GRB_DoubleAttr_X, vars, nVars);
        sol->read(x, gurobi.get(GRB_DoubleAttr_ObjVal));

        if (LOWER_BOUND && sol->objval == LOWER_BOUND) {
            output << "***\nGlobal optimum solution found." << endl;
            cout << endl << "***\nGlobal optimum solution found." << endl;
            break;
        }

        if (TIME_LIMIT > 0 && sol->getTime() > TIME_LIMIT) break;

        lastObjval = sol->objval;
        unfixVars_All(sol, gurobi, vars);
    }
    return 0;
}

bool solve(Instance &instance)
{
    srand(RSEED);

    try {
        Solution *sol = solutions[0];
        string outfile = OUTFILE;
        ofstream output((outfile + ".log").c_str(), ios::out);
        ofstream timelog((outfile + ".timelog").c_str(), ios::out);
        cout << setiosflags(ios::fixed);
        output << setiosflags(ios::fixed);
        timelog << setiosflags(ios::fixed) << "Time(s) ; Obj" << endl;

        // importing formulation into gurobi
        GRBEnv env;
        GRBModel gurobi(env, LPFILE);
        GRBVar *vars = gurobi.getVars();

        double objval = INISOLFILE ? sol->objval : DBL_MAX;

        gurobi.set(GRB_IntParam_OutputFlag, 0);

        // setting callbacks
        //cplex.use(infoCallback(env, cplex, vars, sol, objval, timeStart, timeImproved, timelog));

        // updating gurobi parameters
        gurobi.set(GRB_IntParam_MIPFocus, MIP_FOCUS);
        if (THREADS > 0) gurobi.set(GRB_IntParam_Threads, THREADS);
        if (TIME_LIMIT > 0) gurobi.set(GRB_DoubleParam_TimeLimit, TIME_LIMIT);
        if (ITER_LIMIT > 0) gurobi.set(GRB_DoubleParam_IterationLimit, ITER_LIMIT);

        // adicionando os indices de conversao
        int nVars = gurobi.get(GRB_IntAttr_NumVars);
        int nurse, day;
        char shift[10];
        for (int j = 0; j < nVars; j++) {
            string varName = vars[j].get(GRB_StringAttr_VarName);
            if (sscanf(varName.c_str(), "x(N%d,%[a-zA-Z],%d)", &nurse, shift, &day) == 3) {
                for (int s = 0; s < solutions.size(); s++) {
                    solutions[s]->setIdx(nurse, string(shift), day - 1, j);
                }
                sol->setIdx(nurse, string(shift), day - 1, j);
            }
        }

        // adicionando a solucao inicial ao cplex
        if (INISOLFILE) {
            GRBVar *newVars = new GRBVar[nVars];
            double *newValues = new double[nVars];
            for (int i = 0; i < sol->x.size(); i++) {
                //if (sol->idx[i] == -1) continue;
                newVars[i] = vars[sol->idx[i]];
                newValues[i] = sol->x[i];
            }
            gurobi.set(GRB_DoubleAttr_Start, newVars, newValues, (int) sol->x.size());
        }

        // imprimindo no arquivo de saida um pequeno relatorio do resultado final
        {
            output << "##########################################################################################" << endl;
            output << endl;
            output << "########################" << endl;
            output << "# EXECUTION PARAMETERS #" << endl;
            output << "########################" << endl;
            output << endl;
            output << "LPFILE  = " << LPFILE << endl;
            output << "DATFILE = " << DATFILE << endl;
            output << "OUTFILE = " << OUTFILE << ".*" << endl;
            if (INISOLFILE) output << "INISOLFILE = " << INISOLFILE << endl;
            output << endl;
            output << "METHOD        = " << METHOD << endl;
            output << "K             = " << K << endl;
            output << "RANDOM_CHOICE = " << RANDOM_CHOICE << endl;
            output << "NURSE_CHOICE  = " << NURSE_CHOICE << endl;
            output << "N_DAYS         = " << N_DAYS << endl;
            output << "DAYS_STEP     = " << DAYS_STEP << endl;
            output << "MAX_NDAYS     = " << MAX_NDAYS << endl;
            output << "MIP_EMPHASIS  = " << MIP_EMPHASIS << endl;
            output << "THREADS       = " << THREADS << endl;
            if (TIME_LIMIT) output << "TIME_LIMIT   = " << TIME_LIMIT << endl;
            if (LOCAL_TIME_LIMIT) output << "LTIME_LIMIT  = " << LOCAL_TIME_LIMIT << endl;
            if (MAX_NODES) output << "MAX_NODES    = " << MAX_NODES << endl;
            if (ITER_LIMIT) output << "ITER_LIMIT   = " << ITER_LIMIT << endl;
            output << endl;
        }

        int ini_objval = sol->objval;

        timelog << setprecision(0) << sol->getTime() << ";" << sol->objval << ";" << setprecision(2) << sol->getGap(LOWER_BOUND) << endl;

        if (METHOD == -1);//validateSolution(cplex, model, vars, sol, objval, timeStart, timeImproved, output);
        else if (METHOD == ORIG_PROBLEM);//runOrigProblem(cplex, vars, sol, objval, timeStart, timeImproved, output);
        else if (METHOD == FIX_SHIFTS);//runFixShifts(cplex, model, vars, sol, objval, timeStart, timeImproved, output, K, false);
        else if (METHOD == FIX_DAYS);//runFixDays(cplex, model, vars, sol, objval, timeStart, timeImproved, output, K, false);
        else if (METHOD == VND);//runVND(cplex, model, vars, sol, objval, timeStart, timeImproved, output);
        else if (METHOD == DECOMP)
            runDecompositionHeuristic(instance, gurobi, vars, nVars, sol, objval, output, timelog);

        timelog << setprecision(0) << sol->getTime() << ";" << sol->objval << ";" << setprecision(2) << sol->getGap(LOWER_BOUND) << endl;

        output.close();
        timelog.close();

        sol->save((string(OUTFILE) + ".sol").c_str());
    }
    catch (GRBException &e) {
        cerr << "Gurobi exception caught: " << e.getMessage() << endl;
        return 0;
    }
    catch (...) {
        cerr << "Unknown exception caught" << endl;
        return 0;
    }

    return 1;
}

#endif

// --------------------------------------------------------------------------

int main(int argc, char *argv[])
{
    // reading arguments
    int index = 0;
    while (index < argc - 1) {
        string arg = string(argv[++index]);

        if (arg == "-dat")
            DATFILE = argv[++index];
        else if (arg == "-lp")
            LPFILE = argv[++index];
        else if (arg == "-out")
            OUTFILE = argv[++index];
        else if (arg == "-ini")
            INISOLFILE = argv[++index];

        else if (arg == "-method")
            METHOD = atoi(argv[++index]);
        else if (arg == "-k")
            K = atoi(argv[++index]);
        else if (arg == "-rseed")
            RSEED = atoi(argv[++index]);
        else if (arg == "-random_choice")
            RANDOM_CHOICE = atoi(argv[++index]);
        else if (arg == "-nurse_choice")
            NURSE_CHOICE = atoi(argv[++index]);
        else if (arg == "-n_days")
            N_DAYS = atoi(argv[++index]);
        else if (arg == "-days_step")
            DAYS_STEP = atoi(argv[++index]);
        else if (arg == "-max_ndays")
            MAX_NDAYS = atoi(argv[++index]);
        else if (arg == "-n_nurses")
            N_NURSES = atoi(argv[++index]);
        else if (arg == "-nurses_step")
            NURSES_STEP = atoi(argv[++index]);
        else if (arg == "-max_nnurses")
            MAX_NNURSES = atoi(argv[++index]);
        else if (arg == "-mip_emphasis")
            MIP_EMPHASIS = atoi(argv[++index]);
        else if (arg == "-threads")
            THREADS = atoi(argv[++index]);
        else if (arg == "-time_limit")
            TIME_LIMIT = atoi(argv[++index]);
        else if (arg == "-ltime_limit")
            LOCAL_TIME_LIMIT = atoi(argv[++index]);
        else if (arg == "-iter_limit")
            ITER_LIMIT = atoi(argv[++index]);
        else if (arg == "-max_nodes")
            MAX_NODES = atoi(argv[++index]);
        else if (arg == "-obj_cut")
            OBJ_CUT = atoi(argv[++index]);
        else if (arg == "-cplex_cuts")
            CPLEX_CUTS = atoi(argv[++index]);
        else if (arg == "-callbacks")
            CALLBACKS = atoi(argv[++index]);
        else if (arg == "-lb")
            LOWER_BOUND = atoi(argv[++index]);

        else if (arg == "--validate-only")
            METHOD = -1;

        else {
            usage(argv[0]);
            cout << "ERROR: Invalid parameter: " << arg << endl << endl;
            exit(EXIT_FAILURE);
        }
    }

    // checking whether all required parameters are present
    if (!DATFILE || !LPFILE || (METHOD >= 0 && !OUTFILE)) {
        usage(argv[0]);
        cout << "ERROR: You must inform the parameters -dat, -lp and -out" << endl << endl;
        exit(EXIT_FAILURE);
    }
    else if (!INISOLFILE && METHOD != ORIG_PROBLEM) {
        usage(argv[0]);
        cout << "ERROR: You must inform the initial solution for the chosen method." << endl << endl;
        exit(EXIT_FAILURE);
    }

    // reading initial solution(s)
    Instance instance(DATFILE);
    if (INISOLFILE) {
        ifstream solfile(INISOLFILE, ifstream::in);
        if (!solfile.good()) {
            cout << "Invalid solution file." << endl << endl;
            exit(EXIT_FAILURE);
        }
        for (int i = 0; i < SOL_COUNT; i++) {
            solutions.push_back(new Solution(instance));
            solutions[i]->read(solfile);
        }
        solfile.close();
    }
    else {
        solutions.push_back(new Solution(instance));
    }

    // sorting pool of solutions
    sort(solutions.begin(), solutions.end(), Solution::compare);

    // executing the solver
    if (solve(instance)) {
        // sorting updated pool of solutions
        sort(solutions.begin(), solutions.end(), Solution::compare);

        cout << endl;
        cout << endl;
        cout << "Objective value: " << solutions[0]->objval << endl;
        cout << "Total time (s): " << solutions[0]->getTime() << endl;
    }
    else {
        cout << endl;
        cout << "The solver has thrown a runtime error.";
    }

    return EXIT_SUCCESS;
}

static void usage(const char *progname)
{
    cout << endl;
    cout << "Usage: " << progname << " <arguments> [options]" << endl;
    cout << endl;
    cout << "Program arguments:" << endl;
    cout << "    -dat <file.dat> : nurse problem" << endl;
    cout << "    -lp <file.lp>   : nurse problem (in cplex lp format)" << endl;
    cout << "    -ini <file.sol> : provide an initial solution" << endl;
    cout << "    -out <file.txt> : prefix of files where the solutions will be saved" << endl;
    cout << endl;
    cout << "Optional parameters (example):" << endl;
    //cout << "    -method <int>       : algorithm to be used (default: " << METHOD << "):" << endl;
    //cout << "                              " << ORIG_PROBLEM << "   original problem" << endl;
    //cout << "                              " << DECOMP << "   latest decomposition-based heuristic" << endl;
    //cout << "                              " << VND << "   variable neighborhood descent" << endl;
    cout << "    -rseed <int>        : random seed (default: " << RSEED << ")" << endl;
    cout << "    -mip_emphasis <int> : mip emphasis (only cplex; default: " << MIP_EMPHASIS << ")." << endl;
    cout << "                             0   Balance optimality and feasibility" << endl;
    cout << "                             1   Emphasize feasibility over optimality" << endl;
    cout << "                             2   Emphasize optimality over feasibility" << endl;
    cout << "                             3   Emphasize moving best bound" << endl;
    cout << "                             4   Emphasize hidden feasible solutions" << endl;
    cout << "    -mip_focus <int>    : mip focus (only gurobi; default: " << MIP_FOCUS << ")." << endl;
    cout << "                             0   Balance optimality and feasibility" << endl;
    cout << "                             1   Focus on feasibility over optimality" << endl;
    cout << "                             2   Focus on optimality over feasibility" << endl;
    cout << "                             3   Focus on improving the best bound" << endl;
    cout << "    -threads <int>      : number of thread the solver may use (default: " << THREADS << ")." << endl;
    cout << "    -time_limit <int>   : runtime limit (default: unlimited)." << endl;
    cout << "    -ltime_limit <int>  : solver runtime limit to solve each subproblem  (default: unlimited)." << endl;
    cout << "    -iter_limit <int>   : maximum number of iterations (default: unlimited)." << endl;
    //cout << "    -obj_cut 0          : if set to 0, the program won't use the objective cut; otherwise it will." << endl;
    //cout << "    -cplex_cuts 1       : if set to 0, the program won't use any of the Cplex cuts; otherwise it will." << endl;
    //cout << "    -callbacks 0        : if set to 0, the program won't use any Callback function; otherwise it will." << endl;
    cout << "    -n_days <int>       : initial number of days in a \"time-based\" sub-problem (default: " << N_DAYS << ")." << endl;
    cout << "    -days_step <int>    : value of the increment on n_days for the \"time-based\" decomposition. (default: " << DAYS_STEP << ")." << endl;
    cout << "    -max_ndays <int>    : maximum value for parameter ndays (default: " << MAX_NDAYS << ")." << endl;
    cout << "    -n_nurses <int>     : initial number of nurses in a \"nurse-based\" sub-problem. (default: " << N_NURSES << ")." << endl;
    cout << "    -nurses_step <int>  : value of the increment on n_nurses for the \"nurse-based\" decomposition (default: " << NURSES_STEP << ")." << endl;
    cout << "    -max_nnurses <int>  : maximum value for parameter n_nurses (default: " << MAX_NNURSES << ")." << endl;
    cout << "    -lb <int>           : value of the best known lower bound or global optimum (default: " << LOWER_BOUND << ")." << endl;
    cout << endl;
    cout << "Special parameters:" << endl;
    cout << "    --validate-only     : if this parameter is set, the program will only evaluate the solution." << endl;
    cout << endl << flush;
}
