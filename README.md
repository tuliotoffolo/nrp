# NRP Solver

## Nurse Rostering Problem Solver

Written by Túlio Toffolo (instance reader by Haroldo Santos).  
Note that the current code-base employs the formulation from <a href="doi.org/10.1007/s10479-014-1594-6" target="_blank">Santos, Toffolo, Gomes and Ribas (2016)</a>.

More information: <a href="https://benchmark.gent.cs.kuleuven.be/nrp" target="_blank">https://benchmark.gent.cs.kuleuven.be/nrp</a>

Please address all contributions, suggestions, and inquiries to the project administrator.

## Latest improvements

The current version includes generates subproblems by employing three decomposition approaches, based on the *Time*, *Nurse* and *Shift* indexes of decision variables **x**. The subproblems are solved in a random order or following certain priorities. This resulted in quicker convergence than previous algorithms, in addition to improvements over some best known solutions.

## Getting Started

The source code includes the decomposition-based heuristic for the NRP.  
Note that specific input files are required to execute the solver (available at <a href="https://benchmark.gent.cs.kuleuven.be/nrp" target="_blank">https://benchmark.gent.cs.kuleuven.be/nrp</a>).

Usage examples:

- ./nurse -dat models/long01.dat -lp models/long01.lp.gz -ini initial_solutions/long01.sol -out output/long01  
- ./nurse -dat models/long01.dat -lp models/long01.lp.gz -ini initial_solutions/long01.sol -out output/long01 -lb 197

## Requirements

CPLEX 12.7 or Gurobi 7.5

## Arguments

```
Usage: ./nurse <arguments> [options]

Program arguments:
    -dat <file.dat> : nurse problem
    -lp <file.lp>   : nurse problem (in cplex lp format)
    -ini <file.sol> : provide an initial solution
    -out <file.txt> : prefix of files where the solutions will be saved

Optional parameters (example):
    -rseed <int>        : random seed (default: 9)
    -mip_emphasis <int> : mip emphasis (only cplex; default: 2).
                             0   Balance optimality and feasibility
                             1   Emphasize feasibility over optimality
                             2   Emphasize optimality over feasibility
                             3   Emphasize moving best bound
                             4   Emphasize hidden feasible solutions
    -mip_focus <int>    : mip focus (only gurobi; default: 1).
                             0   Balance optimality and feasibility
                             1   Focus on feasibility over optimality
                             2   Focus on optimality over feasibility
                             3   Focus on improving the best bound
    -threads <int>      : number of thread the solver may use (default: 0).
    -time_limit <int>   : runtime limit (default: unlimited).
    -ltime_limit <int>  : solver runtime limit to solve each subproblem  (default: unlimited).
    -iter_limit <int>   : maximum number of iterations (default: unlimited).
    -n_days <int>       : initial number of days in a "time-based" sub-problem (default: 4).
    -days_step <int>    : value of the increment on n_days for the "time-based" decomposition. (default: 3).
    -max_ndays <int>    : maximum value for parameter ndays (default: 30).
    -n_nurses <int>     : initial number of nurses in a "nurse-based" sub-problem. (default: 2).
    -nurses_step <int>  : value of the increment on n_nurses for the "nurse-based" decomposition (default: 4).
    -max_nnurses <int>  : maximum value for parameter n_nurses (default: 5).".
    -lb <int>           : value of the best known lower bound or global optimum (default: 0).

Special parameters:
    --validate-only     : if this parameter is set, the program will only evaluate the solution.

```

## Questions?

If you have any questions, please feel free to contact me.

Thanks!
