cmake_minimum_required(VERSION 3.8)
project(nurse)

add_definitions(-DILOUSESTL)
add_definitions(-DCPLEX)

set(CMAKE_CXX_STANDARD 11)

# cplex folder and system
set(CPLEX_HOME "/opt/cplex127")
set(CPLEX_SYSTEM "x86-64_osx")
set(GUROBI_HOME "/opt/gurobi751/mac64")

# library directories
include_directories("${CPLEX_HOME}/cplex/include")
include_directories("${CPLEX_HOME}/concert/include")
include_directories("${GUROBI_HOME}/include")
link_directories("${CPLEX_HOME}/cplex/lib/${CPLEX_SYSTEM}/static_pic")
link_directories("${CPLEX_HOME}/concert/lib/${CPLEX_SYSTEM}/static_pic")
link_directories("${GUROBI_HOME}/lib")

set(SOURCE_FILES
    src/decomposition.cpp
    src/decomposition.h
    src/instance.cpp
    src/instance.h
    src/main.cpp
    src/solution.cpp
    src/solution.h)

add_executable(nurse ${SOURCE_FILES})

target_link_libraries(nurse cplex)
target_link_libraries(nurse concert)
target_link_libraries(nurse ilocplex)
target_link_libraries(nurse gurobi_c++)
target_link_libraries(nurse gurobi75)
