#include "interface.hpp"
#include "differential_equation.hpp"
#include <string>
#include <iostream>

int main(int argc, char* argv[])
{
    auto problem = Interface::initProblem(
        std::stod(argv[1]), std::stod(argv[2]),
        std::stod(argv[3]), std::stod(argv[4]),
        std::stod(argv[5]), std::stod(argv[6]),
        std::stod(argv[7]), std::stoull(argv[8]));

    auto res1 = RungeKuttaMethod::solve(problem);
    auto res2 = ExactMethod::solve(problem); 

    Interface::printRes(std::cout, res2, res1);
}