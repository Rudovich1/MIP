#include "interface.hpp"
#include "differential_equation.hpp"
#include <string>
#include <iostream>
#include <fstream>

int main(int argc, char* argv[])
{   
    std::ofstream out("test.txt");

    auto problem = Interface::initProblem(
        std::stod(argv[1]), std::stod(argv[2]),
        std::stod(argv[3]), std::stod(argv[4]),
        std::stod(argv[5]), std::stod(argv[6]),
        std::stod(argv[7]), std::stoull(argv[8]));

    auto res1 = RungeKuttaMethod::solve(problem);
    auto res2 = ExactMethod::solve(problem); 

    out << "m = " << argv[1] << '\n'
        << "p = " << argv[2] << '\n'
        << "k = " << argv[3] << '\n'
        << "step = " << argv[4] << '\n'
        << "x_0 = " << argv[5] << '\n'
        << "y_0 = " << argv[6] << '\n'
        << "dy_0 = " << argv[7] << '\n'
        << "num_steps = " << argv[8] << '\n';

    Interface::printRes(out, res2, res1);
}