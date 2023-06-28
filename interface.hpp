#pragma once

#include "differential_equation.hpp"
#include <iostream>
#include <iomanip>

struct Interface
{
    static Problem initProblem(double m, double p, double k, double step, double start_x, double start_y, double start_dy, size_t num_steps);

    static void printRes(std::ostream& out, const std::vector<Iteration>& res1, const std::vector<Iteration>& res2);
};