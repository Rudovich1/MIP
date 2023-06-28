#include "interface.hpp"

Problem Interface::initProblem(double m, double p, double k, double step, double start_x, double start_y, double start_dy, size_t num_steps)
{
    return Problem(m, p, k, step, start_x, start_y, start_dy, num_steps);
}

void Interface::printRes(std::ostream& out, const std::vector<Iteration>& res1, const std::vector<Iteration>& res2)
{
    out << std::setw(12) << std::left << "step" << " | " 
        << std::setw(12) << std::left << "x" << " | "
        << std::setw(12) << std::left << "y_exact" << " | "
        << std::setw(12) << std::left << "y_r_k" << " | "
        << std::setw(12) << std::left << "dy_exact" << " | "
        << std::setw(12) << std::left << "dy_r_k" << " | \n";

    out << std::setprecision(8) << std::fixed;

    for (size_t i = 0; i < res1.size(); ++i)
    {
        out << std::setw(12) << std::left << i << " | " 
            << std::setw(12) << std::left << res1[i].x_ << " | "
            << std::setw(12) << std::left << res1[i].y_ << " | "
            << std::setw(12) << std::left << res2[i].y_ << " | "
            << std::setw(12) << std::left << res1[i].dy_ << " | "
            << std::setw(12) << std::left << res2[i].dy_ << " | \n";
    }
}