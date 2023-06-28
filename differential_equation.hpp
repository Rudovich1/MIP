#pragma once

#include <vector>
#include <functional>
#include <initializer_list>
#include <cmath>

struct Iteration
{
    double x_;
    double y_;
    double dy_;

    Iteration(){}
    Iteration(double x, double y, double dy): x_(x), y_(y), dy_(dy){}

    void operator=(std::initializer_list<double> lst);
};

struct Problem
{
    double m_;
    double p_;
    double k_;
    double step_;
    Iteration start_;
    size_t num_steps_;

    Problem(double m, double p, double k, double step, double start_x, double start_y, double start_dy, size_t num_steps):
        m_(m), p_(p), k_(k), step_(step), start_(start_x, start_y, start_dy), num_steps_(num_steps){} 
};

struct RungeKuttaMethod
{
    static std::vector<Iteration> solve(Problem problem);
};

struct ExactMethod
{
    static std::vector<Iteration> solve(Problem problem);
};