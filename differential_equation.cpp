#include "differential_equation.hpp"

void Iteration::operator=(std::initializer_list<double> lst)
{
    x_ = *lst.begin();
    y_ = *(lst.begin() + 1);
    dy_ = *(lst.begin() + 2);
}

std::vector<Iteration> RungeKuttaMethod::solve(Problem problem)
{
    double m = problem.m_;
    double p = problem.p_;
    double k = problem.k_;
    double h = problem.step_;

    std::function<double(Iteration)> fun1 = [](Iteration iteration){return iteration.dy_;};
    std::function<double(Iteration)> fun2 = [m, p, k](Iteration iteration){return 1/m * (-p * iteration.dy_ - k * iteration.y_);};
    std::function<double(double, double, double, double)> del_fun = [](double A1, double A2, double A3, double A4){return 1./6. * (A1 + 2*A2 + 2*A3 + A4);};

    double K1, K2, K3, K4;
    double L1, L2, L3, L4;
    std::vector<Iteration> res(problem.num_steps_ + 1);
    res[0] = problem.start_;
    Iteration temp_iter = problem.start_;

    for (size_t i = 0; i < problem.num_steps_; ++i)
    {
        temp_iter = res[i];
        
        K1 = h * fun1(temp_iter);
        L1 = h * fun2(temp_iter);
        
        temp_iter.x_ = res[i].x_ + h/2;
        temp_iter.y_ = res[i].y_ + K1/2;
        temp_iter.dy_ = res[i].dy_ + L1/2;

        K2 = h * fun1(temp_iter);
        L2 = h * fun2(temp_iter);

        temp_iter.y_ = res[i].y_ + K2/2;
        temp_iter.dy_ = res[i].dy_ + L2/2;

        K3 = h * fun1(temp_iter);
        L3 = h * fun2(temp_iter);

        temp_iter.x_ = res[i].x_ + h;
        temp_iter.y_ = res[i].y_ + K3;
        temp_iter.dy_ = res[i].dy_ + L3;

        K4 = h * fun1(temp_iter);
        L4 = h * fun2(temp_iter);

        res[i + 1] = {res[i].x_ + h, res[i].y_ + del_fun(K1, K2, K3, K4), res[i].dy_ + del_fun(L1, L2, L3, L4)};
    }

    return res;
}


std::vector<Iteration> ExactMethod::solve(Problem problem)
{
    double eps = 1e-18;

    double m = problem.m_;
    double p = problem.p_;
    double k = problem.k_;
    double x = problem.start_.x_;
    double y = problem.start_.y_;
    double dy = problem.start_.dy_;
    std::function<double(double)> fun1, fun2;

    double D = p*p-4.*m*k;

    if (std::abs(D) < eps)
    {
        double lambda = -p / (2*m);
        double C1 = (y * (1. - lambda * x) - dy * x) / (std::exp(x * lambda));
        double C2 = (dy - y * lambda) / (std::exp(x * lambda));
        fun1 = [C1, C2, lambda](double x){return C1 * std::exp(x * lambda) + C2 * x * std::exp(x * lambda);};
        fun2 = [C1, C2, lambda](double x){return std::exp(x * lambda) * (C1 * lambda + C2 * (1 + x * lambda));};
    }
    else if (D > 0)
    {
        double lambda1 = (-p + std::sqrt(D)) / (2*m);
        double lambda2 = (-p - std::sqrt(D)) / (2*m);
        double C1 = (dy - lambda2 * y) / (std::exp(x * lambda1) * (lambda1 - lambda2));
        double C2 = (dy - lambda1 * y) / (std::exp(x * lambda2) * (lambda2 - lambda1));
        fun1 = [C1, C2, lambda1, lambda2](double x){return C1 * std::exp(x * lambda1) + C2 * std::exp(x * lambda2);}; 
        fun2 = [C1, C2, lambda1, lambda2](double x){return C1 * lambda1 * std::exp(x * lambda1) + C2 * lambda2 * std::exp(x * lambda2);};
    }
    else
    {
        double alpha = -p / (2 * m);
        double beta1 = std::sqrt(-D) / (2 * m);
        double beta2 = -std::sqrt(-D) / (2 * m);
        double A = std::exp(x * alpha) * std::sin(beta1 * x);
        double B = std::exp(x * alpha) * std::cos(beta2 * x);
        double _A = std::exp(x * alpha) * (alpha * std::sin(beta1 * x) + beta1 * std::cos(beta1 * x));
        double _B = std::exp(x * alpha) * (alpha * std::cos(beta2 * x) - beta2 * std::sin(beta2 * x));
        double C1 = (B * dy - _B * y) / (B * _A - _B * A);
        double C2 = (A * dy - _A * y) / (A * _B - _A * B);
        fun1 = [C1, C2, alpha, beta1, beta2](double x){return C1 * std::exp(x * alpha) * std::sin(beta1 * x) + C2 * std::exp(x * alpha) * std::cos(beta2 * x);};
        fun2 = [C1, C2, alpha, beta1, beta2](double x){
            return C1 * std::exp(x * alpha) * (alpha * std::sin(beta1 * x) + beta1 * std::cos(beta1 * x)) 
                 + C2 * std::exp(x * alpha) * (alpha * std::cos(beta2 * x) - beta2 * std::sin(beta2 * x));};
    }
    
    std::vector<Iteration> res(problem.num_steps_ + 1);
    res[0] = problem.start_;

    for (size_t i = 1; i < problem.num_steps_ + 1; ++i)
    {
        x += problem.step_;
        res[i] = {x, fun1(x), fun2(x)};
    }

    return res;
}