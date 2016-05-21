#include "nlopt.hpp"
#include <iostream>
#include <string>

#include <cassert>

double f(unsigned int n, const double* x, double* grad, void*)
{
    static int it = 0;

    double a = x[0];
    double b = x[1];
    double c = x[2];

    grad[0] = 2*a*b + 2*a*c + b*b + c*c + 1;
    grad[1] = a*a + c*c + 2*b*a + 2*b*c + 1;
    grad[2] = a*a + b*b + 2*c*a + 2*c*b + 1;

    std::cout << it++ << std::endl;

    std::cout << "J = " << a*a*b + a*a*c + b*b*a + b*b*c + c*c*a + c*c*b + a + b + c << std::endl;
    std::cout << "x = " << x[0] << ", " << x[1] << ", " << x[2] << std::endl;
    std::cout << "grad = " << grad[0] << ", " << grad[1] << ", " << grad[2] << std::endl;

    return a*a*b + a*a*c + b*b*a + b*b*c + c*c*a + c*c*b + a + b + c;

}

double g(unsigned int n, const double* x, double* grad, void*)
{
    static int it = 0;

    double a = x[0];
    double b = x[1];
    double c = x[2];

    grad[0] = 2*(a - 3)*std::exp((a - 3)*(a - 3) - b*b + c*c) - 2*(a - 3)*std::exp(-(a - 3)*(a - 3) + b*b - c*c + 3);
    grad[1] = -2*b*std::exp((a - 3)*(a - 3) - b*b + c*c) + 2*b*std::exp(-(a - 3)*(a - 3) + b*b - c*c + 3);
    grad[2] = 2*c*std::exp((a - 3)*(a - 3) - b*b + c*c) - 2*c*std::exp(-(a - 3)*(a - 3) + b*b - c*c + 3);

    std::cout << it++ << std::endl;

    std::cout << "J = " << std::exp((a - 3)*(a - 3) - b*b + c*c) + std::exp(-(a - 3)*(a - 3) + b*b - c*c + 3) << std::endl;
    std::cout << "x = " << x[0] << ", " << x[1] << ", " << x[2] << std::endl;
    std::cout << "grad = " << grad[0] << ", " << grad[1] << ", " << grad[2] << std::endl;

    return std::exp((a - 3)*(a - 3) - b*b + c*c) + std::exp(-(a - 3)*(a - 3) + b*b - c*c + 3);

}

int main(int argc, char *argv[])
{
    if(argc != 2)
     return -1;
  
    std::string algorithm = argv[1];
  
    nlopt::opt optimizator;
    int size = 3;
    
    if(algorithm == "lbfgs")
     optimizator = nlopt::opt(nlopt::LD_LBFGS, size);
    else if (algorithm == "mma")
     optimizator = nlopt::opt(nlopt::LD_MMA, size);
    
    
    //optimizator.set_min_objective(g, NULL);
    optimizator.set_min_objective(f, NULL);
    optimizator.set_xtol_rel(1e-8);
    optimizator.set_ftol_rel(1e-12);
    optimizator.set_maxeval(5);

    std::vector<double> x(size, 10);
    double J;
    optimizator.optimize(x, J);

    std::cout << "x = " << x[0] << ", " << x[1] << ", " << x[2] << std::endl;
    std::cout << "J = " << J << std::endl;
    
   
    return 0;

}
