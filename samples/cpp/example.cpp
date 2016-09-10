#include <iostream>
#include <vector>
#include <cmath>
#include "nlopt.hpp"
using std::vector;
using std::cout;
using std::cerr;
using std::endl;

class ObjectiveFunc
{
public:
    ObjectiveFunc() : count(0) {}
    double operator()(const vector<double> &x, vector<double> &grad)
    {
        ++count;
        if (!grad.empty()) {
            grad[0] = 0.0;
            grad[1] = 0.5 / sqrt(x[1]);
        }
        return sqrt(x[1]);
    }
    static double wrap(const vector<double> &x, vector<double> &grad, void *data)
    {
        return (*reinterpret_cast<ObjectiveFunc*>(data))(x, grad);
    }
public:
    int count;
};

class ConstraintFunc
{
public:
    ConstraintFunc(double _a, double _b) : count(0), a(_a), b(_b) {}
    double operator()(const vector<double> &x, vector<double> &grad)
    {
        if (!grad.empty()) {
            grad[0] = 3 * a * (a*x[0] + b) * (a*x[0] + b);
            grad[1] = -1.0;
        }
        return ((a*x[0] + b) * (a*x[0] + b) * (a*x[0] + b) - x[1]);
    }
    static double wrap(const vector<double> &x, vector<double> &grad, void *data)
    {
        return (*reinterpret_cast<ConstraintFunc*>(data))(x, grad);
    }
public:
    int count;
    double a, b;
};

int main()
{
    // algorithm and dimensionality
    nlopt::opt opt(nlopt::LD_MMA, 2);

    // lower bounds
    vector<double> lb(2);
    lb[0] = -HUGE_VAL;
    lb[1] = 0;
    opt.set_lower_bounds(lb);

    ObjectiveFunc myfunc;
    opt.set_min_objective(ObjectiveFunc::wrap, &myfunc);
    opt.set_xtol_rel(1e-4);

    ConstraintFunc myconst1(2,0), myconst2(-1,1);
    opt.add_inequality_constraint(ConstraintFunc::wrap, &myconst1, 1e-8);
    opt.add_inequality_constraint(ConstraintFunc::wrap, &myconst2, 1e-8);

    // some initial guess
    vector<double> x(2);
    x[0] = 1.234;
    x[1] = 5.678;

    double minf;  // the minimum objective value, upon return
    nlopt::result result = opt.optimize(x, minf);
    if (result < 0)
        cerr << "nlopt failed" << endl;
    else {
        cout << "found minimum after "
            << myfunc.count << " evaluations" << endl;
        cout << "found minimum at f("
            << x[0] << "," << x[1] << ") = " << minf << endl;
    }

    return 0;
}
