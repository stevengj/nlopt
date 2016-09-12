#include <iostream>
using namespace std;

#include "nlopt-util.h"
#include "global.h"
#include "tools.h"
#include "linalg.h"
#include "testfun.h"

int main()
{
  int dim = 2;

  Pdom Dom = Domain_Rosenbrock;
  Pobj Obj = Objective_Rosenbrock;
  Pgrad Grad = Gradient_Rosenbrock;

  nlopt_stopping stop = {0};
  stop.n = dim;
  stop.maxtime = 1;

  GlobalParams params;
  params.stop = &stop;
  params.det_pnts = 2*dim+1; params.rnd_pnts = 0;
  params.eps_cl = 0.1; params.rshift = 0.3;
  params.mu = 1.0E-4;

  TBox D(dim);
  Dom(D);
  Global Problem(D, Obj, Grad, params);
  RVector dummyvec(dim);
  Problem.Search(-1, dummyvec);

  cout << "Minimizers found\n";
  Problem.DispMinimizers();

  return 0;
}
