/*
   Multi Dimensional Global Search.

   Author: Steinn Gudmundsson
   Email: steinng@hotmail.com

   This program is supplied without any warranty whatsoever.

   NB The RNGs seed should be initialized using some timer
*/

#include <iostream>
#include <time.h>

#include <iterator>
#include <algorithm>
#include <stack>

#include "stogo_config.h"
#include "global.h"
#include "local.h"

// Timer stuff
time_t   StartTime;

double MacEpsilon ;
int FC=0, GC=0 ;

int stogo_verbose = 0; /* set to nonzero for verbose output */

Global::Global(RTBox D, Pobj o, Pgrad g, GlobalParams P): Domain(D) {

  dim=Domain.GetDim();
  Objective=o;
  Gradient=g;

  // Initialize parameters
  maxtime=P.maxtime;
  maxeval = P.maxeval;
  numeval = 0;
  eps_cl=P.eps_cl; mu=P.mu; rshift=P.rshift;
  det_pnts=P.det_pnts; rnd_pnts=P.rnd_pnts;
  fbound=DBL_MAX;
}

#if 0 // not necessary; default copy is sufficient 
Global& Global::operator=(const Global &G) {
  // Copy the problem info and parameter settings
  Domain=G.Domain; Objective=G.Objective;  Gradient=G.Gradient;
  maxtime=G.maxtime;
  maxeval = G.maxeval;
  numeval = G.numeval;
  eps_cl=G.eps_cl; mu=G.mu; rshift=G.rshift;
  det_pnts=G.det_pnts; rnd_pnts=G.rnd_pnts;
  return *this;
}
#endif

void Global::FillRegular(RTBox SampleBox, RTBox box) {
  // Generation of regular sampling points
  double w;
  int i, flag, dir;
  Trial tmpTrial(dim);
  RVector m(dim), x(dim);

  if (det_pnts>0) {
    // Add midpoint
    box.Midpoint(m) ;
    tmpTrial.xvals=m ; tmpTrial.objval=DBL_MAX ;
    SampleBox.AddTrial(tmpTrial) ;
    // Add the rest
    i=1 ; flag=1 ; dir=0 ;
    x=m ; w=box.Width(dir) ;
    while (i<det_pnts) {
      x(dir)=m(dir)+flag*rshift*w ;
      tmpTrial.xvals=x ; 
      SampleBox.AddTrial(tmpTrial) ;
      flag=-flag;
      if (flag==1 && dir<dim) {
	x(dir)=m(dir) ;
	dir++ ;
	w=box.Width(dir) ;
      }
      i++ ;
    }
  }
}

void Global::FillRandom(RTBox SampleBox, RTBox box) {
  // Generation of stochastic sampling points
  Trial tmpTrial(dim);

  tmpTrial.objval=DBL_MAX;
  for (int i=1 ; i<=rnd_pnts ; i++) {
    for (int dir=0 ; dir<dim ; dir++)
      tmpTrial.xvals(dir) =
	box.lb(dir)+(box.ub(dir)-box.lb(dir))*(double(rand())/RAND_MAX) ;
    SampleBox.AddTrial(tmpTrial) ;
  }
}

double Global::NewtonTest(RTBox box, int axis, RCRVector x_av, int *noutside) {
  // Perform the Newton test

  int info,nout=0;
  Trial tmpTrial(dim);
  TBox SampleBox(dim) ;
  double maxgrad=0 ;

  // Create sampling points
  FillRegular(SampleBox, box);
  FillRandom(SampleBox, box);

  // Perform the actual sampling
  while ( !SampleBox.EmptyBox() ) {
    SampleBox.RemoveTrial(tmpTrial) ;
    info = local(tmpTrial, box, Domain, eps_cl, &maxgrad, *this,
		 axis, x_av) ;
    // What should we do when info=LS_Unstable?
    if (info == LS_Out)
      nout++;
    if (info == LS_New ) {
      box.AddTrial(tmpTrial) ;

      if (tmpTrial.objval<=fbound+mu && tmpTrial.objval<=box.fmin+mu) {
	if (stogo_verbose) {
	  cout << "Found a candidate, x=" << tmpTrial.xvals;
	  cout << " F=" <<tmpTrial.objval << " FC=" << FC << endl;
	}
	SolSet.push_back(tmpTrial);
      }
#ifdef GS_DEBUG
      cout << "Found a stationary point, X= " << tmpTrial.xvals;
      cout <<" objval=" << tmpTrial.objval << endl;
#endif
    }

    if (!InTime())
      break;
  }
  *noutside=nout;
  return maxgrad;
}

void Global::ReduceOrSubdivide(RTBox box, int axis, RCRVector x_av) {
  TBox B1(dim), B2(dim);
  Trial tmpTrial(dim);
  double maxgrad;
  int ns,nout;

  // Monotonicity test has not been implemented yet
  maxgrad=NewtonTest(box, axis, x_av, &nout);
  ns=box.NStationary() ;
  if (ns==0) {
    // All iterates outside
    // NB result=Intersection(B,boundary(Domain))
    Garbage.push(box) ;
  }
  else
    if (ns==1 && nout==0) {
      // All iterates converge to same point
      Garbage.push(box) ;
    }
    else
      if ( (ns>1) && (box.LowerBound(maxgrad)>fbound) ) {
	// Several stationary points found and lower bound > fbound
	Garbage.push(box) ;
      }
      else {
	// Subdivision
	B1.ClearBox() ; B2.ClearBox() ;
	box.split(B1,B2) ;
	CandSet.push(B1) ; CandSet.push(B2) ;
      }

  // Update fbound
  if (box.fmin < fbound) {
    fbound=box.fmin ;
#ifdef GS_DEBUG
    cout <<"*** Improving fbound, fbound=" << fbound << endl;
#endif
  }
}

void Global::Search(int axis, RCRVector x_av){
  Trial tmpTrial(dim) ;
  TBox box(dim), B1(dim), B2(dim);
  RVector m(dim), x(dim);
  int inner_iter, outer_iter;

  MacEpsilon=eps(); // Get machine precision
  if (det_pnts>2*dim+1) {
    det_pnts=2*dim+1;
    if (stogo_verbose)
      cout << "Warning: Reducing det_pnts to " << det_pnts << endl;
  }

  // Initialize timer
  time(&StartTime);

  // Clear priority_queues
  while (!Garbage.empty())
    Garbage.pop();
  while (!CandSet.empty())
    CandSet.pop();

  box=Domain;
  CandSet.push(box);
  int done=0 ; outer_iter=0 ;

  while (!done) {
    outer_iter++ ;

    // Inner loop
    inner_iter=0 ;
    while (!CandSet.empty()) {
      inner_iter++ ;
      // Get best box from Candidate set
      box=CandSet.top() ; CandSet.pop() ;

#ifdef GS_DEBUG
      cout << "Iteration..." << inner_iter << " #CS=" << CandSet.size()+1 ;
      cout << " Processing " << box.NStationary() << " trials in the box " <<box;
#endif
      ReduceOrSubdivide(box, axis, x_av);

      if (!InTime()) {
        done=TRUE;
	if (stogo_verbose)
	  cout << "The program has run out of time or function evaluations\n";
        break;
      }

    } // inner while-loop
    if (stogo_verbose)
      cout << endl << "*** Inner loop completed ***" << endl ;
    
    // Reduce SolSet if necessary
    SolSet.erase(remove_if(SolSet.begin(), SolSet.end(),
			   TrialGT(fbound+mu)),SolSet.end());
    if (InTime()) {
      if (stogo_verbose) {
	cout << "Current set of minimizers (" << SolSet.size() << ")" << endl ;
	DispMinimizers() ;
      }

      while (!Garbage.empty()) {
        box=Garbage.top() ;
        Garbage.pop() ;
        // Split box
        B1.ClearBox() ; B2.ClearBox() ;
        box.split(B1,B2) ;
        // Add boxes to Candidate set
        CandSet.push(B1) ; CandSet.push(B2) ;
      }
    }
  } // Outer while-loop

  if (stogo_verbose) {
    cout << "Number of outer iterations : " << outer_iter << endl;
    cout << "Number of unexplored boxes : " << CandSet.size() << endl;
    cout << "Number of boxes in garbage : " << Garbage.size() << endl;
    cout << "Number of elements in SolSet : " << SolSet.size() << endl;
    cout << "Number of function evaluations : " << FC << endl;
    cout << "Number of gradient evaluations : " << GC << endl;
  }

  if (axis != -1) {
    // Return minimizer when doing the AV method
    tmpTrial=SolSet.back();
    x_av(axis)=tmpTrial.xvals(0);
  }
}

/************* Various utility functions ****************/
double Global::GetTime()
{
 time_t ctime; time(&ctime);
 return difftime(ctime,StartTime);
}

bool Global::InTime()
{
 return (maxtime <= 0.0 || GetTime()<maxtime) && (!maxeval || numeval<maxeval);
}

double Global::GetMinValue() {
  return fbound;
}

void Global::SetMinValue(double new_fb) {
  fbound=new_fb;
}

void Global::SetDomain(RTBox box) {
  Domain=box;
}

void Global::GetDomain(RTBox box) {
  box=Domain;
}

void Global::DispMinimizers() {
  copy(SolSet.begin(), SolSet.end(), ostream_iterator<Trial>(cout));
}

double Global::OneMinimizer(RCRVector x) {
  if (NoMinimizers()) return 0.0;
  for (int i=0;i<x.GetLength();i++) x(i) = SolSet.front().xvals(i);
  return SolSet.front().objval;
}

bool Global::NoMinimizers() {
  return SolSet.empty();
}

void Global::ClearSolSet() {
  SolSet.erase(SolSet.begin(), SolSet.end()) ;
}

void Global::AddPoint(RCRVector x, double f) {
  Trial T(dim);
  T.xvals=x; T.objval=f;
  Domain.AddTrial(T);
  SolSet.push_back(T);
}
