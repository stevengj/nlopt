% Usage: [xopt, fopt, retcode] = nlopt_optimize(opt, xinit)
%
% Optimizes (minimizes or maximizes) a nonlinear function under
% nonlinear constraints from the starting guess xinit, where the
% objective, constraints, stopping criteria, and other options are 
% specified in the structure opt described below.  A variety of local
% and global optimization algorithms can be used, as specified by the 
% opt.algorithm parameter described below.  Returns the optimum
% function value fopt, the location xopt of the optimum, and a
% return code retcode described below (> 0 on success).
%
% This function is a front-end for the external routine nlopt_optimize
% in the free NLopt nonlinear-optimization library, which is a wrapper
% around a number of free/open-source optimization subroutines.  More
% details can be found on the NLopt web page (ab-initio.mit.edu/nlopt)
% and also under 'man nlopt_minimize' on Unix.
