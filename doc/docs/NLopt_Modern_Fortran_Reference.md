# NLopt Modern Fortran Reference

NLopt is written in C and the C NLopt programming interface (API), as described in the [NLopt Reference](NLopt_Reference.md)., is interoperable with Modern Fortran by using appropriate interfaces.

However, we also provide another Modern Fortran (from here on just Fortran) module in file nlopt.f90, that wraps a more natural Fortran interface around the NLopt API, which may be more convenient for Fortran programmers used to using derived types.

The main distinctions of the Modern Fortran API are:

* Use of an `nlopt` module.
* Use of a Fortran derived-type `opt`, with constructors, destructors, and member functions.
* An optional exception tracker.

The main purpose of this section is to document the syntax and unique features of the Modern Fortran API; for more details on the underlying features, please refer to the C documention in the [NLopt Reference](NLopt_Reference.md).

## Compiling and linking your program to NLopt

An NLopt program in Fortran should first import the `nlopt` module:

```Fortran
use nlopt
```

On Unix you would normally link your program exactly as for the C API, with a command like:

*`compiler`*` `*`...source/object` `files...`*` -lnlopt -lm -o myprogram`

where *compiler* is `gfortran` or whatever is appropriate for your machine/system.

A distinction in Fortran is that the path to the .mod files should also be passed.

## The `opt` Fortran derived-type

The Fortran API revolves around an object with the derived-type `opt`. Via member functions of this object, all of the parameters of the optimization are specified.

The algorithm and dimension parameters of the object are immutable (cannot be changed without constructing a new object), but you can query them for a given object by the subroutines:

```
call myopt%get_algorithm(algorithm)
call myopt%get_dimension(n)
```

You can get a Fortran-style string description of the algorithm via the function:
```Fortran
myopt%get_algorithm_name()
```
that returns a 

When you are finished with your object, automatic deallocation will occur through a so-called finalizer, that calls the `nlopt_destroy` function.


## Objective function

The objective function is specified by calling one of the methods:

```Fortran

```

depending on whether one wishes to minimize or maximize the objective function `f`, respectively. The function `f` should be of the form:

function f(n,x,grad,f_data)
    integer, intent(in) :: n
    real(wp), intent(in) :: x(n)
    real(wp), intent(out), optional :: grad(n)
    type(c_ptr), value :: f_data
end function

The return value should be the value of the function at the point `x`, where `x` is a vector of length `n` of the optimization parameters (the same dimension as passed to the constructor).

In addition, if the argument `grad` is present, then `grad` is a vector of length `n` which should (upon return) be set to the gradient of the function with respect to the optimization parameters at `x`. That is, `grad(i)` should upon return contain the partial derivative $\partial f / \partial x_i$, for $0 \leq i < n$, if `grad` is present. Not all of the optimization algorithms (below) use the gradient information: for algorithms listed as "derivative-free," the `grad` argument will always be empty and need never be computed. (For algorithms that do use gradient information, however, `grad` may still be empty for some calls.)

## Bound constraints

The [bound constraints](NLopt_Reference#Bound_constraints.md) can be specified by calling the methods:

```Fortran
real(wp) :: lb(n), ub(n)

call myopt%set_lower_bounds(lb [, ires])
call myopt%set_upper_bounds(ub [, ires])
```
where `lb` and `ub` are arrays of length `n` (the same as the dimensions passed to the constructor `opt` constructor), and `ires` is an optional integer [return value]() (positive on success). For convenience, these are overloaded with subroutines that take a single number as arguments, in order to set the lower/upper bounds for all optimization parameters to a single constant:

```Fortran
real(wp) :: lb, ub

call myopt%set_lower_bounds(lb [, ires])
call myopt%set_upper_bounds(ub [, ires])
```

To retrieve the values of the lower/upper bounds, you can call one of:

```Fortran
call myopt%get_lower_bounds(lb [, ires])
call myopt%get_upper_bounds(ub [, ires])
```

To specify an unbounded dimension, you can use +-`huge(lb)` in Fortran to specify $\pm \infty$, where `huge` is a Fortran intrinsic function.

## Nonlinear constraints

Just as for [nonlinear constraints in C](NLopt_Reference#Nonlinear_constraints.md), you can specify nonlinear inequality and equality constraints by the methods:

```Fortran

call myopt%add_inequality_constraint(fc,tol)
call myopt%add_equality_constraint(h,tol)
```

where the arguments `fc` and `h` have the same from as the objective function above.


To remove all of the inequality and/or equality constraints from a given problem, you can call the following methods:

```Fortran
call myopt%remove_inequality_constraints([ires])
call myopt%remove_equality_constraints([ires])
```

### Vector-valued constraints

Just as for [nonlinear constraints in C](NLopt_Reference#Vector-valued_constraints.md), you can specify nonlinear inequality and equality constraints by the methods:

```Fortran
call myopt%add_inequality_mconstraint(m)
call myopt%add_equality_mconstraint(m)
```

Here, `m` is the dimensionality of the constraint results and `tol` is a length-`m` array of the tolerances in each constraint dimension. The constraint subroutine `c` must be of the form:

```Fortran
subroutine c(m,result,n,x,grad,f_data)
    
end subroutine
```

(You can add multiple vector-valued constraints and/or scalar constraints in the same problem.)

Stopping criteria
-----------------

As explained in the [C API Reference](NLopt_Reference#Stopping_criteria.md) and the [Introduction](NLopt_Introduction#Termination_conditions.md)), you have multiple options for different stopping criteria that you can specify. (Unspecified stopping criteria are disabled; i.e., they have innocuous defaults.)

For each stopping criteria, there are (at least) two method: a `set` method to specify the stopping criterion, and a `get` method to retrieve the current value for that criterion. The meanings of each criterion are exactly the same as in the C API.

```Fortran
call myopt%set_stopval()
call myopt%get_stopval(stopval)
```
Stop when an objective value of at least `stopval` is found.

```Fortran
```
Set relative tolerance on function value.

```Fortran
```
Set absolute tolerance on function value.

```Fortran
```
Set relative tolerance on optimization parameters.

```Fortran
```
Set absolute tolerances on optimization parameters. The `tol` input must be an array of length `n` (the dimension specified in the `opt` constructor). This subroutine with a version that accepts a single number `tol1` in order to set the same tolerance for all optimization parameters.

```Fortran
```
Stop when the number of function evaluations exceeds the integer `maxeval`. (Zero or negative for no limit).

```Fortran
call myopt%get_numevals(nevals)
```
Request the number of evaluations.

```Fortran
```
Stop when the optimization time (in seconds) exceeds `maxtime` (double precision). (Zero or negative for no limit.)

### Forced termination



Performing the optimization
---------------------------

Once all the desired optimization parameters have been specified in a given object `myopt`, you can perform the optimization by calling:

```Fortran
real(wp) :: x(n), optf
call myopt%optimize(x, optf [,ires])
```
On input, `x` is an array of length `n` (the dimension of the problem from the `opt` constructor) giving an initial guess for the optimization parameters. Upon succesful return, `x` contains the optimized values of the optimization parameters and `optf` contains the corresponding value of the objective function.

The return value (see below) is positive on success, indicating the reason for termination. On failure (negative return codes), it throws an exception (see [Exceptions](#Exceptions.md), below).

You can also call the following methods to retrieve the `optf` value from the last `optimize` call, and the return value (including negative/failure return values) from the last `optimize` call:

```Fortran
real(wp) :: optf
integer :: ires

optf = myopt%last_optimum_value()
ires = myopt%last_optimize_result()
```

### Return values

The possible return values are the same as the [return values in the C API](NLopt_Reference#Return_values.md), except that the `NLOPT_` prefix is replaced with the `nlopt::` namespace. That is, `NLOPT_SUCCESS` becomes `nlopt::SUCCESS`, etcetera.

## Exceptions

Local/subsidiary optimization algorithm
---------------------------------------

Some of the algorithms, especially MLSL and AUGLAG, use a different optimization algorithm as a subroutine, typically for local optimization. You can change the local search algorithm and its tolerances by calling:

```
type(opt) :: local_opt

! ... set local_opt settings ...

call myopt%set_local_optimizer(local_opt)
```

Here, `local_opt` is another object of `type(opt)` whose parameters are used to determine the local search algorithm, its stopping criteria, and other algorithm parameters. (However, the objective function, bounds, and nonlinear-constraint parameters of `local_opt` are ignored.) The dimension `n` of `local_opt` must match that of `myopt`.

This function makes a copy of the `local_opt` object, so you can freely destroy your original `local_opt` afterwards.


Initial step size
-----------------

Just as in the C API, you can [get and set the initial step sizes](NLopt_Reference#Initial_step_size.md) for derivative-free optimization algorithms. The Modern Fortran equivalents of the C functions are the following methods:

```Fortran
real(wp) :: x(n), dx(n), dx1
call myopt%set_initial_step(dx)
call myopt%set_initial_step(dx1)
call myopt%get_initial_step(x,dx)
```

Here, `dx` is an array of the (nonzero) initial steps for each dimension. For convenience, you can also pass a single number `dx1` to `set_initial_step` if you wish to use the same initial steps for all dimensions. A call to `get_initial_step()` sets `dx` to the initial step that will be used for a starting guess of `x` in a call to `myopt%optimize(x,optf)`.

Stochastic population
---------------------

Just as in the C API, you can [get and set the initial population](NLopt_Reference#Stochastic_population.md) for stochastic optimization algorithms, by the methods:

```
call myopt%set_population(ipop [,ires])
call myopt%get_population(ipop)
```
where `ipop` is an integer and `ires` is an optional integer return value (positive on success).

(An `ipop` of zero implies that the heuristic default will be used.)

Pseudorandom numbers
--------------------

For stochastic optimization algorithms, we use pseudorandom numbers generated by the [Mersenne Twister](https://en.wikipedia.org/wiki/Mersenne_twister) algorithm, based on code from Makoto Matsumoto. By default, the [seed](https://en.wikipedia.org/wiki/Random_seed) for the random numbers is generated from the system time, so that you will get a different sequence of pseudorandom numbers each time you run your program. If you want to use a "deterministic" sequence of pseudorandom numbers, i.e. the same sequence from run to run, you can set the seed by calling:

```Fortran
integer :: iseed
call srand(iseed)
```
To reset the seed based on the system time, you can call:

```Fortran
call srand_time()
```


(Normally, you don't need to call this as it is called automatically. However, it might be useful if you want to "re-randomize" the pseudorandom numbers after calling `nlopt::srand` to set a deterministic seed.)

Vector storage for limited-memory quasi-Newton algorithms
---------------------------------------------------------

Just as in the C API, you can get and set the [number *M* of stored vectors](NLopt_Reference#Vector_storage_for_limited-memory_quasi-Newton_algorithms.md) for limited-memory quasi-Newton algorithms, via the methods:

```Fortran
call myopt%set_vector_storage(M [,ires])
call myopt%get_vector_storage(M)
```
(The default is *M*=0, in which case NLopt uses a heuristic nonzero value.)

Version number
--------------

To determine the version number of NLopt at runtime, you can call:

```Fortran
integer :: major, minor, bugfix
call version(major,minor,bugfix)
```

For example, NLopt version 3.1.4 would return `major=3`, `minor=1`, and `bugfix=4`. You can also retrieve these three values individually by calling:

```Fortran
major = version_major()
minor = version_minor()
bugfix = version_bugfix()
```


[Category:NLopt](index.md)
