---
# NLopt Java Reference
---

The NLopt includes an interface callable from the [Java programming language](https://en.wikipedia.org/wiki/Java_(programming_language)).

The main purpose of this section is to document the syntax and unique features of the Java API; for more detail on the underlying features, please refer to the C documentation in the [NLopt Reference](NLopt_Reference.md).

[TOC]

Using the NLopt Java API
------------------------

To use NLopt in Java, your Java program should usually include the line:

```java
import nlopt.*;
```

which imports the complete `nlopt` package, or alternatively, individual imports of every class you intend to use (which some Java IDEs will maintain automatically for you). Java also allows using the package with explicit namespacing, e.g., `nlopt.Opt opt = new nlopt.Opt(nlopt.Algorithm.LD_MMA, 2);`, but this is typically not recommended.

In addition, your Java program *must* ensure that the JNI library `nloptjni` is loaded before using the `nlopt` package, because Java will otherwise be unable to find the native methods of the NLopt Java binding and throw a runtime error. This can be done with the line:

```java
System.loadLibrary("nloptjni");
```

Simple programs will typically call this at the beginning of `main`. In more complex applications, the most suitable spot is a constructor or a static initializer in the application class interfacing to NLopt. The application may also already have a wrapper around `System.loadLibrary` that sets library search paths and/or extracts libraries from JARs. In that case, that wrapper can also be used for `nloptjni`.

The `nlopt.Opt` class
---------------------

The NLopt API revolves around an object of type `nlopt.Opt`, analogous to (and internally wrapping) `nlopt::opt` in C++. Via methods of this object, all of the parameters of the optimization are specified (dimensions, algorithm, stopping criteria, constraints, objective function, etcetera), and then one finally calls the `Opt.optimize` method in order to perform the optimization. The object should normally be created via the constructor:

```java
Opt opt = new Opt(algorithm, n);
```

given an `algorithm` (see [NLopt Algorithms](NLopt_Algorithms.md) for possible values) and the dimensionality of the problem (`n`, the number of optimization parameters). (Writing just `Opt` assumes that `nlopt.Opt` has been imported with `import nlopt.*;` or `import nlopt.Opt;` in the import section of the class. Otherwise, you have to write `nlopt.Opt` explicitly.) The code snippets below assume that `opt` is an instance of `nlopt.Opt`, constructed as above.

Whereas the C algorithms are specified by `nlopt_algorithm` constants of the form `NLOPT_LD_MMA`, `NLOPT_LN_COBYLA`, etcetera, the Java `Algorithm` values are of the form `nlopt.Algorithm.LD_MMA`, `nlopt.Algorithm.LN_COBYLA`, etcetera (with the `NLOPT_` prefix replaced by the `nlopt.Algorithm.` enum). (Again, the `nlopt.` package can be omitted if `nlopt.Algorithm` has been imported. It is also possible to import the individual enum entries using `import static`, e.g., `import static nlopt.Algorithm.*;` or `import static nlopt.Algorithm.LD_MMA`, allowing to write, e.g., just `LD_MMA` in the code below.)

There are also a copy constructor `nlopt.Opt(Opt)` to make a copy of a given object (equivalent to `nlopt_copy` in the C API).

If there is an error in the constructor (or copy constructor, or assignment), a runtime exception such as `IllegalArgumentException` or `NullPointerException`, or an `OutOfMemoryError` is thrown.

The algorithm and dimension parameters of the object are immutable (cannot be changed without constructing a new object), but you can query them for a given object by the methods:

```java
opt.getAlgorithm()
opt.getDimension()
```

You can get a string description of the algorithm via:

```java
opt.getAlgorithmName()
```

Objective function
------------------

The objective function is specified by calling one of the methods:

```java
opt.setMinObjective(f)
opt.setMaxObjective(f)
```

depending on whether one wishes to minimize or maximize the objective function `f`, respectively. The function `f` must implement the interface `nlopt.Opt.Func`, which specifies a method `double apply(double[] x, double[] gradient)`. This can be done in 3 ways:

1. explicitly, by declaring a named or anonymous class implementing the interface. (In Java versions prior to 1.8, this was the only way.)
2. as an explicit lambda, of the form
   ```java
   (x, grad) -> {
     if (grad != null) {
       ...set grad to gradient, in-place...
     }
     return ...value of f(x)...;
   }
   ```
3. as a static method reference of the form `MyClass::f` where `MyClass` contains a method of the form:
   ```java
   private static double f(double[] x, double[] grad) {
     if (grad != null) {
       ...set grad to gradient, in-place...
     }
     return ...value of f(x)...;
   }
   ```
   Note that, if the reference `MyClass::f` is within the same class `MyClass`, the method `f` can and should be `private`. Otherwise, it needs a higher visibility, e.g., `public`. Also note that `MyClass::f` is just a shortcut for the lambda `(x, grad) -> MyClass.f(x,grad)`, which is itself a shortcut for `(x, grad) -> {return MyClass.f(x,grad);}`.

The return value should be the value of the function at the point `x`, where `x` is a `double[]` array of length `n` of the optimization parameters (the same as the dimension passed to the constructor).

In addition, if the argument `grad` is not null, i.e. `grad != null`, then `grad` is a `double[]` array of length `n` which should (upon return) be set to the gradient of the function with respect to the optimization parameters at `x`. That is, `grad[i]` should upon return contain the partial derivative $\partial f / \partial x_i$, for $0 \leq i < n$, if `grad` is non-null. Not all of the optimization algorithms (below) use the gradient information: for algorithms listed as "derivative-free," the `grad` argument will always be null and need never be computed. (For algorithms that do use gradient information, however, `grad` may still be null for some calls.)

Note that `grad` must be modified *in-place* by your function `f`. Generally, this means using indexing operations `grad[...] = ...;` to overwrite the contents of `grad`, as described below.

### Assigning results in-place

Your objective and constraint functions must overwrite the contents of the `grad` (gradient) argument in-place (although of course you can allocate whatever additional storage you might need, in addition to overwriting `grad`). However, typical Java assignment operations do *not* do this. For example:

```java
grad = Arrays.stream(x).map(t -> 2*t).toArray();
```

might seem like the gradient of the function `sum(x*x)`, but it will *not work* with NLopt because this expression actually allocates a *new* array to store `2*x` and re-assigns `grad` to point to it, rather than overwriting the old contents of `grad`. Instead, you should either explicitly copy the local array to `grad` using `System.arraycopy`:

```java
double[] mygrad = Arrays.stream(x).map(t -> 2*t).toArray();
System.arraycopy(mygrad, 0, grad, 0, grad.length);
```

or set `grad` in place to begin with, e.g.:

```java
Arrays.setAll(grad, i -> 2*x[i]);
```

or simply:

```java
int n = x.length;
for (int i = 0; i < n; i++) {
  grad[i] = 2*x[i];
}
```

Bound constraints
-----------------

The [bound constraints](NLopt_Reference.md#bound-constraints) can be specified by calling the methods:

```java
opt.setLowerBounds(lb);
opt.setUpperBounds(ub);
```

where `lb` and `ub` are `DoubleVector` instances of length *n* (the same as the dimension passed to the `nlopt.Opt` constructor). For convenience, these are overloaded with functions that take a single number as arguments, in order to set the lower/upper bounds for all optimization parameters to a single constant.

To retrieve the values of the lower/upper bounds, you can call one of:

```java
opt.getLowerBounds();
opt.getUpperBounds();
```

both of which return `DoubleVector` instances.

To specify an unbounded dimension, you can use `Double.POSITIVE_INFINITY` or `Double.NEGATIVE_INFINITY` in Java to specify $\pm\infty$, respectively.

Nonlinear constraints
---------------------

Just as for [nonlinear constraints in C](NLopt_Reference.md#nonlinear-constraints), you can specify nonlinear inequality and equality constraints by the methods:

```java
opt.addInequalityConstraint(fc, tol);
opt.addEqualityConstraint(h, tol);
```

where the arguments `fc` and `h` have the same form as the objective function above. The (optional) `tol` arguments specify a tolerance in judging feasibility for the purposes of stopping the optimization, as in C (defaulting to zero if they are omitted).

To remove all of the inequality and/or equality constraints from a given problem, you can call the following methods:

```java
opt.removeInequalityConstraints();
opt.removeEqualityConstraints();
```


### Vector-valued constraints

Just as for [nonlinear constraints in C](NLopt_Reference.md#vector-valued-constraints), you can specify vector-valued nonlinear inequality and equality constraints by the methods

```java
opt.addInequalityMconstraint(c, tol)
opt.addEqualityMconstraint(c, tol)
```

Here, `tol` is a `DoubleVector` of the tolerances in each constraint dimension; the dimensionality *m* of the constraint is determined by `tol.size`. The constraint function `c` must implement the interface `nlopt.Opt.MFunc`, which specifies a method `double[] apply(double[] x, double[] gradient);`. It can be implemented in the same 3 ways as the `nlopt.Opt.Func` interface, but now the lambda must be of the form:

```java
(x, grad) -> {
  if (grad != null) {
    ...set grad to gradient, in-place...
  }
  double[] result = new double[m];
  result[0] = ...value of c_0(x)...
  result[1] = ...value of c_1(x)...
  return result;
}
```

It should return a `double[]` array whose length equals the dimensionality *m* of the constraint (same as the length of `tol` above) and containing the constraint results at the point `x` (a `double[]` array whose length *n* is the same as the dimension passed to the constructor).

In addition, if the argument `grad` is not null, i.e. `grad != null`, then `grad` is a `double[]` array of length `m*n` which should (upon return) be set in-place ([see above](#assigning-results-in-place)) to the Jacobian (i.e., the matrix of gradient rows, in row-major order) of the constraints with respect to the optimization parameters at `x`. That is, `grad[i*n+j]` should upon return contain the partial derivative $\partial c_i / \partial x_j$, for $0 \leq j < n$, if `grad` is non-null. Not all of the optimization algorithms (below) use the gradient information: for algorithms listed as "derivative-free," the `grad` argument will always be null and need never be computed. (For algorithms that do use gradient information, however, `grad` may still be null for some calls.)

An inequality constraint corresponds to $c_i \le 0$ for $0 \le i < m$, and an equality constraint corresponds to $c_i = 0$, in both cases with tolerance `tol[i]` for purposes of termination criteria.

(You can add multiple vector-valued constraints and/or scalar constraints in the same problem.)

Stopping criteria
-----------------

As explained in the [C API Reference](NLopt_Reference.md#stopping-criteria) and the [Introduction](NLopt_Introduction.md#termination-conditions)), you have multiple options for different stopping criteria that you can specify. (Unspecified stopping criteria are disabled; i.e., they have innocuous defaults.)

For each stopping criteria, there are (at least) two methods: a `set` method to specify the stopping criterion, and a `get` method to retrieve the current value for that criterion. The meanings of each criterion are exactly the same as in the C API.

```java
opt.setStopval(stopval);
opt.getStopval();
```

Stop when an objective value of at least `stopval` is found.

```java
opt.setFtolRel(tol);
opt.getFtolRel();
```

Set relative tolerance on function value.

```java
opt.setFtolAbs(tol);
opt.getFtolAbs();
```

Set absolute tolerance on function value.

```java
opt.setXtolRel(tol);
opt.getXtolRel();
```

Set relative tolerance on optimization parameters.

```java
opt.setXtolAbs(tol);
opt.getXtolAbs();
```

Set absolute tolerances on optimization parameters. The `tol` input must be a `DoubleVector` of length `n` (the dimension specified in the `nlopt.Opt` constructor); alternatively, you can pass a single number in order to set the same tolerance for all optimization parameters. `getXtolAbs()` returns the tolerances as a `DoubleVector`.

```java
opt.setXWeights(w);
opt.getXWeights();
```

Set the weights used when the computing L₁ norm for the `xtolRel` stopping criterion above.

```java
opt.setMaxeval(maxeval);
opt.getMaxeval();
```

Stop when the number of function evaluations exceeds `maxeval`. (0 or negative for no limit.)

```java
opt.setMaxtime(maxtime);
opt.getMaxtime();
```

Stop when the optimization time (in seconds) exceeds `maxtime`. (0 or negative for no limit.)

```java
opt.getNumevals();
```

Request the number of evaluations.

### Forced termination

In certain cases, the caller may wish to *force* the optimization to halt, for some reason unknown to NLopt. For example, if the user presses Ctrl-C, or there is an error of some sort in the objective function. You can do this by raising *any* exception inside your objective/constraint functions: the optimization will be halted gracefully, and the same exception will be raised to the caller. See [Exceptions](#exceptions), below. The Java equivalent of `nlopt_forced_stop` from the [C API](NLopt_Reference.md#forced-termination) is to throw an `nlopt.ForcedStopException`.

Algorithm-specific parameters
-----------------------------

Certain NLopt optimization algorithms allow you to specify additional parameters by calling
```java
opt.setParam("name", val);
opt.hasParam("name");
opt.getParam("name", defaultval);
opt.numParams();
opt.nthParam(n);
```
where the string `"name"` is the name of an algorithm-specific parameter and `val` is the value you are setting the parameter to. These functions are equivalent to the [C API](NLopt_Reference.md#algorithm-specific-parameters) functions of the corresponding names.

Performing the optimization
---------------------------

Once all of the desired optimization parameters have been specified in a given object `opt`, you can perform the optimization by calling:

```java
DoubleVector xopt = opt.optimize(x);
```

On input, `x` is a `DoubleVector` of length `n` (the dimension of the problem from the `nlopt.Opt` constructor) giving an initial guess for the optimization parameters. The return value `xopt` is a `DoubleVector` containing the optimized values of the optimization parameters.

You can call the following methods to retrieve the optimized objective function value from the last `optimize` call, and also the return code (including negative/failure return values) from the last `optimize` call:

```java
double opt_val = opt.lastOptimumValue();
Result result = opt.lastOptimizeResult();
```

The return code (see below) is positive on success, indicating the reason for termination. On failure (negative return codes), by default, `optimize()` throws an exception (see [Exceptions](#exceptions), below).

### Return values

The possible return values are the same as the [return values in the C API](NLopt_Reference.md#return-values), except that the `NLOPT_` prefix is replaced with the `nlopt.Result` enum. That is, `NLOPT_SUCCESS` becomes `Result.SUCCESS`, etcetera.

Exceptions
----------

If exceptions are enabled (the default), the [Error codes (negative return values)](NLopt_Reference.md#error-codes-negative-return-values) in the C API are replaced in the Java API by thrown exceptions. The following exceptions are thrown by the various routines:

```
RuntimeException
```

Generic failure, equivalent to `NLOPT_FAILURE`. Note that, since other, more specific exceptions will typically be subclasses of `RuntimeException`, this should be caught *last*.

```
IllegalArgumentException
```

Invalid arguments (e.g. lower bounds are bigger than upper bounds, an unknown algorithm was specified, etcetera), equivalent to `NLOPT_INVALID_ARGS`.

```
OutOfMemoryError
```

Ran out of memory (a memory allocation failed), equivalent to `NLOPT_OUT_OF_MEMORY`.

`nlopt.RoundoffLimitedException` (subclass of `RuntimeException`)
Halted because roundoff errors limited progress, equivalent to `NLOPT_ROUNDOFF_LIMITED`.

`nlopt.ForcedStopException` (subclass of `RuntimeException`)
Halted because of a [forced termination](#forced-termination): the user called `opt.forceStop()` from the user’s objective function or threw an `nlopt.ForcedStop` exception. Equivalent to `NLOPT_FORCED_STOP`.

Whether this behavior is enabled or whether `nlopt.Opt.optimize` just returns the error code as is is controlled by the `enableExceptions` flag in `nlopt.Opt`, which can be set and retrieved with the methods below.

```java
opt.setExceptionsEnabled(enable)
opt.getExceptionsEnabled()
```

The default is `true`, i.e., to throw an exception. When setting `opt.setExceptionsEnabled(false)`, it is the caller's responsibility to *manually* check `opt.lastOptimizeResult()`. While that makes the `false` setting more error-prone, it has the advantage that the best point found (which can be quite good even in some error cases) can still be returned through the return value of `optimize`, so is not lost, whereas if exceptions are enabled through `opt.setExceptionsEnabled(true)`, the exception prevents the best point from being returned.

If your objective/constraint functions throw *any* (runtime) exception during the execution of `opt.optimize`, it will be caught by NLopt and the optimization will be halted gracefully, and `opt.optimize` will re-throw the *same* exception to its caller. (Note that the Java compiler will not allow you to throw a checked exception from your callbacks, only a runtime exception.) For Java, the exception *will always* be rethrown, even if exceptions are otherwise disabled (`opt.setExceptionsEnabled(false)`).

Local/subsidiary optimization algorithm
---------------------------------------

Some of the algorithms, especially MLSL and AUGLAG, use a different optimization algorithm as a subroutine, typically for local optimization. You can change the local search algorithm and its tolerances by calling:

```java
opt.setLocalOptimizer(localOpt);
```

Here, `localOpt` is another `nlopt.Opt` object whose parameters are used to determine the local search algorithm, its stopping criteria, and other algorithm parameters. (However, the objective function, bounds, and nonlinear-constraint parameters of `localOpt` are ignored.) The dimension `n` of `localOpt` must match that of `opt`.

This function makes a copy of the `localOpt` object, so you can freely change your original `localOpt` afterwards without affecting `opt`.

Initial step size
-----------------

Just as in the C API, you can [get and set the initial step sizes](NLopt_Reference.md#initial-step-size) for derivative-free optimization algorithms. The Java equivalents of the C functions are the following methods:

```java
opt.setInitialStep(dx);
DoubleVector dx = opt.getInitialStep(x);
```

Here, `dx` is a `DoubleVector` of the (nonzero) initial steps for each dimension, or a single number if you wish to use the same initial steps for all dimensions. `opt.getInitialStep(x)` returns the initial step that will be used for a starting guess of `x` in `opt.optimize(x)`.

Stochastic population
---------------------

Just as in the C API, you can [get and set the initial population](NLopt_Reference.md#stochastic-population) for stochastic optimization algorithms, by the methods:

```java
opt.setPopulation(pop);
opt.getPopulation();
```

(A `pop` of zero implies that the heuristic default will be used.)

Pseudorandom numbers
--------------------

For stochastic optimization algorithms, we use pseudorandom numbers generated by the [Mersenne Twister](https://en.wikipedia.org/wiki/Mersenne_twister) algorithm, based on code from Makoto Matsumoto. By default, the [seed](https://en.wikipedia.org/wiki/Random_seed) for the random numbers is generated from the system time, so that you will get a different sequence of pseudorandom numbers each time you run your program. If you want to use a "deterministic" sequence of pseudorandom numbers, i.e. the same sequence from run to run, you can set the seed by calling:

```java
NLopt.srand(seed);
```

where `seed` is an integer. To reset the seed based on the system time, you can call:

```java
NLopt.srandTime();
```

(Normally, you don't need to call this as it is called automatically. However, it might be useful if you want to "re-randomize" the pseudorandom numbers after calling `NLopt.srand` to set a deterministic seed.)

Vector storage for limited-memory quasi-Newton algorithms
---------------------------------------------------------

Just as in the C API, you can get and set the [number *M* of stored vectors](NLopt_Reference.md#vector-storage-for-limited-memory-quasi-newton-algorithms) for limited-memory quasi-Newton algorithms, via the methods:

```java
opt.setVectorStorage(M);
opt.getVectorStorage();
```

(The default is *M*=0, in which case NLopt uses a heuristic nonzero value.)

Version number
--------------

To determine the version number of NLopt at runtime, you can call:

```java
int major = NLopt.versionMajor();
int minor = NLopt.versionMinor();
int bugfix = NLopt.versionBugfix();
```

For example, NLopt version 3.1.4 would return `major=3`, `minor=1`, and `bugfix=4`.


