Example in Modern Fortran
-------------------------

Using Modern Fortran, an equivalent of the C example above would be as follows. First, we would write our functions and package them in a module as:

```Fortran
module myopt_mod

    use, intrinsic :: iso_c_binding, only: c_int, c_double
    use nlopt_dt, only: nlopt_user_func ! the abstract type we want to extend
    implicit none
    private

    public :: func, constraint

    type, extends(nlopt_user_func) :: func
        ! no data needed here
    contains
        procedure, public :: eval => eval_func
    end type

    type, extends(nlopt_user_func) :: constraint
        real(c_double) :: a, b
    contains
        procedure, public :: eval => eval_constraint
    end type

contains

    real(c_double) function eval_func(this,n,x,grad)
        class(func), intent(in) :: this
        integer(c_int), intent(in) :: n
        real(c_double), intent(in) :: x(n)
        real(c_double), intent(out), optional :: grad

        if (present(grad)) then
            grad(1) = 0.0_c_double
            grad(2) = 0.5_c_double/sqrt(x(2))
        end if
        eval_func = sqrt(x(2))
    end function

    real(c_double) function eval_constraint(this,n,x,grad)
        class(constraint), intent(in) :: this
        integer(c_int), intent(in) :: n
        real(c_double), intent(in) :: x(n)
        real(c_double), intent(out), optional :: grad

        if (present(grad)) then
            grad(1) = 3._c_double*this%a*(this%a*x(1) + this%b)**2
            grad(2) = -1.0_c_double
        end if
        eval_constraint = (this%a*x(1) + this%b)**3 - x(2)
    end function

end module
```
Notice that the derived-types `func` and `constraint` are now extensions of the abstract type `nlopt_user_func` from the module `nlopt_dt`. The objective function `eval_func` and constraint `eval_constraint` are then used to overload the deferred procedure `eval` as required by the definition of the abstract type. Any data needed by the function can be stored inside the extended type (e. g. `a` and `b` in `constraint`) along with procedures for initialization if desired.

The gradient of the functions are returned through an optional argument `grad`. In a call from C with a null pointer for the optional argument it will be treated as absent. Whether `grad` is present can be checked by using the Fortran intrinsic function `present(a)` that returns the value true if `a` is present and false otherwise. (This is one of the improvements compared to the old Fortran interface where a logical value `need_gradient` had to be passed separately).

Then, to run the optimization, we can use the following Fortran program:

```Fortran
program main
    
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    use nlopt, only: opt, version
    use nlopt_algorithms, only: NLOPT_LD_MMA
    use myopt_mod, only: func, constraint
    implicit none

    integer, parameter :: ip = c_int    ! integer precision
    integer, parameter :: wp = c_double ! working precision

    integer(ip), parameter :: n = 2
    real(wp) :: lb(n), x(n), optf
    type(opt) :: myopt
    type(func) :: myfunc
    type(constraint) :: myc1, myc2
    integer(ip) :: ires 
    
    call version(major,minor,bugfix)
    write(*,"(NLopt version I0.(I0).(I0)") major, minor, bugfix

    myopt = opt(a=NLOPT_LD_MMA,n=n)

    lb = [-huge(lb),0.0_wp]
    call myopt%set_lower_bounds(lb)

    call myopt%set_min_objective(myfunc)

    myc1 = constraint(2.0_wp,0.0_wp)
    call myopt%add_inequality_constraint(myc1,tol=1.e-8_wp)

    myc2 = constraint(-1._wp,1._wp)
    call myopt%add_inequality_constraint(myc2,tol=1.e-8_wp)

    call myopt%set_xtol_rel(tol=1.e-4_wp)

    x = [1.234_wp,5.678_wp]
    call myopt%optimize(x,optf,ires)

    if (ires < 0) then
        write(*,*) 'NLopt failed!'
    else
        write(*,*) 'Found mininum at ', x
        write(*,*) ' val = ', optf
    end if

end program
```
A few things to note here:
* The

### Example using the "raw" interface

The Modern Fortran interface demonstrated above relies on the derived-type `opt` to encapsulate almost all C aspects of NLopt (only function data has to be unpacked in the routines). Using the `nlopt_interface` module however, it is possible to completely mirror a C-style implementation. (In fact, one could even solve an optimization problem with function in C and a driver program written in Fortran.) A demonstration of calling the "raw" C interface from Fortran is given as:


```Fortran
module myopt_mod

    use, intrinsic :: iso_c_binding
    implicit none
    private

    public :: myfunc, myconstraint, mydata

    type :: mydata
        real(c_double) :: a, b
    end type

contains

    real(c_double) function myfunc(n,x,grad,data)
        integer(c_int), intent(in) :: n
        real(c_double), intent(in) :: x(n)
        real(c_double), intent(out), optional :: grad
        type(c_ptr), intent(in), value :: data ! not used in this example
        if (present(grad)) then
            grad(1) = 0.0_c_double
            grad(2) = 0.5_c_double/sqrt(x(2))
        end if
        myfunc = sqrt(x(2))
    end function

    real(c_double) function myconstraint(n,x,grad,c_data)
        integer(c_int), intent(in) :: n
        real(c_double), intent(in) :: x(n)
        real(c_double), intent(out), optional :: grad
        type(c_ptr), intent(in), value :: c_data
        type(mydata), pointer :: data
        call c_f_pointer(c_data,data) ! unpack data from c_ptr
        if (present(grad)) then
            grad(1) = 3._c_double*data%a*(data%a*x(1) + data%b)**2
            grad(2) = -1.0_c_double
        end if
        myconstraint = (data%a*x(1) + data%b)**3 - x(2)
    end function

end module
```

```Fortran
program main

    use, intrinsic :: iso_c_binding
    use nlopt_interface
    use nlopt_algorithms, only: NLOPT_LD_MMA
    use myopt_mod, only: myfunc, myconstraint, mydata
    implicit none

    integer(c_int), parameter :: n = 2
    real(c_double) :: lb(n), x(n), optf
    type(c_ptr) :: opt
    type(c_funptr) :: c_func, c_constraint  ! function and constraint pointers
    type(mydata) :: d1, d2                  ! constraint data

    opt = nlopt_create(NLOPT_LD_MMA,n)

    lb = [-huge(lb),0._c_double]
    ires = nlopt_set_lower_bounds(opt,lb)

    c_func = c_funloc(myfunc)               ! c-pointer to objective function
    ires = nlopt_set_min_objective(opt,c_func,c_null_ptr)

    c_constraint = c_funloc(myconstraint)   ! c-pointer to constraint function
    d1 = data(2._c_double,0._c_double)
    ires = nlopt_add_inequality_constraint(opt,c_constraint,c_loc(d1),1.e-8_c_double)
    d2 = data(-1._c_double,0._c_double)
    ires = nlopt_add_inequality_constraint(opt,c_constraint,c_loc(d2),1.e-8_c_double)

    ires = nlopt_set_xtol_rel(opt,1.e-4_c_double)

    x = [1.234_c_double,5.678_c_double] ! initial guess
    ires = nlopt_optimize(opt,x,optf)   ! optf contains minimum objective value upon return

    if (ires < 0) then
        write(*,*) 'NLopt failed!'
    else
        write(*,*) 'found min at ', x
        write(*,*) 'opt val = ', optf
    end if

    call nlopt_destroy(opt)

end program
```

There are a few things to note here:
* Just like in the C-API most `nlopt_` functions return a integer result flag.
* The "`nlopt_opt`" variable `opt` is now declared as a `type(c_ptr)` (and not an `integer*8` like in the old Fortran interface).
* Data needed by the user functions should be first converted to a `type(c_ptr)` using the intrinsic function `c_loc()`. For functions that do not require any data a `c_null_ptr` can be passed instead.
* The objective function and constraint function are converted to C function pointers using the intrinsic function `c_funloc(f)` that returns a `type(c_funptr)` object.