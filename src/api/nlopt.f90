module nlopt_enums

    implicit none
    public

    ! Naming conventions:
    !
    !   NLOPT_{G/L}{D/N}_*
    !   = global/local derivative/no-derivative optimization,
    !   respectively
    !
    !   *_RAND algorithms involve some randomization.

    !   *_NOSCAL algorithms are *not* scaled to a unit hypercube
    !   (i.e. they are sensitive to the units of x)
    
    !
    ! nlopt_algorithm
    ! 
    integer, parameter :: NLOPT_GN_DIRECT = 0
    integer, parameter :: NLOPT_GN_DIRECT_L = 1
    integer, parameter :: NLOPT_GN_DIRECT_L_RAND = 2
    integer, parameter :: NLOPT_GN_DIRECT_NOSCAL = 3
    integer, parameter :: NLOPT_GN_DIRECT_L_NOSCAL = 4
    integer, parameter :: NLOPT_GN_DIRECT_L_RAND_NOSCAL = 5

    integer, parameter :: NLOPT_GN_ORIG_DIRECT = 6
    integer, parameter :: NLOPT_GN_ORIG_DIRECT_L = 7

    integer, parameter :: NLOPT_GD_STOGO = 8
    integer, parameter :: NLOPT_GD_STOGO_RAND = 9

    integer, parameter :: NLOPT_LD_LBFGS_NOCEDAL = 10

    integer, parameter :: NLOPT_LD_LBFGS = 11

    integer, parameter :: NLOPT_LN_PRAXIS = 12

    integer, parameter :: NLOPT_LD_VAR1 = 13
    integer, parameter :: NLOPT_LD_VAR2 = 14

    integer, parameter :: NLOPT_LD_TNEWTON = 15
    integer, parameter :: NLOPT_LD_TNEWTON_RESTART = 16
    integer, parameter :: NLOPT_LD_TNEWTON_PRECOND = 17
    integer, parameter :: NLOPT_LD_TNEWTON_PRECOND_RESTART = 18

    integer, parameter :: NLOPT_GN_CRS2_LM = 19

    integer, parameter :: NLOPT_GN_MLSL = 20
    integer, parameter :: NLOPT_GD_MLSL = 21
    integer, parameter :: NLOPT_GN_MLSL_LDS = 22
    integer, parameter :: NLOPT_GD_MLSL_LDS = 23

    integer, parameter :: NLOPT_LD_MMA = 24

    integer, parameter :: NLOPT_LN_COBYLA = 25

    integer, parameter :: NLOPT_LN_NEWUOA = 26
    integer, parameter :: NLOPT_LN_NEWUOA_BOUND = 27

    integer, parameter :: NLOPT_LN_NELDERMEAD = 28
    integer, parameter :: NLOPT_LN_SBPLX = 29

    integer, parameter :: NLOPT_LN_AUGLAG = 30
    integer, parameter :: NLOPT_LD_AUGLAG = 31
    integer, parameter :: NLOPT_LN_AUGLAG_EQ = 32
    integer, parameter :: NLOPT_LD_AUGLAG_EQ = 33

    integer, parameter :: NLOPT_LN_BOBYQA = 34

    integer, parameter :: NLOPT_GN_ISRES = 35

    ! new variants that require local_optimizer to be set,
    ! not with older constants for backwards compatibility
    integer, parameter :: NLOPT_AUGLAG = 36
    integer, parameter :: NLOPT_AUGLAG_EQ = 37
    integer, parameter :: NLOPT_G_MLSL = 38
    integer, parameter :: NLOPT_G_MLSL_LDS = 39

    integer, parameter :: NLOPT_LD_SLSQP = 40

    integer, parameter :: NLOPT_LD_CCSAQ = 41

    integer, parameter :: NLOPT_GN_ESCH = 42

    integer, parameter :: NUM_ALGORITHMS = 43 ! not an algorithm, just the number of them


    !
    ! nlopt_result
    !
    integer, parameter :: NLOPT_FAILURE = -1             ! generic failure code
    integer, parameter :: NLOPT_INVALID_ARGS = -2
    integer, parameter :: NLOPT_OUT_OF_MEMORY = -3
    integer, parameter :: NLOPT_ROUNDOFF_LIMITED = -4
    integer, parameter :: NLOPT_FORCED_STOP = -5
    integer, parameter :: NLOPT_SUCCESS = 1              ! generic success code
    integer, parameter :: NLOPT_STOPVAL_REACHED = 2
    integer, parameter :: NLOPT_FTOL_REACHED = 3
    integer, parameter :: NLOPT_XTOL_REACHED = 4
    integer, parameter :: NLOPT_MAXEVAL_REACHED = 5
    integer, parameter :: NLOPT_MAXTIME_REACHED = 6

end module

module nlopt_interfaces

    use iso_c_binding
    implicit none

    public

    abstract interface
        real(c_double) function nlopt_func(n,x,gradient,func_data) bind(c)
            use iso_c_binding, only: c_double, c_int, c_ptr
            integer(c_int), intent(in), value :: n
            real(c_double), intent(in) :: x(n)
            real(c_double), intent(out), optional :: gradient(n)
            type(c_ptr), value :: func_data
        end function
        subroutine nlopt_mfunc(m,result,n,x,gradient,func_data) bind(c)
            import c_int, c_double, c_ptr
            integer(c_int), intent(in), value :: m
            real(c_double), intent(out) :: result(m)
            integer(c_int), intent(in), value :: n
            real(c_double), intent(in) :: x(n)
            real(c_double), intent(out), optional :: gradient(m,n)
            type(c_ptr), value :: func_data
        end subroutine
        ! A preconditioner, which preconditions v at x to return vpre.
        ! (The meaning of "preconditioning" is algorithm-dependent.)
        subroutine nlopt_precond(n,x,v,vpre,data) bind(c)
            import c_int, c_double, c_ptr
            integer(c_int) :: n
            real(c_double), intent(in) :: x(n)
            real(c_double), intent(in) :: v(n)
            real(c_double), intent(out) :: vpre(n)
            type(c_ptr), value :: data
        end subroutine
    end interface

    interface
        function nlopt_algorithm_name(algorithm) bind(c,name="nlopt_algorithm_name")
            import c_int, c_ptr
            type(c_ptr) :: nlopt_algorithm_name
            integer(c_int), value :: algorithm
        end function

        subroutine nlopt_srand(seed) bind(c,name="nlopt_srand")
            use iso_c_binding, only : c_int
            integer(c_int), value :: seed
        end subroutine
        subroutine nlopt_srand_time() bind(c, name="nlopt_srand_time")
        end subroutine
        subroutine nlopt_version(major,minor,bugfix) bind(c,name="nlopt_version")
            use iso_c_binding, only : c_int
            integer(c_int), intent(out) :: major
            integer(c_int), intent(out) :: minor
            integer(c_int), intent(out) :: bugfix
        end subroutine

        ! the only immutable parameters of an optimization are the algortihm and
        ! the dimension n of the problem, since changing either of these could
        ! have side effects on lots of other parameters
        type(c_ptr) function nlopt_create(algorithm,n) bind(c,name="nlopt_create")
            import c_ptr, c_int
            integer(c_int), value :: algorithm, n
        end function
        subroutine nlopt_destroy(opt) bind(c,name="nlopt_destroy")
            import c_ptr
            type(c_ptr), value :: opt
        end subroutine
        type(c_ptr) function nlopt_copy(opt) bind(c,name="nlopt_copy")
            import c_ptr
            type(c_ptr), intent(in), value :: opt
        end function

        integer(c_int) function nlopt_optimize(opt,x,opt_f) bind(c,name="nlopt_optimize")
            import c_int, c_ptr, c_double
            type(c_ptr), intent(in), value :: opt
            real(c_double), intent(inout) :: x(nlopt_get_dimension(opt))
            real(c_double), intent(inout) :: opt_f
        end function

        integer(c_int) function nlopt_set_min_objective(opt,f,f_data) bind(c,name="nlopt_set_min_objective")
            import c_int, c_ptr, c_funptr
            type(c_ptr), value :: opt
            type(c_funptr), value :: f
            type(c_ptr), value :: f_data
        end function
        integer(c_int) function nlopt_set_max_objective(opt,f,f_data) bind(c,name="nlopt_set_max_objective")
            import c_int, c_ptr, c_funptr
            type(c_ptr), value :: opt
            type(c_funptr), value :: f
            type(c_ptr), value :: f_data
        end function
        integer(c_int) function nlopt_set_precond_min_objective(opt,f,pre,f_data) bind(c,name="nlopt_set_precond_min_objective")
            import c_int, c_ptr, c_funptr
            type(c_ptr), value :: opt
            type(c_funptr), value :: f
            type(c_funptr), value :: pre
            type(c_ptr), value :: f_data
        end function
        integer(c_int) function nlopt_set_precond_max_objective(opt,f,pre,f_data) bind(c,name="nlopt_set_precond_max_objective")
            import c_int, c_ptr, c_funptr
            type(c_ptr), value :: opt
            type(c_funptr), value :: f
            type(c_funptr), value :: pre
            type(c_ptr), value :: f_data
        end function

        pure integer(c_int) function nlopt_get_algorithm(opt) bind(c,name="nlopt_get_algorithm")
            import c_int, c_ptr
            type(c_ptr), value :: opt
        end function
        pure integer(c_int) function nlopt_get_dimension(opt) bind(c,name="nlopt_get_dimension")
            import c_int, c_ptr
            type(c_ptr), value :: opt
        end function


        ! function nlopt_get_errmsg(opt) bind(c,name="nlopt_get_errmsg")
        !     import c_ptr, c_char
        !     type(c_ptr), value :: opt
        ! end function

        !
        ! constraints
        !
        integer(c_int) function nlopt_set_lower_bounds(opt,lb) bind(c,name="nlopt_set_lower_bounds")
            import c_int, c_ptr, c_double
            type(c_ptr), value :: opt
            real(c_double), intent(in) :: lb(nlopt_get_dimension(opt))
        end function
        integer(c_int) function nlopt_set_lower_bounds1(opt,lb) bind(c,name="nlopt_set_lower_bounds1")
            import c_int, c_ptr, c_double
            type(c_ptr), value :: opt
            real(c_double), intent(in), value :: lb
        end function
        integer(c_int) function nlopt_get_lower_bounds(opt,lb) bind(c,name="nlopt_get_lower_bounds")
            import c_int, c_ptr, c_double
            type(c_ptr), intent(in), value :: opt
            real(c_double), intent(out) :: lb(nlopt_get_dimension(opt))
        end function
        integer(c_int) function nlopt_set_upper_bounds(opt,ub) bind(c,name="nlopt_set_upper_bounds")
            import c_int, c_ptr, c_double
            type(c_ptr), value :: opt
            real(c_double), intent(in) :: ub(nlopt_get_dimension(opt))
        end function
        integer(c_int) function nlopt_set_upper_bounds1(opt,ub) bind(c,name="nlopt_set_lower_bounds1")
            import c_int, c_ptr, c_double
            type(c_ptr), intent(in), value :: opt
            real(c_double), intent(in), value :: ub
        end function
        integer(c_int) function nlopt_get_upper_bounds(opt,ub) bind(c,name="nlopt_get_upper_bounds")
            import c_int, c_ptr, c_double
            type(c_ptr), intent(in), value :: opt
            real(c_double), intent(out) :: ub(nlopt_get_dimension(opt))
        end function


        integer(c_int) function nlopt_remove_inequality_constraints(opt) bind(c,name="nlopt_remove_inequality_constraints")
            import c_int, c_ptr
            type(c_ptr), intent(in) :: opt
        end function
        integer(c_int) function nlopt_add_inequality_constraint(opt,fc,fc_data,tol) bind(c,name="nlopt_add_inequality_constraint")
            import c_int, c_ptr, c_funptr, c_double
            type(c_ptr), value :: opt
            type(c_funptr), value :: fc
            type(c_ptr), value :: fc_data
            real(c_double), value :: tol
        end function
        !integer(c_int) nlopt_add_precond_inequality_constraint
        !integer(c_int) nlopt_add_inequality_mconstraint

        integer(c_int) function nlopt_remove_equality_constraints(opt) bind(c,name="nlopt_remove_equality_constraints")
            import c_int, c_ptr
            type(c_ptr), intent(in) :: opt
        end function
        integer(c_int) function nlopt_add_equality_constraint(opt,h,h_data,tol) bind(c,name="nlopt_add_equality_constraint")
            import c_int, c_ptr, c_funptr, c_double
            type(c_ptr), value :: opt
            type(c_funptr), value :: h
            type(c_ptr), value :: h_data
            real(c_double), value :: tol
        end function
        !integer(c_int) nlo_add_precond_equality_constraint
        !integer(c_int) nlo_add_equality_mconstraint


        !
        ! stopping criteria
        !
        integer(c_int) function nlopt_set_stopval(opt,stopval) bind(c,name="nlopt_set_stopval")
            import c_int, c_double, c_ptr
            type(c_ptr), value :: opt
            real(c_double), value :: stopval
        end function
        real(c_double) function nlopt_get_stopval(opt) bind(c,name="nlopt_get_stopval")
            import c_double, c_ptr
            type(c_ptr), value :: opt
        end function

        integer(c_int) function nlopt_set_ftol_rel(opt,tol) bind(c,name="nlopt_set_ftol_rel")
            import c_int, c_ptr, c_double
            type(c_ptr), value :: opt
            real(c_double), value :: tol
        end function
        real(c_double) function nlopt_get_ftol_rel(opt) bind(c,name="nlopt_get_ftol_rel")
            import c_double, c_ptr
            type(c_ptr), value :: opt
        end function
        integer(c_int) function nlopt_set_ftol_abs(opt,tol) bind(c,name="nlopt_set_ftol_abs")
            import c_int, c_ptr, c_double
            type(c_ptr), value :: opt
            real(c_double), value :: tol
        end function
        real(c_double) function nlopt_get_ftol_abs(opt) bind(c,name="nlopt_get_ftol_abs")
            import c_double, c_ptr
            type(c_ptr), value :: opt
        end function

        integer(c_int) function nlopt_set_xtol_rel(opt,tol) bind(c,name="nlopt_set_xtol_rel")
            import c_int, c_ptr, c_double
            type(c_ptr), value :: opt
            real(c_double), value :: tol
        end function
        real(c_double) function nlopt_get_xtol_rel(opt) bind(c,name="nlopt_get_xtol_rel")
            import c_double, c_ptr
            type(c_ptr), value :: opt
        end function
        integer(c_int) function nlopt_set_xtol_abs1(opt,tol) bind(c,name="nlopt_set_xtol_abs1")
            import c_int, c_ptr, c_double
            type(c_ptr), value :: opt
            real(c_double), value :: tol
        end function
        integer(c_int) function nlopt_set_xtol_abs(opt,tol) bind(c,name="nlopt_set_xtol_abs")
            import c_int, c_ptr, c_double
            type(c_ptr), value :: opt
            real(c_double), intent(in) :: tol(nlopt_get_dimension(opt))
        end function
        integer(c_int) function nlopt_get_xtol_abs(opt,tol) bind(c,name="nlopt_get_xtol_abs")
            import c_int, c_ptr, c_double
            type(c_ptr), value :: opt
            real(c_double), intent(out) :: tol(nlopt_get_dimension(opt))
        end function

        integer(c_int) function nlopt_set_maxeval(opt,maxeval) bind(c,name="nlopt_set_maxeval")
            import c_int, c_ptr
            type(c_ptr), value :: opt
            integer(c_int), value :: maxeval
        end function
        integer(c_int) function nlopt_get_maxeval(opt) bind(c,name="nlopt_get_maxeval")
            import c_int, c_ptr
            type(c_ptr), value :: opt
        end function

        integer(c_int) function nlopt_get_numevals(opt) bind(c,name="nlopt_get_numevals")
            import c_int, c_ptr
            type(c_ptr), value :: opt
        end function

        integer(c_int) function nlopt_set_maxtime(opt,maxtime) bind(c,name="nlopt_set_maxtime")
            import c_int, c_double, c_ptr
            type(c_ptr), value :: opt
            real(c_double), value :: maxtime
        end function
        real(c_double) function nlopt_get_maxtime(opt) bind(c,name="nlopt_get_maxtime")
            import c_double, c_ptr
            type(c_ptr), value :: opt
        end function

        integer(c_int) function nlopt_force_stop(opt) bind(c,name="nlopt_force_stop")
            import c_int, c_ptr
            type(c_ptr), value :: opt
        end function
        integer(c_int) function nlopt_set_force_stop(opt,ival) bind(c,name="nlopt_set_force_stop")
            import c_int, c_ptr
            type(c_ptr), value :: opt
            integer(c_int), value :: ival
        end function
        integer(c_int) function nlopt_get_force_stop(opt) bind(c,name="nlopt_get_force_stop")
            import c_int, c_ptr
            type(c_ptr), value :: opt
        end function


        !
        ! more algorithm-specific parameters
        !
        integer(c_int) function nlopt_set_local_optimizer(opt,local_opt) bind(c,name="nlopt_set_local_optimizer")
            import c_int, c_ptr
            type(c_ptr), value :: opt
            type(c_ptr), value :: local_opt
        end function

        integer(c_int) function nlopt_set_population(opt,pop) bind(c,name="nlopt_set_population")
            import c_int, c_ptr
            type(c_ptr), value :: opt
            integer(c_int), value :: pop
        end function
        integer(c_int) function nlopt_get_population(opt) bind(c,name="nlopt_get_population")
            import c_int, c_ptr
            type(c_ptr), value :: opt
        end function

        integer(c_int) function nlopt_set_vector_storage(opt,dim) bind(c,name="nlopt_set_vector_storage")
            import c_int, c_ptr
            type(c_ptr), value :: opt
            integer(c_int), value :: dim
        end function
        integer(c_int) function nlopt_get_vector_storage(opt) bind(c,name="nlopt_get_vecotr_storage")
            import c_int, c_ptr
            type(c_ptr), value :: opt
        end function

        integer(c_int) function nlopt_set_default_initial_step(opt,x) bind(c,name="nlopt_set_default_initial_step")
            import c_int, c_ptr, c_double
            type(c_ptr), value :: opt
            real(c_double), intent(in) :: x(nlopt_get_dimension(opt))
        end function
        integer(c_int) function nlopt_set_initial_step(opt,dx) bind(c,name="nlopt_set_initial_step")
            import c_int, c_ptr, c_double
            type(c_ptr), value :: opt
            real(c_double), intent(in) :: dx(nlopt_get_dimension(opt))
        end function
        integer(c_int) function nlopt_set_initial_step1(opt,dx) bind(c,name="nlopt_set_initial_step1")
            import c_int, c_ptr, c_double
            type(c_ptr), value :: opt
            real(c_double), value :: dx
        end function
        integer(c_int) function nlopt_get_initial_step(opt,x,dx) bind(c,name="nlopt_get_initial_step")
            import c_int, c_ptr, c_double
            type(c_ptr), value :: opt
            real(c_double), intent(in) :: x(nlopt_get_dimension(opt))
            real(c_double), intent(out) :: dx(nlopt_get_dimension(opt))
        end function
    end interface

end module

module nlopt

    use nlopt_enums
    use nlopt_interfaces, srand => nlopt_srand, srand_time => nlopt_srand_time, &
    func => nlopt_func, mfunc => nlopt_mfunc, version => nlopt_version

    implicit none


    ! abstract interface
    !     real(c_double) function vfunc(x,grad,data)
    !         real(c_double), intent(in) :: x(:)
    !         real(c_double), intent(in) :: grad(:)
    !         type(c_ptr), value :: data
    !     end function
    ! end interface

    public :: opt

    type :: opt
        type(c_ptr), private :: o = c_null_ptr
        integer(c_int), private :: last_result = NLOPT_FAILURE
        real(c_double), private :: last_optf = huge(1.0_c_double)
        integer(c_int), private :: forced_stop_reason = NLOPT_FORCED_STOP
    contains
        ! procedure, private :: mythrow

        procedure, public :: optimize

        procedure, public :: last_optimize_result
        procedure, public :: last_optimum_value

        procedure, public :: get_algorithm
        procedure, public :: get_algorithm_name
        procedure, public :: get_dimension

        procedure, public :: set_min_objective
        procedure, public :: set_max_objective

        procedure, public :: remove_inequality_constraints
        procedure, public :: add_inequality_constraint
        procedure, public :: remove_equality_constraints
        procedure, public :: add_equality_constraint

        procedure, public :: set_lower_bounds
        procedure, public :: get_lower_bounds
        procedure, public :: set_upper_bounds
        procedure, public :: get_upper_bounds
        
        ! stopping criteria
        procedure, public :: set_stopval
        procedure, public :: get_stopval
        procedure, public :: set_ftol_rel
        procedure, public :: get_ftol_rel
        procedure, public :: set_ftol_abs
        procedure, public :: get_ftol_abs
        procedure, public :: set_xtol_rel
        procedure, public :: get_xtol_rel
        procedure, public :: set_xtol_abs
        procedure, public :: get_xtol_abs

        procedure, public :: set_maxeval
        procedure, public :: get_maxeval

        procedure, public :: get_numevals

        procedure, public :: get_maxtime
        procedure, public :: set_maxtime

        procedure, public :: get_force_stop
        procedure, public :: set_force_stop
        procedure, public :: force_stop

        ! algorithm-specific parameters
        procedure, public :: set_local_optimizer
        procedure, public :: set_default_initial_step
        procedure, public :: get_initial_step

        procedure, private :: assign_opt
        generic, public :: assignment(=) => assign_opt

        ! Destructor
        final :: destroy
    end type

    public :: srand, srand_time, version

    ! Constructors
    interface opt_init
        module procedure new_opt
        module procedure copy_opt
    end interface

contains

    type(opt) function new_opt(a,n) result(this)
        integer(c_int), intent(in) :: a
        integer(c_int), intent(in) :: n

        this%o = nlopt_create(a,n)
        this%last_result = NLOPT_FAILURE
        this%last_optf = huge(this%last_optf)
        this%forced_stop_reason = NLOPT_FORCED_STOP
    end function

    type(opt) function copy_opt(f) result(this)
        type(opt), intent(in) :: f

        this%o = nlopt_copy(f%o)
        this%last_result = f%last_result
        this%last_optf = f%last_optf
        this%forced_stop_reason = f%forced_stop_reason
    end function

    subroutine assign_opt(lhs,rhs)
        class(opt), intent(inout) :: lhs
        class(opt), intent(in) :: rhs

        call nlopt_destroy(lhs%o)
        lhs%o = nlopt_copy(rhs%o)
        lhs%last_result = rhs%last_result
        lhs%last_optf = rhs%last_optf
        lhs%forced_stop_reason = rhs%forced_stop_reason
    end subroutine


    subroutine destroy(this)
        type(opt), intent(inout) :: this
        call nlopt_destroy(this%o)
    end subroutine

    integer(c_int) function optimize(this,x,opt_f)
        class(opt), intent(inout) :: this
        real(c_double), intent(inout) :: x(get_dimension(this))
        real(c_double), intent(inout) :: opt_f
        integer(c_int) :: ret
        
        this%forced_stop_reason = NLOPT_FORCED_STOP
        ret = nlopt_optimize(this%o,x,opt_f)
        this%last_result = ret
        this%last_optf = opt_f
        ! if (ret == NLOPT_FORCED_STOP) call mythrow(this%forced_stop_reason)
        ! call mythrow(ret)
        optimize = ret
    end function

    ! last_optimize_result
    integer(c_int) function last_optimize_result(this)
        class(opt), intent(in) :: this
        last_optimize_result = this%last_result
    end function

    real(c_double) function last_optimum_value(this)
        class(opt), intent(in) :: this
        last_optimum_value = this%last_optf
    end function


    ! Accessors
    integer(c_int) function get_algorithm(this)
        class(opt), intent(in) :: this

        get_algorithm = nlopt_get_algorithm(this%o)
    end function

    function get_algorithm_name(this) result(name)
        class(opt), intent(in) :: this
        character(len=256), allocatable :: name
        character(len=256), pointer :: local
        type(c_ptr) :: cptr
        cptr = nlopt_algorithm_name(this%get_algorithm())
        call c_f_pointer(cptr,local)
        name = local
    end function

    pure integer(c_int) function get_dimension(this)
        class(opt), intent(in) :: this
        get_dimension = nlopt_get_dimension(this%o)
    end function

    ! Set the objective function

    subroutine set_min_objective(this,f,f_data)
        class(opt), intent(inout) :: this
        procedure(func) :: f
        type(c_ptr), intent(in) :: f_data
        integer(c_int) :: ret
        type(c_funptr) :: cf
        cf = c_funloc(f)
        ret = nlopt_set_min_objective(this%o,cf,f_data)
        ! call mythrow(ret)
    end subroutine
    subroutine set_max_objective(this,f,f_data)
        class(opt), intent(inout) :: this
        procedure(func) :: f
        type(c_ptr), intent(in) :: f_data
        integer(c_int) :: ret
        type(c_funptr) :: cf
        cf = c_funloc(f)
        ret = nlopt_set_max_objective(this%o,cf,f_data)
        ! call mythrow(ret)
    end subroutine

    ! Nonlinear constraints

    subroutine remove_inequality_constraints(this)
        class(opt), intent(inout) :: this
        integer(c_int) :: ret
        ret = nlopt_remove_equality_constraints(this%o)
        ! call mythrow(ret)
    end subroutine
    subroutine add_inequality_constraint(this,f,f_data,tol)
        class(opt), intent(inout) :: this
        procedure(func) :: f
        type(c_ptr), value :: f_data
        real(c_double), optional :: tol
        real(c_double) :: tol_
        integer(c_int) :: ret
        type(c_funptr) :: cf
        tol_ = 0
        if (present(tol)) tol_ = tol

        cf = c_funloc(f)

        ret = nlopt_add_inequality_constraint(this%o,cf,f_data,tol_)
        ! call mythrow(ret)
    end subroutine

    subroutine remove_equality_constraints(this)
        class(opt), intent(inout) :: this
        integer(c_int) :: ret
        ret = nlopt_remove_inequality_constraints(this%o)
        ! call mythrow(ret)
    end subroutine
    subroutine add_equality_constraint(this,f,f_data,tol)
        class(opt), intent(inout) :: this
        procedure(func) :: f
        type(c_ptr), value :: f_data
        real(c_double), optional :: tol
        real(c_double) :: tol_
        integer(c_int) :: ret

        type(c_funptr) :: cf

        tol_ = 0
        if (present(tol)) tol_ = tol


        cf = c_funloc(f)

        ret = nlopt_add_equality_constraint(this%o,cf,f_data,tol_)
        ! call mythrow(ret)
    end subroutine
    ! subroutine add_equality_mconstraint(this,mf,f_data,tol)
    !     class(opt), intent(inout) :: this
    !     procedure(mfunc) :: mf
    !     type(c_ptr), value :: f_data
    !     real(c_double), optional :: tol()
    !     real(c_double) :: tol_

    !     tol_ = 0
    !     if (present(tol)) tol_ = tol

    !     ret = nlopt_add_equality_mconstraint(this%o,f,f_data,tol_)
    !     ! call mythrow(ret)
    ! end subroutine


    subroutine set_lower_bounds(this,lb)
        class(opt), intent(inout) :: this
        real(c_double), intent(in) :: lb(get_dimension(this))
        integer(c_int) :: ret
        ret = nlopt_set_lower_bounds(this%o,lb)
    end subroutine
    subroutine get_lower_bounds(this,lb)
        class(opt), intent(in) :: this
        real(c_double), intent(out) :: lb(get_dimension(this))
        integer(c_int) :: ret
        ret = nlopt_get_lower_bounds(this%o,lb)
    end subroutine
    subroutine set_upper_bounds(this,ub)
        class(opt), intent(inout) :: this
        real(c_double), intent(in) :: ub(get_dimension(this))
        integer(c_int) :: ret
        ret = nlopt_set_upper_bounds(this%o,ub)
    end subroutine
    subroutine get_upper_bounds(this,ub)
        class(opt), intent(in) :: this
        real(c_double), intent(out) :: ub(get_dimension(this))
        integer(c_int) :: ret
        ret = nlopt_get_upper_bounds(this%o,ub)
    end subroutine


    real(c_double) function get_stopval(this)
        class(opt), intent(in) :: this
        get_stopval = nlopt_get_stopval(this%o)
    end function
    subroutine set_stopval(this,stopval)
        class(opt), intent(inout) :: this
        real(c_double), intent(in) :: stopval
        integer(c_int) :: ret
        ret = nlopt_set_stopval(this%o,stopval)
    end subroutine

    real(c_double) function get_ftol_rel(this)
        class(opt), intent(in) :: this
        get_ftol_rel = nlopt_get_ftol_rel(this%o)
    end function
    subroutine set_ftol_rel(this,tol)
        class(opt), intent(inout) :: this
        real(c_double), intent(in) :: tol
        integer(c_int) :: ret
        ret = nlopt_set_ftol_rel(this%o,tol)
    end subroutine

    real(c_double) function get_ftol_abs(this)
        class(opt), intent(in) :: this
        get_ftol_abs = nlopt_get_ftol_abs(this%o)
    end function
    subroutine set_ftol_abs(this,tol)
        class(opt), intent(inout) :: this
        real(c_double), intent(in) :: tol
        integer(c_int) :: ret
        ret = nlopt_set_ftol_abs(this%o,tol)
    end subroutine

    real(c_double) function get_xtol_rel(this)
        class(opt), intent(in) :: this
        get_xtol_rel = nlopt_get_xtol_rel(this%o)
    end function
    subroutine set_xtol_rel(this,tol)
        class(opt), intent(inout) :: this
        real(c_double), intent(in) :: tol
        integer(c_int) :: ret
        ret = nlopt_set_xtol_rel(this%o,tol)
    end subroutine
    subroutine set_xtol_abs(this,tol)
        class(opt), intent(inout) :: this
        real(c_double), intent(in) :: tol(get_dimension(this))
        integer(c_int) :: ret
        ret = nlopt_set_xtol_abs(this%o,tol)
    end subroutine
    subroutine get_xtol_abs(this,tol)
        class(opt), intent(in) :: this
        real(c_double), intent(out) :: tol(get_dimension(this))
        integer(c_int) :: ret
        ret = nlopt_get_xtol_abs(this%o,tol)
    end subroutine


    integer(c_int) function get_maxeval(this)
        class(opt), intent(in) :: this
        ! if (.not. associated(this%o)) 
        get_maxeval = nlopt_get_maxeval(this%o)
    end function
    subroutine set_maxeval(this,maxeval)
        class(opt), intent(inout) :: this
        integer(c_int), intent(in) :: maxeval
        integer(c_int) :: ret
        ret = nlopt_set_maxeval(this%o,maxeval)
        ! call mythrow(ret)
    end subroutine

    integer(c_int) function get_numevals(this)
        class(opt), intent(in) :: this
        get_numevals = nlopt_get_numevals(this%o)
    end function


    real(c_double) function get_maxtime(this)
        class(opt), intent(in) :: this
        get_maxtime = nlopt_get_maxtime(this%o)
    end function
    subroutine set_maxtime(this,maxtime)
        class(opt), intent(inout) :: this
        real(c_double), intent(in) :: maxtime
        integer(c_int) :: ret
        ret = nlopt_set_maxtime(this%o,maxtime)
    end subroutine

    integer(c_int) function get_force_stop(this)
        class(opt), intent(in) :: this
        ! if (.not. associated(this%o)) 
        get_force_stop = nlopt_get_force_stop(this%o)
    end function
    subroutine set_force_stop(this,ival)
        class(opt), intent(inout) :: this
        integer(c_int), intent(in) :: ival
        integer(c_int) :: ret
        ret = nlopt_set_force_stop(this%o,ival)
        ! call mythrow(ret)
    end subroutine
    subroutine force_stop(this)
        class(opt), intent(inout) :: this
        call this%set_force_stop(1_c_int)
    end subroutine


    subroutine set_local_optimizer(this,lo)
        class(opt), intent(inout) :: this
        class(opt), intent(in) :: lo
        integer(c_int) :: ret
        ret = nlopt_set_local_optimizer(this%o,lo%o)
        ! call mythrow(ret)
    end subroutine

    subroutine set_default_initial_step(this,x)
        class(opt), intent(inout) :: this
        real(c_double), intent(in) :: x(get_dimension(this))
        integer(c_int) :: ret

        ret = nlopt_set_default_initial_step(this%o,x)
        ! call mythrow(ret)
    end subroutine

    subroutine get_initial_step(this,x,dx)
        class(opt), intent(in) :: this    
        real(c_double), intent(in) :: x(get_dimension(this))
        real(c_double), intent(out) :: dx(get_dimension(this))
        integer(c_int) :: ret
        ret = nlopt_get_initial_step(this%o,x,dx)
        ! call mythrow(ret)
    end subroutine

    ! subroutine mythrow(ret)
    !     integer, intent(in) :: ret

    !     select case(ret)
    !     case(NLOPT_FAILURE)
    !     case(NLOPT_OUT_OF_MEMORY)
    !     case(NLOPT_INVALID_ARGS)
    !     case(NLOPT_ROUNDOFF_LIMITED)
    !     case(NLOPT_FORCED_STOP)
    !     case default
    ! end subroutine


    integer function version_major() result(major)
        integer :: minor, bugfix
        call version(major,minor,bugfix)
    end function

    integer function version_minor() result(minor)
        integer :: major, bugfix
        call version(major,minor,bugfix)
    end function

    integer function version_bugfix() result(bugfix)
        integer :: major, minor
        call version(major,minor,bugfix)
    end function

    function algorithm_name(a) result(name)
        integer(c_int), intent(in) :: a
        character(len=256), allocatable :: name
        character(len=256), pointer :: local
        type(c_ptr) :: cptr
        cptr = nlopt_algorithm_name(a)
        call c_f_pointer(cptr,local)
        name = local
    end function

end module nlopt

module my_test_problem

    use iso_c_binding
    implicit none

    public :: myfunc, myconstraint

contains

    real(c_double) function myfunc(n,x,grad,data)
        integer(c_int), intent(in), value :: n
        real(c_double), intent(in) :: x(n)
        real(c_double), intent(out), optional :: grad(n)
        type(c_ptr), value :: data

        if (present(grad)) then
            grad(1) = 0.0_c_double
            grad(2) = 0.5_c_double/sqrt(x(2))
        end if
        myfunc = sqrt(x(2))
    end function

    real(c_double) function myconstraint(n,x,gradient,func_data) bind(c)
        integer(c_int), intent(in), value :: n
        real(c_double), intent(in) :: x(n)
        real(c_double), intent(out), optional :: gradient(n)
        type(c_ptr), value :: func_data

        real(c_double), pointer :: d(:)
        real(c_double) :: a, b

        call c_f_pointer(func_data,d,[2])

        a = d(1)
        b = d(2)

        if (present(gradient)) then
            gradient(1) = 3._c_double*a*(a*x(1) + b)**2
            gradient(2) = 0.0_c_double
        end if
        myconstraint = (a*x(1) + b)**3 - x(2)
    end function

end module

 program main

    use my_test_problem
    use nlopt
    use nlopt_interfaces
    use nlopt_enums, only : NLOPT_LD_MMA

    implicit none

    call procedural_example

    call oo_example

contains

    subroutine procedural_example()

        integer :: major, minor, bugfix
        integer, parameter :: n = 2

        type(c_ptr) :: opt, opt_copy, cd1,cd2

        real(8), dimension(n) :: x, lb
        real(8), target :: d1(2), d2(2)
        integer :: ires
        real(8) :: optf

        type(c_funptr) :: c_func, c_constraint

        print *, "========= PROCEDURAL EXAMPLE =========="

        call nlopt_version(major,minor,bugfix)
        print *, "NLopt version ",major,minor,bugfix

        opt = nlopt_create(NLOPT_LD_MMA,n)
        print *, "Algorithm = ", algorithm_name(nlopt_get_algorithm(opt))
        print *, "Dimension = ", nlopt_get_dimension(opt)

        lb = [-huge(1.0_c_double),0.0_c_double]
        ires = nlopt_set_lower_bounds(opt,lb)
        print *, "Lower bounds are = ", lb, " with return status ", ires

        ! Fortran function to C function pointer
        c_func = c_funloc(myfunc)
        ires = nlopt_set_min_objective(opt,c_func,c_null_ptr)


        ! Fortran constraint to C function pointer
        c_constraint = c_funloc(myconstraint)

        d1 = [2.0_c_double,0.0_c_double]
        cd1 = c_loc(d1)
        ires = nlopt_add_inequality_constraint(opt,c_constraint,cd1,1.d-8)
        
        d2 = [-1._c_double, 1.0_c_double]
        cd2 = c_loc(d2)
        ires = nlopt_add_inequality_constraint(opt,c_constraint,cd2,1.d-8)

        ires = nlopt_set_xtol_rel(opt,1.d-4)


        x = [1.234_c_double,5.678_c_double]
        ires = nlopt_optimize(opt,x,optf)

        if (ires < 0) then
            print *, "Nlopt failed!"
        else
            print *, "Found minimum at ", x
            print *, "Minimum value = ", optf
        end if

        ! opt_copy = nlopt_copy(opt)
        call nlopt_destroy(opt)

        ! print *, opt, opt_copy
    end subroutine


    subroutine oo_example()

        integer(c_int) :: major, minor, bugfix
        integer(c_int), parameter :: n = 2

        type(opt) :: myopt
        type(c_ptr) :: cd1, cd2

        real(c_double), dimension(n) :: x, lb
        real(c_double), target :: d1(2), d2(2)
        integer(c_int) :: ires
        real(c_double) :: optf

        type(c_funptr) :: c_func, c_constraint

        print *, "========= OO EXAMPLE =========="

        call version(major,minor,bugfix)
        print *, "NLopt version ",major,minor,bugfix

        myopt = new_opt(NLOPT_LD_MMA,n)
        print *, "Algorithm = ", myopt%get_algorithm_name()
        print *, "Dimension = ", myopt%get_dimension()

        lb = [-huge(1.0_c_double),0.0_c_double]
        call myopt%set_lower_bounds(lb)

        ! Fortran function to C function pointer
        ! c_func = c_funloc(myfunc)
        call myopt%set_min_objective(myfunc,c_null_ptr)

        d1 = [2.0_c_double,0.0_c_double]
        cd1 = c_loc(d1)
        call myopt%add_inequality_constraint(myconstraint,cd1,1.d-8)
        
        d2 = [-1._c_double, 1.0_c_double]
        cd2 = c_loc(d2)
        call myopt%add_inequality_constraint(myconstraint,cd2,1.d-8)

        call myopt%set_xtol_rel(1.d-6)
        
        x = [1.234_c_double,5.678_c_double]
        ires = myopt%optimize(x,optf)

        if (myopt%last_optimize_result() < 0) then
            print *, "Nlopt failed!"
        else
            print *, "Found minimum at ", x
            print *, "Minimum value = ", optf
        end if

        ! automatic finalization
    end subroutine
end program