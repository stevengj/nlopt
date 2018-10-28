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

    integer, parameter :: NUM_ALGORITHMS = 43


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

    abstract interface
        real(c_double) function nlopt_func(n,x,gradient,func_data) bind(c,name="nlopt_func")
            import c_double, c_int, c_ptr
            integer(c_int), intent(in), value :: n
            real(c_double), intent(in) :: x(n)
            real(c_double), intent(out), optional :: gradient(n)
            type(c_ptr), value :: func_data
        end function
        subroutine nlopt_mfunc(m,result,n,x,gradient,func_data) bind(c,name="nlopt_mfunc")
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
        subroutine nlopt_precond(n,x,v,vpre,data) bind(c,name="nlopt_precond")
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
            type(c_ptr) :: nlo_algorithm_name
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

        integer(c_int) function nlo_optimize(opt,x,opt_f) bind(c,name="nlopt_optimize")
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
            import only : c_int, c_ptr, c_funptr
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

        integer(c_int) function nlo_remove_equality_constraints(opt) bind(c,name="nlopt_remove_equality_constraints")
            import c_int, c_ptr
            type(c_ptr), intent(in) :: opt
        end function
        integer(c_int) function nlo_add_equality_constraint(opt,h,h_data,tol) bind(c,name="nlopt_add_equality_constraint")
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
            real(c_double), intent(in) :: tol(nlo_get_dimension(opt))
        end function
        integer(c_int) function nlopt_get_xtol_abs(opt,tol) bind(c,name="nlopt_get_xtol_abs")
            import c_int, c_ptr, c_double
            type(c_ptr), value :: opt
            real(c_double), intent(out) :: tol(nlo_get_dimension(opt))
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
            real(c_double), intent(in) :: x(nlo_get_dimension(opt))
        end function
        integer(c_int) function nlopt_set_initial_step(opt,dx) bind(c,name="nlopt_set_initial_step")
            import c_int, c_ptr, c_double
            type(c_ptr), value :: opt
            real(c_double), intent(in) :: dx(nlo_get_dimension(opt))
        end function
        integer(c_int) function nlopt_set_initial_step1(opt,dx) bind(c,name="nlopt_set_initial_step1")
            import c_int, c_ptr, c_double
            type(c_ptr), value :: opt
            real(c_double), value :: dx
        end function
        integer(c_int) function nlopt_get_initial_step(opt,x,dx) bind(c,name="nlopt_get_initial_step")
            import c_int, c_ptr, c_double
            type(c_ptr), value :: opt
            real(c_double), intent(in) :: x(nlo_get_dimension(opt))
            real(c_double), intent(out) :: dx(nlo_get_dimension(opt))
        end function      
    end interface

module nlopt

    use nlopt_interfaces, srand => nlopt_srand, srand_time => nlopt_srand_time

    implicit none

    type :: opt
        type(c_ptr) :: o
    contains

    end type

    public :: srand, srand_time, version

contains

    function new_opt(a,n) result(this)
        type(algorithm), intent(in) :: a
        integer, intent(in) :: n

    end function

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
    end interface

end module nlopt