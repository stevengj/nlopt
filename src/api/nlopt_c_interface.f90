module nlopt_c_interface

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
            integer(c_int), intent(in), value :: n
            real(c_double), intent(in) :: x(n)
            real(c_double), intent(in) :: v(n)
            real(c_double), intent(out) :: vpre(n)
            type(c_ptr), value :: data
        end subroutine

    end interface


    interface

        type(c_ptr) function nlopt_algorithm_name(algorithm) bind(c,name="nlopt_algorithm_name")
            import c_int, c_ptr
            integer(c_int), value :: algorithm
        end function

        subroutine nlopt_srand(seed) bind(c,name="nlopt_srand")
            import c_int
            integer(c_int), value :: seed
        end subroutine
        subroutine nlopt_srand_time() bind(c, name="nlopt_srand_time")
        end subroutine

        subroutine nlopt_version(major,minor,bugfix) bind(c,name="nlopt_version")
            import c_int
            integer(c_int), intent(out) :: major, minor, bugfix
        end subroutine

        ! the only immutable parameters of an optimization are the algorithm and
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
            type(c_ptr),  value :: opt
            real(c_double), intent(inout) :: x(nlopt_get_dimension(opt))
            real(c_double), intent(inout) :: opt_f
        end function

        integer(c_int) function nlopt_set_min_objective(opt,f,f_data) bind(c,name="nlopt_set_min_objective")
            import c_int, c_ptr, c_funptr
            type(c_ptr), value :: opt
            type(c_funptr), intent(in), value :: f
            type(c_ptr), intent(in), value :: f_data
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


        integer(c_int) function nlopt_get_algorithm(opt) bind(c,name="nlopt_get_algorithm")
            import c_int, c_ptr
            type(c_ptr), intent(in), value :: opt
        end function

        ! This function must be pure, so that we can access the dimension in other procedures.
        pure integer(c_int) function nlopt_get_dimension(opt) bind(c,name="nlopt_get_dimension")
            import c_int, c_ptr
            type(c_ptr), intent(in), value :: opt
        end function

        type(c_ptr) function nlopt_get_errmsg(opt) bind(c,name="nlopt_get_errmsg")
            import c_ptr
            type(c_ptr), value :: opt
        end function

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
        integer(c_int) function nlopt_set_upper_bounds1(opt,ub) bind(c,name="nlopt_set_upper_bounds1")
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
            type(c_ptr), value :: opt
        end function
        integer(c_int) function nlopt_add_inequality_constraint(opt,fc,fc_data,tol) bind(c,name="nlopt_add_inequality_constraint")
            import c_int, c_ptr, c_funptr, c_double
            type(c_ptr), value :: opt
            type(c_funptr), value :: fc
            type(c_ptr), value :: fc_data
            real(c_double), value :: tol
        end function
        integer(c_int) function nlopt_add_precond_inequality_constraint(opt,fc,pre,fc_data,tol) &
                   bind(c,name="nlopt_add_precond_inequality_constraint")
            import c_int, c_ptr, c_funptr, c_double
            type(c_ptr), value :: opt
            type(c_funptr), value :: fc
            type(c_funptr), value :: pre
            type(c_ptr), value :: fc_data
            real(c_double), value :: tol
        end function
        integer(c_int) function nlopt_add_inequality_mconstraint(opt,m,fc,fc_data,tol) &
                   bind(c,name="nlopt_add_inequality_mconstraint")
            import c_int, c_ptr, c_funptr, c_double
            type(c_ptr), value :: opt
            integer(c_int), value :: m
            type(c_funptr), value :: fc
            type(c_ptr), value :: fc_data
            real(c_double), intent(in) :: tol(nlopt_get_dimension(opt))
        end function

        integer(c_int) function nlopt_remove_equality_constraints(opt) bind(c,name="nlopt_remove_equality_constraints")
            import c_int, c_ptr
            type(c_ptr), value :: opt
        end function
        integer(c_int) function nlopt_add_equality_constraint(opt,h,h_data,tol) bind(c,name="nlopt_add_equality_constraint")
            import c_int, c_ptr, c_funptr, c_double
            type(c_ptr), value :: opt
            type(c_funptr), value :: h
            type(c_ptr), value :: h_data
            real(c_double), value :: tol
        end function
        integer(c_int) function nlopt_add_precond_equality_constraint(opt,h,pre,h_data,tol) &
                   bind(c,name="nlopt_add_precond_equality_constraint")
            import c_int, c_ptr, c_funptr, c_double
            type(c_ptr), value :: opt
            type(c_funptr), value :: h
            type(c_funptr), value :: pre
            type(c_ptr), value :: h_data
            real(c_double), value :: tol
        end function
        integer(c_int) function nlopt_add_equality_mconstraint(opt,m,h,h_data,tol) &
                   bind(c,name="nlopt_add_equality_mconstraint")
            import c_int, c_ptr, c_funptr, c_double
            type(c_ptr), value :: opt
            integer(c_int), value :: m
            type(c_funptr), value :: h
            type(c_ptr), value :: h_data
            real(c_double), intent(in) :: tol(nlopt_get_dimension(opt))
        end function

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
            type(c_ptr), intent(in), value :: opt
        end function

        integer(c_int) function nlopt_set_ftol_rel(opt,tol) bind(c,name="nlopt_set_ftol_rel")
            import c_int, c_ptr, c_double
            type(c_ptr), value :: opt
            real(c_double), value :: tol
        end function
        real(c_double) function nlopt_get_ftol_rel(opt) bind(c,name="nlopt_get_ftol_rel")
            import c_double, c_ptr
            type(c_ptr), intent(in), value :: opt
        end function
        integer(c_int) function nlopt_set_ftol_abs(opt,tol) bind(c,name="nlopt_set_ftol_abs")
            import c_int, c_ptr, c_double
            type(c_ptr), value :: opt
            real(c_double), value :: tol
        end function
        real(c_double) function nlopt_get_ftol_abs(opt) bind(c,name="nlopt_get_ftol_abs")
            import c_double, c_ptr
            type(c_ptr), intent(in), value :: opt
        end function

        integer(c_int) function nlopt_set_xtol_rel(opt,tol) bind(c,name="nlopt_set_xtol_rel")
            import c_int, c_ptr, c_double
            type(c_ptr), value :: opt
            real(c_double), value :: tol
        end function
        real(c_double) function nlopt_get_xtol_rel(opt) bind(c,name="nlopt_get_xtol_rel")
            import c_double, c_ptr
            type(c_ptr), intent(in), value :: opt
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
            type(c_ptr), intent(in), value :: opt
            real(c_double), intent(out) :: tol(nlopt_get_dimension(opt))
        end function

        integer(c_int) function nlopt_set_maxeval(opt,maxeval) bind(c,name="nlopt_set_maxeval")
            import c_int, c_ptr
            type(c_ptr), value :: opt
            integer(c_int), value :: maxeval
        end function
        integer(c_int) function nlopt_get_maxeval(opt) bind(c,name="nlopt_get_maxeval")
            import c_int, c_ptr
            type(c_ptr), intent(in), value :: opt
        end function

        integer(c_int) function nlopt_get_numevals(opt) bind(c,name="nlopt_get_numevals")
            import c_int, c_ptr
            type(c_ptr), intent(in), value :: opt
        end function

        integer(c_int) function nlopt_set_maxtime(opt,maxtime) bind(c,name="nlopt_set_maxtime")
            import c_int, c_double, c_ptr
            type(c_ptr), value :: opt
            real(c_double), value :: maxtime
        end function
        real(c_double) function nlopt_get_maxtime(opt) bind(c,name="nlopt_get_maxtime")
            import c_double, c_ptr
            type(c_ptr), intent(in), value :: opt
        end function

        integer(c_int) function nlopt_force_stop(opt) bind(c,name="nlopt_force_stop")
            import c_int, c_ptr
            type(c_ptr), value :: opt
        end function
        integer(c_int) function nlopt_set_force_stop(opt,val) bind(c,name="nlopt_set_force_stop")
            import c_int, c_ptr
            type(c_ptr), value :: opt
            integer(c_int), value :: val
        end function
        integer(c_int) function nlopt_get_force_stop(opt) bind(c,name="nlopt_get_force_stop")
            import c_int, c_ptr
            type(c_ptr), intent(in), value :: opt
        end function

        !
        ! more algorithm-specific parameters
        !

        integer(c_int) function nlopt_set_local_optimizer(opt,local_opt) bind(c,name="nlopt_set_local_optimizer")
            import c_int, c_ptr
            type(c_ptr), value :: opt
            type(c_ptr), intent(in), value :: local_opt
        end function

        integer(c_int) function nlopt_set_population(opt,pop) bind(c,name="nlopt_set_population")
            import c_int, c_ptr
            type(c_ptr), value :: opt
            integer(c_int), value :: pop
        end function
        integer(c_int) function nlopt_get_population(opt) bind(c,name="nlopt_get_population")
            import c_int, c_ptr
            type(c_ptr), intent(in), value :: opt
        end function

        integer(c_int) function nlopt_set_vector_storage(opt,dim) bind(c,name="nlopt_set_vector_storage")
            import c_int, c_ptr
            type(c_ptr), value :: opt
            integer(c_int), value :: dim
        end function
        integer(c_int) function nlopt_get_vector_storage(opt) bind(c,name="nlopt_get_vector_storage")
            import c_int, c_ptr
            type(c_ptr), intent(in), value :: opt
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
            type(c_ptr), intent(in), value :: opt
            real(c_double), intent(in) :: x(nlopt_get_dimension(opt))
            real(c_double), intent(out) :: dx(nlopt_get_dimension(opt))
        end function

    end interface

end module nlopt_c_interface