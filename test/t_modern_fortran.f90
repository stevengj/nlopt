module my_test_problem

    use iso_c_binding
    implicit none
    private

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

    real(c_double) function myconstraint(n,x,gradient,func_data)
        integer(c_int), intent(in), value :: n
        real(c_double), intent(in) :: x(n)
        real(c_double), intent(out), optional :: gradient(n)
        type(c_ptr), value :: func_data

        real(c_double), pointer :: d(:), a, b

        call c_f_pointer(func_data,d,[2]) ! unpack data from C pointer
        a => d(1)
        b => d(2)

        if (present(gradient)) then
            gradient(1) = 3._c_double*a*(a*x(1) + b)**2
            gradient(2) = 0.0_c_double
        end if
        myconstraint = (a*x(1) + b)**3 - x(2)
    end function

end module

module functor_mod

    use iso_c_binding, only: c_int, c_double
    use nlopt, only: nlopt_user_func

    implicit none
    private

    public :: square, cubic

    type, extends(nlopt_user_func) :: square
    contains
        procedure, public :: eval => eval_square
    end type

    type, extends(nlopt_user_func) :: cubic
        real(c_double) :: a, b
        real(c_double) :: leak(10000) = 0
    contains
        procedure, public :: eval => eval_cubic
    end type

contains

    real(c_double) function eval_square(this,n,x,grad)
        class(square), intent(in) :: this
        integer(c_int), intent(in) :: n
        real(c_double), intent(in) :: x(n)
        real(c_double), intent(out), optional :: grad(n)

        if (present(grad)) then
            grad(1) = 0.0_c_double
            grad(2) = 0.5_c_double/sqrt(x(2))
        end if
        eval_square = sqrt(x(2))
    end function

    real(c_double) function eval_cubic(this,n,x,grad)
        class(cubic), intent(in) :: this
        integer(c_int), intent(in) :: n
        real(c_double), intent(in) :: x(n)
        real(c_double), intent(out), optional :: grad(n)

        associate(a => this%a, b=>this%b)
        if (present(grad)) then
            grad(1) = 3._c_double*a*(a*x(1) + b)**2
            grad(2) = 0.0_c_double
        end if
        eval_cubic = (a*x(1) + b)**3 - x(2)
        end associate
    end function
end module functor_mod

program main

    use iso_c_binding
    use my_test_problem, only: myfunc, myconstraint

    use nlopt, only: algorithm_name, version, opt

    use nlopt_c_interface
    
    use nlopt_enum, only : NLOPT_LD_MMA

    use functor_mod, only: square, cubic

    implicit none


    call procedural_example
    call oo_example

contains

    subroutine procedural_example()

        integer(c_int) :: major, minor, bugfix
        integer(c_int), parameter :: n = 2

        type(c_ptr) :: opt, cd1,cd2

        real(c_double), dimension(n) :: x, lb
        real(c_double), target :: d1(2), d2(2)
        integer(c_int) :: ires
        real(c_double) :: optf

        integer :: irepeat

        type(c_funptr) :: c_myfunc, c_myconstraint

        print *, "========= PROCEDURAL EXAMPLE =========="

        call version(major,minor,bugfix)
        print *, "NLopt version ",major,minor,bugfix

        opt = nlopt_create(NLOPT_LD_MMA,n)
        print *, "Algorithm = ", algorithm_name(nlopt_get_algorithm(opt))
        print *, "Dimension = ", nlopt_get_dimension(opt)

        lb = [-huge(1.0_c_double),0.0_c_double]
        ires = nlopt_set_lower_bounds(opt,lb)
        print *, "Lower bounds are = ", lb, " with return status ", ires

        ! Fortran function to C function pointer
        c_myfunc = c_funloc(myfunc)
        ires = nlopt_set_min_objective(opt,c_myfunc,c_null_ptr)

        ! Fortran constraint to C function pointer
        c_myconstraint = c_funloc(myconstraint)


        d1 = [2.0_c_double,0.0_c_double]
        ires = nlopt_add_inequality_constraint(opt,c_myconstraint,c_loc(d1),1.d-8)
        
        d2 = [-1._c_double, 1.0_c_double]
        ires = nlopt_add_inequality_constraint(opt,c_myconstraint,c_loc(d2),1.d-8)

        do irepeat = 1, 80000
            ires = nlopt_remove_inequality_constraints(opt)

            d1 = [2.0_c_double,0.0_c_double]
            ires = nlopt_add_inequality_constraint(opt,c_myconstraint,c_loc(d1),1.d-8)
            
            d2 = [-1._c_double, 1.0_c_double]
            ires = nlopt_add_inequality_constraint(opt,c_myconstraint,c_loc(d2),1.d-8)

        end do
        ! stop ires


        ires = nlopt_set_xtol_rel(opt,1.d-4)
        print*, "set_xtol_rel", ires

        x = [1.234_c_double,5.678_c_double]
        ires = nlopt_optimize(opt,x,optf)
        print *, ires, nlopt_get_errmsg(opt)
        if (ires < 0) then
            print *, "Nlopt failed!"
        else
            print *, "Found minimum at ", x
            print *, "Minimum value = ", optf
        end if

        call nlopt_destroy(opt)
    end subroutine


    subroutine oo_example()

        integer(c_int) :: major, minor, bugfix
        integer(c_int), parameter :: n = 2

        type(opt) :: myopt

        real(c_double), dimension(n) :: x, lb
        real(c_double), target :: d1(2), d2(2)
        integer(c_int) :: ires
        real(c_double) :: optf
        type(square) :: square_func
        type(cubic) :: constraint
        integer :: i

        print *, "========= OO EXAMPLE =========="

        call version(major,minor,bugfix)
        print *, "NLopt version ",major,minor,bugfix

        myopt = opt(a=NLOPT_LD_MMA,n=n)
        print *, "Algorithm = ", myopt%get_algorithm_name()
        print *, "Dimension = ", myopt%get_dimension()

        lb = [-huge(1.0_c_double),0.0_c_double]
        call myopt%set_lower_bounds(lb)
        ! Fortran function to C function pointer
        ! c_func = c_funloc(myfunc)
        ! call myopt%set_min_objective(myfunc,c_null_ptr)

        call myopt%set_min_objective(square_func)


        d1 = [2.0_c_double,0.0_c_double]
        call myopt%add_inequality_constraint(myconstraint,c_loc(d1),tol=1.d-8,ires=ires)
        if (ires < 0) then
            write(*,*) "something went wrong"
            stop myopt%get_errmsg()
        end if

        constraint%a = -1._c_double
        constraint%b = 1.0_c_double
        d2 = [-1._c_double, 1.0_c_double]
        ! call myopt%add_inequality_constraint(constraint,1.d-8,ires)
        call myopt%add_inequality_constraint(myconstraint,c_loc(d2),tol=1.d-8,ires=ires)
        if (ires < 0) then
            print *, ires
            stop myopt%get_errmsg()
        end if

        ! pin down memory leak
        ! do i = 1, 90000
        !     call myopt%remove_inequality_constraints(ires)
        !     ! print *, d1, constraint%a, constraint%b
        !     if (ires < 0) then
        !         print *, ires
        !         stop myopt%get_errmsg()
        !     end if

        !     d1 = [2.0_c_double,0.0_c_double]
        !     call myopt%add_inequality_constraint(myconstraint,c_loc(d1),tol=1.d-8,ires=ires)
        !     if (ires < 0) then
        !         write(*,*) "something went wrong"
        !         stop myopt%get_errmsg()
        !     end if

        !     constraint%a = -1._c_double
        !     constraint%b = 1.0_c_double
        !     ! d2 = [-1._c_double, 1.0_c_double]
        !     call myopt%add_inequality_constraint(constraint,1.d-8,ires)
        !     ! call myopt%add_inequality_constraint(myconstraint,c_loc(d2),tol=1.d-8,ires=ires)
        !     if (ires < 0) then
        !         print *, ires
        !         stop myopt%get_errmsg()
        !     end if
        ! end do

        call myopt%set_xtol_rel(1.d-4)
        
        x = [1.234_c_double,5.678_c_double]
        ires = myopt%optimize(x,optf)

        if (ires < 0) then
            print *, "Nlopt failed!"
        else
            print *, "Found minimum at ", x
            print *, "Minimum value = ", optf
        end if
    end subroutine

end program