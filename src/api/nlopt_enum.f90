module nlopt_enum

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

end module nlopt_enum