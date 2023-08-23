PROGRAM riemann
   USE arbitrary_eos_lagrangian_lambda_module
   IMPLICIT NONE
   INTEGER, PARAMETER           :: NUMBER = KIND(1.d0)
   REAL(KIND=NUMBER), PARAMETER :: gamma_ideal=1.4d0
   REAL(KIND=NUMBER), PARAMETER :: gamma_vdw=1.02d0, a_vdw=1.d0, b_vdw=1.d0
   REAL(KIND=NUMBER)            :: rhol, el, rhor, er
   REAL(KIND=NUMBER)            :: taul, taur
   REAL(KIND=NUMBER)            :: ul, pl, ur, pr
   REAL(KIND=NUMBER)            :: lambda_maxl, lambda_maxr, pstar, tol, t1, t2
   INTEGER                      :: k, n, num_cases, it, unit = 21
   CHARACTER(LEN=3)             :: case_number
   CHARACTER(LEN=11)            :: header, cmd_line_arg
   LOGICAL                      :: OKAY
   LOGICAL                      :: WANT_ITER = .TRUE.

   !===Handle command line argument
   IF (COMMAND_ARGUMENT_COUNT() .GT. 0) THEN
      CALL GET_COMMAND_ARGUMENT(1, cmd_line_arg)

      IF (TRIM(cmd_line_arg) .EQ. 'true') THEN
         WANT_ITER = .TRUE.
      ELSE IF (TRIM(cmd_line_arg) .EQ. 'false') THEN
         WANT_ITER = .FALSE.
      ELSE
         WRITE (*, *) "Invalid command line argument!"
         WRITE (*, *) "Please use 'true', 'false', or leave blank (default is true)."
         STOP
      END IF
   END IF

   !===Read data file
   OPEN (UNIT=unit, FILE='data', FORM='formatted', STATUS='unknown')
   CALL find_string(unit, '===Number of cases', OKAY)

   !===Terminate if unable to find the string
   IF (.NOT. OKAY) THEN
      WRITE (*, *) 'string "===Number of cases" not found.'
      STOP
   END IF

   !===Read in the number of test cases
   READ (unit, *) num_cases

   WRITE(*,*) "Starting main loop"
   !===Main loop for computing the test problems
   DO it = 1, num_cases
      WRITE(*,*) "Iteration: ", it
      WRITE (case_number, '(I3)') it
      header = '===Case '//TRIM(ADJUSTL(case_number))
      CALL find_string(unit, header, OKAY)
      IF (.NOT. OKAY) THEN
         WRITE (*, *) '===The end.'
         STOP
      END IF
      READ (21, *) rhol, rhor, ul, ur, pl, pr
      READ (21, *) tol

      IF (it < 15) THEN 
            el = gamma_law_internal(rhol, pl)
            er = gamma_law_internal(rhor, pr)
      ELSE
            el = van_der_waals_internal(rhol, pl)
            er = van_der_waals_internal(rhor, pr)
      END IF 

      taul = 1.d0/rhol
      taur = 1.d0/rhor

      CALL CPU_TIME(t1)
      DO n = 1, 1 !1000000
         CALL lambda_arbitrary_eos(taul, ul, el, pl, taur, ur, er, pr, tol, WANT_ITER, &
                                   lambda_maxl, lambda_maxr, pstar, k)
      END DO
      CALL CPU_TIME(t2)
      WRITE (*, *) header
      WRITE(*,*) "rhol: ", rhol, ", rhor: ", rhor
      !WRITE(*,'(A,e23.17)') 'CPU ', t2-t1
      WRITE (*, '(2(A,e23.17,x),A,I1)') ' lambda_max=', &
         MAX(ABS(lambda_maxl), ABS(lambda_maxr)), 'pstar=', pstar, 'k=', k

      IF (num_cases == 1) THEN
         WRITE (*, *) 'Warning. These results only hold for data generated from the ideal gas law.'
         WRITE (*, *) 'gamma', gamma_ideal
         WRITE (*, *) 'rhoL', rhol, 'rhostarL', rhostar(pstar, rhol, pl, gamma_ideal), &
            'rhostarR', rhostar(pstar, rhor, pr, gamma_ideal), 'rhor', rhor
         WRITE (*, *) 'uL', ul, 'ustar', ustar(pstar), 'ur', ur
         WRITE (*, *) 'pL', pl, 'pstar', pstar, 'pr', pr
         WRITE (*, *) 'relative Residual', phi(pstar)/MAX(ABS(phi(pl)), ABS(phi(pr)))
      ELSE
         WRITE (*, *) 'relative Residual', phi(pstar)/MAX(ABS(phi(pl)), ABS(phi(pr)))
      END IF
   END DO
   CLOSE (21)
CONTAINS
   !===
   SUBROUTINE find_string(unit_file, string, IS_OKAY)
      IMPLICIT NONE
      INTEGER, PARAMETER                 :: long_max = 128
      INTEGER, INTENT(IN) :: unit_file
      CHARACTER(LEN=*), INTENT(IN) :: string
      CHARACTER(len=long_max)            :: control
      LOGICAL                            :: IS_OKAY
      IS_OKAY = .TRUE.
      REWIND (unit_file)
      DO WHILE (.TRUE.)
         READ (unit_file, '(64A)', ERR=11, END=22) control
         IF (TRIM(ADJUSTL(control)) == string) RETURN
      END DO
11    WRITE (*, *) ' Error in find_string'
      STOP
22    IS_OKAY = .FALSE.
      RETURN
   END SUBROUTINE find_string

   !===Compute the specific internal energy from the covolume EOS
   function gamma_law_internal(rho, p) RESULT(e)
      IMPLICIT NONE
      REAL(KIND=NUMBER) :: rho, p
      REAL(KIND=NUMBER) :: e
      e = p*(1.d0 - b_covolume*rho)/((gamma_ideal - 1.d0)*rho)
   END function gamma_law_internal
   
   !===Compute the specific internal energy from the van der Waals EOS
   function van_der_waals_internal(rho, p) RESULT(e)
      IMPLICIT NONE
      REAL(KIND=NUMBER) :: rho, p
      REAL(KIND=NUMBER) :: e
      e = (p + a_vdw*rho*rho)*(1.d0 - b_vdw*rho)/((gamma_vdw - 1.d0)*rho) - a_vdw*rho
   END function van_der_waals_internal

END PROGRAM riemann
