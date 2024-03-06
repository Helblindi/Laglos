PROGRAM riemann
  USE arbitrary_eos_lambda_module
  IMPLICIT NONE
  INTEGER, PARAMETER:: NUMBER = KIND(1.d0)
  REAL(KIND=NUMBER) :: gamma=1.4d0
  REAL(KIND=NUMBER) :: rhol, el, rhor, er
  REAL(KIND=NUMBER) :: ul, pl, ur, pr
  REAL(KIND=NUMBER) :: lambda_maxl, lambda_maxr, pstar, tol, t1, t2
  INTEGER           :: k, n, nb_case, it, unit=21
  CHARACTER(LEN=3)  :: case_nb
  CHARACTER(LEN=11) :: header
  LOGICAL           :: okay
  LOGICAL(1)        :: no_iter
  
  OPEN(UNIT = unit, FILE = 'data', FORM = 'formatted', STATUS = 'unknown')
  CALL find_string(unit, '===Number of cases', okay)
  IF (.NOT.okay) THEN
     WRITE(*,*) 'string ===Number of cases not found.'
     STOP
  END IF
  READ(unit,*) nb_case
  
  DO it = 1, nb_case
     WRITE(case_nb,'(I3)') it
     header = '===Case '//trim(adjustl(case_nb))
     CALL find_string(unit,header,okay)
     IF (.NOT.okay) THEN
        WRITE(*,*) '===The end.'
        STOP
     END IF
     READ(21,*) rhol, rhor, ul, ur, pl, pr
     READ(21,*) tol
     !===Testing the covolume EOS
     b_covolume = 0.1/MAXVAL([rhol,rhor])
     
     el = gamma_law_internal(rhol,pl,gamma)
     er = gamma_law_internal(rhor,pr,gamma)

     no_iter = .false.

     CALL CPU_TIME(t1)
     DO n = 1, 1 !1000000 
        CALL lambda_arbitrary_eos(rhol,ul,el,pl,rhor,ur,er,pr,tol,no_iter,&
             lambda_maxl,lambda_maxr,pstar,k,b_covolume)
     END DO
     CALL CPU_TIME(t2)
     write(*,*) header
     WRITE(*,'(A,e23.17)') 'CPU ', t2-t1
     WRITE(*,'(2(A,e23.17,x),A,I1)') ' lambda_max=', &
          MAXVAL([abs(lambda_maxl),abs(lambda_maxr)]), 'pstar=', pstar, 'k=', k

     IF (nb_case==1) THEN
        WRITE(*,*) 'gamma', gamma
        WRITE(*,*) 'rhoL', rhol, 'rhostarL', rhostar(pstar,rhol,pl,gamma), &
             'rhostarR', rhostar(pstar,rhor,pr,gamma), 'rhor', rhor
        WRITE(*,*) 'uL', ul, 'ustar', ustar(pstar), 'ur', ur
        WRITE(*,*) 'pL', pl, 'pstar', pstar, 'pr', pr
        WRITE(*,*) 'relative Residual', phi(pstar)/MAXVAL([abs(phi(pl)),abs(phi(pr))]) 
     ELSE
        WRITE(*,*) 'relative Residual', phi(pstar)/MAXVAL([abs(phi(pl)),abs(phi(pr))]) 
     END IF
  END DO
  CLOSE(21)
CONTAINS
  SUBROUTINE find_string(unit, string, okay)
    IMPLICIT NONE
    INTEGER, PARAMETER                 :: long_max=128
    INTEGER,                INTENT(IN) :: unit
    CHARACTER(LEN=*),       INTENT(IN) :: string
    CHARACTER(len=long_max)            :: control
    LOGICAL                            :: okay
    okay = .TRUE.
    REWIND(unit)
    DO WHILE (.TRUE.)
       READ(unit,'(64A)',ERR=11,END=22) control
       IF (trim(adjustl(control))==string) RETURN
    END DO
11  WRITE(*,*) ' Error in find_string'
    STOP
22  okay = .FALSE.
    RETURN
  END SUBROUTINE find_string

  function gamma_law_internal(rho,p,gamma) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND=NUMBER) :: rho, p, gamma
    REAL(KIND=NUMBER) :: vv
    vv = p*(1-b_covolume*rho)/((gamma-1)*rho)
  END function gamma_law_internal


END PROGRAM riemann
