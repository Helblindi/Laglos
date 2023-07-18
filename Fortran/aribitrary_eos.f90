!===Authors: Bennett Clayton, Jean-Luc Guermond, and Bojan Popov, Texas A&M, April 5, 2021
!===Revised September 1, 2022 for Lagrangian coortdinates
MODULE arbitrary_eos_Lagrangian_lambda_module
  IMPLICIT NONE
  PUBLIC               :: Lagrangian_lambda_arbitrary_eos !===Main function
  PUBLIC               :: rhostar, ustar, phi  !===Optional functions. Can be removed
  PRIVATE
  INTEGER, PARAMETER:: NUMBER = KIND(1.d0)
  REAL(KIND=NUMBER), PARAMETER :: zero = 0
  REAL(KIND=NUMBER), PARAMETER :: one = 1
  REAL(KIND=NUMBER), PARAMETER :: half = 0.5d0
  REAL(KIND=NUMBER), PARAMETER :: five_third = 5.d0/3.d0
  REAL(KIND=NUMBER) :: rhol, ul, pl, el
  REAL(KIND=NUMBER) :: rhor, ur, pr, er
  REAL(KIND=NUMBER) :: gammal, al, alphal, capAl, capBl, capCl, expol
  REAL(KIND=NUMBER) :: gammar, ar, alphar, capAr, capBr, capCr, expor
  REAL(KIND=NUMBER) :: p_min, rho_min, gamma_min, alpha_min, capA_min, capB_min, phi_pmin
  REAL(KIND=NUMBER) :: p_max, rho_max, gamma_max, alpha_max, capC_max, expo_max, phi_pmax
  REAL(KIND=NUMBER) :: gamma_lm, expo_lm
  REAL(KIND=NUMBER) :: gamma_uM, expo_uM
  REAL(KIND=NUMBER) :: numerator, vacuum
  CHARACTER(LEN=1)  :: gamma_min_index, gamma_lm_index

CONTAINS
  !===Inputs: in_rhol,in_ul,in_el,in_pl,in_rhor,in_ur,in_er,in_pr,in_tol,no_iter
  !===Outputs: lambda_maxl_out,lambda_maxr_out,pstar,vstar,k
  !===use no_iter=.true. if one does not want iterative line search. Then the parameter ``in_tol'' is useless.
  !===use no_iter=.false. to enable the line search with ``in_tol'' tolerance. Then ``k'' is the numer of iterations used.
  SUBROUTINE Lagrangian_lambda_arbitrary_eos(in_rhol,in_ul,in_el,in_pl,in_rhor,in_ur,in_er,in_pr,in_tol,no_iter,&
       lambda_maxl_out,lambda_maxr_out,pstar,vstar,k, b_covolume)
    IMPLICIT NONE
    REAL(KIND=8), INTENT(IN) :: in_rhol, in_el, in_rhor, in_er, in_tol
    REAL(KIND=8), INTENT(IN), TARGET :: in_ul, in_pl, in_ur, in_pr
    LOGICAL,      INTENT(IN) :: no_iter
    REAL(KIND=8), INTENT(IN) :: b_covolume
    REAL(KIND=8), INTENT(OUT):: lambda_maxl_out, lambda_maxr_out, pstar
    REAL(KIND=8), INTENT(INOUT):: vstar
    INTEGER,      INTENT(OUT):: k
    REAL(KIND=NUMBER)        :: p1, phi1, phi11, p2, phi2, phi22, phi12, phi112, phi221
    LOGICAL                  :: check
    !===Initialization
   !  WRITE(*,*) "Computing lambda max in fortran"
    rhol= in_rhol
    ul = in_ul
    pl = in_pl
    el = in_el
    rhor =in_rhor
    ur = in_ur
    pr = in_pr
    er = in_er
    k = 0
    CALL init(rhol,el,pl,b_covolume,gammal,al,alphal,capAl,capBl,capCl,expol)
    CALL init(rhor,er,pr,b_covolume,gammar,ar,alphar,capAr,capBr,capCr,expor)
    IF (pl.LE.pr) THEN
       p_min     = pl
       rho_min   = rhol
       gamma_min = gammal
       gamma_min_index = 'l'
       alpha_min = alphal
       capA_min  = capAl
       capB_min  = capBl
       p_max     = pr
       rho_max   = rhor
       gamma_max = gammar
       alpha_max = alphar
       capC_max  = capCr
    ELSE
       p_min     = pr
       rho_min   = rhor
       gamma_min = gammar
       gamma_min_index = 'r'
       alpha_min = alphar
       capA_min  = capAr
       capB_min  = capBr
       p_max     = pl
       rho_max   = rhol
       gamma_max = gammal
       alpha_max = alphal
       capC_max  = capCl 
    END IF
    IF (gammal.LE.gammar) THEN
       gamma_lm = gammal
       gamma_lm_index = 'l'
       gamma_uM = gammar 
    ELSE
       gamma_lm = gammar
       gamma_lm_index = 'r'
       gamma_uM = gammal
    END IF
    expo_lm = (gamma_lm-1)/(2*gamma_lm)
    expo_uM = (gamma_uM-1)/(2*gamma_uM)
    expo_max = (gamma_max-1)/(2*gamma_max)
    numerator = alphal+alphar-ur+ul
    vacuum = capCl+capCr+ul-ur
    phi_pmin = capC_max*((p_min/p_max)**expo_max-1)+ur-ul 
    phi_pmax = (p_max-p_min)*SQRT(capA_min/(p_max+capB_min))+ur-ul

    !===Initialize p1 and p2
   !  WRITE(*,*) "F: Initializing p1 and p2"
    CALL initialize_p1_p2(p1,p2)

    IF (no_iter) THEN
       pstar = p2
      !  WRITE(*,*) "F: No iter"
       CALL no_iter_update_lambda(rhol,pl,al,gammal,rhor,pr,ar,gammar,p2,lambda_maxl_out,lambda_maxr_out)
       RETURN
    ELSE
       !===Iterations
      ! WRITE(*,*) "F: iter"
       p1 = MAX(p1,p2-phi(p2)/phi_prime(p2))
       DO WHILE(.TRUE.)
          CALL update_lambda(rhol,pl,al,gammal,rhor,pr,ar,gammar,p1,p2,in_tol,&
               lambda_maxl_out,lambda_maxr_out,check)
          pstar = p2
          IF (check) EXIT !RETURN
          phi1 =  phi(p1)
          phi11 = phi_prime(p1)
          phi2 =  phi(p2)
          phi22 = phi_prime(p2)
          IF (phi1>zero) THEN
             lambda_maxl_out = lambdaz(rhol,pl,al,gammal,p1,-1)
             lambda_maxr_out = lambdaz(rhor,pr,ar,gammar,p1, 1)
             pstar = p1
             EXIT
            !  RETURN
          END IF
          IF (phi2<zero) EXIT!RETURN
          phi12 = (phi2-phi1)/(p2-p1) 
          phi112 = (phi12-phi11)/(p2-p1)
          phi221 = (phi22-phi12)/(p2-p1)
          p1 = p1 - 2*phi1/(phi11 + SQRT(phi11**2 - 4*phi1*phi112))
          p2 = p2 - 2*phi2/(phi22 + SQRT(phi22**2 - 4*phi2*phi221))
          k = k+1
       END DO
    END IF
    ! Madison's contribution to compute vstar
    vstar = ustar(pstar)
   !  WRITE(*,*) "In Fortran, vstar is: ", vstar
  END SUBROUTINE Lagrangian_lambda_arbitrary_eos


  SUBROUTINE init(rho,e,p,b_covolume,gamma,a,alpha,capA,capB,capC,expo)
    IMPLICIT NONE
    REAL(KIND=NUMBER), INTENT(IN)  :: rho, e, p, b_covolume
    REAL(KIND=NUMBER), INTENT(OUT) :: gamma, a, alpha, capA, capB, capC, expo
    REAL(KIND=NUMBER) :: x

    x = 1-b_covolume*rho
    gamma = 1 + p*x/(rho*e)
    
   !  WRITE(*,*) "b_covolume: ", b_covolume
   !  WRITE(*,*) "p: ", p
   !  WRITE(*,*) "x: ", x
   !  WRITE(*,*) "rho: ", rho
   !  WRITE(*,*) "e: ", e
   !  write(*,*) "gamma: ", gamma

    a = SQRT(gamma*p/(rho*x))
    capC = 2*a*x/(gamma-1)
    alpha = cc(gamma)*capC
    capA = 2*x/((gamma+1)*rho)
    capB = p*(gamma-1)/(gamma+1)
    expo = (gamma-1)/(2*gamma)
  CONTAINS
    FUNCTION cc(gamma) RESULT(vv)
      IMPLICIT NONE
      REAL(KIND=NUMBER), INTENT(IN) :: gamma
      REAL(KIND=NUMBER)             :: vv
      IF (gamma.LE.1) THEN
         WRITE(*,*) "BUG: gamma .LE. 1"
         WRITE(*,*) "Gamma = ", gamma
         STOP
      ELSE IF (gamma .LE. five_third) THEN
         vv = one 
      ELSE IF (gamma .LE. 3) THEN
         vv = SQRT((3*gamma+11)/(6*(gamma+1)))
      ELSE
         expo = (4-2*gamma)/(gamma-1)
         vv = SQRT(half+2*3**expo/(gamma-1))
      END IF
    END FUNCTION cc
  END SUBROUTINE init

  SUBROUTINE initialize_p1_p2(p1,p2)
    IMPLICIT NONE
    REAL(KIND=NUMBER), INTENT(OUT) :: p1, p2
    REAL(KIND=NUMBER) :: phat1, phat2, r, rp
    REAL(KIND=NUMBER) :: xl, xr, a, b, c
    IF (vacuum.LE.zero) THEN
       p1=zero
       p2=zero
    ELSE IF (zero.LE.phi_pmin) THEN
       p1 = zero
       phat1 = p_min*(numerator/(alpha_min +alpha_max*(p_min/p_max)**expo_uM))**(1/expo_uM)
       p2 = MIN(p_min,phat1)
    ELSE IF (zero.LE.phi_pmax) THEN  
       p1 = p_min
       rp = p_min/p_max
       r = (rp)**((gamma_uM-gamma_lm)/(2*gamma_lm*gamma_uM))
       IF (gamma_min_index==gamma_lm_index) THEN
          phat1 = p_min*(numerator/(r*alpha_min +alpha_max*(rp)**expo_uM))**(1/expo_uM)
          phat2 = p_min*(numerator/(alpha_min +r*alpha_max*(rp)**expo_lm))**(1/expo_lm)
       ELSE
          phat1 = p_min*(numerator/(alpha_min +alpha_max*(rp)**expo_lm))**(1/expo_lm)
          phat2 = p_min*(numerator/(alpha_min +alpha_max*(rp)**expo_uM))**(1/expo_uM)
       END IF
       p2 = MIN(p_max,phat1,phat2)
    ELSE
       p1 = p_max
       p2 = p_min*(numerator/(alpha_min +alpha_max*(p_min/p_max)**expo_lm))**(1/expo_lm)
       xl = SQRT(capAl/(1+capBl/p_max))
       xr = SQRT(capAr/(1+capBr/p_max))
       a = xl+xr
       b = ur-ul
       c = -pl*xl-pr*xr
       phat2 = ((-b+SQRT(b*b-4*a*c))/(2*a))**2
       p2 = MIN(p2,phat2)
    END IF
  END SUBROUTINE initialize_p1_p2

  SUBROUTINE no_iter_update_lambda(rhol,pl,al,gammal,rhor,pr,ar,gammar,p2,lambda_maxl,lambda_maxr)
    IMPLICIT NONE
    REAL(KIND=NUMBER), INTENT(IN)  :: rhol, pl, al, gammal, rhor, pr, ar, gammar, p2
    REAL(KIND=NUMBER), INTENT(OUT) :: lambda_maxl, lambda_maxr
    REAL(KIND=NUMBER) :: v11, v32, lambda_max
    v11 = lambdaz(rhol,pl,al,gammal,p2,-1)
    v32 = lambdaz(rhor,pr,ar,gammar,p2,1)
    lambda_maxl = MAX(-v11,zero)
    lambda_maxr = MAX(v32,zero)
    lambda_max = MAX(lambda_maxl,lambda_maxr)
  END SUBROUTINE no_iter_update_lambda

  FUNCTION lambdaz(rhoz,pz,az,gammaz,pstar,z) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND=NUMBER), INTENT(IN) :: rhoz,pz,az,gammaz,pstar
    INTEGER,           INTENT(IN) :: z
    REAL(KIND=NUMBER)             :: vv
    vv = z*az*rhoz*SQRT(1+MAX((pstar-pz)/pz,zero)*(gammaz+1)/(2*gammaz))
  END FUNCTION lambdaz
  !===end of code if no iteration

  !=== code below is needed for iterative solver
  SUBROUTINE update_lambda(rhol,pl,al,gammal,rhor,pr,ar,gammar,p1,p2,tol,lambda_maxl,lambda_maxr,check)
    IMPLICIT NONE
    REAL(KIND=NUMBER), INTENT(IN)  :: rhol, pl, al, gammal, rhor, pr, ar, gammar
    REAL(KIND=NUMBER), INTENT(IN)  :: p1, p2, tol
    REAL(KIND=NUMBER), INTENT(OUT) :: lambda_maxl, lambda_maxr
    LOGICAL,           INTENT(OUT) :: check
    REAL(KIND=NUMBER) :: v11, v12, v31, v32, lambda_max, err1, err3
    v11 = lambdaz(rhol,pl,al,gammal,p2,-1)
    v12 = lambdaz(rhol,pl,al,gammal,p1,-1)
    v31 = lambdaz(rhor,pr,ar,gammar,p1,1)
    v32 = lambdaz(rhor,pr,ar,gammar,p2,1)
    lambda_maxl = MAX(-v11,zero)
    lambda_maxr = MAX(v32,zero)
    lambda_max = MAX(lambda_maxl,lambda_maxr)
    err3 =  ABS(v32 - v31)/lambda_max
    err1 =  ABS(v12 - v11)/lambda_max
    IF (MAX(err1,err3).LE.tol) THEN
       check = .TRUE.
    ELSE
       check = .FALSE.
    END IF
  END SUBROUTINE update_lambda

  FUNCTION phi(p) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND=NUMBER), INTENT(IN) :: p
    REAL(KIND=NUMBER)             :: vv
    vv = f(p,pl,capAl,capBl,capCl,expol) + f(p,pr,capAr,capBr,capCr,expor)+ur-ul
  END FUNCTION phi

  FUNCTION f(p,pz,capAz,capBz,capCz,expoz) RESULT(ff)
    REAL(KIND=NUMBER), INTENT(IN) :: p, pz, capAz, capBz, capCz, expoz
    REAL(KIND=NUMBER)             :: ff
    IF (p.LE.pz) THEN
       ff = capCz*((p/pz)**expoz-1)
    ELSE
       ff = (p-pz)*SQRT(capAz/(p+capBz))
    END IF
  END FUNCTION f

  FUNCTION phi_prime(p) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND=NUMBER), INTENT(IN) :: p
    REAL(KIND=NUMBER)             :: vv
    vv = fp(p,pl,capAl,capBl,capCl,expol) + fp(p,pr,capAr,capBr,capCr,expor)
  CONTAINS
    FUNCTION fp(p,pz,capAz,capBz,capCz,expoz) RESULT(ff)
      REAL(KIND=NUMBER), INTENT(IN) :: p, pz, capAz, capBz, capCz, expoz
      REAL(KIND=NUMBER)             :: ff
      IF (p.LE.pz) THEN
         ff = capCz*expoz*(p/pz)**(expoz-1)/pz
      ELSE
         ff = SQRT(capAz/(p+capBz))*(1-(p-pz)/(2*(capBz+p)))
      END IF
    END FUNCTION fp
  END FUNCTION phi_prime

  !===Optional functions to compute rhostar and ustar (not necessary)
  FUNCTION ustar(pstar) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND=NUMBER), INTENT(IN) :: pstar
    REAL(KIND=NUMBER)             :: vv
    vv = half*(ul-f(pstar,pl,capAl,capBl,capCl,expol)+ur+f(pstar,pr,capAr,capBr,capCr,expor))
  END FUNCTION ustar

  FUNCTION rhostar(pstar,rhoz,pz,gammaz,b_covolume) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND=NUMBER), INTENT(IN) :: pstar, rhoz, pz, gammaz, b_covolume
    REAL(KIND=NUMBER)             :: vv
    IF (pstar.LE.pz) THEN
       vv = rhoz /(b_covolume*rhoz+(1-b_covolume*rhoz)*(pz/pstar)**(1/gammaz))
    ELSE
       vv = rhoz*(pstar/pz+(gammaz-1)/(gammaz+1))/&
            (((gammaz-1+2*b_covolume*rhoz)*pstar)/((gammaz+1)*pz) &
            + (gammaz+1-2*b_covolume*rhoz)/(gammaz+1))
    END IF
  END FUNCTION rhostar

END MODULE arbitrary_eos_Lagrangian_lambda_module
