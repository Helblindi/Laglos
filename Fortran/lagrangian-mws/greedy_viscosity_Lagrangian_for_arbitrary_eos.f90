!===Authors: Bennett Clayton, Jean-Luc Guermond, and Bojan Popov, Texas A&M, April 5, 2021
MODULE arbitrary_eos_lagrangian_greedy_lambda_module
  IMPLICIT NONE
  PUBLIC               :: greedy_lambda_arbitrary_eos !===Main function
  PUBLIC               :: rhostar, ustar, phi  !===Optional functions. Can be removed
  REAL(KIND=8), PUBLIC :: b_covolume = 0.0d0   !===Covolume constant, if known
  PRIVATE
  INTEGER, PARAMETER:: NUMBER = KIND(1.d0)
  REAL(KIND=NUMBER), PARAMETER :: zero = 0
  REAL(KIND=NUMBER), PARAMETER :: one = 1
  REAL(KIND=NUMBER), PARAMETER :: half = 0.5d0
  REAL(KIND=NUMBER), PARAMETER :: five_third = 5.d0/3.d0, epsilon=1.d-10
  REAL(KIND=NUMBER) :: rhol, ul, pl, el
  REAL(KIND=NUMBER) :: rhor, ur, pr, er
  REAL(KIND=NUMBER) :: cmin, gamma
  REAL(KIND=NUMBER) :: gammal, al, alphal, capAl, capBl, capCl, expol
  REAL(KIND=NUMBER) :: gammar, ar, alphar, capAr, capBr, capCr, expor
  REAL(KIND=NUMBER) :: p_min, rho_min, gamma_min, alpha_min, capA_min, capB_min, phi_pmin, e_min, c_max, S_max
  REAL(KIND=NUMBER) :: p_max, rho_max, gamma_max, alpha_max, capC_max, expo_max, phi_pmax, e_max, c_min, S_min
  REAL(KIND=NUMBER) :: gamma_lm, expo_lm
  REAL(KIND=NUMBER) :: gamma_uM, expo_uM
  REAL(KIND=NUMBER) :: numerator, vacuum
  CHARACTER(LEN=1)  :: gamma_min_index, gamma_lm_index

CONTAINS

  SUBROUTINE greedy_lambda_arbitrary_eos(in_rhol,in_ul,in_el,in_pl,in_rhor,in_ur,in_er,in_pr,in_tol,no_iter,&
       lambda_max,pstar,k)
    IMPLICIT NONE
    REAL(KIND=8), INTENT(IN) :: in_rhol, in_el, in_rhor, in_er, in_tol
    REAL(KIND=8), INTENT(IN), TARGET :: in_ul, in_pl, in_ur, in_pr
    LOGICAL,      INTENT(IN) :: no_iter
    REAL(KIND=8), INTENT(OUT):: lambda_max, pstar
    INTEGER,      INTENT(OUT):: k
    REAL(KIND=NUMBER)        :: p1, phi1, phi11, p2, phi2, phi22, phi12, phi112, phi221, &
         lambda_maxl_out,lambda_maxr_out
    LOGICAL                  :: check
    !===Initialization
    rhol= in_rhol
    ul = in_ul
    pl = in_pl
    el = in_el
    rhor =in_rhor
    ur = in_ur
    pr = in_pr
    er = in_er
    k = 0
    CALL init(rhol,el,pl,gammal,al,alphal,capAl,capBl,capCl,expol)
    CALL init(rhor,er,pr,gammar,ar,alphar,capAr,capBr,capCr,expor)
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
       e_min = el
       e_max = er
       S_min = el*(1/rhol-b_covolume)**(gammal-1)
       S_max = er*(1/rhor-b_covolume)**(gammar-1)
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
       e_min = er
       e_max = el
       S_max = el*(1/rhol-b_covolume)**(gammal-1)
       S_min = er*(1/rhor-b_covolume)**(gammar-1)
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
    CALL initialize_p1_p2(p1,p2)

    IF (no_iter) THEN
       pstar = p2
       !phi2 =  phi(p2)
       !phi22 = phi_prime(p2)
       !p1 = p2 - phi2/phi22
       !CALL update_lambda(ul,pl,al,gammal,ur,pr,ar,gammar,p1,p2,in_tol,&
       !     lambda_maxl_out,lambda_maxr_out,check)
       CALL no_iter_update_lambda(ul,pl,al,gammal,ur,pr,ar,gammar,p1,p2,lambda_max)
       RETURN
    ELSE
       !===Iterations
       p1 = MAX(p1,p2-phi(p2)/phi_prime(p2))
       DO WHILE(.TRUE.)
          CALL update_lambda(ul,pl,al,gammal,ur,pr,ar,gammar,p1,p2,in_tol,&
               lambda_maxl_out,lambda_maxr_out,check)
          pstar = p2
          IF (check) THEN
             lambda_max = MAX(lambda_maxl_out, lambda_maxr_out)
             RETURN
          END IF
          phi1 =  phi(p1)
          phi11 = phi_prime(p1)
          phi2 =  phi(p2)
          phi22 = phi_prime(p2)
          IF (phi1>zero) THEN
             lambda_maxl_out = lambdaz(rhol,ul,pl,al,gammal,p1,-1)
             lambda_maxr_out = lambdaz(rhor,ur,pr,ar,gammar,p1, 1)
             pstar = p1
             lambda_max = MAX(lambda_maxl_out, lambda_maxr_out)
             RETURN
          END IF
          IF (phi2<zero) THEN
             lambda_max = MAX(lambda_maxl_out, lambda_maxr_out)
             RETURN
          END IF
          phi12 = (phi2-phi1)/(p2-p1) 
          phi112 = (phi12-phi11)/(p2-p1)
          phi221 = (phi22-phi12)/(p2-p1)
          p1 = p1 - 2*phi1/(phi11 + SQRT(phi11**2 - 4*phi1*phi112))
          p2 = p2 - 2*phi2/(phi22 + SQRT(phi22**2 - 4*phi2*phi221))
          k = k+1
       END DO
    END IF
  END SUBROUTINE greedy_lambda_arbitrary_eos

  SUBROUTINE init(rho,e,p,gamma,a,alpha,capA,capB,capC,expo)
    IMPLICIT NONE
    REAL(KIND=NUMBER), INTENT(IN)  :: rho, e, p
    REAL(KIND=NUMBER), INTENT(OUT) :: gamma, a, alpha, capA, capB, capC, expo
    REAL(KIND=NUMBER) :: x
    x = 1-b_covolume*rho
    gamma = 1 + p*x/(rho*e)
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

  SUBROUTINE no_iter_update_lambda(ul,pl,al,gammal,ur,pr,ar,gammar,p1,p2,lambda_max)
    IMPLICIT NONE
    REAL(KIND=NUMBER), INTENT(IN)  :: ul, pl, al, gammal, ur, pr, ar, gammar, p1
    REAL(KIND=NUMBER) :: p2
    REAL(KIND=NUMBER), INTENT(OUT) :: lambda_max 
    REAL(KIND=NUMBER) :: v11, v32, x, y, lambda_maxl, lambda_maxr, &
         tau_l, tau_r, tau_max, tau_min, tau_small, lambda_max_safe, lambda_max_small, E_small, &
         tau_star_min, tau_star_max, ttt
    REAL(KIND=8) :: unl(3),unr(3), Plr(3)
    v11 = lambdaz(rhol,ul,pl,al,gammal,p2,-1)
    v32 = lambdaz(rhor,ur,pr,ar,gammar,p2,1)
    lambda_maxl = MAX(-v11,zero)
    lambda_maxr = MAX(v32,zero)
    lambda_max_safe = MAX(lambda_maxl,lambda_maxr)
    lambda_max_small = epsilon*lambda_max_safe
    E_small = epsilon*MAX(el,er)
    tau_max = 1/MIN(rhol,rhor)
    tau_min = 1/MAX(rhol,rhor)
    IF (0.d0< phi_pmin) THEN !===two expansions
       x = b_covolume + (1/rhol - b_covolume)*(pl/p2)**(1/gammal)
       y = b_covolume + (1/rhor - b_covolume)*(pr/p2)**(1/gammar)
       tau_max = MAX(tau_max, x, y)
    ELSE IF (0.d0< phi_pmax) THEN !===one expansions one shock
       !===Shock
       x = rho_min*(p2/p_min + (gamma_min-1)/(gamma_min+1)) &
            /((p2/p_min)*(gamma_min-1+2*b_covolume*rho_min)/ &
            (gamma_min+1)+(gamma_min+1-2*b_covolume*rho_min)/(gamma_min+1))
       tau_min = MIN(tau_min,1/x)
       !===Expansion
       x = b_covolume + (1/rho_max - b_covolume)*(p_max/p2)**(1/gammal)
       tau_max = MAX(tau_max, x)
    ELSE !===two shocks
       !===Shock
       x = rho_min*(p2/p_min + (gamma_min-1)/(gamma_min+1)) &
            /((p2/p_min)*(gamma_min-1+2*b_covolume*rho_min)/(gamma_min+1)+(gamma_min+1-2*b_covolume*rho_min)/(gamma_min+1))
       y = rho_max*(p2/p_max + (gamma_max-1)/(gamma_max+1)) &
            /((p2/p_max)*(gamma_max-1+2*b_covolume*rho_max)/(gamma_max+1)+(gamma_max+1-2*b_covolume*rho_max)/(gamma_max+1))
       tau_min = MIN(tau_min,1/x,1/y)
    END IF

    tau_small = epsilon*tau_max
    tau_l = 1/rhol
    tau_r = 1/rhor
    x = (ur-ul)/(2*tau_max-tau_l-tau_r+tau_small)
    y = (ul-ur)/(tau_l+tau_r-2*tau_min+tau_small)
    lambda_max = MAX(x,y,lambda_max_small)
    
    !===Now minimum principle on the entropy
    IF (0.d0< phi_pmin) THEN
       tau_star_min =  b_covolume + (1/rho_min - b_covolume)*(p_min/p1)**(1/gamma_min)
       c_min = S_min*(tau_star_min - b_covolume)**(gamma_lm-gamma_min)
       tau_star_max =  b_covolume + (1/rho_max - b_covolume)*(p_max/p1)**(1/gamma_max)
       c_max = S_max*(tau_star_max - b_covolume)**(gamma_lm-gamma_max)  
    ELSE IF (0.d0< phi_pmax) THEN
       c_min = S_min*(tau_min - b_covolume)**(gamma_lm-gamma_min)
       tau_star_max =  b_covolume + (1/rho_max - b_covolume)*(p_max/p1)**(1/gamma_max)
       c_max = S_max*(tau_star_max - b_covolume)**(gamma_lm-gamma_max)  
    ELSE
       c_min = S_min*(tau_min - b_covolume)**(gamma_lm-gamma_min)
       c_max = S_max*(tau_max - b_covolume)**(gamma_lm-gamma_max)
    END IF
    cmin = 0.99*MIN(c_min,c_max)

    Plr(1) = (ur-ul)/2
    Plr(2) = -(pr-pl)/2
    Plr(3) = -(pr*ur-pl*ul)/2
    unl(1) = 1/rhol
    unl(2) = ul
    unl(3) = el + ul**2/2
    unr(1) = 1/rhor
    unr(2) = ur
    unr(3) = er + ur**2/2

    gamma = gamma_min
    ttt = 1.d0/lambda_max
    CALL Newton_secant(unl,unr,Plr,ttt,psi_entrop,psi_entrop_prime,E_small)
    IF (ttt.LE.0.d0) THEN
       !WRITE(*,*) ' BUG ttt<0', ttt
       !STOP
       ttt = 1.d3/lambda_max_safe
    END IF
    lambda_max = MAX(lambda_max,1.d0/ttt)
    IF (lambda_max > 1.01*lambda_max_safe) THEN
       WRITE(*,*) ' Problem here ', lambda_max, lambda_max_safe
       write(*,*) rhol, rhor, ul, ur, pl, pr
       STOP
    END IF
    !write(*,*) 'lambda_max, lambda_max_safe', lambda_max, lambda_max_safe
    
  END SUBROUTINE no_iter_update_lambda

  FUNCTION lambdaz(rhoz,uz,pz,az,gammaz,pstar,z) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND=NUMBER), INTENT(IN) :: rhoz,uz,pz,az,gammaz,pstar
    INTEGER,           INTENT(IN) :: z
    REAL(KIND=NUMBER)             :: vv
    vv = z*rhoz*az*SQRT(1+MAX((pstar-pz)/pz,zero)*(gammaz+1)/(2*gammaz))
  END FUNCTION lambdaz
  !===end of code if no iteration

  !=== code below is needed for iterative solver
  SUBROUTINE update_lambda(ul,pl,al,gammal,ur,pr,ar,gammar,p1,p2,tol,lambda_maxl,lambda_maxr,check)
    IMPLICIT NONE
    REAL(KIND=NUMBER), INTENT(IN)  :: ul, pl, al, gammal, ur, pr, ar, gammar
    REAL(KIND=NUMBER), INTENT(IN)  :: p1, p2, tol
    REAL(KIND=NUMBER), INTENT(OUT) :: lambda_maxl, lambda_maxr
    LOGICAL,           INTENT(OUT) :: check
    REAL(KIND=NUMBER) :: v11, v12, v31, v32, lambda_max, err1, err3
    v11 = lambdaz(rhol,ul,pl,al,gammal,p2,-1)
    v12 = lambdaz(rhol,ul,pl,al,gammal,p1,-1)
    v31 = lambdaz(rhor,ur,pr,ar,gammar,p1,1)
    v32 = lambdaz(rhor,ur,pr,ar,gammar,p2,1)
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

  FUNCTION rhostar(pstar,rhoz,pz,gammaz) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND=NUMBER), INTENT(IN) :: pstar, rhoz, pz, gammaz
    REAL(KIND=NUMBER)             :: vv
    IF (pstar.LE.pz) THEN
       vv = rhoz /(b_covolume*rhoz+(1-b_covolume*rhoz)*(pz/pstar)**(1/gammaz))
    ELSE
       vv = rhoz*(pstar/pz+(gammaz-1)/(gammaz+1))/&
            (((gammaz-1+2*b_covolume*rhoz)*pstar)/((gammaz+1)*pz) &
            + (gammaz+1-2*b_covolume*rhoz)/(gammaz+1))
    END IF
  END FUNCTION rhostar


  SUBROUTINE Newton_secant(ui, uj, Pij,limiter,psi_func,psi_prime_func,psi_small)
    IMPLICIT NONE
    INTERFACE
       FUNCTION psi_func(u) RESULT(psi)
         IMPLICIT NONE
         REAL(KIND=8), DIMENSION(3), INTENT(IN)  :: u
         REAL(KIND=8)                     :: psi
       END FUNCTION psi_func
       FUNCTION psi_prime_func(Pij,u) RESULT(psi)
         IMPLICIT NONE
         REAL(KIND=8), DIMENSION(3), INTENT(IN)  :: u, Pij
         REAL(KIND=8)                     :: psi
       END FUNCTION psi_prime_func
    END INTERFACE
    REAL(KIND=8), DIMENSION(3), INTENT(IN) :: Pij
    REAL(KIND=8), DIMENSION(3), INTENT(IN) :: ui, uj
    REAL(KIND=8), INTENT(IN)    :: psi_small
    REAL(KIND=8), INTENT(INOUT) :: limiter
    REAL(KIND=8), DIMENSION(3) ::  ul, ur, ulr
    REAL(KIND=8) :: psil, psir, ll, lr, llold, lrold
    LOGICAL :: once

    ul = ui
    ur = uj
    ulr = (ul+ur)/2

    lr = limiter
    ur = ulr + lr*Pij
    psir = psi_func(ur)
    IF (psir.GE.0.d0) THEN 
       !===input limiter is okay
       RETURN
    END IF
    ll = 0.d0
    ul = ulr
    psil = MAX(psi_func(ulr),0.d0)
    once=.TRUE.
    DO WHILE (ABS(psil-psir) .GT. psi_small .OR. once)
       once =.FALSE.
       llold = ll
       lrold = lr
       ll = ll - psil*(lr-ll)/(psir-psil)
       lr = lr - psir/psi_prime_func(Pij,ur)
       IF (ll.GE.lr) THEN
          ll = lr
          EXIT
       END IF
       IF (ll< llold) THEN
          ll = llold
          EXIT
       END IF
       IF (lr > lrold) THEN
          lr = lrold
          EXIT
       END IF
       ul = ulr + ll*Pij
       ur = ulr + lr*Pij
       psil = psi_func(ul)
       psir = psi_func(ur)
    END DO

    IF (psir.GE.0.d0) THEN 
       limiter = lr
    ELSE
       limiter = ll
    END IF

  END SUBROUTINE Newton_secant

  FUNCTION psi_entrop(u) RESULT(psi)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(3), INTENT(IN) :: u
    REAL(KIND=8)                                 :: psi
    psi = u(3) - u(2)**2/2 - cmin*(u(1)-b_covolume)**(1-gamma)
  END FUNCTION psi_entrop

  FUNCTION psi_entrop_prime(Pij,u) RESULT(psi)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(3), INTENT(IN) :: u, Pij
    REAL(KIND=8)                           :: psi
    psi =(gamma-1)*cmin*Pij(1)*(u(1)-b_covolume)**(-gamma)&
         -u(2)*Pij(2)+Pij(3)
  END FUNCTION psi_entrop_prime

END MODULE arbitrary_eos_lagrangian_greedy_lambda_module
