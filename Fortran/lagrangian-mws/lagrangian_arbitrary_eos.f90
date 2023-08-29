!===Authors: Bennett Clayton, Jean-Luc Guermond, and Bojan Popov, Texas A&M, April 5, 2021
MODULE arbitrary_eos_lagrangian_lambda_module
   IMPLICIT NONE
   PUBLIC               :: lambda_arbitrary_eos !===Main function
   PUBLIC               :: rhostar, ustar, phi  !===Optional functions. Can be removed
   REAL(KIND=8), PUBLIC :: b_covolume = 0.d0   !===Covolume constant, if known
   !=== Van der waals constants, if known
   REAL(KIND=8), PUBLIC :: a_vdw = 1.d0, b_vdw = 1.d0, gamma_vdw = 1.02d0 
   REAL(KIND=8), PARAMETER, PUBLIC :: p_infty = 0.d0      !===Reference pressure
   REAL(KIND=8), PARAMETER, PUBLIC :: q = 0.d0            !===Reference specific internal energy
   PRIVATE
   INTEGER, PARAMETER:: NUMBER = KIND(1.d0)
   REAL(KIND=NUMBER), PARAMETER :: five_third = 5.d0/3.d0
   REAL(KIND=NUMBER) :: taul, ul, pl, el
   REAL(KIND=NUMBER) :: taur, ur, pr, er
   REAL(KIND=NUMBER) :: gammal, al, alphal, capAl, capBl, capCl, expol
   REAL(KIND=NUMBER) :: gammar, ar, alphar, capAr, capBr, capCr, expor
   REAL(KIND=NUMBER) :: p_min, tau_min, gamma_min, alpha_min, capA_min, capB_min, phi_pmin
   REAL(KIND=NUMBER) :: p_max, tau_max, gamma_max, alpha_max, capC_max, expo_max, phi_pmax
   REAL(KIND=NUMBER) :: gamma_lm, expo_lm
   REAL(KIND=NUMBER) :: gamma_uM, expo_uM
   REAL(KIND=NUMBER) :: numerator, vacuum
   CHARACTER(LEN=1)  :: gamma_min_index, gamma_lm_index

CONTAINS

   SUBROUTINE lambda_arbitrary_eos(in_rhol, in_ul, in_el, in_pl, in_rhor, in_ur, in_er, in_pr, in_tol, &
                                   WANT_ITERATION, lambda_maxl_out, lambda_maxr_out, pstar, k, &
                                   b_covolume_in, a_vdw_in, b_vdw_in, gamma_vdw_in)
      use iso_c_binding
      IMPLICIT NONE
      REAL(KIND=8), INTENT(IN) :: in_rhol, in_el, in_rhor, in_er, in_tol
      REAL(KIND=8), INTENT(IN), TARGET :: in_ul, in_pl, in_ur, in_pr
      REAL(KIND=8), INTENT(IN) :: b_covolume_in, a_vdw_in, b_vdw_in, gamma_vdw_in
      LOGICAL(c_bool), INTENT(IN) :: WANT_ITERATION
      REAL(KIND=8), INTENT(OUT):: lambda_maxl_out, lambda_maxr_out, pstar
      INTEGER, INTENT(OUT):: k
      REAL(KIND=NUMBER)        :: p1, phi1, phi11, p2, phi2, phi22, phi12, phi112, phi221
      LOGICAL                  :: check
      !===EOS initialization
      b_covolume = b_covolume_in
      a_vdw = a_vdw_in
      b_vdw = b_vdw_in
      gamma_vdw = gamma_vdw_in

      !===State vars initialiation
      taul = 1.d0/in_rhol
      ul = in_ul
      pl = in_pl
      el = in_el
      taur = 1.d0/in_rhor
      ur = in_ur
      pr = in_pr
      er = in_er
      k = 0
      CALL init(taul, el, pl, gammal, al, alphal, capAl, capBl, capCl, expol)
      CALL init(taur, er, pr, gammar, ar, alphar, capAr, capBr, capCr, expor)
      IF (pl <= pr) THEN
         p_min = pl
         tau_min = taul
         gamma_min = gammal
         gamma_min_index = 'l'
         alpha_min = alphal
         capA_min = capAl
         capB_min = capBl
         p_max = pr
         tau_max = taur
         gamma_max = gammar
         alpha_max = alphar
         capC_max = capCr
      ELSE
         p_min = pr
         tau_min = taur
         gamma_min = gammar
         gamma_min_index = 'r'
         alpha_min = alphar
         capA_min = capAr
         capB_min = capBr
         p_max = pl
         tau_max = taul
         gamma_max = gammal
         alpha_max = alphal
         capC_max = capCl
      END IF
      IF (gammal <= gammar) THEN
         gamma_lm = gammal
         gamma_lm_index = 'l'
         gamma_uM = gammar
      ELSE
         gamma_lm = gammar
         gamma_lm_index = 'r'
         gamma_uM = gammal
      END IF
      expo_lm = (gamma_lm - 1.d0)/(2.d0*gamma_lm)
      expo_uM = (gamma_uM - 1.d0)/(2.d0*gamma_uM)
      expo_max = (gamma_max - 1.d0)/(2.d0*gamma_max)
      numerator = alphal + alphar - ur + ul
      vacuum = capCl + capCr + ul - ur
      phi_pmin = capC_max*(((p_min + p_infty)/(p_max + p_infty))**expo_max - 1.d0) + ur - ul
      phi_pmax = (p_max - p_min)*SQRT(capA_min/(p_max + p_infty + capB_min)) + ur - ul

      !===Initialize p1 and p2 where p1 <= pstar <= p2
      CALL initialize_p1_p2(p1, p2)

      IF (.NOT. WANT_ITERATION) THEN   
         pstar = p2
         CALL no_iter_update_lambda(taul, pl, al, gammal, taur, pr, ar, gammar, p2, lambda_maxl_out, lambda_maxr_out)
         RETURN
      ELSE
         !===Iterations
         p1 = MAX(p1, p2 - phi(p2)/phi_prime(p2))
         DO WHILE (.TRUE.)
            CALL update_lambda(taul, pl, al, gammal, taur, pr, ar, gammar, p1, p2, in_tol, &
                               lambda_maxl_out, lambda_maxr_out, check)
            pstar = p2
            IF (check) RETURN
            phi1 = phi(p1)
            phi11 = phi_prime(p1)
            phi2 = phi(p2)
            phi22 = phi_prime(p2)
            IF (phi1 > 0.d0) THEN
               lambda_maxl_out = lambdaz(taul, pl, al, gammal, p1, -1)
               lambda_maxr_out = lambdaz(taur, pr, ar, gammar, p1, 1)
               pstar = p1
               RETURN
            END IF
            IF (phi2 < 0.d0) RETURN
            phi12 = (phi2 - phi1)/(p2 - p1)
            phi112 = (phi12 - phi11)/(p2 - p1)
            phi221 = (phi22 - phi12)/(p2 - p1)
            p1 = p1 - 2*phi1/(phi11 + SQRT(phi11**2 - 4*phi1*phi112))
            p2 = p2 - 2*phi2/(phi22 + SQRT(phi22**2 - 4*phi2*phi221))
            k = k + 1
         END DO
      END IF
   END SUBROUTINE lambda_arbitrary_eos

   SUBROUTINE init(tau, e, p, gamma, a, alpha, capA, capB, capC, expo)
      IMPLICIT NONE
      REAL(KIND=NUMBER), INTENT(IN)  :: tau, e, p
      REAL(KIND=NUMBER), INTENT(OUT) :: gamma, a, alpha, capA, capB, capC, expo
      REAL(KIND=NUMBER) :: x
      x = tau - b_covolume
      ! write(*,*) "tau: ", tau, ", b_covolume: ", b_covolume, ", x: ", x
      !===local gamma (gamma_Z)
      ! WRITE(*,*) "p: ", p, ", pinf: ", p_infty, ", x: ", x, ", e: ", e, ", q: ", q
      gamma = 1.d0 + (p + p_infty)*x/(e - q - p_infty*x)
      !===local sound speed (a_Z)
      a = tau*SQRT(gamma*(p + p_infty)/x)
      !===other relevant constants
      capC = 2.d0*a*x/(tau*(gamma - 1.d0))
      alpha = cc(gamma)*capC
      capA = 2.d0*x/(gamma + 1.d0)
      capB = (gamma - 1.d0)/(gamma + 1.d0)*(p + p_infty)
      expo = 0.5d0*(gamma - 1.d0)/gamma
   CONTAINS
      FUNCTION cc(g) RESULT(c_of_gamma)
         IMPLICIT NONE
         REAL(KIND=NUMBER), INTENT(IN) :: g
         REAL(KIND=NUMBER)             :: c_of_gamma, expo_of_three
         IF (g .LE. 1.d0) THEN
            WRITE (*, *) "BUG: gamma less than or equal to 1"
            WRITE (*, *) "gamma = ", g
            STOP
         ELSE IF (g .LE. five_third) THEN
            c_of_gamma = 1.d0
         ELSE IF (g .LE. 3.d0) THEN
            c_of_gamma = SQRT((3.d0*g + 11.d0)/(6.d0*(g + 1.d0)))
         ELSE
            expo_of_three = (4.d0 - 2.d0*g)/(g - 1.d0)
            c_of_gamma = SQRT(0.5d0 + 2.d0*3.d0**expo_of_three/(g - 1.d0))
         END IF
      END FUNCTION cc
   END SUBROUTINE init

   SUBROUTINE initialize_p1_p2(p1, p2)
      IMPLICIT NONE
      REAL(KIND=NUMBER), INTENT(OUT) :: p1, p2
      REAL(KIND=NUMBER) :: phat1, phat2, r, p_ratio
      REAL(KIND=NUMBER) :: xl, xr, a, b, c

      IF (vacuum <= 0.d0) THEN
         p1 = 0.d0
         p2 = 0.d0
         RETURN
      END IF
      p_ratio = (p_min + p_infty)/(p_max + p_infty)
      IF (0.d0 <= phi_pmin) THEN
         p1 = 0.d0
         phat1 = p_min*(numerator/(alpha_min + alpha_max*(p_ratio)**expo_uM))**(1.d0/expo_uM)
         p2 = MIN(p_min, phat1)
      ELSE IF (0.d0 <= phi_pmax) THEN
         p1 = p_min + p_infty
         r = (p_ratio)**((gamma_uM - gamma_lm)/(2.d0*gamma_lm*gamma_uM))
         IF (gamma_min_index == gamma_lm_index) THEN
            phat1 = (p_min + p_infty)*(numerator/(r*alpha_min &
                                                  + alpha_max*(p_ratio)**expo_uM))**(1.d0/expo_uM) - p_infty
            phat2 = (p_min + p_infty)*(numerator/(alpha_min &
                                                  + r*alpha_max*(p_ratio)**expo_lm))**(1.d0/expo_lm) - p_infty
         ELSE
            phat1 = (p_min + p_infty)*(numerator/(alpha_min &
                                                  + alpha_max*(p_ratio)**expo_lm))**(1.d0/expo_lm) - p_infty
            phat2 = (p_min + p_infty)*(numerator/(alpha_min &
                                                  + alpha_max*(p_ratio)**expo_uM))**(1.d0/expo_uM) - p_infty
         END IF
         p2 = MIN(p_max, phat1, phat2)
      ELSE
         p1 = p_max
         p2 = p_min*(numerator/(alpha_min &
                                + alpha_max*(p_ratio)**expo_lm))**(1.d0/expo_lm) - p_infty
         xl = SQRT(capAl*(p_max + p_infty)/(p_max + p_infty + capBl))
         xr = SQRT(capAr*(p_max + p_infty)/(p_max + p_infty + capBr))
         a = xl + xr
         b = ur - ul
         c = -pl*xl - pr*xr
         phat2 = ((-b + SQRT(b*b - 4.d0*a*c))/(2.d0*a))**2
         p2 = MIN(p2, phat2)
      END IF
   END SUBROUTINE initialize_p1_p2

   SUBROUTINE no_iter_update_lambda(tau_L, p_L, a_L, gamma_L, tau_R, p_R, a_R, gamma_R, &
                                    p2, lambda_max_L, lambda_max_R)
      IMPLICIT NONE
      REAL(KIND=NUMBER), INTENT(IN)  :: tau_L, p_L, a_L, gamma_L, tau_R, p_R, a_R, gamma_R, p2
      REAL(KIND=NUMBER), INTENT(OUT) :: lambda_max_L, lambda_max_R
      REAL(KIND=NUMBER) :: v11, v32, lambda_max
      v11 = lambdaz(tau_L, p_L, a_L, gamma_L, p2, -1)
      v32 = lambdaz(tau_R, p_R, a_R, gamma_R, p2, 1)
      lambda_max_L = MAX(-v11, 0.d0)
      lambda_max_R = MAX(v32, 0.d0)
      lambda_max = MAX(lambda_max_L, lambda_max_R)
   END SUBROUTINE no_iter_update_lambda

   FUNCTION lambdaz(tauz, pz, az, gammaz, pstar, z) RESULT(vv)
      IMPLICIT NONE
      REAL(KIND=NUMBER), INTENT(IN) :: tauz, pz, az, gammaz, pstar
      INTEGER, INTENT(IN) :: z
      REAL(KIND=NUMBER)             :: vv
      vv = z*az/tauz*SQRT(1.d0 + MAX((pstar - pz)/(pz + p_infty), 0.d0)*(gammaz + 1.d0)/(2.d0*gammaz))
   END FUNCTION lambdaz
   !===end of code if no iteration

   !=== code below is needed for iterative solver
   SUBROUTINE update_lambda(tau_L, p_L, a_L, gamma_L, tau_R, p_R, a_R, gamma_R, &
                            p1, p2, tol, lambda_max_L, lambda_max_R, check)
      IMPLICIT NONE
      REAL(KIND=NUMBER), INTENT(IN)  :: tau_L, p_L, a_L, gamma_L, tau_R, p_R, a_R, gamma_R
      REAL(KIND=NUMBER), INTENT(IN)  :: p1, p2, tol
      REAL(KIND=NUMBER), INTENT(OUT) :: lambda_max_L, lambda_max_R
      LOGICAL, INTENT(OUT) :: check
      REAL(KIND=NUMBER) :: v11, v12, v31, v32, lambda_max, err1, err3
      v11 = lambdaz(tau_L, p_L, a_L, gamma_L, p2, -1)
      v12 = lambdaz(tau_L, p_L, a_L, gamma_L, p1, -1)
      v31 = lambdaz(tau_R, p_R, a_R, gamma_R, p1, 1)
      v32 = lambdaz(tau_R, p_R, a_R, gamma_R, p2, 1)
      lambda_max_L = MAX(-v11, 0.d0)
      lambda_max_R = MAX(v32, 0.d0)
      lambda_max = MAX(lambda_max_L, lambda_max_R)
      err3 = ABS(v32 - v31)/lambda_max
      err1 = ABS(v12 - v11)/lambda_max
      IF (MAX(err1, err3) <= tol) THEN
         check = .TRUE.
      ELSE
         check = .FALSE.
      END IF
   END SUBROUTINE update_lambda

   FUNCTION phi(p) RESULT(vv)
      IMPLICIT NONE
      REAL(KIND=NUMBER), INTENT(IN) :: p
      REAL(KIND=NUMBER)             :: vv
      vv = f(p, pl, capAl, capBl, capCl, expol) + f(p, pr, capAr, capBr, capCr, expor) + ur - ul
   END FUNCTION phi

   FUNCTION f(p, pz, capAz, capBz, capCz, expoz) RESULT(ff)
      IMPLICIT NONE
      REAL(KIND=NUMBER), INTENT(IN) :: p, pz, capAz, capBz, capCz, expoz
      REAL(KIND=NUMBER)             :: ff
      IF (p <= pz) THEN
         ff = capCz*(((p + p_infty)/(pz + p_infty))**expoz - 1.d0)
      ELSE
         ff = (p - pz)*SQRT(capAz/(p + p_infty + capBz))
      END IF
   END FUNCTION f

   FUNCTION phi_prime(p) RESULT(val_phi)
      IMPLICIT NONE
      REAL(KIND=NUMBER), INTENT(IN) :: p
      REAL(KIND=NUMBER)             :: val_phi
      val_phi = f_prime(p, pl, capAl, capBl, capCl, expol) + f_prime(p, pr, capAr, capBr, capCr, expor)
   CONTAINS
      FUNCTION f_prime(p_var, pz, capAz, capBz, capCz, expoz) RESULT(val_f)
         REAL(KIND=NUMBER), INTENT(IN) :: p_var, pz, capAz, capBz, capCz, expoz
         REAL(KIND=NUMBER)             :: val_f
         IF (p_var <= pz) THEN
            val_f = capCz*expoz*((p_var + p_infty)/(pz + p_infty))**(expoz - 1.d0)/(pz + p_infty)
         ELSE
            val_f = SQRT(capAz/(p_var + p_infty + capBz))*(1.d0 - (p_var - pz)/(2.d0*(p_var + p_infty + capBz)))
         END IF
      END FUNCTION f_prime
   END FUNCTION phi_prime

   !===Optional functions to compute rhostar and ustar (not necessary)
   FUNCTION ustar(pstar) RESULT(vv)
      IMPLICIT NONE
      REAL(KIND=NUMBER), INTENT(IN) :: pstar
      REAL(KIND=NUMBER)             :: vv
      vv = 0.5d0*(ul + f(pstar, pl, capAl, capBl, capCl, expol) + ur + f(pstar, pr, capAr, capBr, capCr, expor))
   END FUNCTION ustar

   FUNCTION rhostar(pstar, rhoz, pz, gammaz) RESULT(vv)
      IMPLICIT NONE
      REAL(KIND=NUMBER), INTENT(IN) :: pstar, rhoz, pz, gammaz
      REAL(KIND=NUMBER)             :: vv, denom
      IF (pstar <= pz) THEN
         denom = b_covolume*rhoz + (1.d0 - b_covolume*rhoz)*((pstar+p_infty)/(pz+p_infty))**(1.d0/gammaz)
         vv = rhoz/denom
      ELSE
         vv = rhoz*(pstar/pz + (gammaz - 1)/(gammaz + 1))/ &
              (((gammaz - 1.d0 + 2.d0*b_covolume*rhoz)*pstar)/((gammaz + 1.d0)*pz) &
               + (gammaz + 1.d0 - 2.d0*b_covolume*rhoz)/(gammaz + 1.d0))
      END IF
   END FUNCTION rhostar

END MODULE arbitrary_eos_lagrangian_lambda_module
