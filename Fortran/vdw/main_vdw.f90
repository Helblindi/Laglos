PROGRAM main
  USE vdw
  IMPLICIT NONE
  INTEGER, PARAMETER :: nmax = 1000
  REAL(KIND=8), DIMENSION(nmax) :: xx, rho, v, p
  REAL(KIND=8) :: rho_plus
  REAL(KIND=8) :: xinit, dx, long=.15d0
  REAL(KIND=8) :: in_rhoL, in_rhoR, in_a, in_b, in_gamma
  REAL(KIND=8) :: out_pL, out_pR, out_vL, out_vR
  INTEGER :: n
  rho_plus = 0.35d0
  in_rhoL = 0.10d0
  in_rhoR = 0.39d0
  in_a = 1.d0
  in_b = 1.d0
  in_gamma = 1.02d0  
  CALL initialize_vdw(rho_plus, in_a, in_b, in_gamma, in_rhoL, in_rhoR, out_vL, out_vR, out_pL, out_pR)
  xinit = -long/2.d0
  dx = long/(nmax-1)
  DO n = 1, nmax
     xx(n) = xinit+(n-1)*dx
  END DO
  CALL rho_v_p_vdw(0.d0,nmax,xx,rho,v,p)
  DO n = 1, nmax
     WRITE(10,*) xx(n), rho(n)
     WRITE(11,*) xx(n), v(n)
     WRITE(12,*) xx(n), p(n)
  END DO
END PROGRAM main
