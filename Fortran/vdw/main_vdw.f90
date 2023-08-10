PROGRAM main
  USE vdw
  IMPLICIT NONE
  INTEGER, PARAMETER :: nmax = 1000
  REAL(KIND=8), DIMENSION(nmax) :: xx, rho, v, p
  REAL(KIND=8) :: rho_plus
  REAL(KIND=8) :: xinit, dx, long=.15d0
  REAL(KIND=8) :: in_state(2), in_data(3), out_state(4)
  INTEGER :: n
  rho_plus = 0.35d0
  in_state(1) = 0.10d0
  in_state(2) = 0.39d0
  in_data(1) = 1.d0
  in_data(2) = 1.d0
  in_data(3) = 1.02d0  
  CALL initialize_vdw(rho_plus, in_state, in_data, out_state)
  xinit = -long/2.d0
  dx = long/(nmax-1)
  DO n = 1, nmax
     xx(n) = xinit+(n-1)*dx
  END DO
  CALL rho_v_p_vdw(0.d0,xx,rho,v,p)
  DO n = 1, nmax
     WRITE(10,*) xx(n), rho(n)
     WRITE(11,*) xx(n), v(n)
     WRITE(12,*) xx(n), p(n)
  END DO
END PROGRAM main
