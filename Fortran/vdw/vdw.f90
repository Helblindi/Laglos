MODULE vdw
  PUBLIC :: rho_v_p_vdw, initialize_vdw, pressure_vdw, f_add, f_test
  REAL(KIND=8), PUBLIC :: avdw, bvdw, gamma_vdw
  PUBLIC :: rhominus, c, vZ
  PRIVATE
  REAL(KIND=8) :: rho_plus
  REAL(KIND=8) :: rhoL, vL, pL, SL, xL, rhoR, vR, pR, SR, xR
  REAL(KIND=8) :: rho_minus, p_minus, p_plus, v_minus, v_plus
CONTAINS
  SUBROUTINE f_test()
    use iso_c_binding
   !  write(*,*) "testing in fortran"
  END SUBROUTINE f_test

  SUBROUTINE f_add(xxx,xxy,xxz) 
    use iso_c_binding
    real(c_double), INTENT(in) :: xxx, xxy
    real(c_double), INTENT(out) :: xxz
   !  write(*,*) "adding in fortran"
    xxz = xxx + xxy
  END SUBROUTINE f_add

  SUBROUTINE initialize_vdw(rhop, in_a, in_b, in_gamma, in_rhoL, in_rhoR, out_vL, out_vR, out_pL, out_pR)
   !  IMPLICIT NONE
    use iso_c_binding
    real(KIND=8), INTENT(in) :: rhop
    real(KIND=8), INTENT(in) :: in_a, in_b, in_gamma !in_data
    real(KIND=8), INTENT(in) :: in_rhoL, in_rhoR !in_state
    real(KIND=8), INTENT(out) :: out_vL, out_vR, out_pL, out_pR !out_state
   !  real(c_double), INTENT(in) :: in_state(*), in_data(*)
   !  real(c_double), INTENT(out) :: out_state(*)
   !  write(*,*) "initialize_vdw call"
    rho_plus = rhop
    avdw      = in_a
    bvdw      = in_b
    gamma_vdw = in_gamma

    rhoL = in_rhoL
    rhoR = in_rhoR

    rho_minus=rhominus(rho_plus)
    p_minus = p_pm(rho_minus)
    p_plus  = p_pm(rho_plus)
    SL = (p_minus+avdw*rho_minus**2)*(1/rho_minus-bvdw)**(gamma_vdw)
    SR = (p_plus +avdw*rho_plus**2) *(1/rho_plus -bvdw)**(gamma_vdw)
    pL = SL*(rhoL/(1-bvdw*rhoL))**gamma_vdw - avdw*rhoL**2
    pR = SR*(rhoR/(1-bvdw*rhoR))**gamma_vdw - avdw*rhoR**2
    v_minus = -c(rho_minus,SL)
    v_plus  = -c(rho_plus,SR)
    VL = vZ(v_minus, rho_minus, in_rhoL, SL)
    VR = vZ(v_plus, rho_plus, in_rhoR, SR)
    xL = vL + c(rhoL,SL)
    xR = vR + c(rhoR,SR)
    out_vL = vL
    out_vR = vR
    out_pL = pL
    out_pR = pR
   !  WRITE(*,*) 'avdw, bvdw, gamma_vdw', avdw, bvdw, gamma_vdw
   !  WRITE(*,*) 'rhominus, SL, SR', rho_minus, SL, SR
   !  WRITE(*,*) 'p_minus, p_plus, v_minus, v_plus', p_minus, p_plus, v_minus, v_plus
   !  WRITE(*,*) 'XL, XR', XL, XR
   !  WRITE(*,*) 'rhoL, rhoR', rhoL, rhoR
   !  WRITE(*,*) 'vL, vR', VL, VR
   !  WRITE(*,*) 'pL, pR', pL, pR
  END SUBROUTINE initialize_vdw
  
  FUNCTION pressure_vdw(un) RESULT(vv)
    IMPLICIT NONE 
    REAL(KIND=8), DIMENSION(:,:),       INTENT(IN) :: un
    REAL(KIND=8), DIMENSION(SIZE(un,2))            :: vv
    vv = (gamma_vdw-1.d0)*(un(3,:)-0.5d0*(un(2,:)**2)/un(1,:) &
         + avdw*un(1,:)**2)/(1-bvdw*un(1,:)) - avdw*un(1,:)**2
  END FUNCTION pressure_vdw
  
  FUNCTION rhominus(rho_plus) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND=8) :: rho_plus, vv
    vv = (1.d0-0.5d0*gamma_vdw)/bvdw - rho_plus
  END FUNCTION rhominus
  
  FUNCTION p_pm(rho) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND=8) :: rho, vv
    vv = (1-bvdw*rho)*(rho*(gamma_vdw-2+2*rho)**2 &
       +4*avdw*(gamma_vdw+1)*rho**2)-2*avdw*gamma_vdw*(gamma_vdw+1)*rho**2
    vv = vv/(2*gamma_vdw*(gamma_vdw+1))
  END FUNCTION p_pm

  FUNCTION vZ(v0, rho0, rhoZ, SZ)  RESULT(vv)
    IMPLICIT NONE
    REAL(KIND=8) :: v0, rho0, rhoZ, SZ, vv
    REAL(KIND=8) :: k1, k2, k3, k4, dx, rho
    INTEGER :: n, nmax= 5000
    dx = (rhoZ-rho0)/nmax
    rho = rho0
    vv  = v0
    DO n = 1, nmax     
       k1 = c(rho,Sz)/rho
       k2 = c(rho+dx/2,Sz)/(rho+dx/2)
       k3 = c(rho+dx/2,Sz)/(rho+dx/2)
       k4 = c(rho+dx,Sz)/(rho+dx)
       vv = vv + (dx/6)*(k1+2*k2+2*k3+k4)
       rho = rho+dx
    END DO
  END FUNCTION vZ

  FUNCTION C(rho,S) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND=8) :: rho, S, vv
    vv  = SQRT(gamma_vdW*S*rho**(gamma_vdw-1)/(1-bvdw*rho)**(gamma_vdw+1)-2*avdw*rho)
  END FUNCTION C

  FUNCTION rk_rho(rho,S) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND=8), INTENT(IN) :: rho, S
    REAL(KIND=8)             :: vv
    REAL(KIND=8) :: num, denom
    num=S*gamma_vdw*rho**(gamma_vdw-1)*(gamma_vdw+1)*(1-bvdw*rho)**(-2-gamma_vdw)-6*avdw*rho
    denom=(1-bvdw*rho)**(-1-gamma_vdw)*(S*gamma_vdw*rho**gamma_vdw-2*avdw*rho**2*(1-bvdw*rho)**(gamma_vdw+1))*rho
    vv = 2*SQRT(denom)/num
  END FUNCTION rk_rho

  SUBROUTINE sol_rho(rhoz,Sz,phi,xx,rho)
    IMPLICIT NONE
    INTERFACE
       FUNCTION phi(r,S) RESULT(vv)
         REAL(KIND=8), INTENT(IN) :: r, S
         REAL(KIND=8)             :: vv
       END FUNCTION phi
    END INTERFACE
    REAL(KIND=8), DIMENSION(:)  :: xx, rho
    REAL(KIND=8)                :: rhoz, Sz
    REAL(KIND=8)                :: dx, k1, k2, k3, k4
    INTEGER                     :: n

    dx = xx(1)-0
    k1 = phi(rhoz,Sz)
    k2 = phi(rhoz+dx*k1/2,Sz)
    k3 = phi(rhoz+dx*k2/2,Sz)
    k4 = phi(rhoz+dx*k3,Sz)
    rho(1) = rhoz + (dx/6)*(k1+2*k2+2*k3+k4)
    
    DO n = 1, SIZE(xx)-1
       dx = xx(n+1)-xx(n)
       k1 = phi(rho(n),Sz)
       k2 = phi(rho(n)+dx*k1/2,Sz)
       k3 = phi(rho(n)+dx*k2/2,Sz)
       k4 = phi(rho(n)+dx*k3,Sz)
       rho(n+1) = rho(n) + (dx/6)*(k1+2*k2+2*k3+k4)
    END DO
  END SUBROUTINE sol_rho

  SUBROUTINE rho_v_p_vdw(x0,nmax,xx,rho,v,p) 
    use iso_c_binding
    REAL(c_double), INTENT(in)  :: xx(*)
    REAL(c_double), INTENT(out) :: rho(*), v(*), p(*)
    REAL(c_double), INTENT(in)  :: x0
    INTEGER(c_int), INTENT(in)  :: nmax
    INTEGER                     :: n, n0, nL, nR

    DO n = 1, nmax
       IF (xx(n).GT.xL) THEN
          nL = n
          EXIT
       END IF
    END DO
    DO n = 1, nmax
       IF (xx(n).GT.x0) THEN
          n0 = n
          EXIT
       END IF
    END DO
    DO n = nmax, 1, -1
       IF (xx(n).LT.xR) THEN
          nR = n
          EXIT
       END IF
    END DO
    rho(1:nL-1) = rhoL
    v(1:nL-1)   = vL
    p(1:nL-1)   = pL
    IF (n0-nL>0) THEN
       CALL sol_rho(rho_minus,SL,rk_rho,xx(n0-1:nL:-1),rho(n0-1:nL:-1))
       DO n = nL, n0-1
          v(n) = xx(n) - c(rho(n),SL)
          p(n) = SL*(rho(n)/(1-bvdw*rho(n)))**gamma_vdw-avdw*rho(n)**2
       END DO
    END IF
    rho(nR+1:nmax) = rhoR
    v(nR+1:nmax)   = vR
    p(nR+1:nmax)   = pR
    IF (nR-n0+1>0) THEN
       CALL sol_rho(rho_plus,SR,rk_rho,xx(n0:nR),rho(n0:nR))
       DO n = n0, nR
          v(n) = xx(n) - c(rho(n),SR)
          p(n) = SR*(rho(n)/(1-bvdw*rho(n)))**gamma_vdw-avdw*rho(n)**2
       END DO
    END IF
  END SUBROUTINE rho_v_p_vdw
END MODULE vdw


