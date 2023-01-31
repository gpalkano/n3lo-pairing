MODULE  potentials
USE gausLeg
USE constants
use vn3lo
IMPLICIT NONE
integer, parameter :: npwave = 500
REAL (kind=8) , SAVE ::  small=10d-4, xx


CONTAINS

!----------------------------------------------------------------------------------------------------
SUBROUTINE setup_potential()
IMPLICIT NONE
real(kind=8) :: xin

read(5,*) xin
xx=xin
END SUBROUTINE setup_potential
!----------------------------------------------------------------------------------------------------
real(kind=8) function vn3lo_1s0(k1,k2)
real(kind=8), intent(in) :: k1, k2
vn3lo_1s0=n3lo_1s0(k1,k2)
end function vn3lo_1s0
!------------------------------------------------------------------------
REAL (kind=8)  FUNCTION V_swave(k1,k2)
! V=lam*(lam-1)*hom*mu**2/(cosh(mu*R)**2  -- Poschl-Teller Potential
IMPLICIT NONE
real (kind=8) , INTENT(IN) :: k1,k2
REAL (kind=8) , PARAMETER :: PI = acos(-1.0), lam = 1.93672d0, mu = 0.79591851143533d0, hom = 41.44252, small=1d-4
real (kind=8) :: beta, k,A
A=-lam*(lam-1)*hom*mu**2
IF (abs(k1-k2)<small) THEN 
  IF (k1==0) THEN
  V_swave=A*pi**2/(12*mu**3)
  ELSE
  V_swave=A*(mu-k1*pi/sinh(k1*pi/mu))/(2*mu**2*k1**2)
  END IF
ELSE
  IF (k1*k2==0) THEN
  k=k1+k2
  beta=pi*k/(2*mu)
  V_swave=A*pi*(beta/tanh(beta)-1)/(2*mu**2*k*sinh(beta))
  ELSE
  V_swave=A*pi*(abs(k1-k2)/sinh(abs(k1-k2)*pi/(2*mu)) - (k1+k2)/sinh((k1+k2)*pi/(2*mu)))/(4*mu**2*k1*k2)
  END IF
END IF
END FUNCTION V_swave
!-------------------------------------------------------------------------------------
real (kind=8) FUNCTION V_pwave(k1,k2)
IMPLICIT NONE
real(kind=8), INTENT(IN) :: k1, k2
real(kind=8) :: a, b, c, small
real(kind=8) :: hom, vs, v0, mu_v, lambda, S, r1, r2, r3, jvj1, jvj2, jvj3, vstef
real(kind=8), DIMENSION(1:npwave) :: x1, w1, x2, w2, x3, w3
integer :: ik1, ik2, i

mu_v = 0.799591851143533d0
lambda = 1.93672d0
hom = 41.44252d0
vstef = 3.9911608d0
vs = (-1)*lambda*(lambda-1)*hom*mu_v**2
!vs = (-1)*vstef*hom*mu_v**2
v0 = xx*vs

small = 0.00001d0
a = 1d0   !this is where 1/cosh^2(mu*r) becomes half 
b = 10d0  !the other two are found with trial and error: looking for convergence & coparing with Mathematica's NIntegrate
c = 10000d0

CALL gauleg(0d0, a, x1, w1, npwave)
CALL gauleg(a, b, x2, w2, npwave)
CALL gauleg(b, c, x3, w3, npwave)

S = 0d0
DO i=1,npwave
        r1 = x1(i)
        r2 = x2(i)
        r3 = x3(i)
       
        IF  ((k1.le.small).or.(k2.le.small))  THEN   
        !Do not change this test to k1*k2.le.small because terms with 
        !small but finite k1 and k2 might be redirected here
                jvj1 = 0d0
                jvj2 = 0d0
                jvj3 = 0d0
        ELSE
                jvj1 = r1**2 * (sin(k1*r1)/(k1*r1) - cos(k1*r1))/(k1*r1) * v0/(cosh(mu_v*r1)**2)* &
                        (sin(k2*r1)/(k2*r1) - cos(k2*r1))/(k2*r1)
                jvj2 = r2**2 * (sin(k1*r2)/(k1*r2) - cos(k1*r2))/(k1*r2) * v0/(cosh(mu_v*r2)**2)* &
                        (sin(k2*r2)/(k2*r2) - cos(k2*r2))/(k2*r2)
                jvj3 = r1**2 * (sin(k1*r3)/(k1*r3) - cos(k1*r3))/(k1*r3) * v0/(cosh(mu_v*r3)**2)* &
                        (sin(k2*r3)/(k2*r3) - cos(k2*r3))/(k2*r3)
        END IF
        S = S + w1(i) * jvj1 + w2(i) * jvj2 + w3(i) * jvj3
END DO
V_pwave= S

END FUNCTION V_pwave

END MODULE
