PROGRAM potential_ex
USE gausLeg
IMPLICIT NONE
integer, parameter :: n=100, Nk=100
real(kind=8) :: a, b, c, small
real(kind=8), DIMENSION(0:Nk,0:Nk) :: v 
real(kind=8), DIMENSION(1:n) :: Delta1, Delta2, Delta3
real(kind=8) :: hom, vs, v0, mu_v, lambda, dk, k1, k2, S, r1, r2, r3, jvj1, jvj2, jvj3
real(kind=8), DIMENSION(1:n) :: x1, w1, x2, w2, x3, w3
integer :: ik1, ik2, i
mu_v = 0.799591851143533d0
lambda = 1.93672d0
hom = 41.44252d0
vs = (-1)*lambda*(lambda-1)*hom*mu_v**2
v0 = 3*vs
dk=0.1d0
small = 0.00001d0
a = 1d0   !this is where 1/cosh^2(mu*r) becomes half 
b = 10d0  !the other two are found with trial and error: looking for convergence & coparing with Mathematica's NIntegrate
c = 10000d0

CALL gauleg(0d0, a, x1, w1, n)
CALL gauleg(a, b, x2, w2, n)
CALL gauleg(b, c, x3, w3, n)

OPEN(UNIT=11, FILE="v0k.dat", STATUS="REPLACE")
OPEN(UNIT=12, FILE="v1k.dat", STATUS="REPLACE")
OPEN(UNIT=13, FILE="vkk.dat", STATUS="REPLACE")

DO ik1=0,Nk
        k1 = ik1*dk
        DO ik2=0,Nk
                k2=ik2*dk
                S = 0d0
                DO i=1,n
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
                v(ik1,ik2) = S
                write(*,*) k1, k2, v(ik1,ik2)
        END DO
        write(11,*) 0, k1, v(0,ik1)
        write(12,*) k1, 20*dk, v(ik1,20)
        write(13,*) k1, k1, v(ik1,ik1)
END DO
CLOSE(UNIT=11,STATUS="KEEP")
CLOSE(UNIT=12,STATUS="KEEP")
CLOSE(UNIT=13,STATUS="KEEP")

END PROGRAM potential_ex
