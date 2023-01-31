PROGRAM makepot

USE constants
USE potentials
USE BCS_eq

IMPLICIT NONE
integer, parameter :: n=500
!real (kind=8) , DIMENSION(1:F) :: C_Matr               !for debugging use        
real(kind=8) :: kf, Density, mu_g, mu_prev, rho, step, a, b, c, TOL_for_Ridder, TOL_for_Iter, ef
real(kind=8) :: mu_0, delmu, eta, df0in, dfout, kf_now, ef_now
integer :: mindex, j, imin, read_potential

read(5,*) eta
read(5,*) mu_prev       !This is the initial guess for mu, e.g., the solution for kf*a=-10
read(5,*) read_potential

kf = eta
ef = 0.5d0 * hom * kf**2d0
rho = ((kf)**(3d0)) / (2d0*3d0 * pi**2d0)      !extra factor of 1/2 due to polarization
a = 50d0 * kf/2d0
b = 200d0 * kf/2d0
c = 100000d0 * kf/2d0
delmu=0.5d0

TOL_for_Ridder = 10**(-4d0)
TOL_for_Iter = 10**(-12d0)

WRITE(*,*) "kf", kf
WRITE(*,*) "ef", ef
WRITE(*,*) "rho", rho
WRITE(*,*) "step", step
WRITE(*,*) "a", a
WRITE(*,*) "b", b
WRITE(*,*) "c", c
WRITE(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
WRITE(*,*) "We are solving for kf =", kf, "and density rho =", rho
WRITE(*,*) "That gives kf*a =", kf * a_par 

!read or make the potential matrix elements
if (read_potential.eq.1) then
        CALL read_potential_matrix(a,b,c,n,kf)
else if (read_potential.eq.0) then
        CALL make_potential_matrix(a,b,c,n,kf)
END IF
END PROGRAM makepot
