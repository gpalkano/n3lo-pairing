PROGRAM Main_TL_2parts

USE constants
USE potentials
USE BCS_eq

IMPLICIT NONE
integer, parameter :: n=1000
!real (kind=8) , DIMENSION(1:F) :: C_Matr               !for debugging use        
real(kind=8) :: kf, Density, mu_g, mu_prev, rho, step, a, b, c, TOL_for_Ridder, TOL_for_Iter, ef
real(kind=8) :: mu_0, delmu, eta, mineq, dfout, kf_now, ef_now
integer :: mindex, j, imin, read_potential
real(kind=8) :: mus(10)

read(5,*) eta
read(5,*) mu_prev       !This is the initial guess for mu, e.g., the solution for kf*a=-10
read(5,*) read_potential

kf = eta
ef = 0.5d0 * hom * kf**2d0
rho = ((kf)**(3d0)) / (3d0 * pi**2d0)      !extra factor of 1/2 due to polarization
!a = 50d0 * kf/2d0
!b = 200d0 * kf/2d0
!c = 100000d0 * kf/2d0
a=5d0
b=10d0
c=100d0
delmu=0.4d0

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

mus = (/ 20d0, 25d0, 30d0, 35d0, 40d0, 45d0, 50d0, 55d0, 0.1d0, 0.2d0 /)
mu_0=mu_prev
OPEN(UNIT=1,FILE="new_df_dfIef_mu_muIef_rho_kf_Ikf_mineq.dat",STATUS='REPLACE')
do mindex=1,10
        !mu_g=mu_0 + delmu*(mindex-1)
        mu_g = mus(mindex)
        write(*,*) 'distance of k0 from k*', kf-sqrt(2d0*mu_g/hom)
        CALL EVALUATION_TL_KHODEL(Density, mu_g, dfout, mineq, n, step, a, b, c, TOL_for_Ridder, TOL_for_Iter)
        kf_now=(Density*3d0*Pi**2d0)**(1d0/3d0)
        ef_now = 0.5 * hom * kf_now**2d0
        write(*,*) mu_g, Density, "grab", dfout, kf_now
        write(1,*) dfout, dfout/ef_now, mu_g, mu_g/ef_now, Density, kf_now, 1d0/kf_now, mineq
end do
CLOSE(UNIT=1,STATUS='KEEP')
END PROGRAM Main_TL_2parts
