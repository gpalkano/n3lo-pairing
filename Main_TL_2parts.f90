PROGRAM Main_TL_2parts

USE constants
USE potentials
USE BCS_eq

IMPLICIT NONE
integer, parameter :: n=500, NfromFile=10
real (kind=8) , ALLOCATABLE :: D1(:), D2(:), D3(:)
real (kind=8) , ALLOCATABLE :: Delta_Final(:), Exc(:)
!real (kind=8) , DIMENSION(1:F) :: C_Matr               !for debugging use        
real(kind=8) :: kf, kfmax, df,Density, mu_s, M1, M2, rho, step, a, b, c, TOL_for_Ridder, TOL_for_Iter, ef
real(kind=8) :: K, xi, minE, Energy, Kin, Pot, Gap, Gap_Benchmark, mu_prev, v2k, u2k, eta
integer :: i, j, imin, read_potential, ikf, nkf
real(kind=8), ALLOCATABLE:: x1(:), x2(:), x3(:), w1(:), w2(:), w3(:)

read(5,*) eta
read(5,*) mu_prev       !This is the initial guess for mu, e.g., the solution for kf*a=-10
read(5,*) read_potential

ALLOCATE(D1(1:n),D2(1:n),D3(1:n))
ALLOCATE(Delta_Final(0:F), exc(0:F))

kfmax = eta
kfmin = 1d0
nkf = 5
a=5d0
b=10d0
c=100d0
if (read_potential.eq.1) then
        CALL read_potential_matrix(a,b,c,n,kf)
else if (read_potential.eq.0) then
        CALL make_potential_matrix(a,b,c,n,kf)
END IF
do ikf=1,nkf
        kf = kfmin + (kfmax-kfmin)*ikf/nkf
        ef = 0.5d0 * hom * kf**2d0
        M1 = 0.9d0 * ef
        M2 = 1.1d0 * ef
        rho = ((kf)**(3d0)) / (3d0 * pi**2d0)      
        step = DBLE(FAR * kf / F)
        
        TOL_for_Ridder = 10**(-4d0)
        TOL_for_Iter = 10**(-12d0)
        
        
        write(*,*) ''
        write(*,*) "point", ikf
        WRITE(*,*) "kf:", kf, "v(kf,kf):", vf
        WRITE(*,*) "rho:", rho
        WRITE(*,*) "ef:", ef, '(and initial guess for mu)'
        WRITE(*,*) "a", a
        WRITE(*,*) "b", b
        WRITE(*,*) "c", c
        
        CALL remake_potential_matrix(a,b,c,n,kf)
        
        !Solve the BCS gap equations
        CALL Root_OF_Gap(Density, D1, D2, D3, df, mu_s, M1, M2, rho, n, step, a, b, c,&
                TOL_for_Iter, TOL_for_Ridder)
        WRITE(*,*) 'kf, df, mu_s', kf, df, mu_s
        
        !!print results
        !OPEN(UNIT=1,FILE="k_gap.dat",STATUS="REPLACE")
        !OPEN(UNIT=2,FILE="k_k2_exc.dat",STATUS="REPLACE")
        !OPEN(UNIT=3,FILE="k_v2k_vkuk_vkouk.dat",STATUS="REPLACE")
        !
        !minE = 100d0 !something large
        !DO i=0,F
        !        K = DBLE(i*step)
        !        xi = 0.5 * hom * K**2d0 - mu_s
        !        Exc(i) = sqrt(xi**2d0 + Delta_Final(i)**2d0)
        !        IF(Exc(i)<minE) THEN
        !                minE = Exc(i)
        !                imin = i
        !        END IF
        !        WRITE(1,*) K, Delta_Final(i)
        !        WRITE(2,*) K, K**2, exc(i)
        !        v2k = 0.5*(1-xi/exc(i))
        !        u2k = dabs(1d0-v2k)
        !        WRITE(3,*) K, v2k, sqrt(v2k*u2k), sqrt(v2k/u2k)
        !END DO
        !Gap = minE
        !
        !CLOSE(UNIT=1)
        !CLOSE(UNIT=2)
        !CLOSE(UNIT=3)
        !
        !
        !WRITE(*,*) "D @ minExc =", minE,'@',imin, "where D = ", Delta_Final(imin), 'and k =', imin * step
        !WRITE(*,*) "mu_s =", mu_s
        !WRITE(*,*) "Benchmarking: Delta =,", Gap_Benchmark
        !
        !CALL Evaluation_Of_the_Energy(Energy, Kin, Pot, D1, D2, D3, mu_s, n, step, a, b, c)
        !!This returns energy density, i.e., E/V
        !
        !write(*,*) 'mu =', mu_s
        !WRITE(*,*) "Gap =", Gap
        !write(*,*) "**"
        !write(*,*) 'mu_s/ef =', mu_s/ef
        !write(*,*) 'D/Ef = ', Gap/ef
        !WRITE(*,*) "Energy / N =", Energy / rho, "(K/N, V/N =,", Kin/rho, Pot/rho, ")"      
        !WRITE(*,*) "Energy / N /3ef/5 =", Energy / rho / (3d0 * ef / 5d0)      
        !WRITE(*,*) "Energy / N /3ef/5 / mu / ef =", Energy / rho / (3d0 * ef / 5d0) / (mu_s / ef)     
end do
END PROGRAM Main_TL_2parts
