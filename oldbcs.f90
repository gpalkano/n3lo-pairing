MODULE BCS_eq
USE potentials
!USE Polynomials
!USE Math_Functions
IMPLICIT NONE
REAL (kind=8) , PARAMETER :: hom=41.44252d0
REAL (kind=8) :: a_par=-18.5d0, re_3D=2.7d0
REAL (kind=8), ALLOCATABLE :: V11(:,:), V12(:,:), V13(:,:), V22(:,:), V23(:,:), V33(:,:)
!REAL (kind=8) :: k_F=(3*rho_0*(PI**2d0))**(1d0/3d0), e_F=0.5*hom*(3*rho_0*PI**2)**(2.0/3.0), a_par=-18.5d0, re_3D=2.7d0
INTEGER  ::  F=1000, FAR=10

CONTAINS
SUBROUTINE Initialization()
END SUBROUTINE Initialization
!-------------------------------------------------------------------------
SUBROUTINE make_potential_matrix(a,b,c,n)
!Calculate and store the potential matrix elements at the integration points
!of Gauss-Legendre. It is Vab(ij) = Vba(j,i)
IMPLICIT NONE
integer, INTENT(IN) :: n
real(kind=8), INTENT(IN) :: a, b, c
real(kind=8) :: xin
real(kind=8), DIMENSION(1:n) :: x1, w1, x2, w2, x3, w3
integer i, j

CALL setup_potential()
write(*,*) "For potential:"
write(*,*) "x,a,b,c,n"
write(*,*) x,a,b,c,n

CALL gauleg(0d0, a, x1, w1, n)
CALL gauleg(a, b, x2, w2, n)
CALL gauleg(b, c, x3, w3, n)

ALLOCATE(V11(1:n,1:n), V12(1:n,1:n), V13(1:n,1:n))
ALLOCATE(V22(1:n,1:n), V23(1:n,1:n))
ALLOCATE(V33(1:n,1:n))

write(11,*) "x, a, b, c, n"
write(11,*) x,a,b,c,n
write(12,*) "x, a, b, c, n"
write(12,*) x,a,b,c,n
write(13,*) "x, a, b, c, n"
write(13,*) x,a,b,c,n
write(22,*) "x, a, b, c, n"
write(22,*) x,a,b,c,n
write(23,*) "x, a, b, c, n"
write(23,*) x,a,b,c,n
write(33,*) "x, a, b, c, n"
write(33,*) x,a,b,c,n

DO i=1,n
        write(*,*) i
        DO j=1,n
                V11(i,j) = V_swave(x1(i),x1(j)) 
                V12(i,j) = V_swave(x1(i),x2(j)) 
                V13(i,j) = V_swave(x1(i),x3(j)) 
                V22(i,j) = V_swave(x2(i),x2(j)) 
                V23(i,j) = V_swave(x2(i),x3(j)) 
                V33(i,j) = V_swave(x3(i),x3(j))

                write(11,*) i, j, V11(i,j)
                write(12,*) i, j, V12(i,j)
                write(13,*) i, j, V13(i,j)
                write(22,*) i, j, V22(i,j)
                write(23,*) i, j, V23(i,j)
                write(33,*) i, j, V33(i,j)
        END DO
END DO

END SUBROUTINE make_potential_matrix
!------------------------------------------------------------
SUBROUTINE read_potential_matrix(a,b,c,n)
!Read and store the potential matrix elements at the integration points
!of Gauss-Legendre. It is Vab(ij) = Vba(j,i)
!This is reading from the V*.dat files and assuming that they were printed
!by make_potential_matrix
IMPLICIT NONE
integer, INTENT(IN) :: n
real(kind=8), INTENT(IN) :: a, b, c
real(kind=8), DIMENSION(1:n) :: x1, w1, x2, w2, x3, w3
real(kind=8) :: x0, xin, ain, bin, cin, nin 
integer i, j, bucketi, bucketj
character(len=50) :: lab
logical :: cond

CALL setup_potential()  !this reads v0/vs: the potential strength 
x0=x                    !for comparison with the read potential
                        !(it cannot be altered if reading the potential)

ALLOCATE(V11(1:n,1:n), V12(1:n,1:n), V13(1:n,1:n))
ALLOCATE(V22(1:n,1:n), V23(1:n,1:n))
ALLOCATE(V33(1:n,1:n))

OPEN(UNIT=11,FILE="V11.dat",STATUS="OLD")
OPEN(UNIT=12,FILE="V12.dat",STATUS="OLD")
OPEN(UNIT=13,FILE="V13.dat",STATUS="OLD")
OPEN(UNIT=22,FILE="V22.dat",STATUS="OLD")
OPEN(UNIT=23,FILE="V23.dat",STATUS="OLD")
OPEN(UNIT=33,FILE="V33.dat",STATUS="OLD")


read(11,*) lab 
read(12,*) lab 
read(13,*) lab 
read(22,*) lab 
read(23,*) lab 
read(33,*) lab 

cond=.FALSE.
read(11,*) xin, ain, bin, cin, nin
cond = cond .OR. dabs(x0-xin).ge.small 
cond = cond .OR. dabs(ain-a).ge.small .OR. dabs(bin-b).ge.small .OR. dabs(cin-c).ge.small
cond = cond .OR. nin.ne.n
read(12,*) xin, ain, bin, cin, nin 
cond = cond .OR. dabs(x0-xin).ge.small 
cond = cond .OR. dabs(ain-a).ge.small .OR. dabs(bin-b).ge.small .OR. dabs(cin-c).ge.small
cond = cond .OR. nin.ne.n
read(13,*) xin, ain, bin, cin, nin 
cond = cond .OR. dabs(x0-xin).ge.small 
cond = cond .OR. dabs(ain-a).ge.small .OR. dabs(bin-b).ge.small .OR. dabs(cin-c).ge.small
cond = cond .OR. nin.ne.n
read(22,*) xin, ain, bin, cin, nin 
cond = cond .OR. dabs(x0-xin).ge.small 
cond = cond .OR. dabs(ain-a).ge.small .OR. dabs(bin-b).ge.small .OR. dabs(cin-c).ge.small
cond = cond .OR. nin.ne.n
read(23,*) xin, ain, bin, cin, nin 
cond = cond .OR. dabs(x0-xin).ge.small 
cond = cond .OR. dabs(ain-a).ge.small .OR. dabs(bin-b).ge.small .OR. dabs(cin-c).ge.small
cond = cond .OR. nin.ne.n
read(33,*) xin, ain, bin, cin, nin 
cond = cond .OR. dabs(x0-xin).ge.small 
cond = cond .OR. dabs(ain-a).ge.small .OR. dabs(bin-b).ge.small .OR. dabs(cin-c).ge.small
cond = cond .OR. nin.ne.n

IF (cond) THEN
        write(*,*) "ATTENTION: you read the wrong potential"
        PAUSE
ENDIF 

write(*,*) "For the potential:"
write(*,*) "x,a,b,c,n"
write(*,*) x,a,b,c,n


DO i=1,n
        write(*,*) i
        DO j=1,n
                read(11,*) bucketi, bucketj, V11(i,j)
                read(12,*) bucketi, bucketj, V12(i,j)
                read(13,*) bucketi, bucketj, V13(i,j)
                read(22,*) bucketi, bucketj, V22(i,j)
                read(23,*) bucketi, bucketj, V23(i,j)
                read(33,*) bucketi, bucketj, V33(i,j)
        END DO
END DO

CLOSE(UNIT=11,STATUS="KEEP")
CLOSE(UNIT=12,STATUS="KEEP")
CLOSE(UNIT=13,STATUS="KEEP")
CLOSE(UNIT=22,STATUS="KEEP")
CLOSE(UNIT=23,STATUS="KEEP")
CLOSE(UNIT=33,STATUS="KEEP")

END SUBROUTINE read_potential_matrix
!------------------------------------------------------------------------
!SUBROUTINE make_potential_matrix(a,b,c,n)
!!Calculate and store the potential matrix elements at the integration points
!!of Gauss-Legendre. It is Vab(ij) = Vba(j,i)
!IMPLICIT NONE
!integer, INTENT(IN) :: n
!real(kind=8), INTENT(IN) :: a, b, c
!real(kind=8), DIMENSION(1:n) :: x1, w1, x2, w2, x3, w3
!integer i, j
!
!CALL gauleg(0d0, a, x1, w1, n)
!CALL gauleg(a, b, x2, w2, n)
!CALL gauleg(b, c, x3, w3, n)
!
!ALLOCATE(V11(1:n,1:n), V12(1:n,1:n), V13(1:n,1:n))
!ALLOCATE(V22(1:n,1:n), V23(1:n,1:n))
!ALLOCATE(V33(1:n,1:n))
!
!DO i=1,n
!        write(*,*) i
!        DO j=1,n
!                V11(i,j) = V_pwave(x1(i),x1(j)) 
!                V12(i,j) = V_pwave(x1(i),x2(j)) 
!                V13(i,j) = V_pwave(x1(i),x3(j)) 
!                V22(i,j) = V_pwave(x2(i),x2(j)) 
!                V23(i,j) = V_pwave(x2(i),x3(j)) 
!                V33(i,j) = V_pwave(x3(i),x3(j))
!
!                write(11,*) i, j, V11(i,j)
!                write(12,*) i, j, V12(i,j)
!                write(13,*) i, j, V13(i,j)
!                write(22,*) i, j, V22(i,j)
!                write(23,*) i, j, V23(i,j)
!                write(33,*) i, j, V33(i,j)
!        END DO
!END DO
!
!END SUBROUTINE make_potential_matrix
!!------------------------------------------------------------------------
!SUBROUTINE read_potential_matrix(n)
!!Read and store the potential matrix elements at the integration points
!!of Gauss-Legendre. It is Vab(ij) = Vba(j,i)
!!This is reading from the V*.dat files and assuming that they were printed
!!by make_potential_matrix
!IMPLICIT NONE
!integer, INTENT(IN) :: n
!real(kind=8), DIMENSION(1:n) :: x1, w1, x2, w2, x3, w3
!integer i, j, bucketi, bucketj
!
!ALLOCATE(V11(1:n,1:n), V12(1:n,1:n), V13(1:n,1:n))
!ALLOCATE(V22(1:n,1:n), V23(1:n,1:n))
!ALLOCATE(V33(1:n,1:n))
!
!OPEN(UNIT=11,FILE="V11.dat",STATUS="OLD")
!OPEN(UNIT=12,FILE="V12.dat",STATUS="OLD")
!OPEN(UNIT=13,FILE="V13.dat",STATUS="OLD")
!OPEN(UNIT=22,FILE="V22.dat",STATUS="OLD")
!OPEN(UNIT=23,FILE="V23.dat",STATUS="OLD")
!OPEN(UNIT=33,FILE="V33.dat",STATUS="OLD")
!
!DO i=1,n
!        write(*,*) i
!        DO j=1,n
!                read(11,*) bucketi, bucketj, V11(i,j)
!                read(12,*) bucketi, bucketj, V12(i,j)
!                read(13,*) bucketi, bucketj, V13(i,j)
!                read(22,*) bucketi, bucketj, V22(i,j)
!                read(23,*) bucketi, bucketj, V23(i,j)
!                read(33,*) bucketi, bucketj, V33(i,j)
!        END DO
!END DO
!
!CLOSE(UNIT=11,STATUS="KEEP")
!CLOSE(UNIT=12,STATUS="KEEP")
!CLOSE(UNIT=13,STATUS="KEEP")
!CLOSE(UNIT=22,STATUS="KEEP")
!CLOSE(UNIT=23,STATUS="KEEP")
!CLOSE(UNIT=33,STATUS="KEEP")
!
!END SUBROUTINE read_potential_matrix
!------------------------------------------------------------------------
!SUBROUTINE Introduction()
!!This subroutine when called will return all the relevant introductory physical quantities 
!IMPLICIT NONE
!WRITE(*,*) "================ INTRODUCTION ====================================================="
!WRITE(*,*) "We are calculating the Gap for Neutron Matter. The physical quantities of relevance are:"
!WRITE(*,*) "Density n =", rho_0, "[fm^-1]"
!WRITE(*,*) "which gives e_F =", mu_ch, "[MeV]", "and k_F =", k_F
!WRITE(*,*) "and E/N in TL 3eF/5 =", 3d0 * mu_ch / 5d0
!WRITE(*,*) "Our interaction is through the S_0 channel with" 
!WRITE(*,*) " - scattering length a =", a_par, "[fm]"
!WRITE(*,*) " - effective range r_e =", re_3D, "[fm]"
!WRITE(*,*) "And thus our dimensional parameter k_F*a =", k_F * a_par
!WRITE(*,*) "===========`~\_/~\_/~\_/~\_/~\_/~\_/~\_/~\_/~\_/~\_/~`============================="
!END SUBROUTINE Introduction
!---------------------------------------------------------------------------
SUBROUTINE Root_Of_Gap(Density, Delta_I, Delta_II, Delta_III, Delta_Final, mu, M1, M2, rho, n, step, a, b, c,&
                        TOL_for_Iter,TOL_for_Ridder)
IMPLICIT NONE
integer, INTENT(IN) :: n
real(kind=8), INTENT(IN) :: M1, M2, rho, step, a, b, c, TOL_for_Iter, TOL_for_Ridder
real(kind=8), INTENT(OUT) :: Density, mu
real(kind=8), DIMENSION(1:n),INTENT(INOUT) :: Delta_I, Delta_II, Delta_III
real(kind=8), DIMENSION(0:F),INTENT(INOUT) :: Delta_Final
real(kind=8) :: Density1, Density2, Density3, Density4, mu1, mu2, mu3, mu4, f1, f2, f3, f4
real(kind=8) :: mu4tocheck
integer :: iter

f1 = 1d0
f2 = 1d0
mu1 = M1
mu2 = M2
iter=0
DO 
        iter = iter +1
        
        write(*,*) "Calculating gap for mu1=", mu1
        CALL EVALUATION_TL(Density1, mu1, n, step, a, b, c, TOL_for_Iter)
        f1 = Density1 - rho

        write(*,*) "Calculating gap for mu2=", mu2
        CALL EVALUATION_TL(Density2, mu2, n, step, a, b, c, TOL_for_Iter)
        f2 = Density2 - rho
        
        IF(f1 * f2 < 0d0) THEN
                WRITE(*,*) "Solution bracketed with f1 =", f1, "and f2 =", f2, "after", iter, "iterations"
                EXIT
        END IF
        WRITE(*,*) "(pre-bracketing) Ridder iteration:", iter, "f1 =", f1, "and f2 =", f2, "in [", mu1, mu2, "]"
        mu1 = 0.8d0 * mu1
        mu2 = 1.2d0 * mu2

END DO

mu4tocheck = 100d0
iter = 0
DO 
        iter = iter + 1
        mu3 = (mu1 + mu2) / 2d0
        CALL EVALUATION_TL(Density3, mu3, n, step, a, b, c, TOL_for_Iter)
        f3 = Density3 - rho
        
        IF(f1 - f2 > 0d0) THEN
                mu4 = mu3 + (mu3 - mu1) * f3 / sqrt(f3**2d0 - f1 * f2)
        ELSE IF(f1 - f2 < 0d0) THEN
                mu4 = mu3 - (mu3 - mu1) * f3 / sqrt(f3**2d0 - f1 * f2)
        ELSE
                WRITE(*,*) "Location: Root_Of_Gap"
                WRITE(*,*) "Failed at conditions f1-f2<0 OR f1-f2>0"
                WRITE(*,*) "With f1 =", f1, "and f2 =", f2
        END IF

        
        CALL EVALUATION_TL(Density4, mu4, n, step, a, b, c, TOL_for_Iter)
        f4 = Density4 - rho
        WRITE(*,*) "(post-bracketing) Ridder iteration:", iter, "At rho - rho_0  =", f4
        IF(DABS(mu4 - mu4tocheck) < TOL_for_Ridder) THEN
                EXIT
        END IF

        IF(f1 * f4 < 0d0) THEN
                mu2 =mu4
                Density2 = Density4
                f2 = f4
        ELSE IF(f1 * f4 > 0d0) THEN
                mu1 =mu4
                Density1 = Density4
                f1 = f4
        ELSE
                WRITE(*,*) "Location: Root_Of_Gap"
                WRITE(*,*) "Failed at conditions f1*f4<0 OR f1*f4>0"
                WRITE(*,*) "With f1 =", f1, "and f4 =", f4
        END IF
        mu4tocheck = mu4
END DO
mu = mu4

CALL EVALUATION_TL_for_Delta(Density, Delta_I, Delta_II, Delta_III, Delta_Final, mu, n, step, a, b, c, TOL_for_Iter)
WRITE(*,*) '--->> Solution found with Density =', Density

END SUBROUTINE Root_Of_Gap
!-------------------------------------------------------------------------------------------------------
SUBROUTINE EVALUATION_TL(Density, mu_g, n, step, a, b, c, TOL_for_Iter)
IMPLICIT NONE
integer, INTENT(IN) :: n
real(kind=8), INTENT(IN) :: step, a, b, c, mu_g, TOL_for_Iter 
real(kind=8), INTENT(INOUT) :: Density
real(kind=8), DIMENSION(1:n) :: Delta1B, Delta2B, Delta3B, Delta1, Delta2, Delta3
real(kind=8), DIMENSION(0:F) :: Delta_FB, Delta_F
real(kind=8) :: er
integer :: i, iter

DO i=1,n
        Delta1B(i) = 1d0
        Delta2B(i) = 1d0
        Delta3B(i) = 1d0
        Delta1(i) = 1d0
        Delta2(i) = 1d0
        Delta3(i) = 1d0
END DO

iter = 0
DO
        iter = iter + 1 
        CALL Evaluation_Of_Gap_n(Delta1, Delta2, Delta3, Delta1B, Delta2B, Delta3B, mu_g, n, a, b, c)
        CALL Comparison_3n(Delta1, Delta2, Delta3, Delta1B, Delta2B, Delta3B, n, er)
        WRITE(*,*) "           error of gap iteration ", er, Delta1(10)
        
        !CALL Evaluation_Of_Gap_F(Delta_FB, Delta1B, Delta2B, Delta3B, mu_g, n, step, a, b, c)
        !CALL Evaluation_Of_Gap_F(Delta_F, Delta1, Delta2, Delta3, mu_g, n, step, a, b, c)
        !CALL Comparison_F(Delta_F, Delta_FB, er)
        
        IF (DABS(er) < TOL_for_Iter) THEN
                EXIT
        END IF
        
        IF (iter > 100000) THEN
                WRITE(*,*) "Location: EVALUATION_TL"
                WRITE(*,*) "Time exceeded at iterations"
                WRITE(*,*) "With Delta(0) =", Delta_F(0)
                PAUSE
        END IF
        Delta1B = Delta1
        Delta2B = Delta2
        Delta3B = Delta3
                
END DO

CALL Evaluation_Of_Density(Density, Delta1, Delta2, Delta3, mu_g, n, a, b, c)
                
END SUBROUTINE EVALUATION_TL
!-------------------------------------------------------------------------------------------------
SUBROUTINE EVALUATION_TL_for_Delta(Density, Delta_I, Delta_II, Delta_III, Delta_Final, mu_g, n, step, a, b, c, TOL_for_Iter)
IMPLICIT NONE
integer, INTENT(IN) :: n
real(kind=8), INTENT(IN) :: step, a, b, c, mu_g, TOL_for_Iter 
real(kind=8), INTENT(OUT) :: Density
real(kind=8), DIMENSION(1:n), INTENT(OUT) :: Delta_I, Delta_II, Delta_III
real(kind=8), DIMENSION(0:F), INTENT(OUT) :: Delta_Final
real(kind=8), DIMENSION(1:n) :: Delta1B, Delta2B, Delta3B, Delta1, Delta2, Delta3
real(kind=8), DIMENSION(0:F) :: Delta_FB, Delta_F
real(kind=8) :: er
integer :: i

DO i=1,n
        Delta1B(i) = 1d0
        Delta2B(i) = 1d0
        Delta3B(i) = 1d0
        Delta1(i) = 1d0
        Delta2(i) = 1d0
        Delta3(i) = 1d0
END DO

DO
        CALL Evaluation_Of_Gap_n(Delta1, Delta2, Delta3, Delta1B, Delta2B, Delta3B, mu_g, n, a, b, c)
        CALL Comparison_3n(Delta1, Delta2, Delta3, Delta1B, Delta2B, Delta3B, n, er)
        !CALL Evaluation_Of_Gap_F(Delta_FB, Delta1B, Delta2B, Delta3B, mu_g, n, step, a, b, c)
        !CALL Evaluation_Of_Gap_F(Delta_F, Delta1, Delta2, Delta3, mu_g, n, step, a, b, c)
        !CALL Comparison_F(Delta_F, Delta_FB, er)
        IF (DABS(er) < TOL_for_Iter) THEN
                EXIT
        END IF
        Delta1B = Delta1
        Delta2B = Delta2
        Delta3B = Delta3
END DO

CALL Evaluation_Of_Gap_n(Delta_I, Delta_II, Delta_III, Delta1, Delta2, Delta3, mu_g, n, a, b, c)
CALL Evaluation_Of_Gap_F(Delta_Final, Delta_I, Delta_II, Delta_III, mu_g, n, step, a, b, c)
CALL Evaluation_Of_Density(Density, Delta_I, Delta_II, Delta_III, mu_g, n, a, b, c)
                
END SUBROUTINE EVALUATION_TL_for_Delta
!--------------------------------------------------------------------------
SUBROUTINE Evaluation_Of_Gap_n(Delta1, Delta2, Delta3, Delta1B, Delta2B, Delta3B, mu_g, n, a, b, c)
IMPLICIT NONE
integer, INTENT(IN) :: n
real(kind=8), INTENT(IN) :: a, b, c, mu_g
real(kind=8), DIMENSION(1:n), INTENT(IN) :: Delta1B, Delta2B, Delta3B
real(kind=8), DIMENSION(1:n), INTENT(INOUT) :: Delta1, Delta2, Delta3
real(kind=8) :: S1, S2, S3, P1, P2, P3, K1, K2, K3, xi1, xi2, xi3, VDE_1to1, VDE_1to2, VDE_1to3, VDE_2to1, VDE_2to2, VDE_2to3 
real(kind=8) :: VDE_3to1, VDE_3to2, VDE_3to3 
real(kind=8), DIMENSION(1:n) :: x1, w1, x2, w2, x3, w3
integer i, j
!CALL p_quadrature_rule ( n, x, w ) ! l_quadrature::Laguerre, p_quadrature::Legendre
CALL gauleg(0d0, a, x1, w1, n)
CALL gauleg(a, b, x2, w2, n)
CALL gauleg(b, c, x3, w3, n)
DO i=1,n
        S1 = 0d0
        S2 = 0d0
        S3 = 0d0
        !P1 = a * (x(i) + 1d0) / 2d0
        !P2 = (b - a) * (x(i) + 1d0) / 2d0 + a
        !P3 = (c - b) * (x(i) + 1d0) / 2d0 + b
        P1 = x1(i)
        P2 = x2(i)
        P3 = x3(i)
        DO j=1,n
                !K1 = a * (x(j) + 1d0) / 2d0
                !K2 = (b - a) * (x(j) + 1d0) / 2d0 + a
                !K3 = (c - b) * (x(j) + 1d0) / 2d0 + b
                K1 = x1(j)
                K2 = x2(j)
                K3 = x3(j)
                xi1 = 0.5d0 * hom * (K1**2d0) - mu_g
                xi2 = 0.5d0 * hom * (K2**2d0) - mu_g
                xi3 = 0.5d0 * hom * (K3**2d0) - mu_g
                
                VDE_1to1 = (-1d0 / Pi) * (K1**2d0) * V11(i,j) * Delta1B(j)&
                /sqrt(xi1**2d0 + Delta1B(j)**2d0)
                
                VDE_1to2 = (-1d0 / Pi) * (K1**2d0) * V12(j,i) * Delta1B(j)&
                /sqrt(xi1**2d0 + Delta1B(j)**2d0)
                
                VDE_1to3 = (-1d0 / Pi) * (K1**2d0) * V13(j,i) * Delta1B(j)&
                /sqrt(xi1**2d0 + Delta1B(j)**2d0)
                
                VDE_2to1 = (-1d0 / Pi) * (K2**2d0) * V12(i,j) * Delta2B(j)&
                /sqrt(xi2**2d0 + Delta2B(j)**2d0)
                
                VDE_2to2 = (-1d0 / Pi) * (K2**2d0) * V22(i,j) * Delta2B(j)&
                /sqrt(xi2**2d0 + Delta2B(j)**2d0)
                
                VDE_2to3 = (-1d0 / Pi) * (K2**2d0) * V23(j,i) * Delta2B(j)&
                /sqrt(xi2**2d0 + Delta2B(j)**2d0)
                
                VDE_3to1 = (-1d0 / Pi) * (K3**2d0) * V13(i,j) * Delta3B(j)&
                /sqrt(xi3**2d0 + Delta3B(j)**2d0)
                
                VDE_3to2 = (-1d0 / Pi) * (K3**2d0) * V23(i,j) * Delta3B(j)&
                /sqrt(xi3**2d0 + Delta3B(j)**2d0)
                
                VDE_3to3 = (-1d0 / Pi) * (K3**2d0) * V33(i,j) * Delta3B(j)&
                /sqrt(xi3**2d0 + Delta3B(j)**2d0)
                
                !S1 = S1 + w(j) * ( a * VDE_1to1 /2d0 + (b - a) * VDE_2to1 / 2d0 + (c - b) * VDE_3to1 / 2d0)
                !S2 = S2 + w(j) * ( a * VDE_1to2 /2d0 + (b - a) * VDE_2to2 / 2d0 + (c - b) * VDE_3to2 / 2d0)
                !S3 = S3 + w(j) * ( a * VDE_1to3 /2d0 + (b - a) * VDE_2to3 / 2d0 + (c - b) * VDE_3to3 / 2d0)
                S1 = S1 + w1(j) * VDE_1to1 + w2(j) * VDE_2to1 + w3(j) * VDE_3to1 
                S2 = S2 + w1(j) * VDE_1to2 + w2(j) * VDE_2to2 + w3(j) * VDE_3to2 
                S3 = S3 + w1(j) * VDE_1to3 + w2(j) * VDE_2to3 + w3(j) * VDE_3to3 
        END DO
        Delta1(i) = S1
        Delta2(i) = S2
        Delta3(i) = S3
END DO

END SUBROUTINE Evaluation_Of_Gap_n
!---------------------------------------------------------------------------
SUBROUTINE Evaluation_Of_Gap_F(Delta_F, Delta1B, Delta2B, Delta3B, mu_g, n, step, a, b, c)
IMPLICIT NONE
integer, INTENT(IN) :: n
real(kind=8), INTENT(IN) :: a, b, c, mu_g, step
real(kind=8), DIMENSION(1:n), INTENT(IN) :: Delta1B, Delta2B, Delta3B
real(kind=8), DIMENSION(0:F), INTENT(INOUT) :: Delta_F
real(kind=8) :: S, P, K1, K2, K3, xi1, xi2, xi3, VDE_1toF, VDE_2toF, VDE_3toF
real(kind=8), DIMENSION(1:n) :: x1, w1, x2, w2, x3, w3
integer i, j
!CALL p_quadrature_rule ( n, x, w ) ! l_quadrature::Laguerre, p_quadrature::Legendre
CALL gauleg(0d0, a, x1, w1, n)
CALL gauleg(a, b, x2, w2, n)
CALL gauleg(b, c, x3, w3, n)

DO i=0,F
        S = 0d0
        P = DBLE(i * step)
        DO j=1,n
                !K1 = a * (x(j) + 1d0) / 2d0
                !K2 = (b - a) * (x(j) + 1d0) / 2d0 + a
                !K3 = (c - b) * (x(j) + 1d0) / 2d0 + b
                K1 = x1(j)
                K2 = x2(j)
                K3 = x3(j)
                xi1 = 0.5d0 * hom * (K1**2d0) - mu_g
                xi2 = 0.5d0 * hom * (K2**2d0) - mu_g
                xi3 = 0.5d0 * hom * (K3**2d0) - mu_g
                
                VDE_1toF = (-1d0 / Pi) * (K1**2d0) * V_swave(P,K1) * Delta1B(j)&
                /sqrt(xi1**2d0 + Delta1B(j)**2d0)
                
                VDE_2toF = (-1d0 / Pi) * (K2**2d0) * V_swave(P,K2) * Delta2B(j)&
                /sqrt(xi2**2d0 + Delta2B(j)**2d0)
                
                VDE_3toF = (-1d0 / Pi) * (K3**2d0) * V_swave(P,K3) * Delta3B(j)&
                /sqrt(xi3**2d0 + Delta3B(j)**2d0)
                
                S = S + w1(j) * VDE_1toF + w2(j) * VDE_2toF + w3(j) * VDE_3toF
        END DO
        Delta_F(i) = S
END DO

END SUBROUTINE Evaluation_Of_Gap_F
!---------------------------------------------------------------------------
SUBROUTINE Evaluation_Of_the_Energy(Energy, Kin, Pot, Delta1, Delta2, Delta3, mu_s, n, step, a, b, c)
IMPLICIT NONE
integer, INTENT(IN) :: n
real(kind=8), INTENT(IN) :: a, b, c, mu_s, step
real(kind=8), DIMENSION(1:n), INTENT(IN) :: Delta1, Delta2, Delta3
real(kind=8), INTENT(INOUT) :: Energy, Kin, Pot
real(kind=8) SK, SV, P, K1, K2, K3, eps1, eps2, eps3, VDE_1toEnergy, VDE_2toEnergy, VDE_3toEnergy
real(kind=8) P1, P2, P3, xi1, xi2, xi3, VDE_1to1, VDE_1to2, VDE_1to3, VDE_2to1, VDE_2to2, VDE_2to3 
real(kind=8) VDE_3to1, VDE_3to2, VDE_3to3, VDE_1toDiag,  VDE_2toDiag, VDE_3toDiag, Diagonal
real(kind=8), DIMENSION(1:n) :: x1, w1, x2, w2, x3, w3
real(kind=8), DIMENSION(1:n) :: v2k_1, v2k_2, v2k_3
real(kind=8), DIMENSION(1:n) :: u2k_1, u2k_2, u2k_3
integer i, j
!CALL p_quadrature_rule ( n, x, w ) ! l_quadrature::Laguerre, p_quadrature::Legendre

CALL gauleg(0d0, a, x1, w1, n)
CALL gauleg(a, b, x2, w2, n)
CALL gauleg(b, c, x3, w3, n)



SK= 0d0
Diagonal = 0d0
DO i=1,n
        K1 = x1(i)
        K2 = x2(i)
        K3 = x3(i)
        eps1 = 0.5d0 * hom * (K1**2d0)
        eps2 = 0.5d0 * hom * (K2**2d0)
        eps3 = 0.5d0 * hom * (K3**2d0)
        v2k_1(i) = (1d0/2d0) * (1d0 - (eps1-mu_s)/sqrt((eps1-mu_s)**2d0 + Delta1(i)**2d0))
        u2k_1(i) = 1d0 - v2k_1(i) 
        v2k_2(i) = (1d0/2d0) * (1d0 - (eps2-mu_s)/sqrt((eps2-mu_s)**2d0 + Delta2(i)**2d0))
        u2k_2(i) = 1d0 - v2k_2(i) 
        v2k_3(i) = (1d0/2d0) * (1d0 - (eps3-mu_s)/sqrt((eps3-mu_s)**2d0 + Delta3(i)**2d0))
        u2k_3(i) = 1d0 - v2k_3(i) 
        !WRITE(*,*) v2k_1, u2k_1  
        VDE_1toEnergy = (1d0 /(2d0 * Pi**2d0)) * (K1**2d0) * 2d0 * v2k_1(i) * eps1
        VDE_2toEnergy = (1d0 /(2d0 * Pi**2d0)) * (K2**2d0) * 2d0 * v2k_2(i) * eps2
        VDE_3toEnergy = (1d0 /(2d0 * Pi**2d0)) * (K3**2d0) * 2d0 * v2k_3(i) * eps3
        
        VDE_1toDiag = (2d0 / Pi) * (K1**2d0) * V_swave(K1,K1) * v2k_1(i) !Might need updating for polarized case
        VDE_2toDiag = (2d0 / Pi) * (K2**2d0) * V_swave(K2,K2) * v2k_2(i) 
        VDE_3toDiag = (2d0 / Pi) * (K3**2d0) * V_swave(K3,K3) * v2k_3(i) 
        
        SK = SK + w1(i) * VDE_1toEnergy + w2(i) * VDE_2toEnergy + w3(i) * VDE_3toEnergy
        Diagonal = Diagonal + w1(i) * VDE_1toDiag + w2(i) * VDE_2toDiag + w3(i) * VDE_3toDiag
END DO

SV=0d0
DO i=1,n
        P1 = x1(i)
        P2 = x2(i)
        P3 = x3(i)
        VDE_1toEnergy = 0d0
        VDE_2toEnergy = 0d0
        VDE_3toEnergy = 0d0
        DO j=1,n
                K1 = x1(j)
                K2 = x2(j)
                K3 = x3(j)
                VDE_1to1 = (1d0 / (Pi**3d0)) * (K1**2d0) * (P1**2d0) * V_swave(P1,K1) &
                * sqrt(DABS(u2k_1(i) * v2k_1(i) * u2k_1(j) * v2k_1(j)))
                
                VDE_1to2 = (1d0 / (Pi**3d0)) * (K1**2d0) * (P2**2d0) * V_swave(P2,K1) &
                * sqrt(DABS(u2k_2(i) * v2k_2(i) * u2k_1(j) * v2k_1(j)))
                
                VDE_1to3 = (1d0 / (Pi**3d0)) * (K1**2d0) * (P3**2d0) * V_swave(P3,K1) &
                * sqrt(DABS(u2k_3(i) * v2k_3(i) * u2k_1(j) * v2k_1(j)))

                VDE_2to1 = (1d0 / (Pi**3d0)) * (K2**2d0) * (P1**2d0) * V_swave(P1,K2) &
                * sqrt(DABS(u2k_1(i) * v2k_1(i) * u2k_2(j) * v2k_2(j)))
                
                VDE_2to2 = (1d0 / (Pi**3d0)) * (K2**2d0) * (P2**2d0) * V_swave(P2,K2) &
                * sqrt(DABS(u2k_2(i) * v2k_2(i) * u2k_2(j) * v2k_2(j)))
                
                VDE_2to3 = (1d0 / (Pi**3d0)) * (K2**2d0) * (P3**2d0) * V_swave(P3,K2) &
                * sqrt(DABS(u2k_3(i) * v2k_3(i) * u2k_2(j) * v2k_2(j)))
                
                VDE_3to1 = (1d0 / (Pi**3d0)) * (K3**2d0) * (P1**2d0) * V_swave(P1,K3) &
                * sqrt(DABS(u2k_1(i) * v2k_1(i) * u2k_3(j) * v2k_3(j)))
                
                VDE_3to2 = (1d0 / (Pi**3d0)) * (K3**2d0) * (P2**2d0) * V_swave(P2,K3) &
                * sqrt(DABS(u2k_2(i) * v2k_2(i) * u2k_3(j) * v2k_3(j)))
                
                VDE_3to3 = (1d0 / (Pi**3d0)) * (K3**2d0) * (P3**2d0) * V_swave(P3,K3) &
                * sqrt(DABS(u2k_3(i) * v2k_3(i) * u2k_3(j) * v2k_3(j)))
                
                VDE_1toEnergy = VDE_1toEnergy + w1(j) * VDE_1to1 + w2(j) * VDE_2to1 + w3(j) * VDE_3to1 
                VDE_2toEnergy = VDE_2toEnergy + w1(j) * VDE_1to2 + w2(j) * VDE_2to2 + w3(j) * VDE_3to2 
                VDE_3toEnergy = VDE_3toEnergy + w1(j) * VDE_1to3 + w2(j) * VDE_2to3 + w3(j) * VDE_3to3 
        END DO
        SV = SV + w1(i) * VDE_1toEnergy + w2(i) * VDE_2toEnergy + w3(i) * VDE_3toEnergy
END DO
write(9,*) SK, SV
Kin =  SK
Pot = SV
Energy = SK + SV!+ Diagonal

END SUBROUTINE Evaluation_Of_the_Energy
!---------------------------------------------------------------------------
SUBROUTINE Evaluation_Of_Density(Density, Delta1, Delta2, Delta3, mu_g, n, a, b, c)
IMPLICIT NONE
integer, INTENT(IN) :: n
real(kind=8), INTENT(IN) :: a, b, c, mu_g
real(kind=8), DIMENSION(1:n), INTENT(IN) :: Delta1, Delta2, Delta3
real(kind=8), INTENT(INOUT) :: Density
real(kind=8) S, P, K1, K2, K3, xi1, xi2, xi3, VDE_1toDensity, VDE_2toDensity, VDE_3toDensity
real(kind=8), DIMENSION(1:n) :: x1, w1, x2, w2, x3, w3
integer i
!CALL p_quadrature_rule ( n, x, w ) ! l_quadrature::Laguerre, p_quadrature::Legendre
CALL gauleg(0d0, a, x1, w1, n)
CALL gauleg(a, b, x2, w2, n)
CALL gauleg(b, c, x3, w3, n)

S = 0d0
DO i=1,n
        !K1 = a * (x(i) + 1d0) / 2d0
        !K2 = (b - a) * (x(i) + 1d0) / 2d0 + a
        !K3 = (c - b) * (x(i) + 1d0) / 2d0 + b
        K1 = x1(i)
        K2 = x2(i)
        K3 = x3(i)
        xi1 = 0.5d0 * hom * (K1**2d0) - mu_g
        xi2 = 0.5d0 * hom * (K2**2d0) - mu_g
        xi3 = 0.5d0 * hom * (K3**2d0) - mu_g
        
        VDE_1toDensity = (1d0 / (2d0 * (Pi**2d0))) * (K1**2d0) *&
        (1d0 - xi1 / sqrt(xi1**2d0 + Delta1(i)**2d0))
        
        VDE_2toDensity = (1d0 / (2d0 * (Pi**2d0))) * (K2**2d0) *&
        (1d0 - xi2 / sqrt(xi2**2d0 + Delta2(i)**2d0))
        
        VDE_3toDensity = (1d0 / (2d0 * (Pi**2d0))) * (K3**2d0) *&
        (1d0 - xi3 / sqrt(xi3**2d0 + Delta3(i)**2d0))
        
        S = S + w1(i) * VDE_1toDensity + w2(i) * VDE_2toDensity + w3(i) * VDE_3toDensity
END DO

Density = S

END SUBROUTINE Evaluation_Of_Density
!---------------------------------------------------------------------------
SUBROUTINE Comparison_F(D1,D2,er)

IMPLICIT NONE

real (kind=8), DIMENSION(0:F), INTENT(IN) :: D1,D2
real (kind=8) , INTENT(INOUT) :: er
integer :: i

er=0d0

DO i=0,F
        er = er + DABS((D1(i)-D2(i))/D1(i))   !consider using D1(i) instead of D1(0); It should change much since this is there to scale
END DO

er=er/DBLE(F)

END SUBROUTINE Comparison_F
!------------------------------------------------------------------------------------------------
SUBROUTINE Comparison_3n(D1,D2,D3,D1B,D2B,D3B,n,er)

IMPLICIT NONE
integer, INTENT(IN) :: n
real (kind=8), DIMENSION(1:n), INTENT(IN) :: D1, D2, D3, D1B, D2B, D3B
real (kind=8) , INTENT(INOUT) :: er
integer :: i

er=0d0
DO i=1,n
        er = er + DABS((D1(i)-D1B(i))/D1(1))   
        er = er + DABS((D2(i)-D2B(i))/D2(1))   
        er = er + DABS((D3(i)-D3B(i))/D3(1))   
END DO

er=er/DBLE(3*n)
END SUBROUTINE Comparison_3n
!------------------------------------------------------------------------------------------------
SUBROUTINE Farmer()

IMPLICIT NONE
integer :: i, ns, clock
integer, DIMENSION(:), ALLOCATABLE :: seed

CALL RANDOM_SEED(SIZE=ns)
ALLOCATE(seed(ns))
CALL SYSTEM_CLOCK(COUNT=clock)
seed= clock +37*(/(i-1, i=1,ns)/)
CALL RANDOM_SEED(PUT=seed)
DEALLOCATE(seed)
END SUBROUTINE Farmer
!------------------------------------------------------------------------------------

END MODULE BCS_eq





