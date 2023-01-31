MODULE BCS_eq
USE potentials
!USE Polynomials
!USE Math_Functions
IMPLICIT NONE
REAL (kind=8) , PARAMETER :: hom=41.44252d0
REAL (kind=8) :: a_par=-18.5d0, re_3D=2.7d0, vf
REAL (kind=8), ALLOCATABLE :: V11(:,:), V12(:,:), V13(:,:), V22(:,:), V23(:,:), V33(:,:)
REAL (kind=8), ALLOCATABLE :: U11(:,:), U12(:,:), U13(:,:), U22(:,:), U23(:,:), U33(:,:), phi1(:), phi2(:), phi3(:)
!REAL (kind=8) :: k_F=(3*rho_0*(PI**2d0))**(1d0/3d0), e_F=0.5*hom*(3*rho_0*PI**2)**(2.0/3.0), a_par=-18.5d0, re_3D=2.7d0
INTEGER  ::  F=1000, FAR=10

CONTAINS
SUBROUTINE Initialization()
END SUBROUTINE Initialization
!-------------------------------------------------------------------------
SUBROUTINE make_potential_matrix(a,b,c,n, kf)
!Calculate and store the potential matrix elements at the integration points
!of Gauss-Legendre. It is Vab(ij) = Vba(j,i)
IMPLICIT NONE
integer, INTENT(IN) :: n
real(kind=8), INTENT(IN) :: a, b, c, kf
real(kind=8) :: xin
real(kind=8), DIMENSION(1:n) :: x1, w1, x2, w2, x3, w3
integer i, j

CALL setup_potential()
write(*,*) "For potential:"
write(*,*) "x,a,b,c,n"
write(*,*) xx,a,b,c,n

CALL gauleg(0d0, a, x1, w1, n)
CALL gauleg(a, b, x2, w2, n)
CALL gauleg(b, c, x3, w3, n)

ALLOCATE(V11(1:n,1:n), V12(1:n,1:n), V13(1:n,1:n))
ALLOCATE(V22(1:n,1:n), V23(1:n,1:n))
ALLOCATE(V33(1:n,1:n))
ALLOCATE(U11(1:n,1:n), U12(1:n,1:n), U13(1:n,1:n))
ALLOCATE(U22(1:n,1:n), U23(1:n,1:n))
ALLOCATE(U33(1:n,1:n))
ALLOCATE(phi1(1:n),phi2(1:n), phi3(1:n))

write(1,*) "x, a, b, c, n"
write(1,*) xx,a,b,c,n
write(2,*) "x, a, b, c, n"
write(2,*) xx,a,b,c,n
write(3,*) "x, a, b, c, n"
write(3,*) xx,a,b,c,n

write(11,*) "x, a, b, c, n"
write(11,*) xx,a,b,c,n
write(12,*) "x, a, b, c, n"
write(12,*) xx,a,b,c,n
write(13,*) "x, a, b, c, n"
write(13,*) xx,a,b,c,n
write(22,*) "x, a, b, c, n"
write(22,*) xx,a,b,c,n
write(23,*) "x, a, b, c, n"
write(23,*) xx,a,b,c,n
write(33,*) "x, a, b, c, n"
write(33,*) xx,a,b,c,n

vf=vn3lo_1s0(kf,kf)
do i=1,n
        phi1(i) = vn3lo_1s0(kf,x1(i))/vf
        phi2(i) = vn3lo_1s0(kf,x2(i))/vf
        phi3(i) = vn3lo_1s0(kf,x3(i))/vf

        write(91,*) i, x1(i), phi1(i)
        write(92,*) i, x2(i), phi2(i)
        write(93,*) i, x3(i), phi3(i)
end do
write(*,*) 'done'
DO i=1,n
        DO j=1,n
                V11(i,j) = vn3lo_1s0(x1(i),x1(j)) 
                U11(i,j) = V11(i,j)-vf*phi1(i)*phi1(j)

                V12(i,j) = vn3lo_1s0(x1(i),x2(j)) 
                U12(i,j) = V12(i,j)-vf*phi1(i)*phi2(j)
                
                V13(i,j) = vn3lo_1s0(x1(i),x3(j)) 
                U13(i,j) = V13(i,j)-vf*phi1(i)*phi3(j)
                
                V22(i,j) = vn3lo_1s0(x2(i),x2(j)) 
                U22(i,j) = V22(i,j)-vf*phi2(i)*phi2(j)
                
                V23(i,j) = vn3lo_1s0(x2(i),x3(j)) 
                U23(i,j) = V23(i,j)-vf*phi2(i)*phi3(j)
                
                V33(i,j) = vn3lo_1s0(x3(i),x3(j))
                U33(i,j) = V33(i,j)-vf*phi3(i)*phi3(j)

                write(11,*) i, j, V11(i,j)
                write(111,*) i, j, U11(i,j)
                
                write(12,*) i, j, V12(i,j)
                write(112,*) i, j, U12(i,j)
                
                write(13,*) i, j, V13(i,j)
                write(113,*) i, j, U13(i,j)
                
                write(22,*) i, j, V22(i,j)
                write(122,*) i, j, U22(i,j)
                
                write(23,*) i, j, V23(i,j)
                write(123,*) i, j, U23(i,j)
                
                write(33,*) i, j, V33(i,j)
                write(133,*) i, j, U33(i,j)
        END DO
        write(*,*) i, x1(i), V11(i,i)
        write(1,*) i, phi1(i)
        write(2,*) i, phi2(i)
        write(3,*) i, phi3(i)
END DO



END SUBROUTINE make_potential_matrix
!------------------------------------------------------------
SUBROUTINE remake_potential_matrix(a,b,c,n, kf)
!Calculate and store the potential matrix elements at the integration points
!of Gauss-Legendre. It is Vab(ij) = Vba(j,i)
IMPLICIT NONE
integer, INTENT(IN) :: n
real(kind=8), INTENT(IN) :: a, b, c, kf
real(kind=8) :: xin
real(kind=8), DIMENSION(1:n) :: x1, w1, x2, w2, x3, w3
integer i, j

CALL gauleg(0d0, a, x1, w1, n)
CALL gauleg(a, b, x2, w2, n)
CALL gauleg(b, c, x3, w3, n)


vf=vn3lo_1s0(kf,kf)
do i=1,n
        phi1(i) = vn3lo_1s0(kf,x1(i))/vf
        phi2(i) = vn3lo_1s0(kf,x2(i))/vf
        phi3(i) = vn3lo_1s0(kf,x3(i))/vf
end do
write(*,*) 'done'
DO i=1,n
        DO j=1,n
                V11(i,j) = vn3lo_1s0(x1(i),x1(j)) 
                U11(i,j) = V11(i,j)-vf*phi1(i)*phi1(j)

                V12(i,j) = vn3lo_1s0(x1(i),x2(j)) 
                U12(i,j) = V12(i,j)-vf*phi1(i)*phi2(j)
                
                V13(i,j) = vn3lo_1s0(x1(i),x3(j)) 
                U13(i,j) = V13(i,j)-vf*phi1(i)*phi3(j)
                
                V22(i,j) = vn3lo_1s0(x2(i),x2(j)) 
                U22(i,j) = V22(i,j)-vf*phi2(i)*phi2(j)
                
                V23(i,j) = vn3lo_1s0(x2(i),x3(j)) 
                U23(i,j) = V23(i,j)-vf*phi2(i)*phi3(j)
                
                V33(i,j) = vn3lo_1s0(x3(i),x3(j))
                U33(i,j) = V33(i,j)-vf*phi3(i)*phi3(j)
        END DO
        write(*,*) i, x1(i), V11(i,i)
END DO
END SUBROUTINE remake_potential_matrix
!--------------------------------------------------
SUBROUTINE read_potential_matrix(a,b,c,n,kf)
!Read and store the potential matrix elements at the integration points
!of Gauss-Legendre. It is Vab(ij) = Vba(j,i)
!This is reading from the V*.dat files and assuming that they were printed
!by make_potential_matrix
IMPLICIT NONE
integer, INTENT(IN) :: n
real(kind=8), INTENT(IN) :: a, b, c, kf
real(kind=8), DIMENSION(1:n) :: x1, w1, x2, w2, x3, w3
real(kind=8) :: x0, xin, ain, bin, cin, nin 
integer i, j, bucketi, bucketj
character(len=50) :: lab
logical :: cond

CALL setup_potential()  !this reads v0/vs: the potential strength 
x0=xx                   !for comparison with the read potential
                        !(it cannot be altered if reading the potential)

ALLOCATE(V11(1:n,1:n), V12(1:n,1:n), V13(1:n,1:n))
ALLOCATE(V22(1:n,1:n), V23(1:n,1:n))
ALLOCATE(V33(1:n,1:n))
ALLOCATE(U11(1:n,1:n), U12(1:n,1:n), U13(1:n,1:n))
ALLOCATE(U22(1:n,1:n), U23(1:n,1:n))
ALLOCATE(U33(1:n,1:n))
ALLOCATE(phi1(1:n),phi2(1:n), phi3(1:n))

OPEN(UNIT=1,FILE="phi1.dat",STATUS="OLD")
OPEN(UNIT=2,FILE="phi2.dat",STATUS="OLD")
OPEN(UNIT=3,FILE="phi3.dat",STATUS="OLD")
OPEN(UNIT=11,FILE="V11.dat",STATUS="OLD")
OPEN(UNIT=111,FILE="U11.dat",STATUS="OLD")
OPEN(UNIT=12,FILE="V12.dat",STATUS="OLD")
OPEN(UNIT=112,FILE="U12.dat",STATUS="OLD")
OPEN(UNIT=13,FILE="V13.dat",STATUS="OLD")
OPEN(UNIT=113,FILE="U13.dat",STATUS="OLD")
OPEN(UNIT=22,FILE="V22.dat",STATUS="OLD")
OPEN(UNIT=122,FILE="U22.dat",STATUS="OLD")
OPEN(UNIT=23,FILE="V23.dat",STATUS="OLD")
OPEN(UNIT=123,FILE="U23.dat",STATUS="OLD")
OPEN(UNIT=33,FILE="V33.dat",STATUS="OLD")
OPEN(UNIT=133,FILE="U33.dat",STATUS="OLD")


read(1,*) lab 
read(2,*) lab 
read(3,*) lab 

read(11,*) lab 
read(12,*) lab 
read(13,*) lab 
read(22,*) lab 
read(23,*) lab 
read(33,*) lab 

cond=.FALSE.
read(1,*) xin, ain, bin, cin, nin
cond = cond .OR. dabs(x0-xin).ge.small 
cond = cond .OR. dabs(ain-a).ge.small .OR. dabs(bin-b).ge.small .OR. dabs(cin-c).ge.small
cond = cond .OR. nin.ne.n
read(2,*) xin, ain, bin, cin, nin
cond = cond .OR. dabs(x0-xin).ge.small 
cond = cond .OR. dabs(ain-a).ge.small .OR. dabs(bin-b).ge.small .OR. dabs(cin-c).ge.small
cond = cond .OR. nin.ne.n
read(3,*) xin, ain, bin, cin, nin
cond = cond .OR. dabs(x0-xin).ge.small 
cond = cond .OR. dabs(ain-a).ge.small .OR. dabs(bin-b).ge.small .OR. dabs(cin-c).ge.small
cond = cond .OR. nin.ne.n
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
        !PAUSE
ENDIF 

write(*,*) "For the potential:"
write(*,*) "x,a,b,c,n"
write(*,*) xx,a,b,c,n


vf=V_pwave(kf,kf)
DO i=1,n
        write(*,*) i
        read(1,*) bucketi, phi1(i)
        read(2,*) bucketi, phi2(i)
        read(3,*) bucketi, phi3(i)
        DO j=1,n
                read(11,*) bucketi, bucketj, V11(i,j)
                read(111,*) bucketi, bucketj, U11(i,j)
                read(12,*) bucketi, bucketj, V12(i,j)
                read(112,*) bucketi, bucketj, U12(i,j)
                read(13,*) bucketi, bucketj, V13(i,j)
                read(113,*) bucketi, bucketj, U13(i,j)
                read(22,*) bucketi, bucketj, V22(i,j)
                read(122,*) bucketi, bucketj, U22(i,j)
                read(23,*) bucketi, bucketj, V23(i,j)
                read(123,*) bucketi, bucketj, U23(i,j)
                read(33,*) bucketi, bucketj, V33(i,j)
                read(133,*) bucketi, bucketj, U33(i,j)
        END DO
END DO

CLOSE(UNIT=1,STATUS="KEEP")
CLOSE(UNIT=2,STATUS="KEEP")
CLOSE(UNIT=3,STATUS="KEEP")
CLOSE(UNIT=11,STATUS="KEEP")
CLOSE(UNIT=111,STATUS="KEEP")
CLOSE(UNIT=12,STATUS="KEEP")
CLOSE(UNIT=112,STATUS="KEEP")
CLOSE(UNIT=13,STATUS="KEEP")
CLOSE(UNIT=113,STATUS="KEEP")
CLOSE(UNIT=22,STATUS="KEEP")
CLOSE(UNIT=122,STATUS="KEEP")
CLOSE(UNIT=23,STATUS="KEEP")
CLOSE(UNIT=123,STATUS="KEEP")
CLOSE(UNIT=33,STATUS="KEEP")
CLOSE(UNIT=133,STATUS="KEEP")

if (.true.) then
        V11=V11*x0
        V12=V12*x0
        V13=V13*x0
        V22=V22*x0
        V23=V23*x0
        V33=V33*x0
        
        U11=U11*x0
        U12=U12*x0
        U13=U13*x0
        U22=U22*x0
        U23=U23*x0
        U33=U33*x0
end if
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
SUBROUTINE Root_Of_Gap(Density, D1, D2, D3, df, mu, M1, M2, rho, n, step, a, b, c,&
                        TOL_for_Iter,TOL_for_Ridder)
IMPLICIT NONE
integer, INTENT(IN) :: n
real(kind=8), INTENT(IN) :: M1, M2, rho, step, a, b, c, TOL_for_Iter, TOL_for_Ridder
real(kind=8), INTENT(OUT) :: Density, mu
real(kind=8), INTENT(INOUT) :: df
real(kind=8), DIMENSION(1:n),INTENT(INOUT) :: D1, D2, D3
real(kind=8) :: Density1, Density2, Density3, Density4, mu1, mu2, mu3, mu4, f1, f2, f3, f4
real(kind=8) :: df1, df2, df3, df4
real(kind=8) :: mu4tocheck, mineq
integer :: iter

f1 = 1d0
f2 = 1d0
mu1 = M1
mu2 = M2


write(*,*) "Calculating gap for mu1=", mu1
CALL EVALUATION_TL_KHODEL(Density1, mu1, df1, mineq, n, step, a, b, c, TOL_for_Ridder,TOL_for_Iter)
f1 = Density1 - rho

write(*,*) "Calculating gap for mu2=", mu2
CALL EVALUATION_TL_KHODEL(Density2, mu2, df2, mineq, n, step, a, b, c, TOL_for_Ridder, TOL_for_Iter)
f2 = Density2 - rho

iter=0
DO 
        iter = iter +1
        IF(f1 * f2 < 0d0) THEN
                WRITE(*,*) "Solution bracketed with f1 =", f1, "and f2 =", f2, "after", iter, "iterations"
                EXIT
        END IF
        WRITE(*,*) "(pre-bracketing) Ridder iteration:", iter, "f1 =", f1, "and f2 =", f2, "in [", mu1, mu2, "]"
        
        if (dabs(f1)<dabs(f2)) then
                !mu1 = 0.7d0 * mu1
                mu1 = mu1/10d0 
                write(*,*) "Calculating gap for mu1=", mu1
                CALL EVALUATION_TL_KHODEL(Density1, mu1, df1, mineq, n, step, a, b, c, TOL_for_Ridder, TOL_for_Iter)
                f1 = Density1 - rho
        else
                !mu2 = 1.3d0 * mu2
                mu2 = mu2*10d0
                write(*,*) "Calculating gap for mu2=", mu2
                CALL EVALUATION_TL_KHODEL(Density2, mu2, df2, mineq, n, step, a, b, c, TOL_for_Ridder, TOL_for_Iter)
                f2 = Density2 - rho
        end if
        if (iter > 100) then
                pause
        end if
END DO

mu4tocheck = 1000d0
iter = 0
DO 
        iter = iter + 1
        mu3 = (mu1 + mu2) / 2d0
        CALL EVALUATION_TL_KHODEL(Density3, mu3, df3, mineq,  n, step, a, b, c, TOL_for_Ridder, TOL_for_Iter)
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

        CALL EVALUATION_TL_KHODEL(Density4, mu4, df4, mineq, n, step, a, b, c, TOL_for_Ridder, TOL_for_Iter)
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

CALL EVALUATION_TL_for_Delta_Khodel(Density, D1, D2, D3, df, mu, df4, n, step, a, b, c, TOL_for_Ridder, TOL_for_Iter)
WRITE(*,*) '--->> Solution found with Density =', Density, 'and df, mineq=', df, mineq

END SUBROUTINE Root_Of_Gap
!-------------------------------------------------------------------------------------------------------
SUBROUTINE get_mineq(mineq, chi1, chi2, chi3, df, mu_g, n, a, b, c)
IMPLICIT NONE
integer, INTENT(IN) :: n
real(kind=8), INTENT(IN) :: a, b, c, mu_g, df
real(kind=8), DIMENSION(1:n), INTENT(IN) :: chi1, chi2, chi3
real(kind=8), INTENT(INOUT) :: mineq
real(kind=8) :: xi, eq, buck
real(kind=8), DIMENSION(1:n) :: x1, x2, x3, w1, w2, w3
real(kind=8), allocatable :: ks(:), chis(:)
integer ik
allocate(ks(1:3*n),chis(1:3*n))
CALL gauleg(0d0, a, x1, w1, n)
CALL gauleg(a, b, x2, w2, n)
CALL gauleg(b, c, x3, w3, n)


chis = (/chi1,chi2,chi3/)
ks = (/x1,x2,x3/)
buck=100d0
do ik=1,3*n
        xi = 0.5*hom*ks(ik)**2-mu_g
        eq = sqrt(xi**2+df**2*chis(ik)**2)
        if (eq<buck) then
                buck = eq
        end if
end do
mineq=buck
end subroutine get_mineq
!---------------------------------------------------------------------------
SUBROUTINE get_df(df, chi1, chi2, chi3, mu_g, df_0, n, a, b, c, TOL_for_Ridder)
IMPLICIT NONE
integer, INTENT(IN) :: n
real(kind=8), INTENT(IN) :: mu_g, df_0, a, b, c, TOL_for_Ridder
real(kind=8), INTENT(IN), dimension(1:n) :: chi1, chi2, chi3
real(kind=8), INTENT(OUT) :: df
real(kind=8) :: ivf1, ivf2, ivf3, ivf4, df1, df2, df3, df4, f1, f2, f3, f4
real(kind=8) :: df4tocheck, ddf4
integer :: iter

f1 = 1d0
f2 = 1d0
df1 = 0.8d0*df_0
df2 = 1.2d0*df_0

!write(*,*) "Calculating Ivf for df1=", df1
CALL eval_df(ivf1, chi1, chi2, chi3, df1, mu_g, n, a, b, c)
f1 = ivf1 - (vf)**(-1d0)

!write(*,*) "Calculating Ivf for df2=", df2
CALL eval_df(ivf2, chi1, chi2, chi3, df2, mu_g, n, a, b, c)
f2 = ivf2 - (vf)**(-1d0)

iter=0
DO 
        iter = iter +1
        IF(f1 * f2 < 0d0) THEN
                !WRITE(*,*) "in-df: Solution bracketed with f1 =", f1, "and f2 =", f2, "after", iter, "iterations"
                EXIT
        END IF
        
        if (iter>500) then
                write(*,*) 'in-df: iter>500 -- maybe stuck with f1 =', &
                        f1, "and f2 =", f2, "in [", df1, df2, "]", vf, 1d0/ivf1, 1d0/ivf2
        end if
        !WRITE(*,*) "            in-df: (pre-bracketing) Ridder iteration:", &
        !        iter, "f1 =", f1, "and f2 =", f2, "in [", df1, df2, "]", vf, 1d0/ivf1, 1d0/ivf2
        
        if (dabs(f1)<dabs(f2)) then
                df1 = df1*0.8d0 
                !write(*,*) "in-df: Calculating gap for df1=", df1
                CALL eval_df(ivf1, chi1, chi2, chi3, df1, mu_g, n, a, b, c)
                f1 = ivf1 - (vf)**(-1d0)
        else
                df2 = df2*1.2d0
                !write(*,*) "in-df: Calculating gap for mu2=", df2
                CALL eval_df(ivf2, chi1, chi2, chi3, df2, mu_g, n, a, b, c)
                f2 = ivf2 - (vf)**(-1d0)
        end if
END DO

df4tocheck = 1000d0
iter = 0
DO 
        iter = iter + 1
        df3 = 0.5*(df1 + df2)
        CALL eval_df(ivf3, chi1, chi2, chi3, df3, mu_g, n, a, b, c)
        f3 = ivf3 - (vf)**(-1d0)
        
        IF(f1 - f2 > 0d0) THEN
                df4 = df3 + (df3 - df1) * f3 / sqrt(f3**2d0 - f1 * f2)
        ELSE IF(f1 - f2 < 0d0) THEN
                df4 = df3 - (df3 - df1) * f3 / sqrt(f3**2d0 - f1 * f2)
        ELSE
                WRITE(*,*) "in-df: Location: get_df"
                WRITE(*,*) "in-df: Failed at conditions f1-f2<0 OR f1-f2>0"
                WRITE(*,*) "in-df: With f1 =", f1, "and f2 =", f2
        END IF

        
        CALL eval_df(ivf4, chi1, chi2, chi3, df4, mu_g, n, a, b, c)
        f4 = ivf4 - (vf)**(-1d0)
        ddf4 = dabs(df4-df4tocheck)
        !WRITE(*,*) "in-df: (post-bracketing) Ridder iteration:", iter, "At ivf - ivf_0  =", f4
        !WRITE(*,*) "in-df:                    and ddf =", ddf4
        IF(ddf4< TOL_for_Ridder) THEN
                EXIT
        END IF

        IF(f1 * f4 < 0d0) THEN
                df2 =df4
                ivf2 = ivf4
                f2 = f4
        ELSE IF(f1 * f4 > 0d0) THEN
                df1 =df4
                ivf1 = ivf4
                f1 = f4
        ELSE
                WRITE(*,*) "in-df: Location: get_df at conditions f1*f4<0 OR f1*f4>0: soft error"
                WRITE(*,*) "in-df: With f1 =", f1, "and f4 =", f4
                WRITE(*,*) "in-df: With df1 =", df1, "and df4 =", df4
                !pause
        END IF
        df4tocheck = df4
END DO
df = df4
CALL eval_df(ivf4, chi1, chi2, chi3, df4, mu_g, n, a, b, c)
!write(*,*) 'solution found with diffference', ivf4-vf**(-1d0)
END SUBROUTINE get_df
!------------------------------------------------------------
SUBROUTINE eval_df(ivf, chi1, chi2, chi3, df, mu_g, n, a, b, c)
IMPLICIT NONE
integer, INTENT(IN) :: n
real(kind=8), INTENT(IN) :: a, b, c, mu_g, df
real(kind=8), DIMENSION(1:n), INTENT(IN) :: chi1, chi2, chi3
real(kind=8), INTENT(INOUT) :: ivf
real(kind=8) :: S, K1, K2, K3, xi1, xi2, xi3, VDE1, VDE2, VDE3, VE1, VE2, VE3
real(kind=8), DIMENSION(1:n) :: x1, x2, x3, w1, w2, w3
integer i
CALL gauleg(0d0, a, x1, w1, n)
CALL gauleg(a, b, x2, w2, n)
CALL gauleg(b, c, x3, w3, n)


S = 0d0
DO i=1,n
        K1 = x1(i)
        K2 = x2(i)
        K3 = x3(i)
        xi1 = 0.5d0 * hom * (K1**2d0) - mu_g
        xi2 = 0.5d0 * hom * (K2**2d0) - mu_g
        xi3 = 0.5d0 * hom * (K3**2d0) - mu_g
        
        VE1 = (1d0 / (Pi)) * (K1**2d0) * phi1(i)**2 &
        /sqrt(xi1**2d0 + df**2*chi1(i)**2d0)
        VDE1 = (1d0 / (Pi)) * (K1**2d0) * phi1(i) * (phi1(i)-chi1(i))&
        /sqrt(xi1**2d0 + df**2*chi1(i)**2d0)
        
        VE2 = (1d0 / (Pi)) * (K2**2d0) * phi2(i)**2 &
        /sqrt(xi2**2d0 + df**2*chi2(i)**2d0)
        VDE2 = (1d0 / (Pi)) * (K2**2d0) * phi2(i) * (phi2(i)-chi2(i))&
        /sqrt(xi2**2d0 + df**2*chi2(i)**2d0)
        
        VE3 = (1d0 / (Pi)) * (K3**2d0) * phi3(i)**2 &
        /sqrt(xi3**2d0 + df**2*chi3(i)**2d0)
        VDE3 = (1d0 / (Pi)) * (K3**2d0) * phi3(i) * (phi3(i)-chi3(i))&
        /sqrt(xi3**2d0 + df**2*chi3(i)**2d0)
        
        S = S + w1(i) * (VDE1 - VE1) 
        S = S + w2(i) * (VDE2 - VE2) 
        S = S + w3(i) * (VDE3 - VE3) 
END DO
ivf=S

END SUBROUTINE eval_df
!--------------------------------------------------
SUBROUTINE EVALUATION_TL_KHODEL_old(Density, mu_g, dfout, df0in, n, step, a, b, c, TOL_for_Ridder, TOL_for_Iter)
IMPLICIT NONE
integer, INTENT(IN) :: n
real(kind=8), INTENT(IN) :: step, a, b, c, mu_g, TOL_for_Ridder, TOL_for_Iter, df0in
real(kind=8), INTENT(INOUT) :: Density, dfout
real(kind=8), DIMENSION(1:n) :: chi1, chi1B, chi2, chi2B, chi3, chi3B
real(kind=8), DIMENSION(0:F) :: Delta_FB, Delta_F
real(kind=8) :: er,df, df0
integer :: i, iter

chi1=1d0
chi2=1d0
chi3=1d0
chi1B=1d0
chi2B=1d0
chi3B=1d0

df0=df0in
CALL get_df(df, chi1B, chi2B, chi3B, mu_g, df0, n, a, b, c, TOL_for_Ridder)

iter = 0
DO
        iter = iter + 1 
        CALL Evaluation_Of_Gap_n_khodel(chi1, chi2, chi3, chi1B, chi2B, chi3B, mu_g, df,n, a, b, c)
        CALL Comparison_3n(chi1, chi2, chi3, chi1B, chi2B, chi3B, n, er)
        WRITE(*,*) "           error of gap iteration ", er, chi1(10)
        
        IF (DABS(er) < TOL_for_Iter) THEN
                EXIT
        END IF
        
        IF (iter > 100000) THEN
                WRITE(*,*) "Location: EVALUATION_TL"
                WRITE(*,*) "Time exceeded at iterations"
                WRITE(*,*) "With Delta(0) =", Delta_F(0)
                PAUSE
        END IF
        chi1B = chi1
        chi2B = chi2
        chi3B = chi3
        df0=df       
        CALL get_df(df, chi1, chi2, chi3, mu_g, df0, n, a, b, c, TOL_for_Ridder)
END DO
dfout=df


CALL Evaluation_Of_Density(Density, chi1, chi2, chi3, df, mu_g, n, a, b, c)


END SUBROUTINE EVALUATION_TL_KHODEL_old
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE EVALUATION_TL_KHODEL(Density, mu_g, dfout, mineqout,  n, step, a, b, c, TOL_for_Ridder, TOL_for_Iter)
IMPLICIT NONE
integer, INTENT(IN) :: n
real(kind=8), INTENT(IN) :: step, a, b, c, mu_g, TOL_for_Ridder, TOL_for_Iter
real(kind=8), INTENT(INOUT) :: Density, dfout, mineqout
real(kind=8), DIMENSION(1:n) :: chi1, chi1B, chi2, chi2B, chi3, chi3B
real(kind=8), DIMENSION(0:F) :: Delta_FB, Delta_F
real(kind=8), DIMENSION(1:n) :: x1, x2, x3, w1, w2, w3
real(kind=8) :: er,df, df0, mineq
integer :: i, ichi, idf
CALL gauleg(0d0, a, x1, w1, n)
CALL gauleg(a, b, x2, w2, n)
CALL gauleg(b, c, x3, w3, n)


chi1=1d0
chi2=1d0
chi3=1d0
chi1B=1d0
chi2B=1d0
chi3B=1d0

!CALL get_df(df, chi1B, chi2B, chi3B, mu_g, df0, n, a, b, c, TOL_for_Ridder)
df=0.1d0
df0=100*df
idf=0
do while (abs(df-df0)>10d-8)
        idf = idf + 1
        ichi = 0
        DO
                ichi = ichi + 1 
                CALL Evaluation_Of_Gap_n_khodel(chi1, chi2, chi3, chi1B, chi2B, chi3B, mu_g, df,n, a, b, c)
                CALL Comparison_3n(chi1, chi2, chi3, chi1B, chi2B, chi3B, n, er)
                !WRITE(*,*) "           error of gap iteration ", er, chi1(10)
                
                IF (DABS(er) < TOL_for_Iter) THEN
                        EXIT
                END IF
                
                IF (ichi > 100000) THEN
                        WRITE(*,*) "Location: EVALUATION_TL"
                        WRITE(*,*) "Time exceeded at iterations"
                        WRITE(*,*) "With chi123(10) =", chi1(10), chi2(10), chi3(10)
                        PAUSE
                END IF
                chi1B = chi1
                chi2B = chi2
                chi3B = chi3
        END DO       
        do i=1,n
                write(91,*) i, x1(i), chi1(i)
                write(92,*) i, x2(i), chi2(i)
                write(93,*) i, x3(i), chi3(i)
        end do
        df0=df
        CALL get_df(df, chi1, chi2, chi3, mu_g, df0, n, a, b, c, TOL_for_Ridder)
        write(*,*) idf, 'df-df0', df-df0
end do
dfout=df


CALL Evaluation_Of_Density(Density, chi1, chi2, chi3, df, mu_g, n, a, b, c)
CALL get_mineq(mineq, chi1, chi2, chi3, df, mu_g, n, a, b, c)

write(*,*) '*mineq', mineq

mineqout=mineq
END SUBROUTINE EVALUATION_TL_KHODEL
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE EVALUATION_TL_for_Delta_KHODEL(Density, D1, D2, D3, dfout, mu_g, dfin, n, step, a, b, c, &
                TOL_for_Ridder, TOL_for_Iter)
IMPLICIT NONE
integer, INTENT(IN) :: n
real(kind=8), INTENT(IN) :: step, a, b, c, mu_g, dfin, TOL_for_Ridder, TOL_for_Iter 
real(kind=8), INTENT(INOUT) :: Density, dfout
real(kind=8), INTENT(INOUT), DIMENSION(1:n) :: D1, D2, D3
real(kind=8), DIMENSION(1:n) :: chi1, chi2, chi3, chi1B, chi2B, chi3B
real(kind=8) :: er, df, df0
integer :: i, iter

chi1=1d0
chi2=1d0
chi3=1d0
chi1B=1d0
chi2B=1d0
chi3B=1d0

df0 = dfin
CALL get_df(df, chi1, chi2, chi3, mu_g, df0, n, a, b, c, TOL_for_Ridder)

iter = 0
DO
        iter = iter + 1 
        CALL Evaluation_Of_Gap_n_khodel(chi1, chi2, chi3, chi1B, chi2B, chi3B, mu_g, df,n, a, b, c)
        CALL Comparison_3n(chi1, chi2, chi3, chi1B, chi2B, chi3B, n, er)
        WRITE(*,*) "           error of gap iteration ", er, chi1(10)
        
        IF (DABS(er) < TOL_for_Iter) THEN
                EXIT
        END IF
        
        IF (iter > 100000) THEN
                WRITE(*,*) "Location: EVALUATION_TL"
                WRITE(*,*) "Time exceeded at iterations"
                PAUSE
        END IF
        chi1B = chi1
        chi2B = chi2
        chi3B = chi3
        df0=df
        CALL get_df(df, chi1B, chi2B, chi3B, mu_g, df0, n, a, b, c, TOL_for_Ridder)
END DO
dfout=df

CALL Evaluation_Of_Density(Density, chi1, chi2, chi3, df, mu_g, n, a, b, c)

D1=df*chi1
D2=df*chi2
D3=df*chi3

END SUBROUTINE EVALUATION_TL_for_Delta_KHODEL
!-------------------------------------------------------------------------------------------------
SUBROUTINE Evaluation_Of_Gap_n_khodel(chi1, chi2, chi3, chi1B, chi2B, chi3B, mu_g, df,n, a, b, c)
IMPLICIT NONE
integer, INTENT(IN) :: n
real(kind=8), INTENT(IN) :: a, b, c, mu_g, df
real(kind=8), DIMENSION(1:n), INTENT(IN) :: chi1B, chi2B, chi3B
real(kind=8), DIMENSION(1:n), INTENT(INOUT) :: chi1, chi2, chi3
real(kind=8) :: S1, S2, S3, P1, P2, P3, K1, K2, K3, xi1, xi2, xi3, VDE_1to1, VDE_1to2, VDE_1to3, VDE_2to1, VDE_2to2, VDE_2to3 
real(kind=8) :: VDE_3to1, VDE_3to2, VDE_3to3 
real(kind=8), DIMENSION(1:n) :: x1, x2, x3, w1, w2, w3
integer i, j
CALL gauleg(0d0, a, x1, w1, n)
CALL gauleg(a, b, x2, w2, n)
CALL gauleg(b, c, x3, w3, n)


DO i=1,n
        S1 = 0d0
        S2 = 0d0
        S3 = 0d0
        P1 = x1(i)
        P2 = x2(i)
        P3 = x3(i)
        DO j=1,n
                K1 = x1(j)
                K2 = x2(j)
                K3 = x3(j)
                xi1 = 0.5d0 * hom * (K1**2d0) - mu_g
                xi2 = 0.5d0 * hom * (K2**2d0) - mu_g
                xi3 = 0.5d0 * hom * (K3**2d0) - mu_g
                
                VDE_1to1 = (-1d0 / Pi) * (K1**2d0) * U11(i,j) * chi1B(j)&
                /sqrt(xi1**2d0 + df**2*chi1B(j)**2d0)
                
                VDE_1to2 = (-1d0 / (Pi)) * (K1**2d0) * U12(j,i) * chi1B(j)&
                /sqrt(xi1**2d0 + df**2*chi1B(j)**2d0)
                
                VDE_1to3 = (-1d0 / (Pi)) * (K1**2d0) * U13(j,i) * chi1B(j)&
                /sqrt(xi1**2d0 + df**2*chi1B(j)**2d0)
                
                VDE_2to1 = (-1d0 / (Pi)) * (K2**2d0) * U12(i,j) * chi2B(j)&
                /sqrt(xi2**2d0 + df**2*chi2B(j)**2d0)
                
                VDE_2to2 = (-1d0 / (Pi)) * (K2**2d0) * U22(i,j) * chi2B(j)&
                /sqrt(xi2**2d0 + df**2*chi2B(j)**2d0)
                
                VDE_2to3 = (-1d0 / (Pi)) * (K2**2d0) * U23(j,i) * chi2B(j)&
                /sqrt(xi2**2d0 + df**2*chi2B(j)**2d0)
                
                VDE_3to1 = (-1d0 / (Pi)) * (K3**2d0) * U13(i,j) * chi3B(j)&
                /sqrt(xi3**2d0 + df**2*chi3B(j)**2d0)
                
                VDE_3to2 = (-1d0 / (Pi)) * (K3**2d0) * U23(i,j) * chi3B(j)&
                /sqrt(xi3**2d0 + df**2*chi3B(j)**2d0)
                
                VDE_3to3 = (-1d0 / (Pi)) * (K3**2d0) * U33(i,j) * chi3B(j)&
                /sqrt(xi3**2d0 + df**2*chi3B(j)**2d0)
                
                S1 = S1 + w1(j) * VDE_1to1 + w2(j) * VDE_2to1 + w3(j) * VDE_3to1 
                S2 = S2 + w1(j) * VDE_1to2 + w2(j) * VDE_2to2 + w3(j) * VDE_3to2 
                S3 = S3 + w1(j) * VDE_1to3 + w2(j) * VDE_2to3 + w3(j) * VDE_3to3 
        END DO
        chi1(i) = S1 + phi1(i)
        chi2(i) = S2 + phi2(i)
        chi3(i) = S3 + phi3(i)
END DO


END SUBROUTINE Evaluation_Of_Gap_n_khodel
!--------------------------------------------------
SUBROUTINE Evaluation_Of_Gap_F(Delta_Final, chi1, chi2, chi3, mu_g, df, n, step, a, b, c)
IMPLICIT NONE
integer, INTENT(IN) :: n
real(kind=8), INTENT(IN) :: a, b, c, mu_g, df, step
real(kind=8), DIMENSION(1:n), INTENT(IN) :: chi1, chi2, chi3
real(kind=8), DIMENSION(0:F), INTENT(INOUT) :: Delta_Final
real(kind=8) :: S, P, K1, K2, K3, xi1, xi2, xi3, VDE_1toF, VDE_2toF, VDE_3toF
real(kind=8), DIMENSION(1:n) :: x1, w1, x2, w2, x3, w3
integer i, j
CALL gauleg(0d0, a, x1, w1, n)
CALL gauleg(a, b, x2, w2, n)
CALL gauleg(b, c, x3, w3, n)

DO i=0,F
        S = 0d0
        P = DBLE(i * step)
        DO j=1,n
                K1 = x1(j)
                K2 = x2(j)
                K3 = x3(j)
                
                xi1 = 0.5d0 * hom * (K1**2d0) - mu_g
                xi2 = 0.5d0 * hom * (K2**2d0) - mu_g
                xi3 = 0.5d0 * hom * (K3**2d0) - mu_g
                
                VDE_1toF = (-1d0 / (Pi)) * (K1**2d0) * V_pwave(P,K1) * df * chi1(j)&
                /sqrt(xi1**2d0 + df**2d0*chi1(j)**2d0)
                
                VDE_2toF = (-1d0 / (Pi)) * (K2**2d0) * V_pwave(P,K2) * df * chi2(j)&
                /sqrt(xi2**2d0 + df**2d0*chi2(j)**2d0)
                
                VDE_3toF = (-1d0 / (Pi)) * (K3**2d0) * V_pwave(P,K3) * df * chi3(j)&
                /sqrt(xi3**2d0 + df**2d0*chi3(j)**2d0)
                
                S = S + w1(j) * VDE_1toF + w2(j) * VDE_2toF + w3(j) * VDE_3toF
        END DO
        Delta_Final(i) = S
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
        VDE_1toEnergy = (1d0 /(4d0 * Pi**2d0)) * (K1**2d0) * 2d0 * v2k_1(i) * eps1
        VDE_2toEnergy = (1d0 /(4d0 * Pi**2d0)) * (K2**2d0) * 2d0 * v2k_2(i) * eps2
        VDE_3toEnergy = (1d0 /(4d0 * Pi**2d0)) * (K3**2d0) * 2d0 * v2k_3(i) * eps3
        
        VDE_1toDiag = (2d0 / Pi) * (K1**2d0) * V_pwave(K1,K1) * v2k_1(i) !Might need updating for polarized case
        VDE_2toDiag = (2d0 / Pi) * (K2**2d0) * V_pwave(K2,K2) * v2k_2(i) 
        VDE_3toDiag = (2d0 / Pi) * (K3**2d0) * V_pwave(K3,K3) * v2k_3(i) 
        
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
                VDE_1to1 = (1d0 / (4d0*Pi**3d0)) * (K1**2d0) * (P1**2d0) * V_pwave(P1,K1) &
                * sqrt(DABS(u2k_1(i) * v2k_1(i) * u2k_1(j) * v2k_1(j)))
                
                VDE_1to2 = (1d0 / (4d0*Pi**3d0)) * (K1**2d0) * (P2**2d0) * V_pwave(P2,K1) &
                * sqrt(DABS(u2k_2(i) * v2k_2(i) * u2k_1(j) * v2k_1(j)))
                
                VDE_1to3 = (1d0 / (4d0*Pi**3d0)) * (K1**2d0) * (P3**2d0) * V_pwave(P3,K1) &
                * sqrt(DABS(u2k_3(i) * v2k_3(i) * u2k_1(j) * v2k_1(j)))

                VDE_2to1 = (1d0 / (4d0*Pi**3d0)) * (K2**2d0) * (P1**2d0) * V_pwave(P1,K2) &
                * sqrt(DABS(u2k_1(i) * v2k_1(i) * u2k_2(j) * v2k_2(j)))
                
                VDE_2to2 = (1d0 / (4d0*Pi**3d0)) * (K2**2d0) * (P2**2d0) * V_pwave(P2,K2) &
                * sqrt(DABS(u2k_2(i) * v2k_2(i) * u2k_2(j) * v2k_2(j)))
                
                VDE_2to3 = (1d0 / (4d0*Pi**3d0)) * (K2**2d0) * (P3**2d0) * V_pwave(P3,K2) &
                * sqrt(DABS(u2k_3(i) * v2k_3(i) * u2k_2(j) * v2k_2(j)))
                
                VDE_3to1 = (1d0 / (4d0*Pi**3d0)) * (K3**2d0) * (P1**2d0) * V_pwave(P1,K3) &
                * sqrt(DABS(u2k_1(i) * v2k_1(i) * u2k_3(j) * v2k_3(j)))
                
                VDE_3to2 = (1d0 / (4d0*Pi**3d0)) * (K3**2d0) * (P2**2d0) * V_pwave(P2,K3) &
                * sqrt(DABS(u2k_2(i) * v2k_2(i) * u2k_3(j) * v2k_3(j)))
                
                VDE_3to3 = (1d0 / (4d0*Pi**3d0)) * (K3**2d0) * (P3**2d0) * V_pwave(P3,K3) &
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
SUBROUTINE Evaluation_Of_Density(Density, chi1, chi2, chi3,  df, mu_g, n, a, b, c)
IMPLICIT NONE
integer, INTENT(IN) :: n
real(kind=8), INTENT(IN) :: a, b, c, mu_g, df
real(kind=8), DIMENSION(1:n), INTENT(IN) :: chi1, chi2, chi3
real(kind=8), INTENT(INOUT) :: Density
real(kind=8) S, P, K1, K2, K3, xi1, xi2, xi3, VDE_1toDensity, VDE_2toDensity, VDE_3toDensity
real(kind=8), DIMENSION(1:n) :: x1, x2, x3, w1, w2, w3
integer i
CALL gauleg(0d0, a, x1, w1, n)
CALL gauleg(a, b, x2, w2, n)
CALL gauleg(b, c, x3, w3, n)

S = 0d0
DO i=1,n
        K1 = x1(i)
        K2 = x2(i)
        K3 = x3(i)
        xi1 = 0.5d0 * hom * (K1**2d0) - mu_g
        xi2 = 0.5d0 * hom * (K2**2d0) - mu_g
        xi3 = 0.5d0 * hom * (K3**2d0) - mu_g
        
        VDE_1toDensity = (1d0 / (2d0 * (Pi**2d0))) * (K1**2d0) *&
        (1d0 - xi1 / sqrt(xi1**2d0 + df**2d0*chi1(i)**2d0))
        
        VDE_2toDensity = (1d0 / (2d0 * (Pi**2d0))) * (K2**2d0) *&
        (1d0 - xi2 / sqrt(xi2**2d0 + df**2d0*chi2(i)**2d0))
        
        VDE_3toDensity = (1d0 / (2d0 * (Pi**2d0))) * (K3**2d0) *&
        (1d0 - xi3 / sqrt(xi3**2d0 + df**2d0*chi3(i)**2d0))
        
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
SUBROUTINE Comparison_n(D1,D1B,n,er)

IMPLICIT NONE
integer, INTENT(IN) :: n
real (kind=8), DIMENSION(1:n), INTENT(IN) :: D1, D1B
real (kind=8) , INTENT(INOUT) :: er
integer :: i

er=0d0
DO i=1,n
        er = er + DABS((D1(i)-D1B(i))/D1(1))   
END DO

er=er/DBLE(n)
END SUBROUTINE Comparison_n
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





