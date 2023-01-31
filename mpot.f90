program mpot
use potentials
implicit none 
real(kind=8) :: dd, k1, k2, hbarc
real(kind=8) :: a, b, c
real(kind=8), allocatable :: x1(:), x2(:), x3(:), w1(:), w2(:), w3(:)
real(kind=8) :: k1i, k1j, k2i, k2j, k3i, k3j
real(kind=8) :: v11, v12, v13, v22, v23, v33
integer :: i, j, n

n=1000
a=5d0
b=10d0
c=100d0

hbarc = 197.327d0

allocate(x1(1:n), x2(1:n), x3(1:n))
allocate(w1(1:n), w2(1:n), w3(1:n))
CALL gauleg(0d0, a, x1, w1, n)
CALL gauleg(a, b, x2, w2, n)
CALL gauleg(b, c, x3, w3, n)

do i=1,n
!        do j=1,n
                k1i = x1(i)
!                k1j = x1(j)
!                k2i = x2(i)
!                k2j = x2(j)
!                k3i = x3(i)
!                k3j = x3(j)
!                
!                v11 = vswave(k1i,k1j)
!                v12 = vswave(k1i,k2j)
!                v13 = vswave(k1i,k3j)
!                v22 = vswave(k2i,k2j)
!                v23 = vswave(k2i,k3j)
!                v33 = vswave(k3i,k3j)
!                
!                write(11,*) i, j, v11
!                
!                write(12,*) i, j, v12
!                
!                write(13,*) i, j, v13
!                
!                write(22,*) i, j, v22
!                
!                write(23,*) i, j, v23
!                
!                write(33,*) i, j, v33
!        END DO
        write(*,*) i, k1i*hbarc,k1i, v_swave(k1i,k1i)
end do

end program mpot
