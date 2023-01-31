MODULE Math_Functions 
IMPLICIT NONE

CONTAINS

!______________________________________________________
integer FUNCTION real_D_Kron(x,y)   !=1 if x=y
IMPLICIT NONE 
real(kind=8), INTENT(IN) :: x,y
real(kind=8), parameter :: tolerance = 10d0**(-9)

IF(DABS(x-y) < tolerance) THEN
        real_D_Kron = 1
ELSE
        real_D_Kron =0
END IF
END FUNCTION real_D_Kron
!------------------------------------------------------------
integer FUNCTION D_Kron(i,j)
IMPLICIT NONE
integer, INTENT(IN) :: i, j

IF (i.EQ.j) THEN
        D_Kron = 1
ELSE
        D_Kron = 0
END IF
END FUNCTION D_Kron
!----------------------------------------------------------------------------------------
integer FUNCTION D3(i,j,k)      !Given a vector (i,j,k), D3(i,j,k) =  # of equal numbers
                                !For example D(i,i,k) = 2, D(i,j,i) = 2, D(i,i,i) = 3, D(i,j,k) = 0
IMPLICIT NONE
integer, INTENT(IN) :: i, j, k

D3 = 2 * D_Kron(i,j) / (1 + D_Kron(j,k) * D_Kron(k,i))&
   + 2 * D_Kron(j,k) / (1 + D_Kron(k,i) * D_Kron(i,j))&
   + 2 * D_Kron(k,i) / (1 + D_Kron(i,j) * D_Kron(j,k))
END FUNCTION D3

END MODULE Math_Functions





