PROGRAM PrintV2K
USE constants
USE BCS_eq

IMPLICIT NONE   
real(kind=8), DIMENSION(0:lk) :: v2k, u2k
real(kind=8) :: L_g, SC, N_g 
integer :: ii
CALL Degeneracy(M)
N_g = 850
L_g = (N_g/rho_0)**(1d0/3d0)
SC = 2d0 * PI / L_g
CALL SOL_filler_ToUse(INT(N_g),v2k,u2k)
OPEN(UNIT=1,FILE='v2k850print.txt',STATUS='REPLACE')
DO ii=0,lk
        WRITE(1,*) M(ii,1), v2k(ii)
END DO
CLOSE(UNIT=1,STATUS='KEEP')

N_g = 10200
L_g = (N_g/rho_0)**(1d0/3d0)
SC = 2d0 * PI / L_g
CALL SOL_filler_ToUse(INT(N_g),v2k,u2k)
OPEN(UNIT=1,FILE='v2k10200print.txt',STATUS='REPLACE')
DO ii=0,lk
        WRITE(1,*) M(ii,1), v2k(ii)
END DO
CLOSE(UNIT=1,STATUS='KEEP')

WRITE(*,*) k_F, SC * M(lk,1)

END PROGRAM PrintV2K



