module vn3lo
implicit real(kind=8) (a-h,o-z)
contains

include 'n3lo450.f90'


real(kind=8) function n3lo_1s0(k1,k2)
real(kind=8), intent(in) :: k1, k2
common /crdwrt/ kread,kwrite,kpunch,kda(9)
common /cpot/   v(6),xmev,ymev
common /cstate/ j,heform,sing,trip,coup,endep,label
common /cnn/ inn
logical heform,sing,trip,coup,endep
character*4 label
pi = acos(-1d0)
hbarc = 197.327053d0

xmev = k1*hbarc
ymev = k2*hbarc

kread = 5
kwrite = 6
!kpunch
!kda(9)
j=0     !1S[0]
heform = .false.
sing = .true.
trip = .true.
coup = .true.
!endep
!label
inn = 3         !neutron-neutron
call n3lo450new
n3lo_1s0 = 0.5*pi*v(1)*hbarc**3d0
end function n3lo_1s0
end module vn3lo
