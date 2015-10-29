!-----------------------------------------------------------------------
module c__eprm
  implicit none
!-----------------------------------------------------------------------
!  integer(kind=4),parameter :: istep = 100000
  integer(kind=4),parameter :: istep = 10000000
!-----------------------------------------------------------------------
  integer(kind=4),parameter :: lx = 5000 , lxa1 = lx+4
  integer(kind=4),parameter :: lpn1 = 256
  integer(kind=4),parameter :: npa1 = 7
!  integer(kind=4),parameter :: npa1 = lpn1*lx
  integer(kind=4),parameter :: nsp = 2 , nqa = 2
  integer(kind=4),parameter :: lx1 = 3 , lx2 = lx+2
!-----------------------------------------------------------------------
end module c__eprm
