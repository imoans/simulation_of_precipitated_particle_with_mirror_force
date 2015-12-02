program main

   use mtmod
   implicit none

   integer, parameter      :: no=1000
   real(8), dimension(0:7) :: r
   integer j,k
!   real(8) grnd

   call sgrnd(4357)


   do j=0,no-1
     r(mod(j,8))=grnd()
     if(mod(j,8).eq.7) then
       write(*,'(8(f9.6,'' ''))') (r(k),k=0,7)
     else if(j.eq.no-1) then
       write(*,'(8(f9.6,'' ''))') (r(k),k=0,mod(no-1,8))
     endif
   enddo

end program main
