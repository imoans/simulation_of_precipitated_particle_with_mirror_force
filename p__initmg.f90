!-----------------------------------------------------------------------
module p__initmg
!-----------------------------------------------------------------------
  use c__eprm  ,only: lxa1
  use v__in    ,only: pi,tpi
  implicit none
!
  real(kind=8)   ,save :: delx
!
  private
  public :: initmg,delx
!
  contains
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
  subroutine initmg( bf,pbr,alt,fce,lat,hh )
!-----------------------------------------------------------------------
    real(kind=8),intent(out) ::  bf(3,lxa1),pbr(lxa1)
    real(kind=8),intent(out) ::  alt(lxa1),fce(lxa1),lat(lxa1),hh(lxa1)
!
    real(kind=8),parameter ::  me = 9.10939d-31
    real(kind=8),parameter ::  RE = 6.378d6
    real(kind=8),parameter ::  R0 = RE + 5.d5
    real(kind=8),parameter ::  bb_eq = 1.4398d-7   ! [nT]@L=6
    real(kind=8),parameter ::  L = 6.d0
!    real(kind=8),parameter ::  bb_eq = 4.8594d-7   ! [nT]@L=4
!    real(kind=8),parameter ::  L = 4.d0
    real(kind=8),parameter ::  qq = 1.602d-19
    real(kind=8),parameter ::  vc = 2.99792458d8
!
    real(kind=8)    ::  Brad,Bram,fact1,rr0,cos2_th0,th0,yy0,sm
    real(kind=8)    ::  bb0,wce,slen,ds,th,rr,factp,dth,ss
    integer(kind=4) ::  ii,ind,lxb1,lxb2
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
    factp = 180.d0/pi
    rr0 = R0/RE
!-----------------------------------------------------------------------
    cos2_th0 = rr0/L
    th0 = acos( sqrt(cos2_th0) )
!### th0 : latitude (0 deg at the magnetic equator)
!-----------------------------------------------------------------------
    yy0 = sqrt(3.d0)*sin(th0)
    fact1 = log(yy0+sqrt(1.d0+yy0**2)) + yy0*sqrt(1.d0+yy0**2)
    sm = rr0/2.d0/sqrt(3.d0)/cos2_th0*fact1
!### rr0 : altitude of the base point of the field line (500 km)
!### sm  : length of field line from rr0 to magnetic equator
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
    fact1 = cos(th0)**6
    fact1 = 1.d0/fact1
    Brad = -2.d0*sin(th0)*fact1
    Bram = cos(th0)*fact1
    bb0 = sqrt( Brad**2 + Bram**2 )*bb_eq
    wce = qq*bb0/me
    slen = vc/wce
!-----------------------------------------------------------------------
    rr0 = rr0*RE/slen
    sm  = sm*RE/slen
!-----------------------------------------------------------------------
    ds = 1.d2/slen
    delx = ds
    th = th0
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
    do ii=1,lxa1
!-----------------------------------------------------------------------
      ss = sm + ds*dble(ii)
      dth = ds*cos2_th0/rr0/sqrt(1.d0+3.d0*sin(th)**2)/cos(th)
      dth = dth
      if(ii.eq.1) dth = 0.d0
      th = th + dth
!-----------------------------------------------------------------------
      rr = rr0*cos(th)**2/cos2_th0
!-----------------------------------------------------------------------
      fact1 = cos(th)**6
      fact1 = 1.d0/fact1
!-----------------------------------------------------------------------
      Brad = -2.d0*sin(th)*fact1
      Bram = cos(th)*fact1
      bf(1,ii) = sqrt( Brad**2 + Bram**2 )*bb_eq
      fce(ii) = qq*bf(1,ii)/me/tpi
      lat(ii) = th*factp
!      hh(ii)  = ds*dble(ii)*slen
      hh(ii)  = ds*dble(ii)
      alt(ii) = (rr*slen - RE)*1.d-3
!-----------------------------------------------------------------------
    end do
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
    bb0 = bf(1,1)
    do ii=1,lxa1
      bf(1,ii) = bf(1,ii)/bb0
      write(91,'(I10,5(1PD12.4))') ii,alt(ii),bf(1,ii),fce(ii),lat(ii),hh(ii)
    end do
!-----------------------------------------------------------------------
    do ii=2,lxa1-1
      pbr(ii) = -0.5d0*( bf(1,ii+1)-bf(1,ii-1) )/2.d0/delx
    end do
    pbr(1) = pbr(2)
    pbr(lxa1) = pbr(lxa1-1)
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
    return
!-----------------------------------------------------------------------
  end subroutine initmg
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
end module p__initmg
