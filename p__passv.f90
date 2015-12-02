!-----------------------------------------------------------------------
module p__passv
!-----------------------------------------------------------------------
  use c__eprm,  only: lxa1,npa1
  use v__in,    only: delt,q
  use p__initmg,only: delx
  use CyclotronWithRungeKutta
  implicit none
!
  private
  public :: passv
!
  contains
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
  subroutine passv( xx,vp1,pm1,lor,bf,ef,pbr )
!-----------------------------------------------------------------------
    real(kind=8)   ,intent(inout) ::  vp1(3,npa1)
    real(kind=8)   ,intent(inout) ::  pm1(3,npa1),lor(npa1)
    real(kind=8)   ,intent(in)    ::  xx(npa1)
    real(kind=8)   ,intent(in)    ::  bf(3,lxa1),ef(3,lxa1),pbr(lxa1)
!
    real(kind=8)  ::  aa,raa1,dd1,dd2
    real(kind=8) ::   tbx,tby,tbz,fbx,fby,fbz,fex,fey,fez
    real(kind=8) ::   fpbr,pperp,rlar,tpbr,tpby,tpbz,vay,vaz
    real(kind=8) ::   lort,hdelt,cnst,factb
    real(kind=8) ::   pax,pay,paz,p1x,p1y,p1z,p0x,p0y,p0z
    real(kind=8) ::   half, fact1, fact2, fact3
    integer           ii,ic
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
    half  = 0.5d0
    fact1 = q(2,1)
    fact2 = q(2,2)
    fact3 = fact2/fact1
    hdelt = delt*half
    cnst = fact3*hdelt
!-----------------------------------------------------------------------
    do ii=1,npa1
!-----------------------------------------------------------------------
!      if(pind(ii).eq.0) then  ! 計算領域に粒子が存在するかの判定
!-----------------------------------------------------------------------
        aa = xx(ii) / delx
        raa1 = dint( aa+half )
        ic = raa1
        dd1 = aa+half-raa1
        dd2 = raa1-aa+half
!-----------------------------------------------------------------------
        fex = ef(1,ic)*dd1 + ef(1,ic-1)*dd2
        fey = ef(2,ic)*dd1 + ef(2,ic-1)*dd2
        fez = ef(3,ic)*dd1 + ef(3,ic-1)*dd2
!-----------------------------------------------------------------------
        fbx = bf(1,ic)*dd1 + bf(1,ic-1)*dd2
        fby = bf(2,ic)*dd1 + bf(2,ic-1)*dd2
        fbz = bf(3,ic)*dd1 + bf(3,ic-1)*dd2
        fpbr = pbr(ic)*dd1 + pbr(ic-1)*dd2
!-----------------------------------------------------------------------
        pax = pm1(1,ii) + fex*cnst
        pay = pm1(2,ii) + fey*cnst
        paz = pm1(3,ii) + fez*cnst
!-----------------------------------------------------------------------
        lort = 1.d0 + pax**2 + pay**2 + paz**2
        lort = sqrt(lort)
!-----------------------------------------------------------------------
        vay = pay/lort/q(2,1)
        vaz = paz/lort/q(2,1)
!-----------------------------------------------------------------------
        rlar = sqrt( vay**2 + vaz**2 )
        rlar = rlar/fbx*lort
!-----------------------------------------------------------------------
        tbx = fbx/lort
        tby = fby/lort
        tbz = fbz/lort
        tpbr = fpbr*rlar/lort
!-----------------------------------------------------------------------
        pperp = sqrt( pay**2 + paz**2 )
        tpby =  tpbr*paz/pperp
        tpbz = -tpbr*pay/pperp
!-----------------------------------------------------------------------
        p0x = pax + cnst*( pay*tbz - paz*tby - pperp*tpbr )
        p0y = pay + cnst*( paz*tbx - pax*tbz - pax*tpbz )
        p0z = paz + cnst*( pax*tby - pay*tbx + pax*tpby )
!-----------------------------------------------------------------------
        factb = tbx**2 + (tby+tpby)**2 + (tbz+tpbz)**2
        factb = 2.d0/( 1.d0 + cnst**2*factb )
!-----------------------------------------------------------------------
        pperp = sqrt( p0y**2 + p0z**2 )
        tpby =  tpbr*p0z/pperp
        tpbz = -tpbr*p0y/pperp
!-----------------------------------------------------------------------
        p1x = pax + factb*cnst*( p0y*tbz - p0z*tby - pperp*tpbr )
        p1y = pay + factb*cnst*( p0z*tbx - p0x*tbz - p0x*tpbz )
        p1z = paz + factb*cnst*( p0x*tby - p0y*tbx + p0x*tpby )
!-----------------------------------------------------------------------
        pm1(1,ii) = p1x + cnst*fex
        pm1(2,ii) = p1y + cnst*fey
        pm1(3,ii) = p1z + cnst*fez
!-----------------------------------------------------------------------
        lor(ii) = 1.d0 + pm1(1,ii)**2 + pm1(2,ii)**2 + pm1(3,ii)**2
        lor(ii) = sqrt( lor(ii) )
!-----------------------------------------------------------------------
        vp1(1,ii) = pm1(1,ii)/lor(ii)/q(2,1)
        vp1(2,ii) = pm1(2,ii)/lor(ii)/q(2,1)
        vp1(3,ii) = pm1(3,ii)/lor(ii)/q(2,1)
!-----------------------------------------------------------------------
!      end if

!-----------------------------------------------------------------------
! rungekutta
        !call run(1.7d-2)

        !vp1(1, ii) = v + rungekutta(v, tbx, acceleration(v, t, fex, fbx))
        !vp1(2, ii) = v + rungekutta(v, tby, acceleration(v, t, fey, fby))
        !vp1(2, ii) = v + rungekutta(v, tbz, acceleration(v, t, fez, fbz))
!-----------------------------------------------------------------------
    end do
!-----------------------------------------------------------------------
    return
!-----------------------------------------------------------------------
  end subroutine passv
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
end module p__passv
