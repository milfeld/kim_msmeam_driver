!
! CDDL HEADER START
!
! The contents of this file are subject to the terms of the Common Development
! and Distribution License Version 1.0 (the "License").
!
! You can obtain a copy of the license at
! http://www.opensource.org/licenses/CDDL-1.0.  See the License for the
! specific language governing permissions and limitations under the License.
!
! When distributing Covered Code, include this CDDL HEADER in each file and
! include the License file in a prominent location with the name LICENSE.CDDL.
! If applicable, add the following below this CDDL HEADER, with the fields
! enclosed by brackets "[]" replaced with your own identifying information:
!
! Portions Copyright (c) [yyyy] [name of copyright owner]. All rights reserved.
!
! CDDL HEADER END
!

!
! Copyright (c) 2013--2014, Regents of the University of Minnesota.
! All rights reserved.
!
! Contributors:
!    Ryan S. Elliott
!    Ellad B. Tadmor
!    Valeriu Smirichinski
!    Stephen M. Whalen
!

!
! Portions Copyright (c) 2015, Joshua S. Gibson.
! All rights reserved.
!
! Contributors:
!    Joshua S. Gibson
!    Srinivasan G. Srivilliputhur
!    Michael I. Baskes
!    Ronald E. Miller
!    Angela K. Wilson
!

!****************************************************************************
!**
!**  MODULE MSMEAM
!**
!**  MSMEAM model driver
!**
!**  Reference: [PENDING]
!**
!**  Language: Fortran 2003
!**
!****************************************************************************


#include "KIM_API_status.h"
#define THIS_FILE_NAME __FILE__
#define TRUEFALSE(TRUTH) merge(1,0,(TRUTH))

module MSMEAM

use, intrinsic :: iso_c_binding
use KIM_API_F03
implicit none

save
private
public Compute_Energy_Forces, &
       model_parameters,      &
       destroy,               &
       elecden,               &
       rscrn,                 &
       get_c,                 &
       calc_phi,              &
       rhof,                  &
       drhof,                 &
       dphif,                 &
       frhoi,                 &
       dfrhoi,                &
       get_s,                 &
       dcmij,                 &
       lagrangefxn,           &
       lagrangedrv

type model_parameters
     ! SRINI: Gibson defined them in "subroutine rscrn(rij,meatable,scrn)"
     real(c_double) :: frcut

     ! SRINI: Move from subroutine get_c()
     real(c_double) :: ctabmax

   integer(c_int)                                    :: NBC
   integer(c_int)                                    :: HalfOrFull
   integer(c_int)                                    :: IterOrLoca
   integer(c_int)                                    :: nel
   integer(c_int)                                    :: ntab
   real(c_double)                                    :: rhotabmin
   real(c_double)                                    :: rhotabmax
   real(c_double)                                    :: rtabmax
   real(c_double)                                    :: rrtabmax
   real(c_double)                                    :: model_cutoff
   real(c_double)                                    :: model_cutsq

   ! SRINI:  Radial Screening Factor computed and saved here
   ! frcutRcSq = frcut*rcutsq; inv1minusFrcutRcSq = 1/((1.0-frcut)*rcutsq)
   real(c_double)                                    :: frcutRcSq, inv1minusFrcutRcSq

   ! SRINI: rtabscale = ( meatable%ntab - 1.0_cd ) / meatable%rrtabmax
   real(c_double)                                    :: rtabscale, ctabscale, rhotabscale
   
   real(c_double), pointer, dimension(:)             :: asubs
   real(c_double), pointer, dimension(:)             :: bsubs
   real(c_double), pointer, dimension(:,:)           :: esubs
   real(c_double), pointer, dimension(:,:)           :: res
   real(c_double), pointer, dimension(:,:)           :: embedtab
   real(c_double), pointer, dimension(:,:,:)         :: phitab
   real(c_double), pointer, dimension(:,:,:,:)       :: rhotab
   real(c_double), pointer, dimension(:,:,:,:)       :: scrp
   real(c_double), pointer, dimension(:,:,:,:)       :: scr0
   real(c_double), pointer, dimension(:,:,:,:)       :: scr1
   real(c_double), pointer, dimension(:,:,:,:)       :: scr2
   real(c_double), pointer, dimension(:,:,:,:)       :: scr3
   real(c_double), pointer, dimension(:,:,:,:)       :: scr54
end type model_parameters

! Below are the definitions and values of all model parameters
integer(c_int), parameter :: cd = c_double  ! used for literal constants
integer(c_int), parameter :: DIM = 3        ! dimensionality of space
integer(c_int), parameter :: speccode = 1   ! internal species code

! SRINI: GIBSON had xmcut/xncut as reals.  Taking power with real exponents is slow
!        So, I changed these to integers
integer(c_int), parameter :: xncut=2
integer(c_int), parameter :: xmcut=4

!-------------------------------------------------------------------------------
! Below are the definitions and values of all additional model parameters
!
! Recall that the Fortran 2003 format for declaring parameters is as follows:
!
! integer(c_int), parameter :: parname = value ! This defines an
!                                              ! integer(c_int) parameter
!                                              ! called `parname' with a value
!                                              ! equal to `value' (a number)
!
! real(c_double), parameter :: parname = value ! This defines a real(c_double)
!                                              ! parameter called `parname'
!                                              ! with a value equal to
!                                              ! `value' (a number)
!-------------------------------------------------------------------------------
integer(c_int), parameter ::  natmax=30000  ! Gibson value= 10000
integer(c_int), parameter ::  nelmax=2
integer(c_int), parameter ::  ntable=10001
integer(c_int), parameter ::  ifree=0
integer(c_int), parameter ::  noscr=0
real(c_double), parameter ::  scnres=0.0_cd
real(c_double), parameter ::  legend00=3./5._cd
real(c_double), parameter ::  legend01=10./9._cd
real(c_double), parameter ::  legend02=5./21._cd

contains

!-------------------------------------------------------------------------------
!
!  Calculate lagrange interpolation
!
!-------------------------------------------------------------------------------
real(c_double) function lagrangefxn ( p, yd )
!DEC$ ATTRIBUTES FORCEINLINE :: lagrangefxn
implicit none

integer(c_int), parameter                    :: n = 4
real(c_double), intent(in)                   :: p
real(c_double), intent(in), dimension(1:n)   :: yd
real(c_double)                               :: c1,c2,c3,c4
real(c_double)                               :: pm1,pp1,pm2
! real(c_double)                               :: lazy
! integer(c_int)                               :: j

pp1=p+1.0_cd; pm1=p-1.0_cd; pm2=p-2.0_cd
!c1 = -0.166666667_cd*p*           (p-1.0_cd)*(p-2.0_cd)
!2 =          0.5_cd*  (p+1.0_cd)*(p-1.0_cd)*(p-2.0_cd)
!3 =         -0.5_cd*p*(p+1.0_cd)*           (p-2.0_cd)
!4 =  0.166666667_cd*p*(p+1.0_cd)*(p-1.0_cd)

! SRINI: COMMENT THIS AND INLINE
! lazy=0.0_cd
!do j=1,4
!   lazy=lazy+cof(j)*yd(j)
!enddo
! lagrangefxn=lazy

! SRINI: WILL THIS PEED THE CODE?
!  lagrangefxn = cof(1)*yd(1) + cof(2)*yd(2) + cof(3)*yd(3) + cof(4)*yd(4)    
   lagrangefxn = (yd(4)*pp1 - yd(1)*pm2)*0.166666667_cd    *pm1* p     + &
                 (yd(2)*pm1 - yd(3)*p  )*0.5_cd        *pm2       *pp1

return
end function lagrangefxn

!-------------------------------------------------------------------------------
!
!  Calculate derivative of lagrange interpolation
!
!-------------------------------------------------------------------------------
real(c_double) function lagrangedrv ( xt, yd )
implicit none

integer(c_int), parameter                    :: ione = 1
integer(c_int), parameter                    :: n = 4
real(c_double), intent(in)                   :: xt
real(c_double)                               :: lazy
real(c_double), intent(in), dimension(1:n)   :: yd
real(c_double), parameter, dimension(1:n)    :: xd = &
 & [-1,0,1,2]
real(c_double), parameter, dimension(1:n)    :: q = &
 & [-6,2,-2,6]

! SRINI: qinv = 1/q (in line above)
real(c_double), parameter, dimension(1:n)    :: qinv = &
 & [-1.0/6.0_cd,0.5_cd, -0.5_cd, 1.0/6.0_cd]

real(c_double), dimension(1:n)               :: d,pa,pb,dpa,dpb,dw
integer(c_int)                               :: j

pa(1)=1.0_cd
dpa(1)=0.0_cd

do j=1,n-ione
   d(j)=xt-xd(j)
   pa(j+ione) = pa(j)*d(j)
   dpa(j+ione)= dpa(j)*d(j)+pa(j)
enddo

d(n)=xt-xd(n)

pb(n) = 1.0_cd
dpb(n)= 0.0_cd

do j=n,2,-1
   pb (j-ione)= pb(j)*d(j)
   dpb(j-ione)= dpb(j)*d(j)+pb(j)
enddo

lazy=0.0_cd

do j=1,n
   !   dw(j)= (pa(j)*dpb(j)+dpa(j)*pb(j))/q(j)

    ! SRINI
    dw(j)= (pa(j)*dpb(j)+dpa(j)*pb(j))*qinv(j)

    lazy = lazy+dw(j)*yd(j)
enddo

lagrangedrv=lazy

return
end function lagrangedrv

!-------------------------------------------------------------------------------
!
!  Calculate screening functions
!
!-------------------------------------------------------------------------------
! subroutine get_s(cx,it,jt,kt,meatable,sp,s0,s1,s2,s3,s54)
! GIBSON
! type(model_parameters), intent(in)  :: meatable
! real(c_double), intent(out)       :: sp,s0,s1,s2,s3,s54
subroutine get_s(cx,it,jt,kt,meatable,sParArray)
real(c_double), intent(in)          :: cx
integer(c_int), intent(in)          :: it,jt,kt
! SRINI2
type(model_parameters), pointer, intent(in) :: meatable
real(c_double),dimension(6),intent(out):: sParArray

real(c_double)                      :: sp,s0,s1,s2,s3,s54
real(c_double)                      :: p
! integer(c_int)                      :: k, ctr, in1
integer(c_int)                      :: k
integer(c_int), parameter           :: Np=4
integer(c_int), parameter           :: Pr=Np/2
real(c_double),dimension(1:Np)      :: sparray,s0array,s1array
real(c_double),dimension(1:Np)      :: s2array,s3array,s5array

! SRINI
real(c_double)                      :: ctabmax
ctabmax = meatable%ctabmax

if(cx.ge.ctabmax) then
     ! GIBSON
     !sp=1.0_cd
     !s0=1.0_cd
     !s1=1.0_cd
     !s2=1.0_cd
     !s3=1.0_cd
     !s54=1.0_cd

     ! SRINI2
     sParArray = 1.0_cd
     return
endif

! GIBSON
! p = cx * ( dble(meatable%ntab) - 1.0_cd ) / ctabmax + 1.0_cd

! SRINI
! meatable%rtabscale          = ( dble(meatable%ntab) - 1.0_cd )/meatable%rrtabmax
! meatable%ctabscale          = ( dble(meatable%ntab) - 1.0_cd )/ctabmax
p = cx * meatable%ctabscale + 1.0_cd

k = int(p)

if ( k .gt. meatable%ntab - 1 ) then
   ! GIBSON
   ! sp = 0.0_cd
   ! s0 = 0.0_cd
   ! s1 = 0.0_cd
   ! s2 = 0.0_cd
   ! s3 = 0.0_cd
   ! s54 = 0.0_cd

   ! SRINI2
   !DIR$ vector nontemporal
   !DIR$ simd
   sParArray = 0.0_cd
   return
endif

k = min( k , meatable%ntab - Pr )
k = max( k , Pr )
p = p - dble( k )

!   make sure that p is less than Pr
!   then if cx is out of range, Pr and c = last value of screening
!

p = min( p, dble(Pr) )
      
if( 1 == 1) then
! SRINI2: DELETED ALL SRINI1 code BELOW and replaced with this array code
sparray(:)=meatable%scrp( (k-Pr+1):(k-Pr+Np), it, jt, kt)
s0array(:)=meatable%scr0( (k-Pr+1):(k-Pr+Np), it, jt, kt)
s1array(:)=meatable%scr1( (k-Pr+1):(k-Pr+Np), it, jt, kt)
s2array(:)=meatable%scr2( (k-Pr+1):(k-Pr+Np), it, jt, kt)
s3array(:)=meatable%scr3( (k-Pr+1):(k-Pr+Np), it, jt, kt)
s5array(:)=meatable%scr54((k-Pr+1):(k-Pr+Np), it, jt, kt)

sp=lagrangefxn(p,sparray)                         !31 seconds start
s0=lagrangefxn(p,s0array)
s1=lagrangefxn(p,s1array)
s2=lagrangefxn(p,s2array)
s3=lagrangefxn(p,s3array)
s54=lagrangefxn(p,s5array)
! SRINI2: DELETED ALL SRINI1 code BELOW and replaced with this array code

else

sp=lagrangefxn(p,meatable%scrp(  (k-Pr+1):(k-Pr+Np), it, jt, kt))
s0 =lagrangefxn(p,meatable%scr0( (k-Pr+1):(k-Pr+Np), it, jt, kt))
s1 =lagrangefxn(p,meatable%scr1( (k-Pr+1):(k-Pr+Np), it, jt, kt))
s2 =lagrangefxn(p,meatable%scr2( (k-Pr+1):(k-Pr+Np), it, jt, kt))
s3 =lagrangefxn(p,meatable%scr3( (k-Pr+1):(k-Pr+Np), it, jt, kt))
s54=lagrangefxn(p,meatable%scr54((k-Pr+1):(k-Pr+Np), it, jt, kt))
endif


sp = min(1.0_cd, max( 0.0_cd , sp))
s0 = min(1.0_cd, max( 0.0_cd , s0))
s1 = min(1.0_cd, max( 0.0_cd , s1))
s2 = min(1.0_cd, max( 0.0_cd , s2))
s3 = min(1.0_cd, max( 0.0_cd , s3))
s54= min(1.0_cd, max( 0.0_cd , s54))             !31 seconds  end

sParArray(1) = sp 
sParArray(2) = s0
sParArray(3) = s1
sParArray(4) = s2
sParArray(5) = s3
sParArray(6) = s54

return
end subroutine get_s
!-------------------------------------------------------------------------------
!
!  Calculate background electron density
!
!-------------------------------------------------------------------------------
subroutine elecden(rho0,AA,dang1,dang2,rhobar)
implicit none

real(c_double)                 :: temp
real(c_double),intent(in)      :: rho0,AA
real(c_double),intent(out)     :: dang1,dang2,rhobar

!   rhobar^2 = rho0^2+AA
!   dang is derivative of rhobar wrt rho0 and  AA

! GIBSON
!    temp = rho0**2 + AA

! SRINI
    temp = (rho0*rho0) + AA

if(temp.le.0.0_cd) then
    rhobar=1d-20
    dang1=0.0_cd
    dang2=0.0_cd
  else
    rhobar=sqrt(temp)
    dang1= rho0/rhobar
    dang2= 0.5_cd/rhobar
endif

return
end subroutine elecden
!-------------------------------------------------------------------------------
!
!  Calculate radial screening
!
!-------------------------------------------------------------------------------
subroutine rscrn(rij,meatable,scrn)
!   calculate r cut-off function
implicit none

! GIBSON
!type(model_parameters), intent(in)  :: meatable

! SRINI: ZZ APRIL 20, 2016
type(model_parameters), pointer, intent(in) :: meatable


real(c_double), intent(in)  :: rij
real(c_double), intent(out) :: scrn
real(c_double)              :: fac,fac2,fac3

! SRINI:  Radial Screening Factor computed and saved here
! fac2 = frcutRcSq = frcut*rcutsq; fac3 = inv1minusFrcutRcSq = 1/((1.0-frcut)*rcutsq)
fac2 = meatable%frcutRcSq
fac3 = meatable%inv1minusFrcutRcSq
fac=(rij-fac2)*fac3

!!!!!!!!!!!!!!!!!!!!!!!!!!
! REMOVED GIBSON'S CODE
!!!!!!!!!!!!!!!!!!!!!!!!!!

if(fac .le. 0.0_cd) then
    scrn=1.0_cd
elseif(fac .ge. 1.0_cd) then
    scrn=0.0_cd
else
    scrn=(1.0_cd-(fac**xmcut))**xncut
endif

return
end subroutine rscrn
!-------------------------------------------------------------------------------
!
!  Calculate screening parameter
!
!-------------------------------------------------------------------------------
real(c_double) function get_c(rij,ril,rjl)
implicit none

real(c_double), intent(in) :: rij,ril,rjl
real(c_double) :: dotil,dotjl,argb,argt,Xi,Xj

! SRINI
real(c_double) :: XiMinXjSQ, sconstrij, rijfactor, dXiXj

real(c_double),parameter :: sconst=1.125

get_c=3.0_cd

! ORIG GIBSON'S CODE
! if ((ril.gt.sconst*rij).or.(rjl.gt.sconst*rij)) return
!
!   If l atom is outside of ellipse there is no screening
!   dotx = 2 r_ij dot r_x
!   this eliminates screeners with distances greater than 0.2*rij
!
!!dotil=ril+rij-rjl
!!if(dotil.le.0.027_cd*rij) return
!!dotjl=rij+rjl-ril
!!if(dotjl.le.0.027_cd*rij) return
! ORIG GIBSON'S CODE ABOVE


! SRINI: REDUCE multiplication within if()
sconstrij = sconst*rij
if ((ril .gt. sconstrij) .or. (rjl .gt. sconstrij)) return

!   If l atom is outside of ellipse there is no screening
!   dotx = 2 r_ij dot r_x
!   this eliminates screeners with distances greater than 0.2*rij
rijfactor = 0.027_cd*rij

dotil=ril+rij-rjl
if(dotil .le. rijfactor) return

dotjl=rij+rjl-ril
if(dotjl .le. rijfactor) return
!!!!!! SRINI-- ABOVE

Xi=(ril/rij)
Xj=(rjl/rij)

! ORIG
!argb=1.0_cd-(Xi-Xj)**2
!argt=2.0_cd*(Xi+Xj)-1.0_cd-(Xi-Xj)**2
!get_c=argt/argb

! SRINI
dXiXj     = (Xi-Xj)
XiMinXjSQ = dXiXj*dXiXj
argb = 1.0_cd - XiMinXjSQ
argt = 2.0_cd*(Xi+Xj) - 1.0_cd - XiMinXjSQ
get_c=argt/argb

return
end function get_c
!-------------------------------------------------------------------------------
!
!  Calculate pair potential phi(i) for atom i
!
!-------------------------------------------------------------------------------
real(c_double) function calc_phi(r,it,jt,meatable, pindx)
implicit none

!-- Transferred variables
real(c_double), intent(in)          :: r, pindx
type(model_parameters), intent(in)  :: meatable
integer(c_int), intent(in)          :: it,jt

!-- Local variables
integer(c_int), parameter           :: Np=4
integer(c_int), parameter           :: Pr=Np/2
real(c_double)                      :: p
real(c_double),dimension(1:Np)      :: phiarray
integer(c_int)                      :: k      ! ,ctr

! SRINI
integer(c_int)                      :: kMinusPr

! GIBSON
! p = r * ( meatable%ntab - 1.0_cd ) / meatable%rrtabmax + 1.0_cd

! SRINI
! meatable%rtabscale          = ( dble(meatable%ntab) - 1.0_cd )/meatable%rrtabmax
! meatable%ctabscale          = ( dble(meatable%ntab) - 1.0_cd )/ctabmax
! p = r * meatable%rtabscale + 1.0_cd
! Added as an argument
p = pindx

k = int(p)

if ( k .gt. meatable%ntab - 1 ) then
     calc_phi = 0.0_cd
     return
endif

k = min(k,meatable%ntab - Pr)
k = max(k,Pr)
p = p - dble(k)

!       make sure that p is less than 2.0
!       then if r is out of range, p = 2.0 and r = last value of phitab
!

p = min(p,dble(Pr))

! SRINI
kMinusPr = k - Pr

!do ctr=1,Np
   ! GIBSON
   ! phiarray(ctr)= meatable%phitab( k - Pr + ctr , it, jt )

   ! SRINI
!   phiarray(ctr)= meatable%phitab(kMinusPr + ctr , it, jt)
!enddo

! SRINI3
phiarray(:)=meatable%phitab( (kMinusPr+1):(kMinusPr+Np), it, jt)


calc_phi = lagrangefxn(p,phiarray)

return
end function calc_phi
!-------------------------------------------------------------------------------
!
!  Calculate rhol+ and rhol- at a given r with atom j
!
!-------------------------------------------------------------------------------
real(c_double) function rhof(r,l,ispin,it,meatable, pindx)
implicit none

!-- Transferred variables
type(model_parameters), intent(in)  :: meatable
integer(c_int),intent(in)           :: l
integer(c_int),intent(in)           :: ispin
integer(c_int),intent(in)           :: it
real(c_double),intent(in)           :: r, pindx

!-- Local variables
integer(c_int), parameter           :: Np=4
integer(c_int), parameter           :: Pr=Np/2
real(c_double)                      :: p
real(c_double),dimension(1:Np)      :: rhoarray
integer(c_int)                      :: k      !,ctr

! SRINI
integer(c_int)                      :: kMinusPr

! GIBSON
! p = r * ( meatable%ntab - 1.0_cd ) / meatable%rrtabmax + 1.0_cd

! SRINI
! meatable%rtabscale          = ( dble(meatable%ntab) - 1.0_cd )/meatable%rrtabmax
! meatable%ctabscale          = ( dble(meatable%ntab) - 1.0_cd )/ctabmax
! p = r * meatable%rtabscale + 1.0_cd
p = pindx

k = int(p)

if ( k .gt. meatable%ntab - 1 ) then
     rhof = 0.0_cd
     return
endif

k = min ( k , meatable%ntab - Pr )
k = max ( k , Pr )
p = p - dble(k)

!       make sure that p is less than 2.0
!       then if r is out of range, p = 2.0 and r = last value of phitab
!
p = min ( p , dble(Pr))

! SRINI
kMinusPr = k - Pr

!do ctr=1,Np
   ! GIBSON
   ! rhoarray(ctr)=meatable%rhotab( k - Pr + ctr , l , ispin , it )

   ! SRINI
!   rhoarray(ctr)=meatable%rhotab(kMinusPr + ctr , l , ispin , it )
!enddo


! SRINI3
rhoarray(:)=meatable%rhotab( (kMinusPr+1):(kMinusPr+Np), l, ispin, it)

rhof = lagrangefxn(p,rhoarray)

return
end function rhof
!-------------------------------------------------------------------------------
!
!  Calculate the derivatives of rhol+ and rhol- at a given r with atom j
!
!-------------------------------------------------------------------------------
real(c_double) function drhof(r,l,ispin,it,meatable,pindx)
implicit none

integer(c_int),intent(in)           :: l
integer(c_int),intent(in)           :: ispin
integer(c_int),intent(in)           :: it
real(c_double),intent(in)           :: r, pindx
type(model_parameters), intent(in)  :: meatable

integer(c_int), parameter           :: Np=4
integer(c_int), parameter           :: Pr=Np/2
real(c_double)                      :: p
real(c_double),dimension(1:Np)      :: rhoarray
integer(c_int)                      :: k     ! ,ctr

! SRINI
integer(c_int)                      :: kMinusPr

! GIBSON
! p = r * ( dble( meatable%ntab ) - 1.0_cd ) / meatable%rrtabmax + 1.0_cd

! SRINI
! meatable%rtabscale          = ( dble(meatable%ntab) - 1.0_cd )/meatable%rrtabmax
! meatable%ctabscale          = ( dble(meatable%ntab) - 1.0_cd )/ctabmax
! p = r * meatable%rtabscale + 1.0_cd
p = pindx

k = int(p)

if ( k .gt. meatable%ntab-1) then
     drhof = 0.0_cd
     return
endif

k = min( k , meatable%ntab - Pr )
k = max( k , Pr )
p = p - dble( k )

!       make sure that p is less than 2.0
!       then if r is out of range, p = 2.0 and r = last value of rhotab
!

p = min( p , dble(Pr) )

! SRINI
kMinusPr = k - Pr

!do ctr=1,Np
   ! GIBSON
   ! rhoarray(ctr)=meatable%rhotab( k - Pr + ctr, l, ispin, it )
   ! SRINI
!   rhoarray(ctr)=meatable%rhotab(kMinusPr + ctr, l, ispin, it )
!enddo

! SRINI3
rhoarray(:)=meatable%rhotab( (kMinusPr+1):(kMinusPr+Np), l, ispin, it)

! GIBSON
! drhof = lagrangedrv(p,rhoarray) * ( dble( meatable%ntab ) - 1.0_cd ) / meatable%rrtabmax
! SRINI
drhof = lagrangedrv(p,rhoarray) * meatable%rtabscale 

return
end function drhof
!-------------------------------------------------------------------------------
!
!  Calculate the derivative of the pair potential dphi(r)
!
!-------------------------------------------------------------------------------
real(c_double) function dphif(r,it,jt,meatable, pindx)
implicit none

!-- Transferred variables
real(c_double), intent(in)          :: r, pindx
integer(c_int), intent(in)          :: it,jt
type(model_parameters), intent(in)  :: meatable

!-- Local variables
integer(c_int), parameter           :: Np=4
integer(c_int), parameter           :: Pr=Np/2
real(c_double)                      :: p
real(c_double),dimension(1:Np)      :: phiarray
integer(c_int)                      :: k     ! ,ctr

! SRINI
integer(c_int)                      :: kMinusPr

!
!  do four-point lagrange interpolation
!

! GIBSON
! p = r * ( dble( meatable%ntab ) - 1.0_cd ) / meatable%rrtabmax  + 1.0_cd

! SRINI
! meatable%rtabscale          = ( dble(meatable%ntab) - 1.0_cd )/meatable%rrtabmax
! meatable%ctabscale          = ( dble(meatable%ntab) - 1.0_cd )/ctabmax
! p = r * meatable%rtabscale + 1.0_cd
p = pindx

k = int(p)

if (k .gt. meatable%ntab-1 ) then
     dphif = 0.0_cd
     return
endif

k = min( k , meatable%ntab - Pr )
k = max( k , Pr )
p = p - dble( k )

!       make sure that p is less than 3.0
!       then if r is out of range, p = 3.0 and r = last value of phitab
!

p = min( p , dble(Pr) )

! SRINI
kMinusPr = k - Pr

!do ctr=1,Np
   ! GIBSON
   ! phiarray(ctr)= meatable%phitab( k - Pr + ctr , it, jt )
   ! SRINI
!   phiarray(ctr)= meatable%phitab(kMinusPr + ctr , it, jt )
!enddo


! SRINI3
phiarray(:)=meatable%phitab( (kMinusPr+1):(kMinusPr+Np), it, jt)

! GIBSON
! dphif = lagrangedrv(p,phiarray) * ( dble( meatable%ntab ) - 1.0_cd ) / meatable%rrtabmax

! SRINI
! meatable%rtabscale          = ( dble(meatable%ntab) - 1.0_cd )/meatable%rrtabmax
! meatable%ctabscale          = ( dble(meatable%ntab) - 1.0_cd )/ctabmax
! meatable%rhotabscale        = ( dble(meatable%ntab) - 1.0_cd )/meatable%rhotabmax
 dphif = lagrangedrv(p,phiarray) * meatable%rtabscale

return
end function dphif
!-------------------------------------------------------------------------------
!
!  Calculate embedding function for atom i
!
!-------------------------------------------------------------------------------
real(c_double) function frhoi(rho,it,asub,bsub,esub,meatable)
implicit none

!-- Transferred variables
real(c_double), intent(in)          :: rho,asub,bsub,esub
integer(c_int), intent(in)          :: it
type(model_parameters), intent(in)  :: meatable

!-- Local variables
integer(c_int), parameter           :: Np=4
integer(c_int), parameter           :: Pr=Np/2
real(c_double)                      :: xrho,p
integer(c_int)                      :: k     ! ,ctr
real(c_double),dimension(1:Np)      :: embarray

! SRINI
integer(c_int)                      :: kMinusPr

if ( rho .le. meatable%rhotabmin ) then
     frhoi = esub * ( asub * log( rho ) + bsub ) * rho
else
     xrho = log ( rho / meatable%rhotabmin )

     ! GIBSON
     ! p= xrho * ( dble(meatable%ntab) - 1.0_cd ) / meatable%rhotabmax + 1.0_cd

     ! SRINI
     ! meatable%rtabscale          = ( dble(meatable%ntab) - 1.0_cd )/meatable%rrtabmax
     ! meatable%ctabscale          = ( dble(meatable%ntab) - 1.0_cd )/ctabmax
     ! meatable%rhotabscale        = ( dble(meatable%ntab) - 1.0_cd )/meatable%rhotabmax
     p =  xrho * meatable%rhotabscale + 1.0_cd

     k = int(p)
     k = min ( k, meatable%ntab - Pr )
     k = max ( k , Pr )
     p = p - dble( k )

     p = min( p , dble(Pr) )

     ! SRINI
     kMinusPr = k - Pr

!    do ctr=1,Np
!     ! GIBSON
     ! embarray(ctr)=meatable%embedtab( k - Pr + ctr , it ) 
     ! SRINI
!      embarray(ctr)=meatable%embedtab(kMinusPr + ctr , it ) 
!   enddo

! SRINI3
embarray(:)=meatable%embedtab( (kMinusPr+1):(kMinusPr+Np), it)

frhoi = lagrangefxn(p,embarray)

endif
return
end function frhoi

!-------------------------------------------------------------------------------
!
!  Calculate first derivative of embedding function for atom i
!
!-------------------------------------------------------------------------------
real(c_double) function dfrhoi(rho,it,asub,bsub,esub,meatable)
implicit none

!-- Transferred variables
real(c_double), intent(in)          :: rho,asub,bsub,esub
integer(c_int), intent(in)          :: it
type(model_parameters), intent(in)  :: meatable

!-- Local variables
integer(c_int), parameter           :: Np=4
integer(c_int), parameter           :: Pr=Np/2
real(c_double)                      :: xrhoi,p
integer(c_int)                      :: k    !,ctr
real(c_double),dimension(1:Np)      :: embarray

! SRINI
integer(c_int)                      :: kMinusPr

if( rho .lt. meatable%rhotabmin  ) then
     dfrhoi = esub * ( asub * log( rho ) + bsub + asub )
else
     xrhoi = log ( rho / meatable%rhotabmin  )

     ! GIBSON
     ! p= xrhoi * ( dble( meatable%ntab ) - 1.0_cd ) / meatable%rhotabmax  + 1.0_cd

     ! SRINI
     ! meatable%rtabscale          = ( dble(meatable%ntab) - 1.0_cd )/meatable%rrtabmax
     ! meatable%ctabscale          = ( dble(meatable%ntab) - 1.0_cd )/ctabmax
     ! meatable%rhotabscale        = ( dble(meatable%ntab) - 1.0_cd )/meatable%rhotabmax
     p =  xrhoi * meatable%rhotabscale + 1.0_cd

     k = int(p)
     k = min( k , meatable%ntab - Pr )
     k = max( k , Pr )
     p = p - dble( k )

     ! make sure that p is less than 1.0
     ! then if r is out of range, p = 1.0 and rho = last value of embedtab
     !
     p = min( p , dble(Pr) )

     ! SRINI
     kMinusPr = k - Pr

!     do ctr=1,Np
!     ! GIBSON
!     ! embarray(ctr)=meatable%embedtab( k - Pr + ctr , it )
!     ! SRINI
!     embarray(ctr)=meatable%embedtab( kMinusPr + ctr , it )
!     enddo

      ! SRINI3
      embarray(:)=meatable%embedtab( (kMinusPr+1):(kMinusPr+Np), it)

     ! GIBSON
     ! dfrhoi = lagrangedrv(p,embarray)*( dble( meatable%ntab ) - 1.0_cd ) / meatable%rhotabmax 

     ! SRINI
     ! meatable%rtabscale          = ( dble(meatable%ntab) - 1.0_cd )/meatable%rrtabmax
     ! meatable%ctabscale          = ( dble(meatable%ntab) - 1.0_cd )/ctabmax
     ! meatable%rhotabscale        = ( dble(meatable%ntab) - 1.0_cd )/meatable%rhotabmax
     dfrhoi = lagrangedrv(p,embarray)*meatable%rhotabscale

     ! CONTINUE WITH GIBSON ORIG
     dfrhoi = dfrhoi / max( rho , meatable%rhotabmin )

endif
return
end function dfrhoi
!-------------------------------------------------------------------------------
!
!  Compute derivative of vector components
!
!-------------------------------------------------------------------------------
subroutine dcmij(drcos,rr,rs)

implicit none

!-- Transferred variables
real(c_double), dimension(1:3,1:3), intent(out)   :: drcos
real(c_double), dimension(1:3), intent(in)        :: rr
real(c_double), intent(in)                        :: rs

!-- Local variables
real(c_double) :: r1,r13,tmp1
integer(c_int) :: i,j

!
!     Find direction cosines
!
        r1 = 1.0_cd / rs
        r13 = r1 * r1 * r1
        do i=1,3
           tmp1 = r13 * rr(i)
           do j=1,3
              drcos(i,j) = -tmp1 * rr(j)
           enddo
           drcos(i,i) = drcos(i,i) + r1
        enddo

end subroutine dcmij
!-------------------------------------------------------------------------------
!
! Compute energy and forces on particles from the positions.
!
!-------------------------------------------------------------------------------
integer(c_int) function Compute_Energy_Forces(pkim) bind(c)
implicit none

!-- Transferred variables
type(c_ptr), intent(in) :: pkim

!-- Local variables
integer(c_int)            :: comp_force,comp_energy
integer(c_int)            :: comp_enepot,comp_virial
integer(c_int)            :: comp_ptvirl
integer(c_int)            :: ispin,kk,l,m,k,s,w,ier
integer(c_int)            :: i,j,jj,ll,mm,mc,atom
integer(c_int)            :: numnei,part_ret
integer(c_int)            :: neighOnI

integer(c_int), parameter :: ity=1
integer(c_int), parameter :: jty=1
integer(c_int), parameter :: kty=1
! GIBSON
! integer(c_int), parameter :: neimax=2000

! SRINI: neimax=2000 is too big; array initializations become slow
integer(c_int), parameter :: neimax=120

real(c_double), parameter :: h=1.0e-5_cd

! SRINI
real(c_double), parameter :: hinv=1.0_cd/h
real(c_double), parameter :: halfhinv=0.5_cd/h
integer(c_int), parameter :: ATMI=1, ATMJ=4, ATML=7  ! Array Indices for coord of atoms-I,J,L

real(c_double), parameter :: hh=h*h
real(c_double), parameter :: h2 = 2.0_cd * h
integer(c_int), parameter :: i2jl=1
integer(c_int), parameter :: j2il=2

real(c_double)           :: r_ij,r1,rej,rs,re, rsinv           ! SRINI added rsinv
real(c_double)           :: negRpij,posRpij,negRpji,posRpji
real(c_double)           :: Rsqij,Rsqil,Rsqjl,Rsqim
real(c_double)           :: Rpjl,Rsqji,Rpji,Rpij
real(c_double)           :: Rsqmj,Rpim,Rpmj,Rpil
real(c_double)           :: dxij,dyij,dzij
real(c_double)           :: dxim,dyim,dzim
real(c_double)           :: dxmj,dymj,dzmj
real(c_double)           :: rho1p,rho2p,rho3p,rho54p
real(c_double)           :: rho1m,rho2m,rho3m,rho54m
real(c_double)           :: cp,cm,cgp,cgm
real(c_double)           :: AA,asub,bsub,esub
real(c_double)           :: cx,t1,t2,t3,t54
real(c_double)           :: dang1,dang2
real(c_double)           :: rho0,phi,plocal,rhobar,roz
real(c_double)           :: smax,smin,term1,term2,term3
real(c_double)           :: angular1,angular2
real(c_double)           :: dfrho,phiif,phiij
real(c_double)           :: rcosj,rcosm,rcosl,rcoss,rcosk,rcosw
real(c_double)           :: rml,rmls,rmlsk,rmlskw
real(c_double)           :: drmllm,drmlsslm,drmlskkslm,drmlskwwkslm
real(c_double)           ::  rlm(  3,3),    rslm(  3,3,3),    rkslm(  3,3,3,3),    rwkslm(  3,3,3,3,3)
real(c_double)           :: drlmkk(3,3,3), drslmkk(3,3,3,3), drkslmkk(3,3,3,3,3), drwkslmkk(3,3,3,3,3,3)
!real(c_double)           :: tmp_lm, tmp_slm, tmp_kslm, tmp_wkslm
real(c_double)           :: drho0,drho0s
real(c_double)           :: drho1,drho1g,drho1gh
real(c_double)           :: drho1s,drho1sg,drho1sgh
real(c_double)           :: drho2,drho2s
real(c_double)           :: drho3,drho3s,drho3h,drho3sh
real(c_double)           :: drho5,drho5s
real(c_double)           :: dscrnabkp,dscrnabk0,dscrnabk1
real(c_double)           :: dscrnabk2,dscrnabk3,dscrnabk54
real(c_double)           :: tmp0p,tmp0m
real(c_double)           :: tmp1p,tmp1m,tmp1gp,tmp1gm,tmp1hp,tmp1hm
real(c_double)           :: tmp2p,tmp2m
real(c_double)           :: tmp3p,tmp3m,tmp3hp,tmp3hm
real(c_double)           :: tmp5p,tmp5m
real(c_double)           :: tmp4,tmp8

! SRINI: LOCAL ARRAYS TO SAVE COORDINATES OF ATOMS I/J/L
real(c_double)           :: RsqijHH, RsqilHH, RsqjiHH, RsqjlHH, invRho0, tmpSG01
real(c_double)           :: RsqimHH, RsqmjHH

! SRINI3: For each ij pair: first entry is Rijsq and second entry is Rij
real(c_double),dimension(2,1:neimax)             :: RijNbor

!! real(c_double),dimension(1:3)                      :: coorI, coorJ, coorL
!! real(c_double)           :: screep,scree0,scree1
!! real(c_double)           :: scree2,scree3,scree54
!  real(c_double)           :: spp,sp0,sp1,sp2,sp3,sp54
!  real(c_double)           :: smp,sm0,sm1,sm2,sm3,sm54
!  real(c_double)           :: screenijp,screenij0,screenij1
!  real(c_double)           :: screenij2,screenij3,screenij54
!  real(c_double)           :: sp,s0,s1,s2,s3,s54
!  real(c_double)           :: dpp,dp0,dp1,dp2,dp3,dp54
!  real(c_double)           :: dmp,dm0,dm1,dm2,dm3,dm54

real(c_double),dimension(1:6)                      :: scrnParArray, sParArray, screeArray
real(c_double),dimension(1:6)                      :: sPlusArray, sMinusArray, sFactorArray
real(c_double),dimension(1:9)                      :: coorIJL
real(c_double)                                     :: rcutsq, pINDX, rtabscle


! GIBSON'S Below
real(c_double),dimension(1:3)                      :: Rij,Ril,Rjl,Rim,rcos,rr,Rji
real(c_double),dimension(1:3)                      :: dscrpp,dscrp0,dscrp1
real(c_double),dimension(1:3)                      :: dscrp2,dscrp3,dscrp54
real(c_double),dimension(1:3)                      :: dscrmp,dscrm0,dscrm1
real(c_double),dimension(1:3)                      :: dscrm2,dscrm3,dscrm54
real(c_double),dimension(1:3)                      :: ap,am,agp,agm,ahp,ahm
real(c_double),dimension(1:3)                      :: dfij,dfim,dfji
real(c_double),dimension(1:3,1:3)                  :: drcos,bp,bm
real(c_double),dimension(1:3,1:3,1:3)              :: dp,dm,dhp,dhm
real(c_double),dimension(1:3,1:3,1:3,1:3,1:3)      :: gp,gm
real(c_double),dimension(1:3,1:3,1:3,1:3,1:3)      :: gp2,gm2
real(c_double),dimension(0:4)                      :: rhjp,rhjm,drhjp,drhjm
real(c_double)                                     :: rhjp4
real(c_double)                                     :: rhjm4

integer(c_int),dimension(1:neimax)                 :: jIdxOnI
real(c_double),dimension(1:neimax)                 :: scrnabp,scrnab0,scrnab1
real(c_double),dimension(1:neimax)                 :: scrnab2,scrnab3,scrnab54
real(c_double),dimension(1:2,1:3,1:neimax)         :: dscrnabp,dscrnab0,dscrnab1
real(c_double),dimension(1:2,1:3,1:neimax)         :: dscrnab2,dscrnab3,dscrnab54
real(c_double),dimension(1:3,1:neimax)             :: Rij_xyzN

real(c_double), allocatable, dimension(:,:,:)      :: slocal

!-- KIM variables
integer(c_int), pointer :: N;                   type(c_ptr) :: pN
real(c_double), pointer :: energy;              type(c_ptr) :: penergy
real(c_double), pointer :: coor(:,:);           type(c_ptr) :: pcoor
real(c_double), pointer :: force(:,:);          type(c_ptr) :: pforce
real(c_double), pointer :: enepot(:);           type(c_ptr) :: penepot
real(c_double), pointer :: Rij_list(:,:);       type(c_ptr) :: pRij_list
integer(c_int), pointer :: numContrib;          type(c_ptr) :: pnumContrib
integer(c_int), pointer :: nei1atom(:);         type(c_ptr) :: pnei1atom
integer(c_int), pointer :: particleSpecies(:);  type(c_ptr) :: pparticleSpecies
real(c_double), pointer :: particleVirial(:,:); type(c_ptr) :: pparticleVirial
real(c_double), pointer :: virial(:);           type(c_ptr) :: pvirial
integer(c_int) numberContrib
integer(c_int) idum
integer(c_int) atom_ret

type(model_parameters), pointer :: meatable; type(c_ptr) :: pmeatable

numberContrib = 0 ! initialize

! get model buffer from KIM object
pmeatable = kim_api_get_model_buffer(pkim, ier)
if (ier.lt.KIM_STATUS_OK) then
   idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                               "kim_api_get_model_buffer", ier)
   goto 42
endif
call c_f_pointer(pmeatable, meatable)

! Check to see if we have been asked to compute the forces, energyperpart,
! energy and virial
!
call kim_api_getm_compute(pkim, ier, &
     "energy",         comp_energy,     1, &
     "forces",         comp_force,      1, &
     "particleEnergy", comp_enepot,     1, &
     "particleVirial", comp_ptvirl,     1, &
     "virial",         comp_virial,     1  &
      )
if (ier.lt.KIM_STATUS_OK) then
   idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                               "kim_api_getm_compute", &
                               ier)
   goto 42
endif

! Unpack data from KIM object
!
call kim_api_getm_data(pkim, ier, &
     "numberOfParticles",           pN,              1,                                    &
     "particleSpecies",             pparticleSpecies,1,                                    &
     "coordinates",                 pcoor,           1,                                    &
     "numberContributingParticles", pnumContrib,     TRUEFALSE(meatable%HalfOrFull.eq.1),  &
     "energy",                      penergy,         TRUEFALSE(comp_energy.eq.1),          &
     "forces",                      pforce,          TRUEFALSE(comp_force.eq.1),           &
     "particleEnergy",              penepot,         TRUEFALSE(comp_enepot.eq.1),          &
     "particleVirial",              pparticleVirial, TRUEFALSE(comp_ptvirl.eq.1),          &
     "virial",                      pvirial,         TRUEFALSE(comp_virial.eq.1)           &
     )
if (ier.lt.KIM_STATUS_OK) then
   idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                               "kim_api_getm_data", ier)
   goto 42
endif

call c_f_pointer(pN, N)
call c_f_pointer(pparticleSpecies, particleSpecies, [N])
call c_f_pointer(pcoor, coor, [DIM,N])
if (meatable%HalfOrFull.eq.1)   call c_f_pointer(pnumContrib, numContrib)
if (comp_energy.eq.1)      call c_f_pointer(penergy, energy)
if (comp_force.eq.1)       call c_f_pointer(pforce,force, [DIM,N])
if (comp_enepot.eq.1)      call c_f_pointer(penepot,enepot, [N])
if (comp_ptvirl.eq.1)      call c_f_pointer(pparticleVirial, particleVirial, [6,N])
if (comp_virial.eq.1)      call c_f_pointer(pvirial, virial, [6])
if (meatable%HalfOrFull.eq.1)   then
   if (meatable%NBC.ne.0) then ! non-CLUSTER cases
      numberContrib = numContrib
   else               ! CLUSTER cases
      numberContrib = N
   endif
endif

! Check to be sure that the atom types are correct
!
ier = KIM_STATUS_FAIL ! assume an error

do i = 1,N
   if (particleSpecies(i).ne.speccode) then
      idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                  "Unexpected species type detected", ier)
      goto 42
   endif
enddo

ier = KIM_STATUS_OK ! everything is ok

!  Initialize energies, forces, and virials
!
if (comp_enepot.eq.1)     enepot(1:N)                = 0.0_cd
if (comp_energy.eq.1)     energy                     = 0.0_cd
if (comp_force.eq.1)      force(1:3,1:N)             = 0.0_cd
if (comp_ptvirl.eq.1)     particleVirial(1:6,1:N)    = 0.0_cd
if (comp_virial.eq.1)     virial(1:6)                = 0.0_cd

! Reset iterator if one is being used
!
if (meatable%IterOrLoca.eq.1) then
   ier = kim_api_get_neigh(pkim,0,0,atom_ret,numnei,pnei1atom,pRij_list)
   if (ier.lt.KIM_STATUS_OK) then
      idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                  "kim_api_get_neigh", ier)
      goto 42
   endif
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Loop over particles in the neighbor list
!  to get neighbor information, screening, and screening derivatives,
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

allocate(slocal(1:3,1:3,1:N))
slocal(1:3,1:3,1:N) = 0.0_cd

! SRINI: Save in local variable2
rcutsq   = meatable%model_cutsq
rtabscle = meatable%rtabscale

! LOOP OVER ALL THE ATOMS
do i = 1,N
   !
   ! Initialize arrays
   !
! GIBSON
!   dfij(1:3)                    = 0.0_cd
!   dfji(1:3)                    = 0.0_cd
!   dfim(1:3)                    = 0.0_cd
!   Rij_xyzN(1:3,1:neimax)       = 0.0_cd
!   jIdxOnI(1:neimax)            = 0
!   scrnabp(1:neimax)            = 0.0_cd
!   scrnab0(1:neimax)            = 0.0_cd
!   scrnab1(1:neimax)            = 0.0_cd
!   scrnab2(1:neimax)            = 0.0_cd
!   scrnab3(1:neimax)            = 0.0_cd
!   scrnab54(1:neimax)           = 0.0_cd
!   dscrnabp(1:2,1:3,1:neimax)   = 0.0_cd
!   dscrnab0(1,2:1:3,1:neimax)   = 0.0_cd
!   dscrnab1(1:2,1:3,1:neimax)   = 0.0_cd
!   dscrnab2(1:2,1:3,1:neimax)   = 0.0_cd
!   dscrnab3(1:2,1:3,1:neimax)   = 0.0_cd
!   dscrnab54(1:2,1:3,1:neimax)  = 0.0_cd

! SRINI: Apr 20, 2016.
   dfij                    = 0.0_cd
   dfji                    = 0.0_cd
   dfim                    = 0.0_cd
   Rij_xyzN                = 0.0_cd
   jIdxOnI                 = 0
   scrnabp                 = 0.0_cd
   scrnab0                 = 0.0_cd
   scrnab1                 = 0.0_cd
   scrnab2                 = 0.0_cd
   scrnab3                 = 0.0_cd
   scrnab54                = 0.0_cd
   dscrnabp                = 0.0_cd
   dscrnab0                = 0.0_cd
   dscrnab1                = 0.0_cd
   dscrnab2                = 0.0_cd
   dscrnab3                = 0.0_cd
   dscrnab54               = 0.0_cd



! ORIGINAL CODE
   neighOnI                = 0
   rho0                    = 0.0_cd
   dang1                   = 0.0_cd
   dang2                   = 0.0_cd
   phi                     = 0.0_cd
   plocal                  = 0.0_cd
   rho1p                   = 0.0_cd
   rho2p                    = 0.0_cd
   rho3p                    = 0.0_cd
   rho54p                   = 0.0_cd
   rho1m                    = 0.0_cd
   rho2m                    = 0.0_cd
   rho3m                    = 0.0_cd
   rho54m                   = 0.0_cd
   cp                       = 0.0_cd
   cm                       = 0.0_cd
   cgp                      = 0.0_cd
   cgm                      = 0.0_cd
   angular1                 = 0.0_cd
   angular2                 = 0.0_cd
   rhjp(0:4)                = 0.0_cd
   rhjm(0:4)                = 0.0_cd
   drhjp(0:4)               = 0.0_cd
   drhjm(0:4)               = 0.0_cd
   ap(1:3)                  = 0.0_cd
   am(1:3)                  = 0.0_cd
   agp(1:3)                 = 0.0_cd
   agm(1:3)                 = 0.0_cd
   ahp(1:3)                 = 0.0_cd
   ahm(1:3)                 = 0.0_cd
   bp(1:3,1:3)              = 0.0_cd
   bm(1:3,1:3)              = 0.0_cd
   dp(1:3,1:3,1:3)          = 0.0_cd
   dm(1:3,1:3,1:3)          = 0.0_cd
   dhp(1:3,1:3,1:3)         = 0.0_cd
   dhm(1:3,1:3,1:3)         = 0.0_cd
   gp(1:3,1:3,1:3,1:3,1:3)  = 0.0_cd
   gm(1:3,1:3,1:3,1:3,1:3)  = 0.0_cd
   gp2(1:3,1:3,1:3,1:3,1:3) = 0.0_cd
   gm2(1:3,1:3,1:3,1:3,1:3) = 0.0_cd

   ier = kim_api_get_neigh(pkim,1,i,part_ret,numnei,pnei1atom,pRij_list)
   if (ier.ne.KIM_STATUS_OK) then
      idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                  "kim_api_get_neigh", ier)
      goto 42
   endif

   if (meatable%NBC.ne.3) call c_f_pointer(pnei1atom, nei1atom, [numnei])
   if (meatable%NBC.eq.2) call c_f_pointer(pRij_list, Rij_list, [DIM,numnei])

   ! SRINI: SAVE Coordinates of atom-I in local array
   ! coorI, coorJ, coorL
   coorIJL(ATMI:ATMI+2) = coor(1:3,i)

   if (numnei .gt. 0) then    ! if(numei > 0)

      do jj = 1,numnei

         j = nei1atom(jj)

         ! SRINI: SAVE Coordinates of atom-J in local array
         ! coorI, coorJ, coorL
	 coorIJL(ATMJ:ATMJ+2) = coor(1:3,j)

         ! GIBSON
         ! Rij(:) = -(coor(:,j) - coor(:,i))

         ! SRINI: See sign change also
	 ! SUM(A*B)
         Rij(:) = (coorIJL(ATMI:ATMI+2) - coorIJL(ATMJ:ATMJ+2))
         ! Rsqij = SUM(Rij*Rij)
         Rsqij = Rij(1)*Rij(1) + Rij(2)*Rij(2) + Rij(3)*Rij(3)

         ! GIBSON
         ! Rsqij = dot_product(Rij,Rij)
         ! Rji(:) = -(coor(:,i) - coor(:,j))

         ! SRINI
         Rji(:) = -Rij(:)
         Rsqji  = Rsqij


	 ! GIBSON
         ! Rsqji = dot_product(Rji,Rji)
         ! if (Rsqij .ge. meatable%model_cutsq) cycle

	 ! SRINI
         if (Rsqij .ge. rcutsq) cycle

         neighOnI = neighOnI + 1

         ! GIBSON
         ! Rij_xyzN(:,neighOnI) = -Rij(:)

         ! SRINI: NO NEED the NEGATIVE SIGN, already done
         Rij_xyzN(:,neighOnI) = Rji(:)

	 ! SRINI3: Save Rij^2 for each 'J' nbor of atom-'I'
	 ! SRINI3: For each ij pair: first entry is Rijsq and second entry is Rij
	 RijNbor(1,neighOnI) = Rsqij
	 RijNbor(2,neighOnI) = sqrt(Rsqij)

         jIdxOnI(neighOnI) = j

	 ! GIBSON
 	 ! call rscrn(Rsqij,meatable,screenijp)
         ! call rscrn(Rsqij,meatable,screenij0)
         ! call rscrn(Rsqij,meatable,screenij1)
         ! call rscrn(Rsqij,meatable,screenij2)
         ! call rscrn(Rsqij,meatable,screenij3)
         ! call rscrn(Rsqij,meatable,screenij54)

	 ! SRINI2: Call once and assign computed values to all local variables
         call rscrn(Rsqij,meatable,tmpSG01)
	 scrnParArray(:) = tmpSG01

         do ll = 1,numnei

            if (ll .eq. jj) cycle

            l = nei1atom(ll)

            ! SRINI: SAVE Coordinates of atom-L in local array
            ! coorI, coorJ, coorL
            coorIJL(ATML:ATML+2) = coor(1:3,l)

	    ! GIBSON
	    ! Ril(:) = -(coor(:,l) - coor(:,i))

	    ! SRINI: See sign change also
            ! SUM(A*B)
            Ril(:) = (coorIJL(ATMI:ATMI+2) - coorIJL(ATML:ATML+2))
            ! Rsqil = SUM(Ril*Ril)
            Rsqil = Ril(1)*Ril(1) + Ril(2)*Ril(2) + Ril(3)*Ril(3)

	    ! GIBSON
	    ! Rsqil = dot_product(Ril,Ril)
            ! Rjl(:) = -(coor(:,l) - coor(:,j))
	    ! Rsqjl = dot_product(Rjl,Rjl)
	    ! GIBSON

	    ! SRINI: See sign change also
	    Rjl(:) = (coorIJL(ATMJ:ATMJ+2) - coorIJL(ATML:ATML+2))
	    ! Rsqjl = SUM(Rjl*Rjl)
            Rsqjl = Rjl(1)*Rjl(1) + Rjl(2)*Rjl(2) + Rjl(3)*Rjl(3)

            cx = get_c(Rsqij,Rsqil,Rsqjl)

	    ! SRINI2
            call get_s(cx,ity,jty,kty,meatable, sParArray)
	    scrnParArray(:) = scrnParArray(:) * sParArray(:)
            smax = maxval(scrnParArray)

            !screenijp   = screenijp * sp
            !screenij0   = screenij0 * s0
            !screenij1   = screenij1 * s1
            !screenij2   = screenij2 * s2
            !screenij3   = screenij3 * s3
            !screenij54  = screenij54* s54
	    ! smax=max(screenijp,screenij0,screenij1,screenij2,screenij3,screenij54)

            if (smax .lt. scnres) exit

         enddo ! end loop on ll

	 ! SRINI2: Old code deleted
         scrnabp(neighOnI) = scrnParArray(1)
         scrnab0(neighOnI) = scrnParArray(2)
         scrnab1(neighOnI) = scrnParArray(3)
         scrnab2(neighOnI) = scrnParArray(4)
         scrnab3(neighOnI) = scrnParArray(5)
         scrnab54(neighOnI)= scrnParArray(6)


	 ! Save locally in array
	 screeArray(:) = scrnParArray(:)

         smin=minval(screeArray)
         smax=maxval(screeArray)

	 ! SRINI: Added Paranthesis inside the if()
         if ((smin .gt. (1.0_cd-scnres)) .or. (smax .lt. scnres)) cycle

         ! Compute dscrnab for i->j
         !

         !dscrpp(1:3) = 0.0_cd
         !dscrp0(1:3) = 0.0_cd
         !dscrp1(1:3) = 0.0_cd
         !dscrp2(1:3) = 0.0_cd
         !dscrp3(1:3) = 0.0_cd
         !dscrp54(1:3) = 0.0_cd

         !dscrmp(1:3) = 0.0_cd
         !dscrm0(1:3) = 0.0_cd
         !dscrm1(1:3) = 0.0_cd
         !dscrm2(1:3) = 0.0_cd
         !dscrm3(1:3) = 0.0_cd
         !dscrm54(1:3) = 0.0_cd


   	 ! SRINI:
	 RsqijHH = Rsqij + hh

         do m = 1,3
	    ! GIBSON
            ! Rpij = Rsqij - h2 * Rij(m) + hh

	    ! SRINI
            Rpij = RsqijHH - (h2 * Rij(m))

            call rscrn(Rpij,meatable,posRpij)
            dscrpp(m) = posRpij
            dscrp0(m) = posRpij
            dscrp1(m) = posRpij
            dscrp2(m) = posRpij
            dscrp3(m) = posRpij
            dscrp54(m) = posRpij

	    ! GIBSON
            ! Rpij = Rsqij + h2 * Rij(m) + hh

	    ! SRINI
            Rpij = RsqijHH + (h2 * Rij(m))

            call rscrn(Rpij,meatable,negRpij)
            dscrmp(m)  = negRpij
            dscrm0(m)  = negRpij
            dscrm1(m)  = negRpij
            dscrm2(m)  = negRpij
            dscrm3(m)  = negRpij
            dscrm54(m) = negRpij

         enddo

         do ll = 1,numnei
  
            if (ll .eq. jj) cycle

            l = nei1atom(ll)

	    ! GIBSON
            ! Ril(:) = -(coor(:,l) - coor(:,i))
	    ! Rsqil = dot_product(Ril,Ril)


            ! SRINI: SAVE Coordinates of atom-L in local array
            ! SRINI: See sign change also
            coorIJL(ATML:ATML+2) = coor(1:3,l)


	    ! SRINI: See sign change also
            Ril(:) = (coorIJL(ATMI:ATMI+2) - coorIJL(ATML:ATML+2))
            ! Rsqil = SUM(Ril*Ril)
            Rsqil = Ril(1)*Ril(1) + Ril(2)*Ril(2) + Ril(3)*Ril(3)


	    ! GIBSON
            ! Rjl(:) = -(coor(:,l) - coor(:,j))
	    ! Rsqjl  = dot_product(Rjl,Rjl)

	    ! SRINI: See sign change also
            Rjl(:) = (coorIJL(ATMJ:ATMJ+2) - coorIJL(ATML:ATML+2))
	    ! Rsqjl  = SUM(Rjl*Rjl)
            Rsqjl = Rjl(1)*Rjl(1) + Rjl(2)*Rjl(2) + Rjl(3)*Rjl(3)


	    ! SRINI:
	    RsqijHH = Rsqij + hh
	    RsqilHH = Rsqil + hh

            do m = 1,3
 
	       ! GIBSON
               ! Rpij = Rsqij - h2 * Rij(m) + hh
               ! Rpil = Rsqil - h2 * Ril(m) + hh

	       ! SRINIZ: Apr 21, 2016
               Rpij = RsqijHH - (h2 * Rij(m))
               Rpil = RsqilHH - (h2 * Ril(m))

               cx = get_c(Rpij,Rpil,Rsqjl)

               ! SRINI2
               call get_s(cx,ity,jty,kty,meatable, sParArray)
               !call get_s(cx,ity,jty,kty,meatable,dpp,dp0,dp1,dp2,dp3,dp54)
               dscrpp(m)  = dscrpp(m)*sParArray(1)
               dscrp0(m)  = dscrp0(m)*sParArray(2)
               dscrp1(m)  = dscrp1(m)*sParArray(3)
               dscrp2(m)  = dscrp2(m)*sParArray(4)
               dscrp3(m)  = dscrp3(m)*sParArray(5)
               dscrp54(m) = dscrp54(m)*sParArray(6)

	       ! GIBSON
               ! Rpij = Rsqij + h2 * Rij(m) + hh
               ! Rpil = Rsqil + h2 * Ril(m) + hh

	       ! SRINIZ: Apr 21, 2016
               Rpij = RsqijHH + (h2 * Rij(m))
               Rpil = RsqilHH + (h2 * Ril(m))

               cx = get_c(Rpij,Rpil,Rsqjl)

	       ! SRINI2
               call get_s(cx,ity,jty,kty,meatable, sParArray)
               !call get_s(cx,ity,jty,kty,meatable,dmp,dm0,dm1,dm2,dm3,dm54)
               dscrmp(m)  = dscrmp(m)*sParArray(1)
               dscrm0(m)  = dscrm0(m)*sParArray(2)
               dscrm1(m)  = dscrm1(m)*sParArray(3)
               dscrm2(m)  = dscrm2(m)*sParArray(4)
               dscrm3(m)  = dscrm3(m)*sParArray(5)
               dscrm54(m) = dscrm54(m)*sParArray(6)

            enddo

         enddo ! end loop on ll

         do m = 1,3
	 
	 ! GIBSON
	 !   dscrnabp(i2jl,m,neighOnI) = 0.5_cd*(dscrpp(m)-dscrmp(m))/h
         !   dscrnab0(i2jl,m,neighOnI) = 0.5_cd*(dscrp0(m)-dscrm0(m))/h
         !   dscrnab1(i2jl,m,neighOnI) = 0.5_cd*(dscrp1(m)-dscrm1(m))/h
         !   dscrnab2(i2jl,m,neighOnI) = 0.5_cd*(dscrp2(m)-dscrm2(m))/h
         !   dscrnab3(i2jl,m,neighOnI) = 0.5_cd*(dscrp3(m)-dscrm3(m))/h
         !   dscrnab54(i2jl,m,neighOnI)= 0.5_cd*(dscrp54(m)-dscrm54(m))/h

	 ! SRINI
            dscrnabp(i2jl,m,neighOnI) = (dscrpp(m)-dscrmp(m))*halfhinv
            dscrnab0(i2jl,m,neighOnI) = (dscrp0(m)-dscrm0(m))*halfhinv
            dscrnab1(i2jl,m,neighOnI) = (dscrp1(m)-dscrm1(m))*halfhinv
            dscrnab2(i2jl,m,neighOnI) = (dscrp2(m)-dscrm2(m))*halfhinv
            dscrnab3(i2jl,m,neighOnI) = (dscrp3(m)-dscrm3(m))*halfhinv
            dscrnab54(i2jl,m,neighOnI)= (dscrp54(m)-dscrm54(m))*halfhinv

         enddo

         ! Compute dscrnab for j->i
         !

         !dscrpp(1:3) = 0.0_cd
         !dscrp0(1:3) = 0.0_cd
         !dscrp1(1:3) = 0.0_cd
         !dscrp2(1:3) = 0.0_cd
         !dscrp3(1:3) = 0.0_cd
         !dscrp54(1:3) = 0.0_cd

         !dscrmp(1:3) = 0.0_cd
         !dscrm0(1:3) = 0.0_cd
         !dscrm1(1:3) = 0.0_cd
         !dscrm2(1:3) = 0.0_cd
         !dscrm3(1:3) = 0.0_cd
         !dscrm54(1:3) = 0.0_cd

	 ! SRINI:
	 RsqjiHH = Rsqji + hh


         do m = 1,3

	    ! GIBSON
            ! Rpji = Rsqji - h2 * Rji(m) + hh


	    ! SRINIZ: Apr 21, 2016
            Rpji = RsqjiHH - (h2 * Rji(m))

            call rscrn(Rpji,meatable,posRpji)
            dscrpp(m) = posRpji
            dscrp0(m) = posRpji
            dscrp1(m) = posRpji
            dscrp2(m) = posRpji
            dscrp3(m) = posRpji
            dscrp54(m) = posRpji

	    ! GIBSON
            ! Rpji = Rsqji + h2 * Rji(m) + hh

	    ! SRINIZ: Apr 21, 2016
            Rpji = RsqjiHH + (h2 * Rji(m))

            call rscrn(Rpji,meatable,negRpji)
            dscrmp(m) = negRpji
            dscrm0(m) = negRpji
            dscrm1(m) = negRpji
            dscrm2(m) = negRpji
            dscrm3(m) = negRpji
            dscrm54(m) = negRpji

         enddo

         do ll = 1,numnei

            if (ll .eq. jj) cycle

            l = nei1atom(ll)

	    ! GIBSON
            ! Ril(:) = -(coor(:,l) - coor(:,i))
	    ! Rsqil = dot_product(Ril,Ril)

            ! SRINI: SAVE Coordinates of atom-L in local array
            ! coorI, coorJ, coorL
            coorIJL(ATML:ATML+2) = coor(1:3,l)

	    ! SRINI: See sign change also
            Ril(:) = (coorIJL(ATMI:ATMI+2) - coorIJL(ATML:ATML+2))
            ! Rsqil = SUM(Ril*Ril)
            Rsqil = Ril(1)*Ril(1) + Ril(2)*Ril(2) + Ril(3)*Ril(3)


	    ! GIBSON
            ! Rjl(:) = -(coor(:,l) - coor(:,j))
	    ! Rsqjl = dot_product(Rjl,Rjl)

	    ! SRINI: See sign change also
	    Rjl(:) = (coorIJL(ATMJ:ATMJ+2) - coorIJL(ATML:ATML+2))            
            ! Rsqjl = SUM(Rjl*Rjl)
            Rsqjl = Rjl(1)*Rjl(1) + Rjl(2)*Rjl(2) + Rjl(3)*Rjl(3)


            ! SRINIZ: Apr 21, 2016
            RsqjiHH = Rsqji + hh
            RsqjlHH = Rsqjl + hh

            do m = 1,3

	       ! GIBSON
               ! Rpji = Rsqji - h2 * Rji(m) + hh
               ! Rpjl = Rsqjl - h2 * Rjl(m) + hh

	       ! SRINI
               Rpji = RsqjiHH - (h2 * Rji(m))
               Rpjl = RsqjlHH - (h2 * Rjl(m))

               cx = get_c(Rpji,Rsqil,Rpjl)

	       ! SRINI2
               call get_s(cx,ity,jty,kty,meatable, sParArray)
               !call get_s(cx,ity,jty,kty,meatable,dpp,dp0,dp1,dp2,dp3,dp54)

               dscrpp(m) = dscrpp(m)*sParArray(1)
               dscrp0(m) = dscrp0(m)*sParArray(2)
               dscrp1(m) = dscrp1(m)*sParArray(3)
               dscrp2(m) = dscrp2(m)*sParArray(4)
               dscrp3(m) = dscrp3(m)*sParArray(5)
               dscrp54(m) = dscrp54(m)*sParArray(6)

	       ! GIBSON
	       ! Rpji = Rsqji + h2 * Rji(m) + hh
               ! Rpjl = Rsqjl + h2 * Rjl(m) + hh

	       ! SRINI
               Rpji = RsqjiHH + (h2 * Rji(m))
               Rpjl = RsqjlHH + (h2 * Rjl(m))

               cx = get_c(Rpji,Rsqil,Rpjl)

               ! SRINI2
               call get_s(cx,ity,jty,kty,meatable, sParArray)
               ! call get_s(cx,ity,jty,kty,meatable,dmp,dm0,dm1,dm2,dm3,dm54)

               dscrmp(m) = dscrmp(m)*sParArray(1)
               dscrm0(m) = dscrm0(m)*sParArray(2)
               dscrm1(m) = dscrm1(m)*sParArray(3)
               dscrm2(m) = dscrm2(m)*sParArray(4)
               dscrm3(m) = dscrm3(m)*sParArray(5)
               dscrm54(m) = dscrm54(m)*sParArray(6)

            enddo

         enddo ! end loop on ll

         do m = 1,3

	   ! GIBSON
           ! dscrnabp(j2il,m,neighOnI) = 0.5_cd*(dscrpp(m)-dscrmp(m))/h
           ! dscrnab0(j2il,m,neighOnI) = 0.5_cd*(dscrp0(m)-dscrm0(m))/h
           ! dscrnab1(j2il,m,neighOnI) = 0.5_cd*(dscrp1(m)-dscrm1(m))/h
           ! dscrnab2(j2il,m,neighOnI) = 0.5_cd*(dscrp2(m)-dscrm2(m))/h
           ! dscrnab3(j2il,m,neighOnI) = 0.5_cd*(dscrp3(m)-dscrm3(m))/h
           ! dscrnab54(j2il,m,neighOnI)= 0.5_cd*(dscrp54(m)-dscrm54(m))/h

	    ! SRINI
            dscrnabp(j2il,m,neighOnI) = (dscrpp(m)-dscrmp(m))*halfhinv
            dscrnab0(j2il,m,neighOnI) = (dscrp0(m)-dscrm0(m))*halfhinv
            dscrnab1(j2il,m,neighOnI) = (dscrp1(m)-dscrm1(m))*halfhinv
            dscrnab2(j2il,m,neighOnI) = (dscrp2(m)-dscrm2(m))*halfhinv
            dscrnab3(j2il,m,neighOnI) = (dscrp3(m)-dscrm3(m))*halfhinv
            dscrnab54(j2il,m,neighOnI)= (dscrp54(m)-dscrm54(m))*halfhinv

         enddo

      enddo  ! end loop on jj

   endif     ! if(numnei > 0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Now compute the pair potential and embedding function for each atom i
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   if (neighOnI .gt. 0) then     ! if(neighOnI > 0)

      do jj = 1,neighOnI

         Rij(:) = -Rij_xyzN(:,jj)

	 ! GIBSON
         ! Rsqij = dot_product(Rij,Rij)


	 ! SRINIZ
	 ! Rsqij = SUM(Rij*Rij)
	 !Rsqij = Rij(1)*Rij(1) + Rij(2)*Rij(2) + Rij(3)*Rij(3)
         !r_ij = sqrt(Rsqij)

         ! SRINI3: Save Rij^2 for each 'J' nbor of atom-'I'
         ! SRINI3: For each ij pair: first entry is Rijsq and second entry is Rij
         Rsqij = RijNbor(1,jj)
	 r_ij  = RijNbor(2,jj)

	 ! SRINI2
         screeArray(1) = scrnabp(jj)
         screeArray(2) = scrnab0(jj)
         screeArray(3) = scrnab1(jj)
         screeArray(4) = scrnab2(jj)
         screeArray(5) = scrnab3(jj)
         screeArray(6) = scrnab54(jj)
         smin=minval(screeArray)
         smax=maxval(screeArray)

         if (smax .le. scnres) cycle

         !
         ! compute phi for atom i
         !

         esub=meatable%esubs(ity,jty)

	 ! GIBSON
	 ! phi=phi+0.5_cd*calc_phi(r_ij,ity,jty,meatable)*screep

	 ! SRINI: ADDED a new argument pINDX
	 pINDX = r_IJ*rtabscle + 1.0_cd

	 ! SRINI2
         phi=phi+0.5_cd*calc_phi(r_ij,ity,jty,meatable, pINDX)*screeArray(1)

         !
         ! compute background electron density for atom i
         !

         rej=meatable%res(ity,jty)

         do ll=0,4
            ispin=1
            ! rhjp(ll)=rhof(r_ij,ll,ispin,jty,meatable)
	    ! SRINI: ADDED a new argument pINDX
            rhjp(ll)=rhof(r_ij,ll,ispin,jty,meatable, pINDX)

            ispin=2
            ! rhjm(ll)=rhof(r_ij,ll,ispin,jty,meatable)
	    ! SRINI: ADDED a new argument pINDX
            rhjm(ll)=rhof(r_ij,ll,ispin,jty,meatable, pINDX)
         enddo

         r1 = 1.0_cd / r_ij
         rcos(:)= r1 * Rij(:)

         do m=1,3
            t1    =rcos(m)
            ap(m) =ap(m)+rhjp(1)*t1*screeArray(3)
            am(m) =am(m)+rhjm(1)*t1*screeArray(3)
            agp(m)=agp(m)+rhjp(3)*t1*screeArray(5)
            agm(m)=agm(m)+rhjm(3)*t1*screeArray(5)
            ahp(m)=ahp(m)+rhjp(4)*t1*screeArray(6)
            ahm(m)=ahm(m)+rhjm(4)*t1*screeArray(6)

            do l=1,3
               t2=rcos(m)*rcos(l)
               bp(l,m)=bp(l,m)+rhjp(2)*t2*screeArray(4)
               bm(l,m)=bm(l,m)+rhjm(2)*t2*screeArray(4)
               do s=1,3
                  t3=rcos(m)*rcos(l)*rcos(s)
                  dp(s,l,m)= dp(s,l,m)+rhjp(3)*t3*screeArray(5)
                  dm(s,l,m)= dm(s,l,m)+rhjm(3)*t3*screeArray(5)
                  dhp(s,l,m)=dhp(s,l,m)+rhjp(4)*t3*screeArray(6)
                  dhm(s,l,m)=dhm(s,l,m)+rhjm(4)*t3*screeArray(6)
                  do k=1,3
                     do w=1,3
                        t54=rcos(m)*rcos(l)*rcos(s)*rcos(k)*rcos(w)
                        gp(w,k,s,l,m)=gp(w,k,s,l,m)+rhjp(4)*t54*screeArray(6)
                        gm(w,k,s,l,m)=gm(w,k,s,l,m)+rhjm(4)*t54*screeArray(6)
                     enddo
                  enddo
               enddo
            enddo
         enddo

         gp2=gp*2.0_cd
         gm2=gm*2.0_cd


         cp=cp+rhjp(0)*screeArray(2)
         cm=cm+rhjm(0)*screeArray(2)

         cgp=cgp+rhjp(2)*screeArray(4)
         cgm=cgm+rhjm(2)*screeArray(4)

      enddo  ! loop on jj

   endif     ! if(neighOnI > 0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Now that we know the + and - partial electron densities, calculate embedding
!  energy and the deriviative of embedding energy, the derivative of the
!  pair potential, to obtain the particle energy, total energy, partical forces,
!  and viral forces.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   if (neighOnI .gt. 0) then     ! if(neighOnI > 0)

      ! rho0=cp**2-cm**2
      ! SRINI: DON'T USE POWER INTRINSIC FUNCTION
      rho0 = (cp*cp) - (cm*cm)

      if(rho0.gt.0.0_cd) then     !if(rho0>0.0)

         rho0=sqrt(rho0)
         asub=meatable%asubs(ity)
         bsub=meatable%bsubs(ity)
         esub=meatable%esubs(ity,jty)

         rho1p=0.0_cd
         rho1m=0.0_cd

      ! SRINI: DON'T USE POWER INTRINSIC FUNCTION
      !   rho2p=-(cgp**2)/3.0_cd
      !   rho2m=-(cgm**2)/3.0_cd
         rho2p=-(cgp*cgp)/3.0_cd
         rho2m=-(cgm*cgm)/3.0_cd

         rho3p=0.0_cd
         rho3m=0.0_cd
         rho54p=0.0_cd
         rho54m=0.0_cd
         do m=1,3
            rho3p =rho3p -legend00*agp(m)**2
            rho3m =rho3m -legend00*agm(m)**2
            rho54p=rho54p+legend02*ahp(m)**2
            rho54m=rho54m+legend02*ahm(m)**2
            do l=1,3
               do s=1,3
                  rho54p=rho54p-legend01*dhp(s,l,m)**2
                  rho54m=rho54m-legend01*dhm(s,l,m)**2
               enddo
            enddo
         enddo

         do m=1,3
            rho1p=rho1p+ap(m)**2
            rho1m=rho1m+am(m)**2
            do l=1,3
               rho2p=rho2p+bp(l,m)**2
               rho2m=rho2m+bm(l,m)**2
               do s=1,3
                  rho3p=rho3p+dp(s,l,m)**2
                  rho3m=rho3m+dm(s,l,m)**2
                  do k=1,3
                     do w=1,3
                        rho54p=rho54p+gp(w,k,s,l,m)**2
                        rho54m=rho54m+gm(w,k,s,l,m)**2
                     enddo
                  enddo
               enddo
            enddo
          enddo

!        Eqn. (9a)
!
          AA=(rho1p+rho2p+rho3p+rho54p)-(rho1m+rho2m+rho3m+rho54m)

          call elecden(rho0,AA,angular1,angular2,rhobar)
          dang1=angular1
          dang2=angular2
          plocal=frhoi(rhobar,ity,asub,bsub,esub,meatable)

          if (rhobar.lt.0.0_cd) then
              print *,'**********warning**********'
              print *, 'for atom ',i,' rhobar= ',rhobar,plocal
              goto 42
          endif

      endif     ! if(rho0 > 0.0)

      if(comp_enepot.eq.1) then
         enepot(i) = enepot(i) + plocal + phi
!!         print *, 'rho1p,rho2p,rho3p,rho5p',rho1p,rho2p,rho3p,rho54p
!!         print *, 'rho1m,rho2m,rho3m,rho5m',rho1m,rho2m,rho3m,rho54m
!!         print *, 'rho0,rho1,rho2,rho3,rho5',rho0,rho1p-rho1m,rho2p-rho2m,rho3p-rho3m,rho54p-rho54m
!!         print *, 'i,enepot,plocal,phi,rhobar ', i,enepot(i),plocal,phi,rhobar
      endif

      if(comp_energy.eq.1) then
         energy = energy + plocal + phi
!!         print *,'atom = ',i,'; energy = ',energy
      endif

   endif     ! if(neighOnI > 0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Also, get the particle forces and the viral forces (if requested)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   if (neighOnI .gt. 0) then

      rho0=cp**2-cm**2
      if(rho0.lt.0.0_cd) goto 42
      rho0=sqrt(rho0)
      asub=meatable%asubs(ity)
      bsub=meatable%bsubs(ity)
      esub=meatable%esubs(ity,jty)
      roz=rhobar
      dfrho=dfrhoi(roz,ity,asub,bsub,esub,meatable)

      do jj = 1,neighOnI

         j = jIdxOnI(jj)

         Rij(:) = -Rij_xyzN(:,jj)

	 ! GIBSON
         ! Rsqij = dot_product(Rij,Rij)
	 

	 ! SRINI
         ! Rsqij = SUM(Rij*Rij)
         ! Rsqij = Rij(1)*Rij(1) + Rij(2)*Rij(2) + Rij(3)*Rij(3)
         ! r_ij = sqrt(Rsqij)

         ! SRINI3: Save Rij^2 for each 'J' nbor of atom-'I'
         ! SRINI3: For each ij pair: first entry is Rijsq and second entry is Rij
         Rsqij = RijNbor(1,jj)
         r_ij  = RijNbor(2,jj)


         ! SRINI: ADDED a new argument pINDX
         pINDX = r_IJ*rtabscle + 1.0_cd
         phiij = calc_phi(r_ij,ity,jty,meatable, pINDX)

	 ! SRINI2
         screeArray(1) = scrnabp(jj)
         screeArray(2) = scrnab0(jj)
         screeArray(3) = scrnab1(jj)
         screeArray(4) = scrnab2(jj)
         screeArray(5) = scrnab3(jj)
         screeArray(6) = scrnab54(jj)
         smin=minval(screeArray(:))
         smax=maxval(screeArray(:))

         if (smax .le. scnres) cycle

         !
         !     Find interatomic distances, normalized directional vectors, etc.
         !

         rs=r_ij
         dxij=Rij(1)
         dyij=Rij(2)
         dzij=Rij(3)
         rr(:)=Rij(:)

	 ! GIBSON
         ! rcos(:)=rr(:)/rs

	 ! SRINI: Use inverse and multiply (not divide)
	 rsinv = 1.0_cd/rs
	 rcos(:) = rr(:)*rsinv

         call dcmij(drcos,rr,rs)

	 ! SRINI: ADDED a new argument pINDX
	 pINDX = r_ij*rtabscle + 1.0_cd

!!!!! Kent  rcos terms

   !!   do m = 1,3
   !!   do l = 1,3
   !!                  tmp_lm = rcos(l)*     rcos(m);     rlm(     l,m)  =    tmp_lm
   !!   do s = 1,3
   !!                 tmp_slm = rcos(s)*     tmp_lm;     rslm(    s,l,m) =   tmp_slm
   !!   do k = 1,3
   !!                tmp_kslm = rcos(k)*    tmp_slm;    rkslm(  k,s,l,m) =  tmp_kslm
   !!   do w = 1,3
   !!               tmp_wkslm = rcos(w)*   tmp_kslm;   rwkslm(w,k,s,l,m) = tmp_wkslm
   !!   end do; end do; end do; end do; end do;
!or
        do m = 1,3
        do l = 1,3
                     rlm(     l,m)     = rcos(l)*    rcos(    m)
        do s = 1,3
                     rslm(    s,l,m)   = rcos(s)*   rlm(    l,m)
        do k = 1,3
                     rkslm(  k,s,l,m)  = rcos(k)*  rslm(  s,l,m)
        do w = 1,3
                     rwkslm(w,k,s,l,m) = rcos(w)* rkslm(k,s,l,m)
        end do; end do; end do; end do; end do;

!!!!!

!!!!! Kent  Derivative terms

        do kk=1,3  
        do  m=1,3
        do  l=1,3
                   drlmkk(         l,m,kk)=   drcos(      m,kk) * rcos(l) + rcos(       m) * drcos(l,kk)
        do  s=1,3
                   drslmkk(      s,l,m,kk)=  drlmkk(    l,m,kk) * rcos(s) + rlm(      l,m) * drcos(s,kk)
        do  k=1,3
                   drkslmkk(   k,s,l,m,kk)= drslmkk(  s,l,m,kk) * rcos(k) + rslm(   s,l,m) * drcos(k,kk)
        do  w=1,3
                   drwkslmkk(w,k,s,l,m,kk)=drkslmkk(k,s,l,m,kk) * rcos(w) + rkslm(k,s,l,m) * drcos(w,kk)
        enddo; enddo; enddo; enddo; enddo; enddo; 
!!!!!

         phiif=dphif(r_ij,ity,jty,meatable, pINDX)

         !
         !     Compile force on i and jindx atoms
         !
         !     note the minus sign on j term as the direc. cosine is opposite
         !

         do m=1,3

            dfij(m)=phiif*rcos(m)*scrnabp(jj)-phiij*dscrnabp(i2jl,m,jj)

            if(comp_force.eq.1) then
               force(m,i) = force(m,i) - dfij(m)
            endif

         enddo

         slocal(1,1,i) = slocal(1,1,i) - dfij(1)*dxij
         slocal(2,2,i) = slocal(2,2,i) - dfij(2)*dyij
         slocal(3,3,i) = slocal(3,3,i) - dfij(3)*dzij
         slocal(2,3,i) = slocal(2,3,i) - dfij(2)*dzij
         slocal(1,3,i) = slocal(1,3,i) - dfij(1)*dzij
         slocal(1,2,i) = slocal(1,2,i) - dfij(1)*dyij

         !
         !     Parameters for calculating embedding potentials
         !

         re=meatable%res(ity,jty)

	 ! WARNING: GIBSON used "rs" (but defined it as r_ij). 
	 ! SRINI USES pINDX calculated from r_ij
         do ll=0,4
            ispin=1
            drhjp(ll)=drhof(rs,ll,ispin,jty,meatable, pINDX)
            rhjp(ll)=rhof(rs,ll,ispin,jty,meatable, pINDX)
            ispin=2
            drhjm(ll)=drhof(rs,ll,ispin,jty,meatable, pINDX)
            rhjm(ll)=rhof(rs,ll,ispin,jty,meatable, pINDX)
         enddo

rhjm4=rhjm(4)
rhjp4=rhjp(4)

         tmp0p = (2.0_cd / 3.0_cd) * cgp
         tmp0m = (2.0_cd / 3.0_cd) * cgm

	 ! SRINI
	 invRho0 = 1.0_cd/rho0

         do kk=1,3     !kk
            rcosj = rcos(kk)

            !
            !     Initialize accumulation variables
            !

	    ! GIBSON
            ! drho0 = (cp*drhjp(0) - cm*drhjm(0))/rho0 * screeArray(2) * rcosj
            ! drho0s = (cp*rhjp(0) - cm*rhjm(0))/rho0

	    ! SRINI
	    drho0 = (cp*drhjp(0) - cm*drhjm(0))*invRho0 * screeArray(2) * rcosj
            drho0s = (cp*rhjp(0) - cm*rhjm(0))*invRho0

            drho1 = 0.0_cd
            drho1g = 0.0_cd
            drho1gh = 0.0_cd
            drho1s = 0.0_cd
            drho1sg = 0.0_cd
            drho1sgh = 0.0_cd
            drho2 = -(tmp0p * drhjp(2) - tmp0m * drhjm(2)) * screeArray(4) * rcosj
            drho2s = -rhjp(2) * tmp0p + rhjm(2) * tmp0m
            drho3 = 0.0_cd
            drho3h = 0.0_cd
            drho3s = 0.0_cd
            drho3sh = 0.0_cd
            drho5 = 0.0_cd
            drho5s = 0.0_cd

            do m = 1 , 3     ! do m
               rcosm = rcos(m)
               tmp1p = 2.0_cd * ap(m)
               tmp1m = 2.0_cd * am(m)
               tmp1gp = 2.0_cd * agp(m)
               tmp1gm = 2.0_cd * agm(m)
               tmp1hp = 2.0_cd * ahp(m)
               tmp1hm = 2.0_cd * ahm(m)

               drho1 = drho1 + tmp1p * (drhjp(1) * rcosj * rcosm + rhjp(1) * drcos(m,kk)) * screeArray(3) &
                             - tmp1m * (drhjm(1) * rcosj * rcosm + rhjm(1) * drcos(m,kk)) * screeArray(3)

               drho1g = drho1g + tmp1gp * (drhjp(3) * rcosj * rcosm + rhjp(3) * drcos(m,kk)) * screeArray(5) &
                               - tmp1gm * (drhjm(3) * rcosj * rcosm + rhjm(3) * drcos(m,kk)) * screeArray(5)

               drho1gh = drho1gh + tmp1hp * (drhjp(4) * rcosj * rcosm + rhjp(4) * drcos(m,kk)) * screeArray(6) &
                                 - tmp1hm * (drhjm(4) * rcosj * rcosm + rhjm(4) * drcos(m,kk)) * screeArray(6)

               drho1s = drho1s + ((tmp1p*rhjp(1))-(tmp1m*rhjm(1)))*rcosm
               drho1sg = drho1sg + ((tmp1gp*rhjp(3))-(tmp1gm*rhjm(3)))*rcosm
               drho1sgh = drho1sgh + ((tmp1hp*rhjp(4))-(tmp1hm*rhjm(4)))*rcosm

               do l = 1, 3     ! do l
                  rcosl = rcos(l)
                  rml = rcosm * rcosl
                  drmllm = drcos(m,kk) * rcosl + rcosm * drcos(l,kk)

                  tmp2p = 2.0_cd * bp(l,m)
                  tmp2m = 2.0_cd * bm(l,m)

                  drho2 = drho2 + tmp2p * (drhjp(2) * rcosj * rml + rhjp(2) * drmllm) * screeArray(4) &
                                - tmp2m * (drhjm(2) * rcosj * rml + rhjm(2) * drmllm) * screeArray(4)

                  drho2s = drho2s + ((tmp2p*rhjp(2))-(tmp2m*rhjm(2)))*rml

                  do s = 1, 3     ! do s
                     rcoss = rcos(s)
                     rmls = rml * rcoss
                     tmp3p = 2.0_cd * dp(s,l,m)
                     tmp3m = 2.0_cd * dm(s,l,m)
                     tmp3hp = 2.0_cd * dhp(s,l,m)
                     tmp3hm = 2.0_cd * dhm(s,l,m)
                     drmlsslm = rml * drcos(s,kk) + rcoss * drmllm

                     drho3 = drho3 + tmp3p * (drhjp(3) * rcosj * rmls + rhjp(3) * drmlsslm) * screeArray(5) &
                                   - tmp3m * (drhjm(3) * rcosj * rmls + rhjm(3) * drmlsslm) * screeArray(5)

                     drho3h = drho3h + tmp3hp * (drhjp(4) * rcosj * rmls + rhjp(4) * drmlsslm) * screeArray(6) &
                                     - tmp3hm * (drhjm(4) * rcosj * rmls + rhjm(4) * drmlsslm) * screeArray(6)

                     drho3s = drho3s + tmp3p * rhjp(3) * rmls &
                                     - tmp3m * rhjm(3) * rmls

                     drho3sh = drho3sh + tmp3hp * rhjp(4) * rmls &
                                       - tmp3hm * rhjm(4) * rmls


                     !$omp SIMD collapse(2) reduction(+:drho5) reduction(+:drho5s)
                     do k = 1, 3    ! do k
                     do w = 1, 3    ! do w

                           tmp5p = gp2(w,k,s,l,m)
                           tmp5m = gm2(w,k,s,l,m)

                          !drho5 = drho5 + tmp5p * (drhjp(4) * rcosj * rmlskw            + rhjp(4) * drmlskwwkslm())          * screeArray(6) &
                          !              - tmp5m * (drhjm(4) * rcosj * rmlskw            + rhjm(4) * drmlskwwkslm())          * screeArray(6)

                           drho5 = drho5 + tmp5p * (drhjp(4) * rcosj * rwkslm(w,k,s,l,m) + rhjp(4) * drwkslmkk(w,k,s,l,m,kk)) * screeArray(6) &
                                         - tmp5m * (drhjm(4) * rcosj * rwkslm(w,k,s,l,m) + rhjm(4) * drwkslmkk(w,k,s,l,m,kk)) * screeArray(6)

                          !drho5s = drho5s + ((tmp5p * rhjp(4)) - (tmp5m * rhjm(4))) * rmlskw  !!3.5

                           drho5s = drho5s + ((tmp5p * rhjp(4)) - (tmp5m * rhjm(4))) * rwkslm(w,k,s,l,m)

                     enddo     !do w
                     enddo     !do k


                  enddo     ! do s
               enddo     ! do l
            enddo     ! do m
                    ! do kk
            !
            !     Compile forces for i -> j and j -> i
            !

            term1 =  drho1 &
                  +  drho2 &
                  + (drho3 - legend00*drho1g) &
                  + (drho5 - legend01*drho3h + legend02*drho1gh)

            term2 =  drho1s * dscrnab1(i2jl,kk,jj) &
                  +  drho2s * dscrnab2(i2jl,kk,jj) &
                  + (drho3s - legend00*drho1sg) * dscrnab3(i2jl,kk,jj)  &
                  + (drho5s - legend01*drho3sh  + legend02*drho1sgh) * dscrnab54(i2jl,kk,jj)

            tmp4 = dscrnab0(i2jl,kk,jj)

            term3 =  drho1s * dscrnab1(j2il,kk,jj) &
                  +  drho2s * dscrnab2(j2il,kk,jj) &
                  + (drho3s - legend00*drho1sg) * dscrnab3(j2il,kk,jj)  &
                  + (drho5s - legend01*drho3sh  + legend02*drho1sgh) * dscrnab54(j2il,kk,jj)

            tmp8 = dscrnab0(j2il,kk,jj)

            dfij(kk) = dfrho * ( (drho0 - drho0s * tmp4) * dang1 + (term1 - term2) * dang2 )
            dfji(kk) = dfrho * ( (drho0 + drho0s * tmp8) * dang1 + (term1 + term3) * dang2 )

            if(comp_force.eq.1) then
               force(kk,i) = force(kk,i) - dfij(kk)
               force(kk,j) = force(kk,j) + dfji(kk)
            endif

         enddo     ! do kk

         slocal(1,1,i) = slocal(1,1,i) - dfij(1)*dxij
         slocal(2,2,i) = slocal(2,2,i) - dfij(2)*dyij
         slocal(3,3,i) = slocal(3,3,i) - dfij(3)*dzij
         slocal(2,3,i) = slocal(2,3,i) - dfij(2)*dzij
         slocal(1,3,i) = slocal(1,3,i) - dfij(1)*dzij
         slocal(1,2,i) = slocal(1,2,i) - dfij(1)*dyij

         slocal(1,1,j) = slocal(1,1,j) - dfji(1)*dxij
         slocal(2,2,j) = slocal(2,2,j) - dfji(2)*dyij
         slocal(3,3,j) = slocal(3,3,j) - dfji(3)*dzij
         slocal(2,3,j) = slocal(2,3,j) - dfji(2)*dzij
         slocal(1,3,j) = slocal(1,3,j) - dfji(1)*dzij
         slocal(1,2,j) = slocal(1,2,j) - dfji(1)*dyij

! SRINI DEBUG
!print *,'**********SRINI_DEBUG**********'
!print *, 'neighOnI/neimax= ', neighOnI, neimax
!print *,'**********SRINI_DEBUG**********'

         do mm = 1,neighOnI    ! do mm (neighOnI)

            if (mm .eq. jj) cycle

            mc = jIdxOnI(mm)

            Rim(:) = -Rij_xyzN(:,mm)

	    ! GIBSON
            ! Rsqim = dot_product(Rim,Rim)

	    ! SRINI
            ! Rsqim = SUM(Rim*Rim)
            Rsqim = Rim(1)*Rim(1) + Rim(2)*Rim(2) + Rim(3)*Rim(3)

            dxim = Rim(1)
            dyim = Rim(2)
            dzim = Rim(3)

            Rsqmj = Rsqij + Rsqim - 2.0_cd * (dxij*dxim + dyij*dyim + dzij*dzim)

            dxmj = 2.0 * dxim - dxij
            dymj = 2.0 * dyim - dyij
            dzmj = 2.0 * dzim - dzij

            cx =  get_c(Rsqij,Rsqim,Rsqmj)

            ! SRINI2
            call get_s(cx,ity,jty,kty,meatable, sParArray)
            ! call get_s(cx,ity,jty,kty,meatable,sp,s0,s1,s2,s3,s54)

	    ! SRINI
	    RsqimHH = Rsqim + hh 
	    RsqmjHH = Rsqmj + hh

	    ! SRINI2: Screening Factor Ration Array
	    sFactorArray(1) = halfhinv/sParArray(1) 
	    sFactorArray(2) = halfhinv/sParArray(2) 
	    sFactorArray(3) = halfhinv/sParArray(3) 
	    sFactorArray(4) = halfhinv/sParArray(4) 
	    sFactorArray(5) = halfhinv/sParArray(5) 
	    sFactorArray(6) = halfhinv/sParArray(6) 

            do kk = 1,3     ! do kk

	       ! GIBSON
               ! Rpim = Rsqim + h2 * Rim(kk) + hh
               ! Rpmj = Rsqmj + h2 * (Rim(kk) - Rij(kk)) + hh

	       ! SRINI
               Rpim = RsqimHH + h2 * Rim(kk)
               Rpmj = RsqmjHH + h2 * (Rim(kk) - Rij(kk))


               cx = get_c(Rsqij,Rpim,Rpmj)

               ! SRINI2
               call get_s(cx,ity,jty,kty,meatable, sPlusArray)
               ! call get_s(cx,ity,jty,kty,meatable,spp,sp0,sp1,sp2,sp3,sp54)

	       ! GIBSON
               ! Rpim = Rsqim - h2 * Rim(kk) + hh
               ! Rpmj = Rsqmj - h2 * (Rim(kk) - Rij(kk)) + hh


	       ! SRINI
               Rpim = RsqimHH - h2 * Rim(kk)
               Rpmj = RsqmjHH - h2 * (Rim(kk) - Rij(kk))

               cx = get_c(Rsqij,Rpim,Rpmj)

               ! SRINI2
               call get_s(cx,ity,jty,kty,meatable, sMinusArray)
	       ! call get_s(cx,ity,jty,kty,meatable,smp,sm0,sm1,sm2,sm3,sm54)

	       ! SRINI2: DEFAULT is ZERO for these variables=> eliminate "else" from following if()
	       dscrnabkp=0.0_cd
               dscrnabk0=0.0_cd
               dscrnabk1=0.0_cd
               dscrnabk2=0.0_cd
               dscrnabk3=0.0_cd
               dscrnabk54=0.0_cd

               if((sParArray(1) .gt. scnres) .and. (sParArray(1) .lt. 1.0_cd)) then
                  dscrnabkp = (sPlusArray(1)-sMinusArray(1))*screeArray(1)*sFactorArray(1)
               endif

               if((sParArray(2) .gt. scnres) .and. (sParArray(2) .lt. 1.0_cd)) then
                  dscrnabk0 = (sPlusArray(2)-sMinusArray(2))*screeArray(2)*sFactorArray(2)
               endif

               if((sParArray(3) .gt. scnres) .and. (sParArray(3) .lt. 1.0_cd)) then
                  dscrnabk1 = (sPlusArray(3)-sMinusArray(3))*screeArray(3)*sFactorArray(3)
               endif

               if((sParArray(4) .gt. scnres) .and. (sParArray(4) .lt. 1.0_cd)) then
                  dscrnabk2 = (sPlusArray(4)-sMinusArray(4))*screeArray(4)*sFactorArray(4)
               endif

               if((sParArray(5) .gt. scnres) .and. (sParArray(5) .lt. 1.0_cd)) then
                  dscrnabk3 = (sPlusArray(5)-sMinusArray(5))*screeArray(5)*sFactorArray(5)
               endif

               if((sParArray(6) .gt. scnres) .and. (sParArray(6) .lt. 1.0_cd)) then
                  dscrnabk54 = (sPlusArray(6)-sMinusArray(6))*screeArray(6)*sFactorArray(6)
               endif

	       ! SRINI4: This is not needed 
               ! rcosj = rcos(kk)

               drho0s = (cp*rhjp(0)-cm*rhjm(0))/rho0
               drho1s = 0.0_cd
               drho1sg = 0.0_cd
               drho1sgh = 0.0_cd
               drho2s = -(rhjp(2) * tmp0p - rhjm(2) * tmp0m)
               drho3s = 0.0_cd
               drho3sh = 0.0_cd
               drho5s = 0.0_cd

               do m = 1,3    ! do m
                  rcosm = rcos(m)
                  tmp1p = 2.0_cd * ap(m)
                  tmp1gp = 2.0_cd * agp(m)
                  tmp1hp = 2.0_cd * ahp(m)
                  tmp1m = 2.0_cd * am(m)
                  tmp1gm = 2.0_cd * agm(m)
                  tmp1hm = 2.0_cd * ahm(m)

		  ! SRINI
                  drho1s = drho1s + ((tmp1p * rhjp(1))- (tmp1m *rhjm(1)))*rcosm
                  drho1sg= drho1sg + ((tmp1gp*rhjp(3))- (tmp1gm*rhjm(3)))*rcosm
                  drho1sgh=drho1sgh+ ((tmp1hp*rhjp(4))- (tmp1hm*rhjm(4)))*rcosm

                  do l = 1,3    ! do l
                     rcosl = rcos(l)
                     rml = rcosm * rcosl
                     tmp2p = 2.0_cd * bp(l,m)
                     tmp2m = 2.0_cd * bm(l,m)

		     ! SRINI: CHANGED GIBSON'S CODE TO THIS
                     drho2s= drho2s + ((tmp2p*rhjp(2))- (tmp2m*rhjm(2)))*rml

                     do s = 1,3    ! do s
                        tmp3p = 2.0_cd * dp(s,l,m)
                        tmp3m = 2.0_cd * dm(s,l,m)
                        tmp3hp = 2.0_cd * dhp(s,l,m)
                        tmp3hm = 2.0_cd * dhm(s,l,m)

			! SRINI
                        drho3s = drho3s + ((tmp3p*rhjp(3))- (tmp3m *rhjm(3)))*rslm(s,l,m)     !!4.74
                        drho3sh= drho3sh+ ((tmp3hp*rhjp(4))-(tmp3hm*rhjm(4)))*rslm(s,l,m)     !!3.39

                        ! https://software.intel.com/en-us/forums/intel-fortran-compiler-for-linux-and-mac-os-x/topic/509157
                        !$omp SIMD collapse(2) reduction(+:drho5s)
                        do k = 1,3    ! do k
                           do w = 1,3    ! do w
			      ! SRINI
                              drho5s=drho5s+(gp2(w,k,s,l,m)*rhjp4-gm2(w,k,s,l,m)*rhjm4)*rwkslm(w,k,s,l,m)  !! 63.2 -->39.9

                           enddo    ! do m
                        enddo    ! do l
                     enddo    ! do s
                  enddo    ! do k
               enddo    ! do w

               term2 =  drho1s * dscrnabk1 &
                     +  drho2s * dscrnabk2 &
                     + (drho3s - legend00 * drho1sg) * dscrnabk3 &
                     + (drho5s - legend01 * drho3sh + legend02 * drho1sgh) * dscrnabk54

               dfim(kk) = 0.5 * phiij * dscrnabkp + dfrho * (dang1 * drho0s * dscrnabk0 + dang2 * term2)

               if(comp_force.eq.1) then
                  force(kk,mc) = force(kk,mc) + dfim(kk)
               endif

            enddo    ! do kk

            slocal(1,1,mc) = slocal(1,1,mc) - dfim(1)*dxmj
            slocal(2,2,mc) = slocal(2,2,mc) - dfim(2)*dymj
            slocal(3,3,mc) = slocal(3,3,mc) - dfim(3)*dzmj
            slocal(2,3,mc) = slocal(2,3,mc) - dfim(2)*dzmj
            slocal(1,3,mc) = slocal(1,3,mc) - dfim(1)*dzmj
            slocal(1,2,mc) = slocal(1,2,mc) - dfim(1)*dymj

         enddo    ! do mm (neighOnI)

      enddo     ! do jj (neighOnI)

   endif     ! if(neighOnI > 0)

enddo    ! do i (N all atoms)


if(comp_ptvirl.eq.1) then
   do atom = 1,N 
      particleVirial(1,atom) = particleVirial(1,atom) - 0.5_cd*slocal(1,1,atom)
      particleVirial(2,atom) = particleVirial(2,atom) - 0.5_cd*slocal(2,2,atom)
      particleVirial(3,atom) = particleVirial(3,atom) - 0.5_cd*slocal(3,3,atom)
      particleVirial(4,atom) = particleVirial(4,atom) - 0.5_cd*slocal(2,3,atom)
      particleVirial(5,atom) = particleVirial(5,atom) - 0.5_cd*slocal(1,3,atom)
      particleVirial(6,atom) = particleVirial(6,atom) - 0.5_cd*slocal(1,2,atom)
   enddo
endif

if(comp_virial.eq.1) then
   do atom = 1,N
      virial(1) = virial(1) - 0.5_cd*slocal(1,1,atom)
      virial(2) = virial(2) - 0.5_cd*slocal(2,2,atom)
      virial(3) = virial(3) - 0.5_cd*slocal(3,3,atom)
      virial(4) = virial(4) - 0.5_cd*slocal(2,3,atom)
      virial(5) = virial(5) - 0.5_cd*slocal(1,3,atom)
      virial(6) = virial(6) - 0.5_cd*slocal(1,2,atom)
   enddo
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

deallocate(slocal)

! Everything is great
!
ier = KIM_STATUS_OK
42 continue
Compute_Energy_Forces = ier

return
end function Compute_Energy_Forces

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Model driver destroy routine
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer(c_int) function destroy(pkim) bind(c)
use KIM_API_F03
implicit none

!-- Transferred variables
type(c_ptr), intent(in) :: pkim

!-- Local variables
integer(c_int) :: ier, idum
type(model_parameters), pointer :: meatable; type(c_ptr) :: pmeatable

pmeatable = kim_api_get_model_buffer(pkim, ier)
if (ier.ne.KIM_STATUS_OK) then
  idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
    "kim_api_get_model_buffer", ier)
  destroy = ier
  return
endif
call c_f_pointer(pmeatable, meatable)

deallocate(meatable%asubs)
deallocate(meatable%bsubs)
deallocate(meatable%esubs)
deallocate(meatable%res)
deallocate(meatable%embedtab)
deallocate(meatable%phitab)
deallocate(meatable%rhotab)
deallocate(meatable%scrp)
deallocate(meatable%scr0)
deallocate(meatable%scr1)
deallocate(meatable%scr2)
deallocate(meatable%scr3)
deallocate(meatable%scr54)

deallocate( meatable )

destroy = KIM_STATUS_OK
return

end function destroy

end module MSMEAM

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Model driver initialization routine (REQUIRED)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer(c_int) function model_driver_init(pkim, pparamfile, nmstrlen, &
                                          numparamfiles) bind(c)
use, intrinsic :: iso_c_binding
use MSMEAM
use KIM_API_F03
implicit none
integer(c_int), parameter :: cd = c_double ! used for literal constants

!-- Transferred variables
type(c_ptr),        intent(in) :: pkim
type(c_ptr), value, intent(in) :: pparamfile
integer(c_int),     intent(in) :: nmstrlen
integer(c_int),     intent(in) :: numparamfiles
character(len=nmstrlen), pointer :: paramfile(:)
character (len=80) :: error_message

!-- Local variables
integer(c_int), parameter :: zero=0
integer(c_int), parameter :: one=1
integer(c_int), parameter :: two=2
integer(c_int), parameter :: three=3
integer(c_int), parameter :: four=4
integer(c_int), parameter :: five=5
integer(c_int)            :: i,j,ier,idum,itab,nelem,mil1
integer(c_int)            :: ispin,it,jt,kt,return_error
type(model_parameters), pointer :: meatable

integer(c_int), parameter :: nelmax=2

! SRINI: Gibson defined them in "subroutine rscrn(rij,meatable,scrn)"
real(c_double),parameter :: frcutSr=0.9_cd

! SRINI: Move from subroutine get_c()
real(c_double), parameter           :: ctabmaxSr = 3.0_cd

real(c_double)            :: rx,embedtemp,phitemp

 ! SRINI
real(c_double)            :: rcutsq

!-- KIM variables
real(c_double), pointer :: cutoff; type(c_ptr) :: pcutoff
character(len=KIM_KEY_STRING_LENGTH) :: NBC_Method

! assume all is well
!
return_error = KIM_STATUS_OK

call c_f_pointer(pparamfile, paramfile, [numparamfiles])

! find first null character and write spaces everywhere afterward
do i = 1, numparamfiles
  j = index(paramfile(1),char(0))
  paramfile(i)(j:)=" "
end do

! store function pointers in KIM object
call kim_api_setm_method(pkim, ier, &
     "compute", one, c_funloc(Compute_Energy_Forces), 1, &
     "destroy", one, c_funloc(destroy),               1)
if (ier.lt.KIM_STATUS_OK) then
   idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                               "kim_api_setm_method", ier)
   return_error = ier
   goto 1000
endif


! allocate parameters and register in KIM object
allocate( meatable )
call kim_api_set_model_buffer(pkim, c_loc(meatable), ier)
if (ier.ne.KIM_STATUS_OK) then
  idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
    "kim_api_set_model_buffer", ier)
  return_error = ier
  goto 1000
endif

! Read in model parameters from parameter file
!
open(10,file=paramfile(1),status="old")
read(10,*,iostat=ier,err=100)
read(10,*,iostat=ier,err=100)
read(10,*,iostat=ier,err=100) &
      &  meatable%ntab,meatable%nel,meatable%rtabmax,meatable%rhotabmin,meatable%rhotabmax
read(10,*,iostat=ier,err=100) meatable%model_cutoff

nelem=meatable%nel
mil1=meatable%ntab

allocate(meatable%asubs(one:nelem))
allocate(meatable%bsubs(one:nelem))
allocate(meatable%esubs(one:nelem,one:nelem))
allocate(meatable%res(one:nelem,one:nelem))
allocate(meatable%embedtab(one:mil1,one:nelem))
allocate(meatable%phitab(one:mil1,one:nelem,one:nelem))
allocate(meatable%rhotab(one:mil1,zero:four,one:two,one:nelem))
allocate(meatable%scrp(one:mil1,one:nelem,one:nelem,one:nelem))
allocate(meatable%scr0(one:mil1,one:nelem,one:nelem,one:nelem))
allocate(meatable%scr1(one:mil1,one:nelem,one:nelem,one:nelem))
allocate(meatable%scr2(one:mil1,one:nelem,one:nelem,one:nelem))
allocate(meatable%scr3(one:mil1,one:nelem,one:nelem,one:nelem))
allocate(meatable%scr54(one:mil1,one:nelem,one:nelem,one:nelem))

do i=1,meatable%nel
    read(10,*,iostat=ier,err=100) &
      &  meatable%asubs(i),meatable%bsubs(i),meatable%res(i,i),meatable%esubs(i,i)
    meatable%rrtabmax=meatable%rtabmax*meatable%res(i,i)
         ispin=1
              do itab=1,mil1
                   read(10,*,iostat=ier,err=100) &
                   meatable%rhotab(itab,zero,ispin,i), &
                   meatable%rhotab(itab,one,ispin,i), &
                   meatable%rhotab(itab,two,ispin,i), &
                   meatable%rhotab(itab,three,ispin,i),&
                   meatable%rhotab(itab,four,ispin,i), &
                   embedtemp
                   meatable%embedtab(itab,i)=embedtemp*meatable%esubs(i,i)

                   !meatable%rhotab(itab,one,ispin,i)=0.0_cd
                   !meatable%rhotab(itab,two,ispin,i)=0.0_cd
                   !meatable%rhotab(itab,three,ispin,i)=0.0_cd
                   !meatable%rhotab(itab,four,ispin,i)=0.0_cd

              enddo
         ispin=2
              do itab=1,mil1
                   read(10,*,iostat=ier,err=100) &
                   meatable%rhotab(itab,zero,ispin,i), &
                   meatable%rhotab(itab,one,ispin,i), &
                   meatable%rhotab(itab,two,ispin,i), &
                   meatable%rhotab(itab,three,ispin,i), &
                   meatable%rhotab(itab,four,ispin,i)

                   !meatable%rhotab(itab,one,ispin,i)=0.0_cd
                   !meatable%rhotab(itab,two,ispin,i)=0.0_cd
                   !meatable%rhotab(itab,three,ispin,i)=0.0_cd
                   !meatable%rhotab(itab,four,ispin,i)=0.0_cd

              enddo
              do itab=1,mil1
                   read(10,*,iostat=ier,err=100) phitemp
                   meatable%phitab(itab,i,i)=phitemp*meatable%esubs(i,i)
                   rx=(dble(itab)/(dble(meatable%ntab)-1.0d0)*&
                     &(meatable%rrtabmax/meatable%res(i,i)))

              enddo
              it=i
              jt=i
              kt=i
              do itab=1,mil1
                  read(10,*,iostat=ier,err=100) &
                  meatable%scrp(itab,it,jt,kt), &
                  meatable%scr0(itab,it,jt,kt), &
                  meatable%scr1(itab,it,jt,kt), &
                  meatable%scr2(itab,it,jt,kt), &
                  meatable%scr3(itab,it,jt,kt), &
                  meatable%scr54(itab,it,jt,kt) 

                  !meatable%scrp(itab,it,jt,kt)=1.0_cd
                  !meatable%scr0(itab,it,jt,kt)=1.0_cd
                  !meatable%scr1(itab,it,jt,kt)=1.0_cd
                  !meatable%scr2(itab,it,jt,kt)=1.0_cd
                  !meatable%scr3(itab,it,jt,kt)=1.0_cd
                  !meatable%scr54(itab,it,jt,kt)=1.0_cd

              enddo
enddo
close(10)

! Determine neighbor list boundary condition (NBC)
! and half versus full mode:
! *****************************
! * HalfOrFull = 1 -- Half
! *            = 2 -- Full
! *****************************
!
!
ier = kim_api_get_nbc_method(pkim, NBC_Method)
if (ier.lt.KIM_STATUS_OK) then
   idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                               "kim_api_get_nbc_method", ier)
   goto 1000
endif
if (index(NBC_Method,"NEIGH_PURE_F").eq.1) then
    meatable%NBC = 2
    meatable%HalfOrFull = 2
else
    ier = KIM_STATUS_FAIL
    idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                              "Unknown NBC method.", ier)
    goto 1000
endif


! Determine neighbor list handling mode
!
   !*****************************
   !* IterOrLoca = 1 -- Iterator
   !*            = 2 -- Locator
   !*****************************
   meatable%IterOrLoca = kim_api_get_neigh_mode(pkim, ier)
   if (ier.lt.KIM_STATUS_OK) then
      idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                    "kim_api_get_neigh_mode", ier)
      goto 1000
   endif
   if (meatable%IterOrLoca.ne.1 .and. meatable%IterOrLoca.ne.2) then
      ier = KIM_STATUS_FAIL
      write(error_message,'(a,i1)') &
         'Unsupported IterOrLoca mode = ',meatable%IterOrLoca
      idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                  error_message, ier)
      goto 1000
   endif


goto 200
100 continue
! reading parameters failed
idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                            "Unable to read potential parameters, kimerror = ",KIM_STATUS_FAIL)
return_error = KIM_STATUS_FAIL
goto 1000

200 continue

! store model cutoff in KIM object
pcutoff =  kim_api_get_data(pkim,"cutoff",ier)
if (ier.lt.KIM_STATUS_OK) then
   idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                               "kim_api_get_data", ier)
   return_error = ier
   goto 1000
endif
call c_f_pointer(pcutoff, cutoff)
cutoff = meatable%model_cutoff
rcutsq = cutoff*cutoff
meatable%model_cutsq = rcutsq     ! meatable%model_cutoff**2


! SRINI
meatable%frcut   = frcutSr
meatable%ctabmax = ctabmaxSr

! SRINI:  Radial Screening Factor computed and saved here
! frcutRcSq = frcut*rcutsq; inv1minusFrcutRcSq = 1/((1.0-frcut)*rcutsq)
meatable%frcutRcSq          = (meatable%frcut)*rcutsq
meatable%inv1minusFrcutRcSq = 1.0_cd/((1.0-meatable%frcut)*rcutsq)
meatable%rtabscale          = ( dble(meatable%ntab) - 1.0_cd )/meatable%rrtabmax
meatable%ctabscale          = ( dble(meatable%ntab) - 1.0_cd )/meatable%ctabmax
meatable%rhotabscale        = ( dble(meatable%ntab) - 1.0_cd )/meatable%rhotabmax


1000 continue
model_driver_init = return_error
return

end function model_driver_init
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
