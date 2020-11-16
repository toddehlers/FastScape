subroutine interpolate (x,y,nu,u,val)

! subroutine used to interpolate fields such as initial_topography or uplift
! using a constant value, a bi-linear interpolation or a set of four bi-linear interpolations
! to allow the user to specify discontinuities in the field

! nu sets the type of interpolation:
! nu = 0 means constant value (needs 1 value to be interpolated)
! nu = 1 means bi-linear interpolation (needs 4 values to be interpolated)
! nu = 2 means four bi-linear interpolations (needs 16 values to be interpolated)

implicit none

double precision x,y,xp,yp,val,h1,h2,h3,h4
double precision u(16),t(16),f(16)
integer nu,nt

if (nu.eq.0) then

val=u(1)

elseif (nu.eq.1) then

h1=(1.d0-x)*(1.d0-y)/4.d0
h2=(1.d0+x)*(1.d0-y)/4.d0
h3=(1.d0+x)*(1.d0+y)/4.d0
h4=(1.d0-x)*(1.d0+y)/4.d0
val=h1*u(1)+h2*u(2)+h3*u(3)+h4*u(4)

elseif (nu.eq.2) then

  if (x.le.0.d0 .and. y.le.0.d0) then
  xp=x*2.d0+1.d0 ; yp=y*2.d0+1.d0
  h1=(1.d0-xp)*(1.d0-yp)/4.d0
  h2=(1.d0+xp)*(1.d0-yp)/4.d0
  h3=(1.d0+xp)*(1.d0+yp)/4.d0
  h4=(1.d0-xp)*(1.d0+yp)/4.d0
  val=h1*u(1)+h2*u(2)+h3*u(3)+h4*u(4)
  elseif (x.ge.0.d0 .and. y.le.0.d0) then
  xp=x*2.d0-1.d0 ; yp=y*2.d0+1.d0
  h1=(1.d0-xp)*(1.d0-yp)/4.d0
  h2=(1.d0+xp)*(1.d0-yp)/4.d0
  h3=(1.d0+xp)*(1.d0+yp)/4.d0
  h4=(1.d0-xp)*(1.d0+yp)/4.d0
  val=h1*u(5)+h2*u(6)+h3*u(7)+h4*u(8)
  elseif (x.ge.0.d0 .and. y.ge.0.d0) then
  xp=x*2.d0-1.d0 ; yp=y*2.d0-1.d0
  h1=(1.d0-xp)*(1.d0-yp)/4.d0
  h2=(1.d0+xp)*(1.d0-yp)/4.d0
  h3=(1.d0+xp)*(1.d0+yp)/4.d0
  h4=(1.d0-xp)*(1.d0+yp)/4.d0
  val=h1*u(9)+h2*u(10)+h3*u(11)+h4*u(12)
  elseif (x.le.0.d0 .and. y.ge.0.d0) then
  xp=x*2.d0+1.d0 ; yp=y*2.d0-1.d0
  h1=(1.d0-xp)*(1.d0-yp)/4.d0
  h2=(1.d0+xp)*(1.d0-yp)/4.d0
  h3=(1.d0+xp)*(1.d0+yp)/4.d0
  h4=(1.d0-xp)*(1.d0+yp)/4.d0
  val=h1*u(13)+h2*u(14)+h3*u(15)+h4*u(16)
  endif

else

print*,'this value of nu ',nu,' is not implemented'
stop

endif

end subroutine interpolate
