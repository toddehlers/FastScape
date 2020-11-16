subroutine precipitation_function (t,p)

! first half of the precipitation.f90 routine in which the user can add lines read from the input file
! thus is where the user sets the precipitation rate

use omp_lib
use definitions

implicit none

type (topography) t
type (param) p

double precision x,y,h,xl,yl,time,precipitation,val_t,rat,dt
double precision w1,w2,w3,w4,w5,w6,w7,w8,w9,w0
integer i1,i2,i3,i4,i5,i6,i7,i8,i9,i0
integer i,j,ij

xl=p%xl
yl=p%yl
time=p%time
dt=p%dt

if (p%precipitation_n.ne.-1) then

  if (p%precipitation_nt.eq.0) then
  val_t=1.d0
  elseif (p%time.lt.p%precipitation_t(1)) then
  val_t=p%precipitation_f(1)
  elseif (p%time.gt.p%precipitation_t(p%precipitation_nt)) then
  val_t=p%precipitation_f(p%precipitation_nt)
  else
    do i=1,p%precipitation_nt-1
      if ((p%time-p%precipitation_t(i))*(p%time-p%precipitation_t(i+1)).le.0.d0) then
        if (abs(p%precipitation_t(i)-p%precipitation_t(i+1)).lt.tiny(x)) then
        val_t=(p%precipitation_f(i)+p%precipitation_f(i+1))/2.d0
        else
        rat=(p%time-p%precipitation_t(i))/(p%precipitation_t(i)-p%precipitation_t(i+1))
        val_t=p%precipitation_f(i)+rat*(p%precipitation_f(i+1)-p%precipitation_f(i))
        endif
      endif
    enddo
  endif
!$OMP parallel shared(t,p,val_t)
!$OMP do private (i,j,ij,x,y,precipitation) schedule(dynamic,p%ny/p%num_threads)
  do j=1,p%ny
    do i=1,p%nx
    ij=i+(j-1)*p%nx
    x=p%xl*float(i-1)/(p%nx-1)
    y=p%yl*float(j-1)/(p%ny-1)
    call interpolate (x/p%xl*2.d0-1.d0,y/p%yl*2.d0-1.d0, &
                      p%precipitation_n,p%precipitation_v, &
                      precipitation)
    t%discharge(ij)=precipitation*p%dx*p%dy*val_t
    enddo
  enddo
!$OMP end do nowait
!$OMP end parallel

else

!$OMP parallel shared(t,p,time,dt,xl,yl)
!$OMP do private (i,j,ij,x,y,w1,w2,w3,w4,w5,w6,w7,w8,w9,w0,i1,i2,i3,i4,i5,i6,i7,i8,i9,i0) &
!$OMP schedule(dynamic,p%ny/p%num_threads)
  do j=1,p%ny
    do i=1,p%nx
    ij=i+(j-1)*p%nx
    x=p%xl*float(i-1)/(p%nx-1)
    y=p%yl*float(j-1)/(p%ny-1)
    precipitation=1.d0
    t%discharge(ij)=precipitation*p%dx*p%dy
    enddo
  enddo
!$OMP end do nowait
!$OMP end parallel

endif

return

end subroutine precipitation_function
