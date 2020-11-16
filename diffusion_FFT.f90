subroutine diffusion (t,g,p,l)

! to compute hillslope diffusion

! here we solve the diffusion equation using an ADI (Alternating Direction Implicit) scheme
! that is implicit and O(n)
! The solution has been checked against an analytical solution for a point source problem

use omp_lib
use definitions

implicit none

type (topography) t
type (topology) g
type (param) p
type (law) l

logical xcycle,ycycle,soil
character c4*4
double precision, dimension(:), allocatable :: f,diag,sup,inf,res,tint,kd,sprev,hprev
integer i,j,i1,i2,j1,j2,ij,ipj,imj,ijp,ijm,ijk,iter,niter
double precision fact,factxp,factxm,factyp,factym,s0,p0,dh,slope
double precision f1,s1,f2,s2,f3,s3,f4,s4

if (l%kd.lt.tiny(l%kd)) return

xcycle=.false.
ycycle=.false.
write (c4,'(i4)') p%ibc
if (c4(1:1).eq.'0' .and. c4(3:3).eq.'0') ycycle=.true.
if (c4(2:2).eq.'0' .and. c4(4:4).eq.'0') xcycle=.true.

!$OMP parallel shared(g,t)
!$OMP do private (i) schedule(dynamic,p%chunk)
  do i=1,p%nn
  if (g%bc(i)) t%pmask(1,i)=.true.
  enddo
!$OMP end do nowait
!$OMP end parallel

call solve_diffusion_FFT (t%h,p%nx,p%ny,p%dt,l%kd,p%xl,p%yl,p%num_threads)

return

end subroutine diffusion
