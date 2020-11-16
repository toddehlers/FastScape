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
integer i,j,i1,i2,j1,j2,ij,ipj,imj,ijp,ijm,ijk,iter,niter,idir
double precision fact,factxp,factxm,factyp,factym,s0,p0,dh,slope
double precision f1,s1,f2,s2,f3,s3,f4,s4,dir,find_kd
external find_kd

if (l%kd.lt.tiny(l%kd)) return

xcycle=.false.
ycycle=.false.
write (c4,'(i4)') p%ibc
if (c4(1:1).eq.'0' .and. c4(3:3).eq.'0') ycycle=.true.
if (c4(2:2).eq.'0' .and. c4(4:4).eq.'0') xcycle=.true.

soil=.true.
if (l%p0.lt.tiny(l%p0)) soil=.false.

!$OMP parallel shared(g,t)
!$OMP do private (i) schedule(dynamic,p%chunk)
  do i=1,p%nn
  if (g%bc(i)) t%pmask(1,i)=.true.
  enddo
!$OMP end do nowait
!$OMP end parallel

! NO SOIL

if (.not.soil) then

allocate (tint(p%nn),kd(p%nn))

!$OMP parallel shared(l,kd,tint)
!$OMP workshare
tint=t%h
!$OMP end workshare nowait
!$OMP end parallel

!$OMP parallel shared(l,t,p)
!$OMP do private (ij) schedule(dynamic,p%chunk)
  do ij=1,p%nn
  kd(ij)=find_kd(ij,l,t,p)
  enddo
!$OMP end do nowait
!$OMP end parallel

j1=2
j2=p%ny-1
f4=1.d0
s4=0.d0
  if (c4(4:4).ne.'1') then
  f4=0.d0
  s4=-1.d0
  endif
f2=1.d0
s2=0.d0
  if (c4(2:2).ne.'1') then
  f2=0.d0
  s2=-1.d0
  endif
!$OMP parallel shared(p,t,kd,j1,j2,tint,c4,f2,s2,f4,s4) &
!$OMP private(f,diag,sup,inf,res,i,j,ij,ipj,imj,ijm,ijp,factxp,factxm,factyp,factym)
allocate (f(p%nx),diag(p%nx),sup(p%nx),inf(p%nx),res(p%nx))
!$OMP do schedule(dynamic,p%ny/p%num_threads)
do j=j1,j2
  do i=2,p%nx-1
  ij=(j-1)*p%nx+i
  ipj=(j-1)*p%nx+i+1
  if (i.eq.p%nx) ipj=(j-1)*p%nx+1
  imj=(j-1)*p%nx+i-1
  if (i.eq.1) imj=(j-1)*p%nx+p%nx
  ijp=(j)*p%nx+i
  if (j.eq.p%ny) ijp=i
  ijm=(j-2)*p%nx+i
  if (j.eq.1) ijm=(p%ny-1)*p%nx+i
  factxp=(kd(ipj)+kd(ij))/2.d0*(p%dt/2.d0)/p%dx**2
  factxm=(kd(imj)+kd(ij))/2.d0*(p%dt/2.d0)/p%dx**2
  factyp=(kd(ijp)+kd(ij))/2.d0*(p%dt/2.d0)/p%dy**2
  factym=(kd(ijm)+kd(ij))/2.d0*(p%dt/2.d0)/p%dy**2
  diag(i)=1.d0+factxp+factxm
  sup(i)=-factxp
  inf(i)=-factxm
  f(i)=tint(ij)+factyp*tint(ijp)-(factyp+factym)*tint(ij)+factym*tint(ijm)
  enddo
diag(1)=1.d0
sup(1)=s4
f(1)=tint((j-1)*p%nx+1)*f4
diag(p%nx)=1.d0
inf(p%nx)=s2
f(p%nx)=tint((j-1)*p%nx+p%nx)*f2
call tridag (inf,diag,sup,f,res,p%nx)
  do i=1,p%nx
  ij=(j-1)*p%nx+i
    if (.not.t%pmask(1,ij)) then
    tint(ij)=res(i)
    endif
  enddo
enddo
!$OMP end do nowait
deallocate (f,diag,sup,inf,res)
!$OMP end parallel

i1=2
i2=p%nx-1
f1=1.d0
s1=0.d0
  if (c4(1:1).ne.'1') then
  f1=0.d0
  s1=-1.d0
  endif
f3=1.d0
s3=0.d0
  if (c4(3:3).ne.'1') then
  f3=0.d0
  s3=-1.d0
  endif
!$OMP parallel shared(p,t,kd,i1,i2,tint,c4,f1,s1,f3,s3)  &
!$OMP private (f,diag,sup,inf,res,i,j,ij,ijm,ijp,imj,ipj,factxp,factxm,factyp,factym)
allocate (f(p%ny),diag(p%ny),sup(p%ny),inf(p%ny),res(p%ny))
!$OMP do schedule(dynamic,p%nx/p%num_threads)
do i=i1,i2
  do j=2,p%ny-1
  ij=(j-1)*p%nx+i
  ipj=(j-1)*p%nx+i+1
  if (i.eq.p%nx) ipj=(j-1)*p%nx+1
  imj=(j-1)*p%nx+i-1
  if (i.eq.1) imj=(j-1)*p%nx+p%nx
  ijp=(j)*p%nx+i
  if (j.eq.p%ny) ijp=i
  ijm=(j-2)*p%nx+i
  if (j.eq.1) ijm=(p%ny-1)*p%nx+i
  factxp=(kd(ipj)+kd(ij))/2.d0*(p%dt/2.d0)/p%dx**2
  factxm=(kd(imj)+kd(ij))/2.d0*(p%dt/2.d0)/p%dx**2
  factyp=(kd(ijp)+kd(ij))/2.d0*(p%dt/2.d0)/p%dy**2
  factym=(kd(ijm)+kd(ij))/2.d0*(p%dt/2.d0)/p%dy**2
  diag(j)=1.+factyp+factym
  sup(j)=-factyp
  inf(j)=-factym
  f(j)=tint(ij)+factxp*tint(ipj)-(factxp+factxm)*tint(ij)+factxm*tint(imj)
  enddo
diag(1)=1.d0
sup(1)=s1
f(1)=tint(i)*f1
diag(p%ny)=1.d0
inf(p%ny)=s3
f(p%ny)=tint((p%ny-1)*p%nx+i)*f3
call tridag (inf,diag,sup,f,res,p%ny)
  do j=1,p%ny
  ij=(j-1)*p%nx+i
    if (.not.t%pmask(1,ij)) then
    tint(ij)=res(j)
    endif
  enddo
enddo
!$OMP end do nowait
deallocate (f,diag,sup,inf,res)
!$OMP end parallel

!$OMP parallel shared(l,kd,tint)
!$OMP workshare
t%h=tint
!$OMP end workshare nowait
!$OMP end parallel

deallocate (kd,tint)

else

! WITH SOIL

s0=l%s0
p0=l%p0
!$OMP parallel shared(t,p,p0,s0)
!$OMP workshare
t%s=t%s+p%dt*p0*exp(-t%s/s0)
!$OMP end workshare nowait
!$OMP end parallel

!$OMP parallel shared(l,t,p)
!$OMP do private (ij) schedule(dynamic,p%chunk)
  do ij=1,p%nn
  kd(ij)=find_kd(ij,l,t,p)
  enddo
!$OMP end do nowait
!$OMP end parallel

!$OMP parallel shared(t,p)
!$OMP do private (ij) schedule(dynamic,p%chunk)
  do ij=1,p%nn
    if (t%pmask(1,ij)) then
    t%h(ij)=max(p%sea_level,t%h(ij)-t%s(ij))
    t%s(ij)=0.d0
    endif
  enddo
!$OMP end do nowait
!$OMP end parallel  

niter=3
allocate (hprev(p%nn),sprev(p%nn))

!$OMP parallel shared(t,p,hprev,sprev)
!$OMP do private (i) schedule(dynamic,p%chunk)
  do i=1,p%nn
  hprev(i)=t%h(i)
  sprev(i)=t%s(i)
  enddo
!$OMP end do nowait
!$OMP end parallel

do iter=1,niter

!$OMP parallel shared(tint,t)
!$OMP workshare
tint=t%h
!$OMP end workshare nowait
!$OMP end parallel

! Threshold slope

if (l%slope.gt.tiny(l%slope)) then
!$OMP parallel shared(t,p,g,kd,l)
!$OMP do private (ij,ijk,slope) schedule(dynamic,p%chunk)
  do ij=1,p%nn
  ijk=g%receiver(ij)
    if (ijk.ne.ij .and. t%length(ij).gt.0.d0) then
    slope=max(0.d0,(t%h(ij)-t%h(ijk))/t%length(ij))
    kd(ij)=min(0.1d0,abs(l%kd/(1.d0-(slope/0.58d0)**2)))
    endif
  enddo
!$OMP end do nowait
!$OMP end parallel
endif

! adjust kd to only transport soil

j1=2
j2=p%ny-1
i1=2
i2=p%nx-1
!$OMP parallel shared(kd,p,t,l,j1,j2,i1,i2) &
!$OMP private(i,j,ij,ijm,imj,ijp,ipj,factxp,factxm,factyp,factym,dh)
!$OMP do  schedule(dynamic,p%ny/p%num_threads)
do j=j1,j2
  do i=i1,i2
  ij=(j-1)*p%nx+i
  ipj=(j-1)*p%nx+i+1
  if (i.eq.p%nx) ipj=(j-1)*p%nx+1
  imj=(j-1)*p%nx+i-1
  if (i.eq.1) imj=(j-1)*p%nx+p%nx
  ijp=(j)*p%nx+i
  if (j.eq.p%ny) ijp=i
  ijm=(j-2)*p%nx+i
  if (j.eq.1) ijm=(p%ny-1)*p%nx+i
  factxp=l%kd*(p%dt)/p%dx**2
  factxm=l%kd*(p%dt)/p%dx**2
  factyp=l%kd*(p%dt)/p%dy**2
  factym=l%kd*(p%dt)/p%dy**2
  dh=factyp*t%h(ijp)+factxp*t%h(ipj)-(factyp+factym+factxp+factxm)*t%h(ij)+factym*t%h(ijm)+factxm*t%h(imj)
  if (dh.lt.-t%s(ij)) kd(ij)=min(kd(ij),-kd(ij)*t%s(ij)/dh)
  enddo
enddo
!$OMP end do nowait
!$OMP end parallel

j1=2
j2=p%ny-1
!$OMP parallel shared(p,t,kd,factxp,factxm,factyp,factym,j1,j2,tint) private(f,diag,sup,inf,res)
allocate (f(p%nx),diag(p%nx),sup(p%nx),inf(p%nx),res(p%nx))
!$OMP do private(i,j,ij,ijm,ijp,ijk,slope) &
!$OMP schedule(dynamic,p%ny/p%num_threads)
do j=j1,j2
  do i=2,p%nx-1
  ij=(j-1)*p%nx+i
  ipj=(j-1)*p%nx+i+1
  if (i.eq.p%nx) ipj=(j-1)*p%nx+1
  imj=(j-1)*p%nx+i-1
  if (i.eq.1) imj=(j-1)*p%nx+p%nx
  ijp=(j)*p%nx+i
  if (j.eq.p%ny) ijp=i
  ijm=(j-2)*p%nx+i
  if (j.eq.1) ijm=(p%ny-1)*p%nx+i
  factxp=(kd(ipj)+kd(ij))/2.d0*(p%dt/2.d0)/p%dx**2
  factxm=(kd(imj)+kd(ij))/2.d0*(p%dt/2.d0)/p%dx**2
  factyp=(kd(ijp)+kd(ij))/2.d0*(p%dt/2.d0)/p%dy**2
  factym=(kd(ijm)+kd(ij))/2.d0*(p%dt/2.d0)/p%dy**2
  diag(i)=1.+factxp+factxm
  sup(i)=-factxp
  inf(i)=-factxm
  f(i)=hprev(ij)+factyp*hprev(ijp)-(factyp+factym)*hprev(ij)+factym*hprev(ijm)
  enddo
diag(1)=1.d0
sup(1)=0.d0
f(1)=t%h((j-1)*p%nx+1)
diag(p%nx)=1.d0
inf(p%nx)=0.d0
f(p%nx)=t%h((j-1)*p%nx+p%nx)
call tridag (inf,diag,sup,f,res,p%nx)
  do i=1,p%nx
  ij=(j-1)*p%nx+i
    if (.not.t%pmask(1,ij)) then
    slope=0.d0
      if (.not.g%bc(ij).and.t%h(ij).gt.p%sea_level) then
      ijk=g%receiver(ij)
      if (ijk.ne.ij .and. t%length(ij).gt.0.d0) slope=atan2(t%h(ij)-t%h(ijk),t%length(ij))
      endif
    t%s(ij)=max(0.d0,sprev(ij)+(res(i)-t%h(ij))*cos(slope))
    tint(ij)=res(i)
    endif
  enddo
enddo
!$OMP end do nowait
deallocate (f,diag,sup,inf,res)
!$OMP end parallel

i1=2
i2=p%nx-1
!$OMP parallel shared(p,t,kd,factxp,factxm,factyp,factym,j1,j2,tint) private (f,diag,sup,inf,res)
allocate (f(p%nx),diag(p%nx),sup(p%nx),inf(p%nx),res(p%nx))
!$OMP do private(i,j,ij,ijm,ijp,slope) schedule(dynamic,p%nx/p%num_threads)
do i=i1,i2
  do j=2,p%ny-1
  ij=(j-1)*p%nx+i
  ipj=(j-1)*p%nx+i+1
  if (i.eq.p%nx) ipj=(j-1)*p%nx+1
  imj=(j-1)*p%nx+i-1
  if (i.eq.1) imj=(j-1)*p%nx+p%nx
  ijp=(j)*p%nx+i
  if (j.eq.p%ny) ijp=i
  ijm=(j-2)*p%nx+i
  if (j.eq.1) ijm=(p%ny-1)*p%nx+i
  factxp=(kd(ipj)+kd(ij))/2.d0*(p%dt/2.d0)/p%dx**2
  factxm=(kd(imj)+kd(ij))/2.d0*(p%dt/2.d0)/p%dx**2
  factyp=(kd(ijp)+kd(ij))/2.d0*(p%dt/2.d0)/p%dy**2
  factym=(kd(ijm)+kd(ij))/2.d0*(p%dt/2.d0)/p%dy**2
  diag(j)=1.+factyp+factym
  sup(j)=-factyp
  inf(j)=-factym
  f(j)=tint(ij)+factxp*tint(ipj)-(factxp+factxm)*tint(ij)+factxm*tint(imj)
  enddo
diag(1)=1.d0
sup(1)=0.d0
f(1)=tint(i)
diag(p%ny)=1.d0
inf(p%ny)=0.d0
f(p%ny)=tint((p%ny-1)*p%nx+i)
call tridag (inf,diag,sup,f,res,p%ny)
  do j=1,p%ny
  ij=(j-1)*p%nx+i
    if (.not.t%pmask(1,ij)) then
    slope=0.d0
      if (.not.g%bc(ij).and.t%h(ij).gt.p%sea_level) then
      ijk=g%receiver(ij)
      if (ijk.ne.ij .and. t%length(ij).gt.0.d0) slope=atan2(t%h(ij)-t%h(ijk),t%length(ij))
      endif
    t%s(ij)=max(0.d0,sprev(ij)+(res(j)-tint(ij))*cos(slope))
    t%h(ij)=res(j)
    endif
  enddo
enddo
!$OMP end do nowait
deallocate (f,diag,sup,inf,res)
!$OMP end parallel

enddo

deallocate (hprev,sprev)
deallocate (tint,kd)

endif

return

end subroutine diffusion
