subroutine local_minima (t,g,p)

! subroutine to remove the local minima from the landscape by simply searching through
! the nodes and finding sills that connect minima catchments to one of their neighbours
! note that we need to go through this procedure iteratively because the receiving
! catchment needs to be connected to a base level node to avoid two or more catchment
! forming a loop that does not connect to the base level

use omp_lib
use definitions

implicit none

type (topography) t
type (topology) g
type (param) p

integer ij,icatch,jcatch,ijk,i,j,ii,jj,nc,niter,nij
double precision, dimension(:), allocatable :: hsill
integer, dimension(:), allocatable :: catch,ijl
logical xcycle,ycycle
character c4*4

xcycle=.false.
ycycle=.false.
write (c4,'(i4)') p%ibc
if (c4(1:1).eq.'1' .and. c4(3:3).eq.'1') ycycle=.true.
if (c4(2:2).eq.'1' .and. c4(4:4).eq.'1') xcycle=.true.

allocate (hsill(p%nn),catch(p%nn))

!$OMP parallel shared(hsill,catch)
!$OMP workshare
hsill=1.d6
catch=-1
!$OMP end workshare nowait
!$OMP end parallel

niter=0
111 continue

!!!!!allocate (ijl(p%nn))

!!!!!!$OMP parallel shared(p,j,ijl)
!!!!!!$OMP do private(i,j,ij,icatch)
!!!!!  do j=2,p%ny-1
!!!!!    do i=2,p%nx-1
!!!!!    ij=i+(j-1)*p%nx
!!!!!    icatch=g%catchment(ij)
!!!!!    ijl(ij)=0
!!!!!    if (icatch.lt.0) ijl(ij)=ij
!!!!!    enddo
!!!!!  enddo
!!!!!!$OMP end do nowait
!!!!!!$OMP end parallel

!!!!!!$OMP parallel shared(g,p,t,hsill,catch,ijl)
!!!!!!$OMP do private(i,ij,icatch,jcatch,jj,ii,ijk)
!!!!!    do i=1,p%nn
!!!!!    if (ijl(i).ne.0) then
!!!!!    ij=ijl(i)
!!!!!    icatch=g%catchment(ij)
!!!!!      do jj=-1,1
!!!!!        do ii=-1,1
!!!!!        ijk=ij+ii+jj*p%nx
!!!!!        jcatch=g%catchment(ijk)
!!!!!          if (jcatch.ne.-icatch .and. jcatch.gt.0 .and. ij.ne.ijk .and. t%h(ijk).lt. hsill(-icatch)) then
!!!!!            hsill(-icatch)=t%h(ijk)
!!!!!            g%receiver(-icatch)=ijk
!!!!!            catch(-icatch)=jcatch
!!!!!            t%length(-icatch)=-1.d0
!!!!!          endif
!!!!!        enddo
!!!!!      enddo
!!!!!      endif
!!!!!    enddo
!!!!!!$OMP end do nowait
!!!!!!$OMP end parallel

!!!!!deallocate (ijl)

!$OMP parallel shared(g,p,t,hsill)
!$OMP do private(i,j,ij,icatch,jcatch,jj,ii,ijk)
  do j=2,p%ny-1
    do i=2,p%nx-1
    ij=i+(j-1)*p%nx
    icatch=g%catchment(ij)
      if (icatch.lt.0) then
        do jj=-1,1
          do ii=-1,1
          ijk=ij+ii+jj*p%nx
          jcatch=g%catchment(ijk)
            if (jcatch.ne.-icatch .and. jcatch.gt.0 .and. ij.ne.ijk .and. t%h(ijk).lt. hsill(-icatch)) then
            hsill(-icatch)=t%h(ijk)
            g%receiver(-icatch)=ijk
            catch(-icatch)=jcatch
            t%length(-icatch)=-1.d0
            endif
          enddo
        enddo
      endif
    enddo
  enddo
!$OMP end do nowait
!$OMP end parallel

if (ycycle) then

j=1
!$OMP parallel shared(g,p,t,j,hsill)
!$OMP do private(i,ij,icatch,jcatch,jj,ii,ijk)
    do i=2,p%nx-1
    ij=i+(j-1)*p%nx
    icatch=g%catchment(ij)
      if (icatch.lt.0) then
        do jj=-1,1
          do ii=-1,1
          ijk=ij+ii+jj*p%nx
          if (jj.eq.-1) ijk=ij+ii+(jj+p%ny)*p%nx
          jcatch=g%catchment(ijk)
            if (jcatch.ne.-icatch .and. jcatch.gt.0 .and. ij.ne.ijk) then
              if (t%h(ijk).lt. hsill(-icatch)) then
              hsill(-icatch)=t%h(ijk)
              g%receiver(-icatch)=ijk
              catch(-icatch)=jcatch
              t%length(-icatch)=-1.d0
              endif
            endif
          enddo
        enddo
      endif
    enddo
!$OMP end do nowait
!$OMP end parallel
j=p%ny
!$OMP parallel shared(g,p,t,j,hsill)
!$OMP do private(i,ij,icatch,jcatch,jj,ii,ijk)
    do i=2,p%nx-1
    ij=i+(j-1)*p%nx
    icatch=g%catchment(ij)
      if (icatch.lt.0) then
        do jj=-1,1
          do ii=-1,1
          ijk=ij+ii+jj*p%nx
          if (jj.eq.1) ijk=ij+ii+(jj-p%ny)*p%nx
          jcatch=g%catchment(ijk)
            if (jcatch.ne.-icatch .and. jcatch.gt.0 .and. ij.ne.ijk) then
              if (t%h(ijk).lt. hsill(-icatch)) then
              hsill(-icatch)=t%h(ijk)
              g%receiver(-icatch)=ijk
              catch(-icatch)=jcatch
              t%length(-icatch)=-1.d0
              endif
            endif
          enddo
        enddo
      endif
    enddo
!$OMP end do nowait
!$OMP end parallel

else

j=1
!$OMP parallel shared(g,p,t,j,hsill)
!$OMP do private(i,ij,icatch,jcatch,jj,ii,ijk)
    do i=2,p%nx-1
    ij=i+(j-1)*p%nx
    icatch=g%catchment(ij)
      if (icatch.lt.0) then
        do jj=0,1
          do ii=-1,1
          ijk=ij+ii+jj*p%nx
          jcatch=g%catchment(ijk)
            if (jcatch.ne.-icatch .and. jcatch.gt.0 .and. ij.ne.ijk) then
              if (t%h(ijk).lt. hsill(-icatch)) then
              hsill(-icatch)=t%h(ijk)
              g%receiver(-icatch)=ijk
              catch(-icatch)=jcatch
              t%length(-icatch)=-1.d0
              endif
            endif
          enddo
        enddo
      endif
    enddo
!$OMP end do nowait
!$OMP end parallel
j=p%ny
!$OMP parallel shared(g,p,t,j,hsill)
!$OMP do private(i,ij,icatch,jcatch,jj,ii,ijk)
    do i=2,p%nx-1
    ij=i+(j-1)*p%nx
    icatch=g%catchment(ij)
      if (icatch.lt.0) then
        do jj=-1,0
          do ii=-1,1
          ijk=ij+ii+jj*p%nx
          jcatch=g%catchment(ijk)
            if (jcatch.ne.-icatch .and. jcatch.gt.0 .and. ij.ne.ijk) then
              if (t%h(ijk).lt. hsill(-icatch)) then
              hsill(-icatch)=t%h(ijk)
              g%receiver(-icatch)=ijk
              catch(-icatch)=jcatch
              t%length(-icatch)=-1.d0
              endif
            endif
          enddo
        enddo
      endif
    enddo
!$OMP end do nowait
!$OMP end parallel

endif

if (xcycle) then

!$OMP parallel shared(g,p,t,j,hsill)
!$OMP do private(i,ij,icatch,jcatch,jj,ii,ijk)
  do j=2,p%ny-1
  i=1
  ij=i+(j-1)*p%nx
  icatch=g%catchment(ij)
    if (icatch.lt.0) then
      do jj=-1,1
        do ii=-1,1
        ijk=ij+ii+jj*p%nx
        if (ii.eq.-1) ijk=ij+ii+p%nx+jj*p%nx
        jcatch=g%catchment(ijk)
          if (jcatch.ne.-icatch .and. jcatch.gt.0 .and. ij.ne.ijk) then
            if (t%h(ijk).lt. hsill(-icatch)) then
            hsill(-icatch)=t%h(ijk)
            g%receiver(-icatch)=ijk
            catch(-icatch)=jcatch
            t%length(-icatch)=-1.d0
            endif
          endif
        enddo
      enddo
    endif
  enddo
!$OMP end do nowait
!$OMP end parallel
!$OMP parallel shared(g,p,t,j,hsill)
!$OMP do private(i,ij,icatch,jcatch,jj,ii,ijk)
  do j=2,p%ny-1
  i=p%nx
  ij=i+(j-1)*p%nx
  icatch=g%catchment(ij)
    if (icatch.lt.0) then
      do jj=-1,1
        do ii=-1,1
        ijk=ij+ii+jj*p%nx
        if (ii.eq.1) ijk=ij+ii-p%nx+jj*p%nx
        jcatch=g%catchment(ijk)
          if (jcatch.ne.-icatch .and. jcatch.gt.0 .and. ij.ne.ijk) then
            if (t%h(ijk).lt. hsill(-icatch)) then
            hsill(-icatch)=t%h(ijk)
            g%receiver(-icatch)=ijk
            catch(-icatch)=jcatch
            t%length(-icatch)=-1.d0
            endif
          endif
        enddo
      enddo
    endif
  enddo
!$OMP end do nowait
!$OMP end parallel

else

!$OMP parallel shared(g,p,t,j,hsill)
!$OMP do private(i,ij,icatch,jcatch,jj,ii,ijk)
  do j=2,p%ny-1
  i=1
  ij=i+(j-1)*p%nx
  icatch=g%catchment(ij)
    if (icatch.lt.0) then
      do jj=-1,1
        do ii=0,1
        ijk=ij+ii+jj*p%nx
        jcatch=g%catchment(ijk)
          if (jcatch.ne.-icatch .and. jcatch.gt.0 .and. ij.ne.ijk) then
            if (t%h(ijk).lt. hsill(-icatch)) then
            hsill(-icatch)=t%h(ijk)
            g%receiver(-icatch)=ijk
            catch(-icatch)=jcatch
            t%length(-icatch)=-1.d0
            endif
          endif
        enddo
      enddo
    endif
  enddo
!$OMP end do nowait
!$OMP end parallel
!$OMP parallel shared(g,p,t,j,hsill)
!$OMP do private(i,ij,icatch,jcatch,jj,ii,ijk)
  do j=2,p%ny-1
  i=p%nx
  ij=i+(j-1)*p%nx
  icatch=g%catchment(ij)
    if (icatch.lt.0) then
      do jj=-1,1
        do ii=-1,0
        ijk=ij+ii+jj*p%nx
        jcatch=g%catchment(ijk)
          if (jcatch.ne.-icatch .and. jcatch.gt.0 .and. ij.ne.ijk) then
            if (t%h(ijk).lt. hsill(-icatch)) then
            hsill(-icatch)=t%h(ijk)
            g%receiver(-icatch)=ijk
            catch(-icatch)=jcatch
            t%length(-icatch)=-1.d0
            endif
          endif
        enddo
      enddo
    endif
  enddo
!$OMP end do nowait
!$OMP end parallel

endif

if (ycycle) then

i=1
j=1
  ij=i+(j-1)*p%nx
  icatch=g%catchment(ij)
    if (icatch.lt.0) then
      do jj=-1,1
        do ii=0,1
        ijk=ij+ii+jj*p%nx
        if (jj.eq.-1) ijk=ij+ii+(jj+p%ny)*p%nx
        jcatch=g%catchment(ijk)
          if (jcatch.ne.-icatch .and. jcatch.gt.0 .and. ij.ne.ijk) then
            if (t%h(ijk).lt. hsill(-icatch)) then
            hsill(-icatch)=t%h(ijk)
            g%receiver(-icatch)=ijk
            catch(-icatch)=jcatch
            t%length(-icatch)=-1.d0
            endif
          endif
        enddo
      enddo
    endif

i=p%nx
j=1
  ij=i+(j-1)*p%nx
  icatch=g%catchment(ij)
    if (icatch.lt.0) then
      do jj=-1,1
        do ii=-1,0
        ijk=ij+ii+jj*p%nx
        if (jj.eq.-1) ijk=ij+ii+(jj+p%ny)*p%nx
        jcatch=g%catchment(ijk)
          if (jcatch.ne.-icatch .and. jcatch.gt.0 .and. ij.ne.ijk) then
            if (t%h(ijk).lt. hsill(-icatch)) then
            hsill(-icatch)=t%h(ijk)
            g%receiver(-icatch)=ijk
            catch(-icatch)=jcatch
            t%length(-icatch)=-1.d0
            endif
          endif
        enddo
      enddo
    endif

i=1
j=p%ny
  ij=i+(j-1)*p%nx
  icatch=g%catchment(ij)
    if (icatch.lt.0) then
      do jj=-1,1
        do ii=0,1
        ijk=ij+ii+jj*p%nx
        if (jj.eq.1) ijk=ij+ii+(jj-p%ny)*p%nx
        jcatch=g%catchment(ijk)
          if (jcatch.ne.-icatch .and. jcatch.gt.0 .and. ij.ne.ijk) then
            if (t%h(ijk).lt. hsill(-icatch)) then
            hsill(-icatch)=t%h(ijk)
            g%receiver(-icatch)=ijk
            catch(-icatch)=jcatch
            t%length(-icatch)=-1.d0
            endif
          endif
        enddo
      enddo
    endif

i=p%nx
j=p%ny
  ij=i+(j-1)*p%nx
  icatch=g%catchment(ij)
    if (icatch.lt.0) then
      do jj=-1,1
        do ii=-1,0
        ijk=ij+ii+jj*p%nx
        if (jj.eq.1) ijk=ij+ii+(jj-p%ny)*p%nx
        jcatch=g%catchment(ijk)
          if (jcatch.ne.-icatch .and. jcatch.gt.0 .and. ij.ne.ijk) then
            if (t%h(ijk).lt. hsill(-icatch)) then
            hsill(-icatch)=t%h(ijk)
            g%receiver(-icatch)=ijk
            catch(-icatch)=jcatch
            t%length(-icatch)=-1.d0
            endif
          endif
        enddo
      enddo
    endif

elseif (xcycle) then

i=1
j=1
  ij=i+(j-1)*p%nx
  icatch=g%catchment(ij)
    if (icatch.lt.0) then
      do jj=0,1
        do ii=-1,1
        ijk=ij+ii+jj*p%nx
        if (ii.eq.-1) ijk=ij+ii+p%nx+jj*p%ny
        jcatch=g%catchment(ijk)
          if (jcatch.ne.-icatch .and. jcatch.gt.0 .and. ij.ne.ijk) then
            if (t%h(ijk).lt. hsill(-icatch)) then
            hsill(-icatch)=t%h(ijk)
            g%receiver(-icatch)=ijk
            catch(-icatch)=jcatch
            t%length(-icatch)=-1.d0
            endif
          endif
        enddo
      enddo
    endif

i=p%nx
j=1
  ij=i+(j-1)*p%nx
  icatch=g%catchment(ij)
    if (icatch.lt.0) then
      do jj=0,1
        do ii=-1,1
        ijk=ij+ii+jj*p%nx
        if (ii.eq.1) ijk=ij+ii-p%nx+jj*p%ny
        jcatch=g%catchment(ijk)
          if (jcatch.ne.-icatch .and. jcatch.gt.0 .and. ij.ne.ijk) then
            if (t%h(ijk).lt. hsill(-icatch)) then
            hsill(-icatch)=t%h(ijk)
            g%receiver(-icatch)=ijk
            catch(-icatch)=jcatch
            t%length(-icatch)=-1.d0
            endif
          endif
        enddo
      enddo
    endif

i=1
j=p%ny
  ij=i+(j-1)*p%nx
  icatch=g%catchment(ij)
    if (icatch.lt.0) then
      do jj=-1,0
        do ii=-1,1
        ijk=ij+ii+jj*p%nx
        if (ii.eq.-1) ijk=ij+ii+p%nx+jj*p%ny
        jcatch=g%catchment(ijk)
          if (jcatch.ne.-icatch .and. jcatch.gt.0 .and. ij.ne.ijk) then
            if (t%h(ijk).lt. hsill(-icatch)) then
            hsill(-icatch)=t%h(ijk)
            g%receiver(-icatch)=ijk
            catch(-icatch)=jcatch
            t%length(-icatch)=-1.d0
            endif
          endif
        enddo
      enddo
    endif

i=p%nx
j=p%ny
  ij=i+(j-1)*p%nx
  icatch=g%catchment(ij)
    if (icatch.lt.0) then
      do jj=-1,0
        do ii=-1,1
        ijk=ij+ii+jj*p%nx
        if (ii.eq.1) ijk=ij+ii-p%nx+jj*p%ny
        jcatch=g%catchment(ijk)
          if (jcatch.ne.-icatch .and. jcatch.gt.0 .and. ij.ne.ijk) then
            if (t%h(ijk).lt. hsill(-icatch)) then
            hsill(-icatch)=t%h(ijk)
            g%receiver(-icatch)=ijk
            catch(-icatch)=jcatch
            t%length(-icatch)=-1.d0
            endif
          endif
        enddo
      enddo
    endif

else

i=1
j=1
  ij=i+(j-1)*p%nx
  icatch=g%catchment(ij)
    if (icatch.lt.0) then
      do jj=0,1
        do ii=0,1
        ijk=ij+ii+jj*p%nx
        jcatch=g%catchment(ijk)
          if (jcatch.ne.-icatch .and. jcatch.gt.0 .and. ij.ne.ijk) then
            if (t%h(ijk).lt. hsill(-icatch)) then
            hsill(-icatch)=t%h(ijk)
            g%receiver(-icatch)=ijk
            catch(-icatch)=jcatch
            t%length(-icatch)=-1.d0
            endif
          endif
        enddo
      enddo
    endif

i=p%nx
j=1
  ij=i+(j-1)*p%nx
  icatch=g%catchment(ij)
    if (icatch.lt.0) then
      do jj=0,1
        do ii=-1,0
        ijk=ij+ii+jj*p%nx
        jcatch=g%catchment(ijk)
          if (jcatch.ne.-icatch .and. jcatch.gt.0 .and. ij.ne.ijk) then
            if (t%h(ijk).lt. hsill(-icatch)) then
            hsill(-icatch)=t%h(ijk)
            g%receiver(-icatch)=ijk
            catch(-icatch)=jcatch
            t%length(-icatch)=-1.d0
            endif
          endif
        enddo
      enddo
    endif

i=1
j=p%ny
  ij=i+(j-1)*p%nx
  icatch=g%catchment(ij)
    if (icatch.lt.0) then
      do jj=-1,0
        do ii=0,1
        ijk=ij+ii+jj*p%nx
        jcatch=g%catchment(ijk)
          if (jcatch.ne.-icatch .and. jcatch.gt.0 .and. ij.ne.ijk) then
            if (t%h(ijk).lt. hsill(-icatch)) then
            hsill(-icatch)=t%h(ijk)
            g%receiver(-icatch)=ijk
            catch(-icatch)=jcatch
            t%length(-icatch)=-1.d0
            endif
          endif
        enddo
      enddo
    endif

i=p%nx
j=p%ny
  ij=i+(j-1)*p%nx
  icatch=g%catchment(ij)
    if (icatch.lt.0) then
      do jj=-1,0
        do ii=-1,0
        ijk=ij+ii+jj*p%nx
        jcatch=g%catchment(ijk)
          if (jcatch.ne.-icatch .and. jcatch.gt.0 .and. ij.ne.ijk) then
            if (t%h(ijk).lt. hsill(-icatch)) then
            hsill(-icatch)=t%h(ijk)
            g%receiver(-icatch)=ijk
            catch(-icatch)=jcatch
            t%length(-icatch)=-1.d0
            endif
          endif
        enddo
      enddo
    endif

endif

1111 continue

!$OMP parallel shared(g,p,catch)
!$OMP do private(ij)
  do ij=1,p%nn
  if (catch(ij).gt.0.d0) g%catchment(ij)=catch(ij)
  enddo
!!!where (catch.gt.0) g%catchment=catch
!$OMP end do nowait
!$OMP end parallel

  do ij=1,p%nn
  icatch=g%catchment(ij)
  if (icatch.lt.0) then
    if (g%catchment(-icatch).gt.0) then
    g%catchment(ij)=g%catchment(-icatch)
    endif
  endif
  enddo

!$OMP parallel shared(nc,g)
!$OMP workshare
nc=count(g%catchment.lt.0)
!$OMP end workshare nowait
!$OMP end parallel

niter=niter+1
if (p%debug) print*,'Iteration ',niter,' Number of local minima :',nc
if (nc.gt.0) goto 111

deallocate (hsill,catch)

return

end subroutine local_minima
