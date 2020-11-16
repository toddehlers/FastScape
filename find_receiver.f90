subroutine find_receiver (g,p,t)

! subroutine to find the receiver (steepest slope) node to each node
! computes also the length between each node and its receiver
! beware that I have used a very tedious, yet efficient way to treat the boundary conditions
! that can be of several types around the domain boundaries (reflective, periodic or base level)

use omp_lib
use definitions

implicit none

type (topography) t
type (topology) g
type (param) p

integer i,j,ij,ii,jj,ijk
double precision slopemax,slopep,lengthp
logical xcycle,ycycle
character c4*4

xcycle=.false.
ycycle=.false.
write (c4,'(i4)') p%ibc
if (c4(1:1).eq.'0' .and. c4(3:3).eq.'0') ycycle=.true.
if (c4(2:2).eq.'0' .and. c4(4:4).eq.'0') xcycle=.true.

!$OMP parallel num_threads(p%num_threads) shared(g,p,t) private(i,j,ij,slopemax,jj,ii,ijk,slopep,lengthp)
!$OMP do schedule(dynamic,p%ny/p%num_threads)
  do j=2,p%ny-1
    do i=2,p%nx-1
    ij=i+(j-1)*p%nx
      if (.not.g%bc(ij)) then
      slopemax=0.d0
        do jj=-1,1
          do ii=-1,1
          ijk=ij+ii+jj*p%nx
          slopep=t%h(ij)-t%h(ijk)
            if (slopep.gt.0.d0) then
            lengthp=sqrt((ii*p%dx)**2+(jj*p%dy)**2)
            slopep=slopep/lengthp
              if (slopep.gt.slopemax) then
              slopemax=slopep
              t%length(ij)=lengthp
              g%receiver(ij)=ijk
              endif
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
!$OMP parallel num_threads(p%num_threads) shared(g,p,t,j) private(i,ij,slopemax,jj,ii,ijk,slopep,lengthp)
!$OMP do schedule(dynamic,p%nx/p%num_threads)
    do i=2,p%nx-1
    ij=i+(j-1)*p%nx
      if (.not.g%bc(ij)) then
      slopemax=0.d0
        do jj=-1,1
          do ii=-1,1
          ijk=ij+ii+jj*p%nx
          if (jj.eq.-1) ijk=ij+ii+(jj+p%ny)*p%nx
          slopep=t%h(ij)-t%h(ijk)
            if (slopep.gt.0.d0) then
            lengthp=sqrt((ii*p%dx)**2+(jj*p%dy)**2)
            slopep=slopep/lengthp
              if (slopep.gt.slopemax) then
              slopemax=slopep
              t%length(ij)=lengthp
              g%receiver(ij)=ijk
              endif
            endif
          enddo
        enddo
      endif
    enddo
!$OMP end do nowait
!$OMP end parallel
j=p%ny
!$OMP parallel num_threads(p%num_threads) shared(g,p,t,j) private(i,ij,slopemax,jj,ii,ijk,slopep,lengthp)
!$OMP do schedule(dynamic,p%nx/p%num_threads)
    do i=2,p%nx-1
    ij=i+(j-1)*p%nx
      if (.not.g%bc(ij)) then
      slopemax=0.d0
        do jj=-1,1
          do ii=-1,1
          ijk=ij+ii+jj*p%nx
          if (jj.eq.1) ijk=ij+ii+(jj-p%ny)*p%nx
          slopep=t%h(ij)-t%h(ijk)
            if (slopep.gt.0.d0) then
            lengthp=sqrt((ii*p%dx)**2+(jj*p%dy)**2)
            slopep=slopep/lengthp
              if (slopep.gt.slopemax) then
              slopemax=slopep
              t%length(ij)=lengthp
              g%receiver(ij)=ijk
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
!$OMP parallel shared(g,p,t,j)
!$OMP do private(i,ij,slopemax,jj,ii,ijk,slopep,lengthp) &
!$OMP schedule(dynamic,p%nx/p%num_threads)
    do i=2,p%nx-1
    ij=i+(j-1)*p%nx
      if (.not.g%bc(ij)) then
      slopemax=0.d0
        do jj=0,1
          do ii=-1,1
          ijk=ij+ii+jj*p%nx
          slopep=t%h(ij)-t%h(ijk)
            if (slopep.gt.0.d0) then
            lengthp=sqrt((ii*p%dx)**2+(jj*p%dy)**2)
            slopep=slopep/lengthp
              if (slopep.gt.slopemax) then
              slopemax=slopep
              t%length(ij)=lengthp
              g%receiver(ij)=ijk
              endif
            endif
          enddo
        enddo
      endif
    enddo
!$OMP end do nowait
!$OMP end parallel
j=p%ny
!$OMP parallel shared(g,p,t,j)
!$OMP do private(i,ij,slopemax,jj,ii,ijk,slopep,lengthp) &
!$OMP schedule(dynamic,p%nx/p%num_threads)
    do i=2,p%nx-1
    ij=i+(j-1)*p%nx
      if (.not.g%bc(ij)) then
      slopemax=0.d0
        do jj=-1,0
          do ii=-1,1
          ijk=ij+ii+jj*p%nx
          slopep=t%h(ij)-t%h(ijk)
            if (slopep.gt.0.d0) then
            lengthp=sqrt((ii*p%dx)**2+(jj*p%dy)**2)
            slopep=slopep/lengthp
              if (slopep.gt.slopemax) then
              slopemax=slopep
              t%length(ij)=lengthp
              g%receiver(ij)=ijk
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

!$OMP parallel shared(g,p,t)
!$OMP do private(i,j,ij,slopemax,jj,ii,ijk,slopep,lengthp) &
!$OMP schedule(dynamic,p%ny/p%num_threads)
  do j=2,p%ny-1
  i=1
  ij=i+(j-1)*p%nx
    if (.not.g%bc(ij)) then
    slopemax=0.d0
      do jj=-1,1
        do ii=-1,1
        ijk=ij+ii+jj*p%nx
        if (ii.eq.-1) ijk=ij+ii+p%nx+jj*p%nx
        slopep=t%h(ij)-t%h(ijk)
          if (slopep.gt.0.d0) then
          lengthp=sqrt((ii*p%dx)**2+(jj*p%dy)**2)
          slopep=slopep/lengthp
            if (slopep.gt.slopemax) then
            slopemax=slopep
            t%length(ij)=lengthp
            g%receiver(ij)=ijk
            endif
          endif
        enddo
      enddo
    endif
  enddo
!$OMP end do nowait
!$OMP end parallel
!$OMP parallel shared(g,p,t)
!$OMP do private(i,j,ij,slopemax,jj,ii,ijk,slopep,lengthp) &
!$OMP schedule(dynamic,p%ny/p%num_threads)
  do j=2,p%ny-1
  i=p%nx
  ij=i+(j-1)*p%nx
    if (.not.g%bc(ij)) then
    slopemax=0.d0
      do jj=-1,1
        do ii=-1,1
        ijk=ij+ii+jj*p%nx
        if (ii.eq.1) ijk=ij+ii-p%nx+jj*p%nx
        slopep=t%h(ij)-t%h(ijk)
          if (slopep.gt.0.d0) then
          lengthp=sqrt((ii*p%dx)**2+(jj*p%dy)**2)
          slopep=slopep/lengthp
            if (slopep.gt.slopemax) then
            slopemax=slopep
            t%length(ij)=lengthp
            g%receiver(ij)=ijk
            endif
          endif
        enddo
      enddo
    endif
  enddo
!$OMP end do nowait
!$OMP end parallel

else

!$OMP parallel shared(g,p,t)
!$OMP do private(i,j,ij,slopemax,jj,ii,ijk,slopep,lengthp) &
!$OMP schedule(dynamic,p%ny/p%num_threads)
  do j=2,p%ny-1
  i=1
  ij=i+(j-1)*p%nx
    if (.not.g%bc(ij)) then
    slopemax=0.d0
      do jj=-1,1
        do ii=0,1
        ijk=ij+ii+jj*p%nx
        slopep=t%h(ij)-t%h(ijk)
          if (slopep.gt.0.d0) then
          lengthp=sqrt((ii*p%dx)**2+(jj*p%dy)**2)
          slopep=slopep/lengthp
            if (slopep.gt.slopemax) then
            slopemax=slopep
            t%length(ij)=lengthp
            g%receiver(ij)=ijk
            endif
          endif
        enddo
      enddo
    endif
  enddo
!$OMP end do nowait
!$OMP end parallel
!$OMP parallel shared(g,p,t)
!$OMP do private(i,j,ij,slopemax,jj,ii,ijk,slopep,lengthp) &
!$OMP schedule(dynamic,p%ny/p%num_threads)
  do j=2,p%ny-1
  i=p%nx
  ij=i+(j-1)*p%nx
    if (.not.g%bc(ij)) then
    slopemax=0.d0
      do jj=-1,1
        do ii=-1,0
        ijk=ij+ii+jj*p%nx
        slopep=t%h(ij)-t%h(ijk)
          if (slopep.gt.0.d0) then
          lengthp=sqrt((ii*p%dx)**2+(jj*p%dy)**2)
          slopep=slopep/lengthp
            if (slopep.gt.slopemax) then
            slopemax=slopep
            t%length(ij)=lengthp
            g%receiver(ij)=ijk
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
    if (.not.g%bc(ij)) then
    slopemax=0.d0
      do jj=-1,1
        do ii=0,1
        ijk=ij+ii+jj*p%nx
        if (jj.eq.-1) ijk=ij+ii+(jj+p%ny)*p%nx
        slopep=t%h(ij)-t%h(ijk)
          if (slopep.gt.0.d0) then
          lengthp=sqrt((ii*p%dx)**2+(jj*p%dy)**2)
          slopep=slopep/lengthp
            if (slopep.gt.slopemax) then
            slopemax=slopep
            t%length(ij)=lengthp
            g%receiver(ij)=ijk
            endif
          endif
        enddo
      enddo
    endif

i=p%nx
j=1
  ij=i+(j-1)*p%nx
    if (.not.g%bc(ij)) then
    slopemax=0.d0
      do jj=-1,1
        do ii=-1,0
        ijk=ij+ii+jj*p%nx
        if (jj.eq.-1) ijk=ij+ii+(jj+p%ny)*p%nx
        slopep=t%h(ij)-t%h(ijk)
          if (slopep.gt.0.d0) then
          lengthp=sqrt((ii*p%dx)**2+(jj*p%dy)**2)
          slopep=slopep/lengthp
            if (slopep.gt.slopemax) then
            slopemax=slopep
            t%length(ij)=lengthp
            g%receiver(ij)=ijk
            endif
          endif
        enddo
      enddo
    endif

i=1
j=p%ny
  ij=i+(j-1)*p%nx
    if (.not.g%bc(ij)) then
    slopemax=0.d0
      do jj=-1,1
        do ii=0,1
        ijk=ij+ii+jj*p%nx
        if (jj.eq.1) ijk=ij+ii+(jj-p%ny)*p%nx
        slopep=t%h(ij)-t%h(ijk)
          if (slopep.gt.0.d0) then
          lengthp=sqrt((ii*p%dx)**2+(jj*p%dy)**2)
          slopep=slopep/lengthp
            if (slopep.gt.slopemax) then
            slopemax=slopep
            t%length(ij)=lengthp
            g%receiver(ij)=ijk
            endif
          endif
        enddo
      enddo
    endif

i=p%nx
j=p%ny
  ij=i+(j-1)*p%nx
    if (.not.g%bc(ij)) then
    slopemax=0.d0
      do jj=-1,1
        do ii=-1,0
        ijk=ij+ii+jj*p%nx
        if (jj.eq.1) ijk=ij+ii+(jj-p%ny)*p%nx
        slopep=t%h(ij)-t%h(ijk)
          if (slopep.gt.0.d0) then
          lengthp=sqrt((ii*p%dx)**2+(jj*p%dy)**2)
          slopep=slopep/lengthp
            if (slopep.gt.slopemax) then
            slopemax=slopep
            t%length(ij)=lengthp
            g%receiver(ij)=ijk
            endif
          endif
        enddo
      enddo
    endif

elseif (xcycle) then

i=1
j=1
  ij=i+(j-1)*p%nx
    if (.not.g%bc(ij)) then
    slopemax=0.d0
      do jj=0,1
        do ii=-1,1
        ijk=ij+ii+jj*p%nx
        if (ii.eq.-1) ijk=ij+ii+p%nx+jj*p%ny
        slopep=t%h(ij)-t%h(ijk)
          if (slopep.gt.0.d0) then
          lengthp=sqrt((ii*p%dx)**2+(jj*p%dy)**2)
          slopep=slopep/lengthp
            if (slopep.gt.slopemax) then
            slopemax=slopep
            t%length(ij)=lengthp
            g%receiver(ij)=ijk
            endif
          endif
        enddo
      enddo
    endif

i=p%nx
j=1
  ij=i+(j-1)*p%nx
    if (.not.g%bc(ij)) then
    slopemax=0.d0
      do jj=0,1
        do ii=-1,1
        ijk=ij+ii+jj*p%nx
        if (ii.eq.1) ijk=ij+ii-p%nx+jj*p%ny
        slopep=t%h(ij)-t%h(ijk)
          if (slopep.gt.0.d0) then
          lengthp=sqrt((ii*p%dx)**2+(jj*p%dy)**2)
          slopep=slopep/lengthp
            if (slopep.gt.slopemax) then
            slopemax=slopep
            t%length(ij)=lengthp
            g%receiver(ij)=ijk
            endif
          endif
        enddo
      enddo
    endif

i=1
j=p%ny
  ij=i+(j-1)*p%nx
    if (.not.g%bc(ij)) then
    slopemax=0.d0
      do jj=-1,0
        do ii=-1,1
        ijk=ij+ii+jj*p%nx
        if (ii.eq.-1) ijk=ij+ii+p%nx+jj*p%ny
        slopep=t%h(ij)-t%h(ijk)
          if (slopep.gt.0.d0) then
          lengthp=sqrt((ii*p%dx)**2+(jj*p%dy)**2)
          slopep=slopep/lengthp
            if (slopep.gt.slopemax) then
            slopemax=slopep
            t%length(ij)=lengthp
            g%receiver(ij)=ijk
            endif
          endif
        enddo
      enddo
    endif

i=p%nx
j=p%ny
  ij=i+(j-1)*p%nx
    if (.not.g%bc(ij)) then
    slopemax=0.d0
      do jj=-1,0
        do ii=-1,1
        ijk=ij+ii+jj*p%nx
        if (ii.eq.1) ijk=ij+ii-p%nx+jj*p%ny
        slopep=t%h(ij)-t%h(ijk)
          if (slopep.gt.0.d0) then
          lengthp=sqrt((ii*p%dx)**2+(jj*p%dy)**2)
          slopep=slopep/lengthp
            if (slopep.gt.slopemax) then
            slopemax=slopep
            t%length(ij)=lengthp
            g%receiver(ij)=ijk
            endif
          endif
        enddo
      enddo
    endif

else

i=1
j=1
  ij=i+(j-1)*p%nx
    if (.not.g%bc(ij)) then
    slopemax=0.d0
      do jj=0,1
        do ii=0,1
        ijk=ij+ii+jj*p%nx
        slopep=t%h(ij)-t%h(ijk)
          if (slopep.gt.0.d0) then
          lengthp=sqrt((ii*p%dx)**2+(jj*p%dy)**2)
          slopep=slopep/lengthp
            if (slopep.gt.slopemax) then
            slopemax=slopep
            t%length(ij)=lengthp
            g%receiver(ij)=ijk
            endif
          endif
        enddo
      enddo
    endif

i=p%nx
j=1
  ij=i+(j-1)*p%nx
    if (.not.g%bc(ij)) then
    slopemax=0.d0
      do jj=0,1
        do ii=-1,0
        ijk=ij+ii+jj*p%nx
        slopep=t%h(ij)-t%h(ijk)
          if (slopep.gt.0.d0) then
          lengthp=sqrt((ii*p%dx)**2+(jj*p%dy)**2)
          slopep=slopep/lengthp
            if (slopep.gt.slopemax) then
            slopemax=slopep
            t%length(ij)=lengthp
            g%receiver(ij)=ijk
            endif
          endif
        enddo
      enddo
    endif

i=1
j=p%ny
  ij=i+(j-1)*p%nx
    if (.not.g%bc(ij)) then
    slopemax=0.d0
      do jj=-1,0
        do ii=0,1
        ijk=ij+ii+jj*p%nx
        slopep=t%h(ij)-t%h(ijk)
          if (slopep.gt.0.d0) then
          lengthp=sqrt((ii*p%dx)**2+(jj*p%dy)**2)
          slopep=slopep/lengthp
            if (slopep.gt.slopemax) then
            slopemax=slopep
            t%length(ij)=lengthp
            g%receiver(ij)=ijk
            endif
          endif
        enddo
      enddo
    endif

i=p%nx
j=p%ny
  ij=i+(j-1)*p%nx
    if (.not.g%bc(ij)) then
    slopemax=0.d0
      do jj=-1,0
        do ii=-1,0
        ijk=ij+ii+jj*p%nx
        slopep=t%h(ij)-t%h(ijk)
          if (slopep.gt.0.d0) then
          lengthp=sqrt((ii*p%dx)**2+(jj*p%dy)**2)
          slopep=slopep/lengthp
            if (slopep.gt.slopemax) then
            slopemax=slopep
            t%length(ij)=lengthp
            g%receiver(ij)=ijk
            endif
          endif
        enddo
      enddo
    endif

endif

return

end subroutine find_receiver
