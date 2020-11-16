subroutine set_bc (g,p,t)

! subroutine to set the boundary conditions
! it reads information from the ibc flag and transforms it into a bc logical array that
! is used later in the program to define base level
! the user could potentially change this routine to add more general boundary conditions

use omp_lib
use definitions

implicit none

type (topology) g
type (param) p
type (topography) t

integer i,j,ij
character c4*4

c4='0000'
write (c4,'(i4)') p%ibc

!$OMP parallel shared(g)
!$OMP workshare
g%bc=.false.
!$OMP end workshare nowait
!$OMP end parallel

!$OMP parallel shared(g,p,c4)
!$OMP do private (i,j,ij) schedule(dynamic,p%ny/p%num_threads)
  do j=1,p%ny
    do i=1,p%nx
    ij=i+(j-1)*p%nx
    if (j.eq.1 .and. c4(1:1).eq.'1') g%bc(ij)=.true.
    if (i.eq.p%nx .and. c4(2:2).eq.'1') g%bc(ij)=.true.
    if (j.eq.p%ny .and. c4(3:3).eq.'1') g%bc(ij)=.true.
    if (i.eq.1 .and. c4(4:4).eq.'1') g%bc(ij)=.true.
    enddo
  enddo
!$OMP end do nowait
!$OMP end parallel

!$OMP parallel shared(g,p,c4)
!$OMP do private (i,j,ij) schedule(dynamic,p%ny/p%num_threads)
  do j=1,p%ny
    do i=1,p%nx
    ij=i+(j-1)*p%nx
    if (i.eq.1 .and. j.eq.1 .and. c4(4:4).eq.'2' .and. c4(1:1).eq.'2') g%bc(ij)=.true.
    if (i.eq.p%nx .and. j.eq.1 .and. c4(1:1).eq.'2' .and. c4(2:2).eq.'2') g%bc(ij)=.true.
    if (i.eq.p%nx .and. j.eq.p%ny .and. c4(2:2).eq.'2' .and. c4(3:3).eq.'2') g%bc(ij)=.true.
    if (i.eq.1 .and. j.eq.p%ny .and. c4(3:3).eq.'2' .and. c4(4:4).eq.'2') g%bc(ij)=.true.
    enddo
  enddo
!$OMP end do nowait
!$OMP end parallel

!$OMP parallel shared(g,p,c4)
!$OMP do private (i,j,ij) schedule(dynamic,p%ny/p%num_threads)
  do j=1,p%ny
    do i=1,p%nx
    ij=i+(j-1)*p%nx
    if (i.eq.p%nx/2+1 .and. j.eq.1 .and. c4(1:1).eq.'3') g%bc(ij)=.true.
    if (i.eq.p%nx .and. j.eq.p%ny/2+1 .and. c4(2:2).eq.'3') g%bc(ij)=.true.
    if (i.eq.p%nx/2+1 .and. j.eq.p%ny .and. c4(3:3).eq.'3') g%bc(ij)=.true.
    if (i.eq.1 .and. j.eq.p%ny/2+1 .and. c4(4:4).eq.'3') g%bc(ij)=.true.
    enddo
  enddo
!$OMP end do nowait
!$OMP end parallel


if (p%restart.lt.0) where (t%h.lt.p%base_level) g%bc=.true.

end subroutine set_bc
