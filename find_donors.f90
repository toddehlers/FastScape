subroutine find_donors (g,p)

! subroutine to find donor list from receiver list
! needs to create a working array to avoid any operations over O(N)

use omp_lib
use definitions

implicit none

type (topology) g
type (param) p

integer, dimension(:),allocatable::nwork
integer ij,ijk

!$OMP parallel shared(g)
!$OMP workshare
g%ndon=0
!$OMP end workshare nowait
!$OMP end parallel

!$OMP parallel shared(p,g)
!$OMP do private (ij,ijk) schedule(dynamic,p%chunk)
  do ij=1,p%nn
  ijk=g%receiver(ij)
  if (ijk.ne.ij) g%ndon(ijk)=g%ndon(ijk)+1
  enddo
!$OMP end do nowait
!$OMP end parallel

g%ndon(p%nn+1)=p%nn+1

  do ij=p%nn,1,-1
  g%ndon(ij)=g%ndon(ij+1)-g%ndon(ij)
  enddo

allocate (nwork(p%nn))
!$OMP parallel shared(nwork)
!$OMP workshare
nwork=0
!$OMP end workshare nowait
!$OMP end parallel

!!$OMP parallel shared(p,nwork,g)
!!$OMP do private (ij,ijk) schedule(dynamic,p%chunk)
  do ij=1,p%nn
  ijk=g%receiver(ij)
    if (ijk.ne.ij) then
    g%donors(g%ndon(ijk)+nwork(ijk))=ij
    nwork(ijk)=nwork(ijk)+1
    endif
  enddo
!!$OMP end do nowait
!!$OMP end parallel

deallocate (nwork)

return

end subroutine find_donors
