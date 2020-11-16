subroutine Interface_Pecube (t,p,istep)

use omp_lib
use definitions

implicit none

type (param) p
type (topography) t

double precision, dimension(:,:), allocatable :: Pecube,temp

character c5*5
integer nc5,i,j,ij,istep

if (istep.lt.10) then
nc5=1
write (c5(1:nc5),'(i1)') istep
elseif (istep.lt.100) then
nc5=2
write (c5(1:nc5),'(i2)') istep
elseif (istep.lt.1000) then
nc5=3
write (c5(1:nc5),'(i3)') istep
elseif (istep.lt.10000) then
nc5=4
write (c5(1:nc5),'(i4)') istep
elseif (istep.lt.100000) then
nc5=5
write (c5(1:nc5),'(i5)') istep
else
stop 'too many time steps in Interface_Pecube'
endif

allocate (Pecube(p%Pecube_nx,p%Pecube_ny),temp(p%nx,p%ny))

open (7,file=p%run//'/topo'//c5(1:nc5),status='unknown')
  do j=1,p%ny
    do i=1,p%nx
    ij=i+(j-1)*p%nx
    temp(i,j)=t%h(ij)
    enddo
  enddo
call InterpolLinearBox(temp,p%nx,p%ny,Pecube,p%Pecube_nx,p%Pecube_ny, &
                       p%Pecube_xmin/p%xl,p%Pecube_xmax/p%xl,p%Pecube_ymin/p%yl,p%Pecube_ymax/p%yl)
  do j=1,p%Pecube_ny
    do i=1,p%Pecube_nx
    write (7,*) Pecube(i,j)
    enddo
  enddo
close (7)

open (7,file=p%run//'/temp'//c5(1:nc5),status='unknown')
  do i=1,p%Pecube_nx*p%Pecube_ny
  write (7,*) 0.d0
  enddo
close (7)

open (7,file=p%run//'/uplift'//c5(1:nc5),status='unknown')
  do j=1,p%ny
    do i=1,p%nx
    ij=i+(j-1)*p%nx
    temp(i,j)=t%u(ij)
    enddo
  enddo
call InterpolLinearBox(temp,p%nx,p%ny,Pecube,p%Pecube_nx,p%Pecube_ny, &
                       p%Pecube_xmin/p%xl,p%Pecube_xmax/p%xl,p%Pecube_ymin/p%yl,p%Pecube_ymax/p%yl)
  do j=1,p%Pecube_ny
    do i=1,p%Pecube_nx
    write (7,*) Pecube(i,j)*1.d3
    enddo
  enddo
close (7)

deallocate (Pecube)

istep=istep+1

end subroutine Interface_Pecube
