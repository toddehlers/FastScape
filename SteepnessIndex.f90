subroutine SteepnessIndex (t,p,g,l,x)

use omp_lib
use definitions

implicit none

type (topography) t
type (topology) g
type (param) p
type (law) l

double precision x(p%nn)
double precision, dimension(:), allocatable :: slope
integer i,j,ij,ia,ib,ic,id,ie,if,ig,ih,ii
double precision dzdx,dzdy,con

allocate (slope(p%nn))

con=90.d0/atan(1.d0)

slope=0.d0
  do j=2,p%ny-1
    do i=2,p%nx-1
    ij=i+(j-1)*p%nx
    ia=ij+p%nx-1
    ib=ia+1
    ic=ib+1
    id=ij-1
    ie=ij
    if=ij+1
    ig=ij-p%nx-1
    ih=ig+1
    ii=ih+1
    dzdx=((t%h(ic)+2.d0*t%h(if)+t%h(ii))-(t%h(ia)+2.d0*t%h(id)+t%h(ig)))/8.d0/p%dx
    dzdy=((t%h(ig)+2.d0*t%h(ih)+t%h(ii))-(t%h(ia)+2.d0*t%h(ib)+t%h(ic)))/8.d0/p%dy
    slope(ij)=dzdx**2+dzdy**2
    if (slope(ij).gt.tiny(slope(ij))) slope(ij)=atan(sqrt(slope(ij)))*con
    enddo
  enddo

  do ij=1,p%nn
  x(ij)=slope(ij)*t%discharge(ij)**(l%m/l%n)
  enddo

deallocate (slope)

end subroutine SteepnessIndex
