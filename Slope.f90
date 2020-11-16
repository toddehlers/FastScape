subroutine Slope (t,p,g,l,x,mean,var,skew,kurt)

use omp_lib
use definitions

implicit none

type (topography) t
type (topology) g
type (param) p
type (law) l

double precision x(p%nn)
integer i,j,ij,ia,ib,ic,id,ie,if,ig,ih,ii
double precision dzdx,dzdy,con
double precision mean,var,skew,kurt

con=45.d0/atan(1.d0)

x=0.d0
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
    x(ij)=dzdx**2+dzdy**2
    if (x(ij).gt.tiny(x(ij))) x(ij)=atan(sqrt(x(ij)))*con
    enddo
  enddo

mean=sum(x)/(p%nx-2)/(p%ny-2)
var=sqrt(sum((x-mean)**2)/(p%nx-2)/(p%ny-2))
skew=sum((x-mean)*(x-mean)*(x-mean))/(p%nx-2)/(p%ny-2)/var**3
kurt=sum((x-mean)**4)/(p%nx-2)/(p%ny-2)/var**4

end subroutine Slope
