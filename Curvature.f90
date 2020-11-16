subroutine Curvature (t,p,g,l,x)

use omp_lib
use definitions

implicit none

type (topography) t
type (topology) g
type (param) p
type (law) l

double precision x(p%nn)
integer i,j,ij,i1,i2,i3,i4,i5,i6,i7,i8,i9
double precision con,a,b,c,d,e,f

con=90.d0/atan(1.d0)

x=0.d0
  do j=2,p%ny-1
    do i=2,p%nx-1
    ij=i+(j-1)*p%nx
    i1=ij+p%nx-1
    i2=i1+1
    i3=i2+1
    i4=ij-1
    i5=ij
    i6=ij+1
    i7=ij-p%nx-1
    i8=i7+1
    i9=i8+1
    a=(t%h(i1)+t%h(i3)+t%h(i4)+t%h(i6)+t%h(i7)+t%h(i9))/p%dx/p%dx/12.d0-(t%h(i2)+t%h(i5)+t%h(i8))/p%dx/p%dx/6.d0
    b=(t%h(i1)+t%h(i2)+t%h(i3)+t%h(i7)+t%h(i8)+t%h(i9))/p%dy/p%dy/12.d0-(t%h(i4)+t%h(i5)+t%h(i6))/p%dy/p%dy/6.d0
    c=(t%h(i3)+t%h(i7)-t%h(i1)-t%h(i9))/p%dx/p%dy/4.d0
    d=(t%h(i3)+t%h(i6)+t%h(i9)-t%h(i1)-t%h(i4)-t%h(i7))/p%dx/6.d0
    e=(t%h(i1)+t%h(i2)+t%h(i3)-t%h(i7)-t%h(i8)-t%h(i9))/p%dy/6.d0
    f=(2.d0*(t%h(i2)+t%h(i4)+t%h(i6)+t%h(i8))-(t%h(i1)+t%h(i3)+t%h(i7)+t%h(i9))+5.d0*t%h(i5))/9.d0
    x(ij)=1.d0+d**2+e**2
    if (x(ij).gt.tiny(x(ij))) x(ij)=(a*(1.d0+e**2)+b*(1.d0+d**2)-c*d*e)/(x(ij)**(3.d0/2.d0))
    enddo
  enddo

end subroutine Curvature
