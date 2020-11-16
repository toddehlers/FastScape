subroutine Concavity (t,p,g,l,x)

! to compute concavity (exponent in Area-Slope relationship) along river profiles

use omp_lib
use definitions

implicit none

type (topography) t
type (topology) g
type (param) p
type (law) l

double precision x(p%nn)
double precision, dimension(:), allocatable :: area,slope,sumx,sumy,sumxx,sumxy
double precision, dimension(:), allocatable :: n
double precision a,b,con,dzdx,dzdy,loga,logs
integer i,j,ij,ijk,ia,ib,ic,id,ie,if,ig,ih,ii

if (p%debug) print*,'calculating Concavity'

! p%stream_threshold_area is the minimum drainage area for a pixel to be considered as a stream

allocate (area(p%nn),slope(p%nn),sumx(p%nn),sumy(p%nn),sumxx(p%nn),sumxy(p%nn),n(p%nn))

!$OMP parallel shared(p,area)
!$OMP workshare
area=p%dx*p%dy
!$OMP end workshare nowait
!$OMP end parallel

! first calculate the drainage area of each point

!$OMP parallel shared(t,g,p)
!$OMP do private(i,j,ij,ijk)
  do j=1,p%num_threads
    do i=g%nstackp,1,-1
    ij=g%stack((j-1)*g%nstackp+i)
      if (ij.ne.0) then
      ijk=g%receiver(ij)
      if (ijk.ne.ij) area(ijk)=area(ijk)+area(ij)
      endif
    enddo
  enddo
!$OMP end do nowait
!$OMP end parallel

! then the slope

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
    if (slope(ij).gt.tiny(x(ij))) slope(ij)=sqrt(slope(ij))
!    if (slope(ij).gt.tiny(x(ij))) slope(ij)=atan(sqrt(slope(ij)))
    enddo
  enddo

! now the regression coefficients

sumx=0.d0
sumy=0.d0
sumxx=0.d0
sumxy=0.d0
n=0.d0
  do j=1,p%num_threads
    do i=g%nstackp,1,-1
    ij=g%stack((j-1)*g%nstackp+i)
      if (ij.ne.0) then
      ijk=g%receiver(ij)
        if (ijk.ne.ij) then
          if (area(ij).gt.tiny(area(ij)).and.slope(ij).gt.tiny(slope(ij))) then
          loga=log10(area(ij))
          logs=log10(slope(ij))
          sumx(ij)=sumx(ij)+loga
          sumy(ij)=sumy(ij)+logs
          sumxx(ij)=sumxx(ij)+loga*loga
          sumxy(ij)=sumxy(ij)+loga*logs
          n(ij)=n(ij)+1.d0
          sumx(ijk)=sumx(ijk)+sumx(ij)
          sumy(ijk)=sumy(ijk)+sumy(ij)
          sumxx(ijk)=sumxx(ijk)+sumxx(ij)
          sumxy(ijk)=sumxy(ijk)+sumxy(ij)
          n(ijk)=n(ijk)+n(ij)
          endif
        endif
      endif
    enddo
  enddo

x=0.d0
  do ij=1,p%nn
  if (n(ij).gt.100.and.n(ij)*sumxx(ij)-sumx(ij)**2.gt.0.d0.and.area(ij).gt.p%stream_threshold_area) &
     x(ij)=-(n(ij)*sumxy(ij)-sumx(ij)*sumy(ij))/(n(ij)*sumxx(ij)-sumx(ij)**2)
  enddo

deallocate (area,slope,sumx,sumy,sumxx,sumxy,n)

if (p%debug) print*,'Concavity calculated'

return

end subroutine Concavity
