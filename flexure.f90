subroutine flexure (p,t,g)

! first half of the flexure.f90 routine that can be complemented by adding lines
! in the input file

! This routine computes the flexural response to the load caused by topography changes

use omp_lib
use definitions

implicit none

type (param) p
type (topography) t
type (topology) g

real*8 hx,drho,elt,dflex,d,xk,xloc,yloc,hh,hh0,h1,h2,h3,h4,r,s
real*8, dimension(:,:), allocatable :: flex,work,fl
real*8 xmean,xfmin,xfmax,dxf,ddxf,ymean,yfmin,yfmax,dyf,ddyf,fmean
real*8 pihx,pihy,pi,fi,fj,tij,xp,yp,hisotot,xflex,yflex,rx,ry,dh
real*8, dimension(:), allocatable :: junk

character c4*4

integer nflex,i,j,k,ii,jj,iref,iflexmin,iflexmax,jflexmin,jflexmax
integer ixflex,iyflex,iflexfirst,iflexlast,jflexfirst,jflexlast
integer ij,iflex,jflex
integer, dimension(:,:), allocatable :: nfl

real*8 x,y,xl,yl,time,h,dx,dy,dt
real*8 w1,w2,w3,w4,w5,w6,w7,w8,w9,w0
integer i1,i2,i3,i4,i5,i6,i7,i8,i9,i0,nx,ny
real*8 drhoe,dist

nflex=512

iref=1

write (c4,'(i4)') p%ibc
if (p%ibc.lt.10) c4(1:3)='000'
if (p%ibc.lt.100) c4(1:2)='00'
if (p%ibc.lt.1000) c4(1:1)='0'

ixflex=1
if (c4(2:2).eq.'0' .and. c4(4:4).eq.'0' .and. p%meanflex) ixflex=0
iyflex=1
if (c4(1:1).eq.'0' .and. c4(3:3).eq.'0' .and. p%meanflex) iyflex=0

hx=3.d0*max(p%xl,p%yl)
drho=p%rhocflex*9.81d0
elt=p%thickflex
dflex=p%ym/12.d0/(1.d0-p%pratio**2)
d=dflex*elt**3
xk=p%rhoaflex*9.81d0

xmean=p%xl/2.
xfmin=xmean-hx/2.
xfmax=xmean+hx/2.
dxf=xfmax-xfmin
ddxf=dxf/(nflex-1)
ymean=p%yl/2.
yfmin=ymean-hx/2.
yfmax=ymean+hx/2.
dyf=yfmax-yfmin
ddyf=dyf/(nflex-1)

allocate (flex(nflex,nflex),work(nflex,nflex),junk(nflex))
!$OMP parallel shared(flex)
!$OMP workshare
flex=0.
!$OMP end workshare nowait
!$OMP end parallel

iflexmin=1+(0.-xfmin)/ddxf
iflexmax=1+(p%xl-xfmin)/ddxf
jflexmin=1+(0.-yfmin)/ddyf
jflexmax=1+(p%yl-yfmin)/ddyf

time=p%time
xl=p%xl
yl=p%yl
nx=p%nx
ny=p%ny
dx=p%dx
dy=p%dy
dt=p%dt

!$OMP parallel shared(flex,ddxf,ddyf,drho,xfmin,yfmin,nx,ny,xl,yl,dx,dy,dt,time) &
!$OMP private(i,j,x,y,h,i1,i2,i3,i4,i5,i6,i7,i8,i9,i0,w1,w2,w3,w4,w5,w6,w7,w8,w9,w0)
!$OMP do schedule(dynamic,nflex/p%num_threads)
do j=1,nflex
y=(j-1)*ddyf+yfmin
do i=1,nflex
x=(i-1)*ddxf+xfmin
h=0.d0
!user supplied dynamic topography start
!user supplied dynamic topography end
flex(i,j)=h*drho*ddxf*ddyf
enddo
enddo
!$OMP end do nowait
!$OMP end parallel

if (ddxf.lt.dx) then

!$OMP parallel shared(flex,ddxf,ddyf,drho,p,xfmin,yfmin,jflexmin,jflexmax,iflexmin,iflexmax,iref) &
!$OMP private(i,j,k,dist,xloc,yloc,ii,jj,r,s,h1,h2,h3,h4,hh,hh0,drhoe)
!$OMP do schedule(dynamic)
do j=jflexmin,jflexmax
yloc=((j-1)*ddyf+yfmin)
jj=1+(p%ny-1)*yloc/p%yl
jj=min(jj,p%ny-1)
jj=max(1,jj)
do i=iflexmin,iflexmax
xloc=((i-1)*ddxf+xfmin)
ii=1+(p%nx-1)*xloc/p%xl
ii=min(ii,p%nx-1)
ii=max(1,ii)
r=(xloc-(ii-1)*p%dx)/p%dx*2.-1.
s=(yloc-(jj-1)*p%dy)/p%dy*2.-1.
h1=t%h(ii+(jj-1)*p%nx)
h2=t%h(ii+1+(jj-1)*p%nx)
h3=t%h(ii+(jj)*p%nx)
h4=t%h(ii+1+(jj)*p%nx)
hh=((1.-r)*(1.-s)*h1+(1.+r)*(1.-s)*h2+(1.-r)*(1.+s)*h3+(1.+r)*(1.+s)*h4)/4.
h1=t%hi_iso(ii+(jj-1)*p%nx)
h2=t%hi_iso(ii+1+(jj-1)*p%nx)
h3=t%hi_iso(ii+(jj)*p%nx)
h4=t%hi_iso(ii+1+(jj)*p%nx)
hh0=iref*((1.-r)*(1.-s)*h1+(1.+r)*(1.-s)*h2+(1.-r)*(1.+s)*h3+(1.+r)*(1.+s)*h4)/4.
drhoe=drho
  do k=1,p%granite_n
  dist=sqrt((xloc-p%granite_x(k))**2/p%granite_rx(k)**2+(yloc-p%granite_y(k))**2/p%granite_ry(k)**2)
  if (dist.lt.1.d0.and. &
      hh0-hh.ge.p%granite_top(k).and.hh0-hh.le.p%granite_bottom(k)) drhoe=drho+p%granite_drho(k)*9.81d0
  enddo
flex(i,j)=flex(i,j)-(hh-hh0)*drhoe*ddxf*ddyf
enddo
enddo
!$OMP end do nowait
!$OMP end parallel

else

allocate (nfl(nflex,nflex),fl(nflex,nflex))
nfl=0
fl=0.d0
!$OMP parallel shared(fl,nfl,p,t,ddxf,ddyf,xfmin,yfmin,iref) private(i,j,xloc,yloc,ii,jj)
!$OMP do schedule(dynamic,p%ny/p%num_threads)
  do j=1,p%ny
  yloc=p%dy*(j-1)+ddyf/2.d0
  jj=1+(yloc-yfmin)/ddyf
    do i=1,p%nx
    xloc=p%dx*(i-1)+ddxf/2.d0
    ii=1+(xloc-xfmin)/ddxf
    fl(ii,jj)=fl(ii,jj)+t%h(i+(j-1)*p%nx)-iref*t%hi_iso(i+(j-1)*p%nx)
    nfl(ii,jj)=nfl(ii,jj)+1
    enddo
  enddo
!$OMP end do nowait
!$OMP end parallel

!$OMP parallel shared(flex,fl,nfl,ddxf,ddyf,drho,p,xfmin,yfmin,jflexmin,jflexmax,iflexmin,iflexmax) &
!$OMP private(i,j,xloc,yloc,ii,jj,r,s,h1,h2,h3,h4,hh,hh0,drhoe,k,dist)
!$OMP do schedule(dynamic)
do j=jflexmin,jflexmax
yloc=((j-1)*ddyf+yfmin)
jj=1+(p%ny-1)*yloc/p%yl
jj=min(jj,p%ny-1)
jj=max(1,jj)
do i=iflexmin,iflexmax
xloc=((i-1)*ddxf+xfmin)
ii=1+(p%nx-1)*xloc/p%xl
ii=min(ii,p%nx-1)
ii=max(1,ii)
r=(xloc-(ii-1)*p%dx)/p%dx*2.-1.
s=(yloc-(jj-1)*p%dy)/p%dy*2.-1.
h1=t%h(ii+(jj-1)*p%nx)
h2=t%h(ii+1+(jj-1)*p%nx)
h3=t%h(ii+(jj)*p%nx)
h4=t%h(ii+1+(jj)*p%nx)
hh=((1.-r)*(1.-s)*h1+(1.+r)*(1.-s)*h2+(1.-r)*(1.+s)*h3+(1.+r)*(1.+s)*h4)/4.
h1=t%hi(ii+(jj-1)*p%nx)
h2=t%hi(ii+1+(jj-1)*p%nx)
h3=t%hi(ii+(jj)*p%nx)
h4=t%hi(ii+1+(jj)*p%nx)
hh0=iref*((1.-r)*(1.-s)*h1+(1.+r)*(1.-s)*h2+(1.-r)*(1.+s)*h3+(1.+r)*(1.+s)*h4)/4.
drhoe=drho
  do k=1,p%granite_n
  dist=sqrt((xloc-p%granite_x(k))**2/p%granite_rx(k)**2+(yloc-p%granite_y(k))**2/p%granite_ry(k)**2)
  if (dist.lt.1.d0.and. &
      hh0-hh.ge.p%granite_top(k).and.hh0-hh.le.p%granite_bottom(k)) drhoe=drho+p%granite_drho(k)*9.81d0
  enddo
if (nfl(i,j).ne.0) flex(i,j)=flex(i,j)-fl(i,j)/nfl(i,j)*drhoe*ddxf*ddyf
enddo
enddo
!$OMP end do nowait
!$OMP end parallel

deallocate (nfl,fl)


endif

if (ixflex.eq.0) then
!$OMP parallel shared(flex,nflex,iflexmin,iflexmax) private(j,fmean)
!$OMP do schedule(dynamic,nflex/p%num_threads)
  do j=1,nflex
  fmean=sum(flex(iflexmin:iflexmax,j))/(iflexmax-iflexmin)
  flex(1:nflex,j)=fmean
  enddo
!$OMP end do nowait
!$OMP end parallel
endif

if (iyflex.eq.0) then
!$OMP parallel shared(flex,nflex,jflexmin,jflexmax) private(i,fmean)
!$OMP do schedule(dynamic,nflex/p%num_threads)
  do i=1,nflex
  fmean=sum(flex(i,jflexmin:jflexmax))/(jflexmax-jflexmin)
  flex(i,1:nflex)=fmean
  enddo
!$OMP end do nowait
!$OMP end parallel
endif

!$OMP parallel shared(flex,nflex) private(j)
!$OMP do schedule(dynamic,nflex/p%num_threads)
do j=1,nflex
call sinft (flex(:,j),nflex)
enddo
!$OMP end do nowait
!$OMP end parallel

!$OMP parallel shared(work,flex)
!$OMP workshare
work=transpose(flex)
!$OMP end workshare nowait
!$OMP end parallel

!$OMP parallel shared(work,nflex) private(i)
!$OMP do schedule(dynamic,nflex/p%num_threads)
do i=1,nflex
call sinft (work(:,i),nflex)
enddo
!$OMP end do nowait
!$OMP end parallel

!$OMP parallel shared(work,hx)
!$OMP workshare
work=work*4./hx/hx
!$OMP end workshare nowait
!$OMP end parallel

pi=3.141592654
pihx=pi/hx
!$OMP parallel shared(work,xk,pihx,d,nflex,ixflex,iyflex) private(i,j,fi,fj,tij)
!$OMP do schedule(dynamic,nflex/p%num_threads)
do j=1,nflex
fj=(j*pihx)**2*iyflex
  do i=1,nflex
  fi=(i*pihx)**2*ixflex
  tij=d/xk*(fi**2+2.*fi*fj+fj**2)+1.d0
  work(j,i)=work(j,i)/xk/tij
  enddo
enddo
!$OMP end do nowait
!$OMP end parallel

!$OMP parallel shared(work,nflex) private(i)
!$OMP do schedule(dynamic,nflex/p%num_threads)
do i=1,nflex
call sinft (work(:,i),nflex)
enddo
!$OMP end do nowait
!$OMP end parallel

!$OMP parallel shared(work,flex)
!$OMP workshare
flex=transpose(work)
!$OMP end workshare nowait
!$OMP end parallel

!$OMP parallel shared(flex,nflex) private(j)
!$OMP do schedule(dynamic,nflex/p%num_threads)
do j=1,nflex
call sinft (flex(:,j),nflex)
enddo
!$OMP end do nowait
!$OMP end parallel

if (ixflex.eq.0) then
iflexfirst=1+(nflex-1)*(0.-xfmin)/dxf
!$OMP parallel shared(flex,nflex,iflexfirst) private(j)
!$OMP do schedule(dynamic,nflex/p%num_threads)
  do j=1,nflex
  flex(iflexfirst,j)=flex(iflexfirst+1,j)
  enddo
!$OMP end do nowait
!$OMP end parallel
iflexlast=1+(nflex-1)*(p%xl-xfmin)/dxf
!$OMP parallel shared(flex,nflex,iflexlast) private(j)
!$OMP do schedule(dynamic,nflex/p%num_threads)
  do j=1,nflex
  flex(iflexlast,j)=flex(iflexlast-1,j)
  enddo
!$OMP end do nowait
!$OMP end parallel
endif

if (iyflex.eq.0) then
jflexfirst=1+(nflex-1)*(0.-yfmin)/dyf
!$OMP parallel shared(flex,nflex,jflexfirst) private(i)
!$OMP do schedule(dynamic,nflex/p%num_threads)
  do i=1,nflex
  flex(i,jflexfirst)=flex(i,jflexfirst+1)
  enddo
!$OMP end do nowait
!$OMP end parallel
jflexlast=1+(nflex-1)*(p%yl-yfmin)/dyf
!$OMP parallel shared(flex,nflex,jflexlast) private(i)
!$OMP do schedule(dynamic,nflex/p%num_threads)
  do i=1,nflex
  flex(i,jflexlast)=flex(i,jflexlast-1)
  enddo
!$OMP end do nowait
!$OMP end parallel
endif

!$OMP parallel shared(t,p,g,dxf,dyf,xfmin,yfmin,ddxf,ddyf,nflex) private(i,j,ij,xp,yp,iflex,jflex,xflex,rx,yflex,ry,hisotot,dh)
!$OMP do schedule(dynamic,p%ny/p%num_threads)
      do j=1,p%ny
        do i=1,p%nx
        ij=i+(j-1)*p%nx
          if (.not.g%bc(ij)) then
          xp=p%xl*float(i-1)/(p%nx-1)
          yp=p%yl*float(j-1)/(p%ny-1)
          iflex=1+(nflex-1)*(xp-xfmin)/dxf
          jflex=1+(nflex-1)*(yp-yfmin)/dyf
          xflex=(iflex-1)*ddxf+xfmin
          rx=(xp-xflex)/ddxf*2.d0-1.d0
          yflex=(jflex-1)*ddyf+yfmin
          ry=(yp-yflex)/ddyf*2.d0-1.d0
          hisotot=(flex(iflex,jflex)*(1.d0-rx)*(1.d0-ry)/4.d0 &
                 +flex(iflex+1,jflex)*(1.d0+rx)*(1.d0-ry)/4.d0 &
                 +flex(iflex+1,jflex+1)*(1.d0+rx)*(1.d0+ry)/4.d0 &
                 +flex(iflex,jflex+1)*(1.d0-rx)*(1.d0+ry)/4.d0)
          dh=hisotot-t%hiso(ij)
          t%hiso(ij)=hisotot
          t%h(ij)=t%h(ij)+dh
          t%hb(ij)=t%hb(ij)+dh
          t%hi(ij)=t%hi(ij)+dh
          t%href(ij)=t%href(ij)+dh
          t%hi_iso(ij)=t%hi_iso(ij)+dh
          t%u(ij)=t%u(ij)+dh/p%dt
          endif
!        t%hi_iso(ij)=t%hi_iso(ij)+dh
        enddo
      enddo
!$OMP end do nowait
!$OMP end parallel

deallocate (flex,work,junk)

end subroutine flexure
