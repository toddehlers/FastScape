subroutine metric (t,p,g,dhdt,fnme,nfnme,flag_scale)

use definitions

implicit none

type (topography) t
type (topology) g
type (param) p

! subroutine to compute landscape metrics depending on the value of the digits of c5; these include
! 10000: average topography in the x- and y-directions
! 01000: topography along a cross-section at the middle of the run in the x- and y-directions
! 00100: slope area relationship but limited to nx/500 by ny/500 nodes picked randomly

character fnme*(*)
integer nfnme,nx,ny,i,j,flag_scale,ii,jj,ij,flag_choice
integer,dimension(:,:),allocatable::receiver
double precision, dimension(:,:),allocatable::h,discharge,length,hb,hi
double precision hmi,hma,hmin,hmax,dx,dy,dt
double precision dhdt(p%nx*p%ny)
double precision slope,slmin,slmax,dzdx,dzdy,con,fltot
double precision fmin,fmax
double precision a,b,c,d,e,f
character cmin*14,cmax*14,c5*5
integer nxskip,nyskip
integer nbin,ibin
double precision, dimension(:),allocatable::fl,al
double precision, dimension(:,:),allocatable::sl,cl

nx=p%nx
ny=p%ny
dt=p%dt
dx=p%dx
dy=p%dy
hmi=p%topo_min
hma=p%topo_max
flag_choice=p%imetric
slmin=p%slope_min
slmax=p%slope_max
allocate (h(nx,ny),discharge(nx,ny),length(nx,ny),hb(nx,ny),hi(nx,ny))
allocate (receiver(nx,ny))

  do j=1,ny
  do i=1,nx
  ij=(j-1)*nx+i
  h(i,j)=t%h(ij)
  hb(i,j)=t%hb(ij)
  hi(i,j)=t%hi(ij)
  discharge(i,j)=t%discharge(ij)
  length(i,j)=t%length(ij)
  receiver(i,j)=g%receiver(ij)
  enddo
  enddo

  if (flag_scale.eq.1) then
  hmin=hmi
  hmax=hma
  else
  hmin=minval(h)
  hmax=maxval(h)
  endif

write (cmin,'(g14.6)') hmin
write (cmax,'(g14.6)') hmax

write (c5,'(i5)') flag_choice

if (c5(1:1).eq.'1') then

open (3,file=fnme(1:nfnme)//'x-av.txt',status='unknown')
do i=1,nx
write (3,'(5g14.6)') minval(h(i,:)),sum(h(i,:))/ny,maxval(h(i,:)),sum(hb(i,:))/ny,sum(hi(i,:))/ny
enddo
close (3)

open (3,file=p%run//'/Metric.R',status='unknown',access='append')
write (3,'(a)') 'pdf("'//fnme(1:nfnme)//'x-av.pdf")'
write (3,'(a)') 'x=read.table("'//fnme(1:nfnme)//'x-av.txt")'
write (3,'(a)') 'plot(x[,3],type="l",main="X-averaged topography - '// &
                'Step='//fnme(nfnme-5:nfnme)//'",'// &
                'ylim=c('//cmin//','//cmax//'),'// &
                'ylab="Height (m)")'
write (3,'(a)') 'lines(x[,2])'
write (3,'(a)') 'lines(x[,1])'
write (3,'(a)') 'lines(x[,4],col="green")'
write (3,'(a)') 'lines(x[,5],col="red")'
write (3,'(a)') 'dev.off()'
write (3,'(a)')
close (3)

open (3,file=fnme(1:nfnme)//'y-av.txt',status='unknown')
do j=1,ny
write (3,'(5g14.6)') minval(h(:,j)),sum(h(:,j))/nx,maxval(h(:,j)),sum(hb(:,j))/nx,sum(hi(:,j))/nx
enddo
close (3)

open (3,file=p%run//'/Metric.R',status='unknown',access='append')
write (3,'(a)') 'pdf("'//fnme(1:nfnme)//'y-av.pdf")'
write (3,'(a)') 'x=read.table("'//fnme(1:nfnme)//'y-av.txt")'
write (3,'(a)') 'plot(x[,3],type="l",main="Y-averaged topography - '// &
                'Step='//fnme(nfnme-5:nfnme)//'",'// &
                'ylim=c('//cmin//','//cmax//'),'// &
                'ylab="Height (m)")'
write (3,'(a)') 'lines(x[,2])'
write (3,'(a)') 'lines(x[,1])'
write (3,'(a)') 'lines(x[,4],col="green")'
write (3,'(a)') 'lines(x[,5],col="red")'
write (3,'(a)') 'dev.off()'
write (3,'(a)')
close (3)

endif

if (c5(2:2).eq.'1') then

open (3,file=fnme(1:nfnme)//'x-md.txt',status='unknown')
do i=1,nx
write (3,'(5g14.6)') h(i,ny/2),hb(i,ny/2),hi(i,ny/2)
enddo
close (3)

open (3,file=p%run//'/Metric.R',status='unknown',access='append')
write (3,'(a)') 'pdf("'//fnme(1:nfnme)//'x-md.pdf")'
write (3,'(a)') 'x=read.table("'//fnme(1:nfnme)//'x-md.txt")'
write (3,'(a)') 'plot(x[,1],type="l",main="X-middle topography - '// &
                'Step='//fnme(nfnme-5:nfnme)//'",'// &
                'ylim=c('//cmin//','//cmax//'),'// &
                'ylab="Height (m)")'
write (3,'(a)') 'lines(x[,2],col="green")'
write (3,'(a)') 'lines(x[,3],col="red")'
write (3,'(a)') 'dev.off()'
write (3,'(a)')
close (3)

open (3,file=fnme(1:nfnme)//'y-md.txt',status='unknown')
do j=1,ny
write (3,'(5g14.6)') h(nx/2,j),hb(nx/2,j),hi(nx/2,j)
enddo
close (3)

open (3,file=p%run//'/Metric.R',status='unknown',access='append')
write (3,'(a)') 'pdf("'//fnme(1:nfnme)//'y-md.pdf")'
write (3,'(a)') 'x=read.table("'//fnme(1:nfnme)//'y-md.txt")'
write (3,'(a)') 'plot(x[,1],type="l",main="Y-middle topography - '// &
                'Step='//fnme(nfnme-5:nfnme)//'",'// &
                'ylim=c('//cmin//','//cmax//'),'// &
                'ylab="Height (m)")'
write (3,'(a)') 'lines(x[,2],col="green")'
write (3,'(a)') 'lines(x[,3],col="red")'
write (3,'(a)') 'dev.off()'
write (3,'(a)')
close (3)

endif

if (c5(3:3).eq.'1') then

open (3,file=fnme(1:nfnme)//'-SA.txt',status='unknown')
nxskip=max(1,nx/500)
nyskip=max(1,ny/500)
do j=1,ny,nyskip
do i=1,nx,nxskip
ij=receiver(i,j)
ii=1+mod(ij-1,nx)
jj=1+int((ij-1)/nx)
  if (ii.ne.i .or. jj.ne.j) then
  slope=(h(i,j)-h(ii,jj))/length(i,j)
  if (slope.gt.0.d0) write (3,'(2g14.6)') discharge(i,j),slope
  endif
enddo
enddo
close (3)

open (3,file=p%run//'/Metric.R',status='unknown',access='append')
write (3,'(a)') 'pdf("'//fnme(1:nfnme)//'-SA.pdf")'
write (3,'(a)') 'x=read.table("'//fnme(1:nfnme)//'-SA.txt")'
write (3,'(a)') 'plot(x[,1],x[,2],main="Slope vs. Area - '// &
                'Step='//fnme(nfnme-5:nfnme)//'",'// &
                'xlab="Area (m^2)",'// &
                'ylab="Slope")'
write (3,'(a)') 'dev.off()'
write (3,'(a)')
close (3)

endif

if (c5(4:4).eq.'1') then
nbin=100
allocate(fl(nbin),al(nbin))
allocate (sl(nx,ny),cl(nx,ny))
sl=0.d0
cl=0.d0
fl=0.d0
al=0.d0
con=45.d0/atan(1.d0)
do j=2,ny-1
do i=2,nx-1
dzdx=((h(i+1,j+1)+2.d0*h(i+1,j)+h(i+1,j-1))-(h(i-1,j+1)+2.d0*h(i-1,j)+h(i-1,j-1)))/8.d0/dx
dzdy=((h(i-1,j-1)+2.d0*h(i,j-1)+h(i+1,j-1))-(h(i-1,j+1)+2.d0*h(i,j+1)+h(i+1,j+1)))/8.d0/dy
sl(i,j)=dzdx**2+dzdy**2
if (sl(i,j).gt.tiny(sl(i,j))) sl(i,j)=atan(sqrt(sl(i,j)))*con
a=(h(i-1,j-1)+h(i+1,j-1)+h(i-1,j)+h(i+1,j)+h(i-1,j+1)+h(i+1,j+1))/p%dx/p%dx/12.d0-(h(i,j-1)+h(i,j)+h(i,j+1))/p%dx/p%dx/6.d0
b=(h(i-1,j-1)+h(i,j-1)+h(i+1,j-1)+h(i-1,j+1)+h(i,j+1)+h(i+1,j+1))/p%dy/p%dy/12.d0-(h(i-1,j)+h(i,j)+h(i+1,j))/p%dy/p%dy/6.d0
c=(h(i+1,j-1)+h(i-1,j+1)-h(i-1,j-1)-h(i+1,j+1))/p%dx/p%dy/4.d0
d=(h(i+1,j-1)+h(i+1,j)+h(i+1,j+1)-h(i-1,j-1)-h(i-1,j)-h(i-1,j+1))/p%dx/6.d0
e=(h(i-1,j-1)+h(i,j-1)+h(i+1,j-1)-h(i-1,j+1)-h(i,j+1)-h(i+1,j+1))/p%dy/6.d0
f=(2.d0*(h(i,j-1)+h(i-1,j)+h(i+1,j)+h(i,j+1))-(h(i-1,j-1)+h(i+1,j-1)+h(i-1,j+1)+h(i+1,j+1))+5.d0*h(i,j))/9.d0
cl(i,j)=1.d0+d**2+e**2
if (abs(cl(i,j)).gt.tiny(cl(i,j))) cl(i,j)=(a*(1.d0+e**2)+b*(1.d0+d**2)-c*d*e)/(cl(i,j)**(3.d0/2.d0))
enddo
enddo
if (slmax.lt.0.d0) then
slmin=minval(sl)
slmax=maxval(sl)
endif
if (slmin.ne.slmax) then
do j=1,ny
do i=1,nx
ibin=1+int((nbin-1)*((sl(i,j)-slmin)/(slmax-slmin)))
ibin=min(ibin,nbin)
ij=(j-1)*nx+i
fl(ibin)=fl(ibin)+dhdt(ij)
al(ibin)=al(ibin)+1.d0
enddo
enddo
  do ibin=1,nbin
  if (al(ibin).gt.0.d0) al(ibin)=fl(ibin)/al(ibin)
  enddo
fltot=sum(fl)
fl=fl/fltot
endif
!open (3,file=fnme(1:nfnme)//'-FS.txt',status='unknown')
!do ibin=1,nbin
!write (3,'(f14.10,",",f14.10,",",f14.10)') slmin+(slmax-slmin)*dfloat(ibin-1)/(nbin-1),fl(ibin),al(ibin)
!enddo
!close (3)

open (3,file=fnme(1:nfnme)//'-HeightHist.txt',status='unknown')
fmin=minval(h,h>0.d0)
fmax=maxval(h,h>0.d0)
if (p%HeightDis_min.gt.-998.d0) fmin=p%HeightDis_min
if (p%HeightDis_max.gt.-998.d0) fmax=p%HeightDis_max
fl=0.d0
do j=1,ny
do i=1,nx
if (h(i,j).gt.0.d0) then
ibin=1+int((nbin-1)*((h(i,j)-fmin)/(fmax-fmin)))
ibin=min(ibin,nbin)
fl(ibin)=fl(ibin)+1.d0
endif
enddo
enddo
do ibin=1,nbin-1
fl(ibin+1)=fl(ibin)+fl(ibin+1)
enddo
fl=fl/fl(nbin)
do ibin=1,nbin
write (3,'(e16.8,",",f14.10)') fmin+(fmax-fmin)*dfloat(ibin-1)/(nbin-1),fl(ibin)
enddo
close (3)

open (3,file=fnme(1:nfnme)//'-SlopeHist.txt',status='unknown')
fmin=minval(sl,h>0.d0)
fmax=maxval(sl,h>0.d0)
if (p%SlopeDis_min.gt.-998.d0) fmin=p%SlopeDis_min
if (p%SlopeDis_max.gt.-998.d0) fmax=p%SlopeDis_max
fl=0.d0
do j=1,ny
do i=1,nx
if (h(i,j).gt.0.d0) then
ibin=1+int((nbin-1)*((sl(i,j)-fmin)/(fmax-fmin)))
ibin=min(ibin,nbin)
fl(ibin)=fl(ibin)+1.d0
endif
enddo
enddo
do ibin=1,nbin-1
fl(ibin+1)=fl(ibin)+fl(ibin+1)
enddo
fl=fl/fl(nbin)
do ibin=1,nbin
write (3,'(e16.8,",",f14.10)') fmin+(fmax-fmin)*dfloat(ibin-1)/(nbin-1),fl(ibin)
enddo
close (3)

open (3,file=fnme(1:nfnme)//'-CurveHist.txt',status='unknown')
fmin=minval(cl,h>0.d0)
fmax=maxval(cl,h>0.d0)
if (p%CurveDis_min.gt.-998.d0) fmin=p%CurveDis_min
if (p%CurveDis_max.gt.-998.d0) fmax=p%CurveDis_max
fl=0.d0
do j=1,ny
do i=1,nx
if (h(i,j).gt.0.d0) then
ibin=1+int((nbin-1)*((cl(i,j)-fmin)/(fmax-fmin)))
ibin=min(ibin,nbin)
ibin=max(ibin,1)
fl(ibin)=fl(ibin)+1.d0
endif
enddo
enddo
do ibin=1,nbin-1
fl(ibin+1)=fl(ibin)+fl(ibin+1)
enddo
fl=fl/fl(nbin)
do ibin=1,nbin
write (3,'(e16.8,",",f14.10)') fmin+(fmax-fmin)*dfloat(ibin-1)/(nbin-1),fl(ibin)
enddo
close (3)

open (3,file=fnme(1:nfnme)//'-DischargeHist.txt',status='unknown')
discharge=log10(discharge)
fmin=minval(discharge,h>0.d0)
fmax=maxval(discharge,h>0.d0)
if (p%DischargeDis_min.gt.-998.d0) fmin=p%DischargeDis_min
if (p%DischargeDis_max.gt.-998.d0) fmax=p%DischargeDis_max
fl=0.d0
do j=1,ny
do i=1,nx
if (h(i,j).gt.0.d0) then
ibin=1+int((nbin-1)*((discharge(i,j)-fmin)/(fmax-fmin)))
ibin=min(ibin,nbin)
fl(ibin)=fl(ibin)+1.d0
endif
enddo
enddo
do ibin=1,nbin-1
fl(ibin+1)=fl(ibin)+fl(ibin+1)
enddo
fl=fl/fl(nbin)
do ibin=1,nbin
write (3,'(e16.8,",",f14.10)') fmin+(fmax-fmin)*dfloat(ibin-1)/(nbin-1),fl(ibin)
enddo
close (3)

deallocate(sl,fl,cl)

!open (3,file=p%run//'/Metric.R',status='unknown',access='append')
!write (3,'(a)') 'pdf("'//fnme(1:nfnme)//'-FS.pdf")'
!write (3,'(a)') 'x=read.table("'//fnme(1:nfnme)//'-FS.txt")'
!write (3,'(a)') 'plot(x[,1],x[,2],main="Flux vs. Slope - '// &
                !'Step='//fnme(nfnme-5:nfnme)//'",'// &
                !'xlab="Slope",'// &
                !'ylab="Flux (m^3/yr)")'
!write (3,'(a)') 'dev.off()'
!write (3,'(a)')
!close (3)

open (3,file=p%run//'/Metric.R',status='unknown',access='append')
write (3,'(a)') 'pdf("'//fnme(1:nfnme)//'-HeightCDF.pdf")'
write (3,'(a)') 'x=read.table("'//fnme(1:nfnme)//'-HeightHist.txt",sep=",")'
write (3,'(a)') 'plot(x[,1],x[,2],main="Cumulative Height Distribution - '// &
                'Step='//fnme(nfnme-5:nfnme)//'",'// &
                'type="l",lwd=3,'// &
                'xlab="Height (m)",'// &
                'ylab="CDF")'
write (3,'(a)') 'dev.off()'
write (3,'(a)')
close (3)

open (3,file=p%run//'/Metric.R',status='unknown',access='append')
write (3,'(a)') 'pdf("'//fnme(1:nfnme)//'-SlopeCDF.pdf")'
write (3,'(a)') 'x=read.table("'//fnme(1:nfnme)//'-SlopeHist.txt",sep=",")'
write (3,'(a)') 'plot(x[,1],x[,2],main="Cumulative Slope Distribution - '// &
                'Step='//fnme(nfnme-5:nfnme)//'",'// &
                'type="l",lwd=3,'// &
                'xlab="Slope (degrees)",'// &
                'ylab="CDF")'
write (3,'(a)') 'dev.off()'
write (3,'(a)')
close (3) 

open (3,file=p%run//'/Metric.R',status='unknown',access='append')
write (3,'(a)') 'pdf("'//fnme(1:nfnme)//'-CurveCDF.pdf")'
write (3,'(a)') 'x=read.table("'//fnme(1:nfnme)//'-CurveHist.txt",sep=",")'
write (3,'(a)') 'plot(x[,1],x[,2],main="Cumulative Curvature Distribution - '// &
                'Step='//fnme(nfnme-5:nfnme)//'",'// &
                'type="l",lwd=3,'// &
                'xlab="Curvature (1/m)",'// &
                'ylab="CDF")'
write (3,'(a)') 'dev.off()'
write (3,'(a)')
close (3) 

open (3,file=p%run//'/Metric.R',status='unknown',access='append')
write (3,'(a)') 'pdf("'//fnme(1:nfnme)//'-DischargeCDF.pdf")'
write (3,'(a)') 'x=read.table("'//fnme(1:nfnme)//'-DischargeHist.txt",sep=",")'
write (3,'(a)') 'plot(x[,1],x[,2],main="Cumulative Discharge  Distribution - '// &
                'Step='//fnme(nfnme-5:nfnme)//'",'// &
                'type="l",lwd=3,'// &
                'xlab="Discharge (m^3/yr)",'// &
                'ylab="CDF")'
write (3,'(a)') 'dev.off()'
write (3,'(a)')
close (3) 

endif

deallocate (h,discharge,length,hb,hi)
deallocate (receiver)

end subroutine metric
