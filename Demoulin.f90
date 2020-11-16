subroutine demoulin (t,g,p,l,sr)

! to compute Alain Demoulin's index based on his two papers on the subject:
! Demoulin, Geomorphology, 2011
! Demoulin, GRL, 2012

use omp_lib
use definitions

implicit none

type (topography) t
type (topology) g
type (param) p
type (law) l

real, dimension(:), allocatable :: bin1,bin2,bin3
real bin1max,bin1min,bin2Min,bin2max,bin3min,bin3max
logical, dimension(:), allocatable :: trunk
double precision, dimension(:), allocatable :: area,ap,rp
double precision ir,ib,r,a,amax,dmax,pi,sr
double precision sumx,sumy,sumxy,sumxx
double precision sumxp,sumyp,sumxyp,sumxxp
integer i,j,k,ij,ijk,minbin,maxbin,kmax,first,ibin,n,np,imouth
integer i1,i2,j1,j2,nmin,nmax

if (p%debug) print*,'calculating Demoulin'

! p%stream_threshold_area is the minimum drainage area for a pixel to be considered as a stream
! p%minimum_area is the minimum area for a drainage basin to be considered for the computation
! of Sr
! minbin and maxbin are the minimum and maximum elevations considered (in m as integers)

!minimum_area=10.d6
minbin=1
maxbin=10000
pi=atan(1.d0)*4.d0

allocate (area(p%nx*p%ny),trunk(p%nx*p%ny))

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

if (minval(area).gt.p%stream_threshold_area) then
print*,' '
print*,'You cant calculate Demoulin s index'
print*,'if the resolution of the model'
print*,'is larger than the maximum area'
print*,'for a stream pixel, i.e.:'
print*,minval(area),'>',p%stream_threshold_area
stop
endif

if (maxval(area).lt.p%minimum_area) then
print*,' '
print*,'You cant calculate Demoulin s index'
print*,'if the maximum drainage area'
print*,'is smaller than the minimum area'
print*,'required to compute R1, i.e.:'
print*,maxval(area),'<',p%minimum_area
stop
endif

! goes through the stack to find the mouth of each basin and compute the
! R1 index for it
first=0
trunk=.false.
n=0
sumx=0.d0;sumy=0.d0;sumxy=0.d0;sumxx=0.d0
!$OMP parallel shared(g,t,p,trunk,n,sumx,sumy,sumxy,sumxx,minbin,maxbin) &
!$OMP private(bin1,bin2,bin3,first,ap,rp,np,sumxp,sumyp,sumxyp,sumxxp)
allocate (bin1(minbin:maxbin),bin2(minbin:maxbin),bin3(minbin:maxbin))
allocate (ap(p%nx*p%ny/p%num_threads),rp(p%nx*p%ny/p%num_threads))
np=0
sumxp=0.d0;sumyp=0.d0;sumxyp=0.d0;sumxxp=0.d0
!$OMP do private(i,j,ij,ijk,ibin,a,ib,ir,r,amax,kmax,imouth,dmax,bin1max,bin1min)
  do j=1,p%num_threads
    do i=1,g%nstackp
    ij=g%stack((j-1)*g%nstackp+i)
      if (ij.ne.0) then
      ibin=int(t%h(ij))
      ibin=min(maxbin,ibin)
      ibin=max(minbin,ibin)
      ijk=g%receiver(ij)
! finds a new catchment
        if (ijk.eq.ij) then
! special treatment for the first catchment
          if (first.eq.0) then
          bin1=0.;bin2=0.;bin3=0.
          a=area(ij)
          trunk(ij)=.true.
          i1=1+mod(ij-1,p%nx)
          j1=1+(ij-1)/p%nx
          dmax=0.d0
          first=1
          else
! compute indices
            if (a.gt.p%minimum_area) then
            nmin=0
              do k=minbin,maxbin
              if (bin1(k).gt.0.5) nmax=k
              enddo
              do k=maxbin,minbin,-1
              if (bin1(k).gt.0.5) nmin=k
              enddo
            if (nmin.eq.nmax) goto 1111
            if (nmin.eq.0) goto 1111
              do k=nmax-1,nmin,-1
              bin1(k)=bin1(k)+bin1(k+1)
              bin2(k)=bin2(k)+bin2(k+1)
              bin3(k)=bin3(k)+bin3(k+1)
              enddo
            bin1max=bin1(nmin)
            if (bin1max.eq.0.d0) goto 1111
            bin2max=bin2(nmin)
            if (bin2max.eq.0.d0) goto 1111
            bin3max=bin3(nmin)
            if (bin3max.eq.0.d0) goto 1111
            bin1(nmin:nmax)=bin1(nmin:nmax)/bin1max
            bin2(nmin:nmax)=bin2(nmin:nmax)/bin2max
            bin3(nmin:nmax)=bin3(nmin:nmax)/bin3max
            ib=sum(bin1(nmin:nmax))-sum(bin2(nmin:nmax))
            ir=sum(bin2(nmin:nmax))-sum(bin3(nmin:nmax))
            if (ib.eq.0.d0) goto 1111
            r=ir/ib
            np=np+1
            ap(np)=log(a/1.d6)
            rp(np)=r/sqrt(4.d0*a/pi/dmax**2)
            sumxp=sumxp+ap(np)
            sumyp=sumyp+rp(np)
            sumxxp=sumxxp+ap(np)*ap(np)
            sumxyp=sumxyp+ap(np)*rp(np)
1111        continue
            endif
          bin1=0.;bin2=0.;bin3=0.
          a=area(ij)
          trunk(ij)=.true.
          i1=1+mod(ij-1,p%nx)
          j1=1+(ij-1)/p%nx
          dmax=0.
          endif
        endif
      bin1(ibin)=bin1(ibin)+1.
      if (area(ij).ge.p%stream_threshold_area) bin2(ibin)=bin2(ibin)+1.
      if (area(ij).ge.p%stream_threshold_area.and.trunk(ij)) bin3(ibin)=bin3(ibin)+1.
      i2=1+mod(ij-1,p%nx)
      j2=1+(ij-1)/p%nx
      dmax=max(dmax,sqrt(((i1-i2)*p%dx)**2+((j1-j2)*p%dy)**2))
        if (trunk(ij)) then
        amax=0.
        kmax=0
          do k=g%ndon(ij),g%ndon(ij+1)-1
            if (area(g%donors(k)).gt.amax) then
            amax=area(g%donors(k))
            kmax=g%donors(k)
            endif
          enddo
          if (kmax.gt.0) then
          if (area(kmax).gt.p%stream_threshold_area) trunk(kmax)=.true.
          endif
        endif
      endif
    enddo
  enddo
!$OMP end do nowait
! adds the bits from the various processors to the sums
n=n+np
sumx=sumx+sumxp
sumy=sumy+sumyp
sumxx=sumxx+sumxxp
sumxy=sumxy+sumxyp
deallocate (bin1,bin2,bin3,ap,rp)
!$OMP end parallel
! computes the linear regression coefficient to obtain sr
sr=(n*sumxy-sumx*sumy)/(n*sumxx-sumx**2)

deallocate (area,trunk)

if (p%debug) print*,'Demoulin calculated'

return

end subroutine demoulin
