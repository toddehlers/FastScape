subroutine fluvial_erosion (t,g,p,l)

! to compute fluvial erosion according to the stream power law equation
! two modes are implemented depending on the value of the slope exponent, n
! ie if it is untiy or not

! note that a third law based on Davy and Lague (2009) paper has been coded
! but is not fully debugged/checked; at this stage it is not documented in the user guide

use omp_lib
use definitions

implicit none

type (topography) t
type (topology) g
type (param) p
type (law) l

double precision kf,m,tauc,ss,fact,factp,dh,n,hn,hnn,tol,depos,deposmax,lf,dhs
double precision x,y,dist,sc,dsc,find_kf,kfp
integer i,j,k,ij,ijk,nk,is,ijp
external find_kf
double precision, dimension(:), allocatable :: scrd

allocate (scrd(p%nn))
call random_number (scrd)

scrd=scrd*2.d0-1.d0

if (l%law.eq.1) then
m=l%m
tauc=l%tauc
sc=tan(l%sc*3.141592654d0/180.d0)
dsc=0.d0
if (l%dsc.gt.0.d0) dsc=l%dsc
t%pmask(1,:)=.false.
!$OMP parallel num_threads(p%num_threads) shared(p,g,t,l,m,tauc) private(j,i,ij,ijk,kf,ss,fact,dh)
!$OMP do schedule(dynamic,1)
  do j=1,p%num_threads
    do i=1,g%nstackp
    ij=g%stack((j-1)*g%nstackp+i)
      if (ij.ne.0) then
        if (.not.g%bc(ij).and.t%h(ij).gt.p%sea_level) then
        ijk=g%receiver(ij)
          if (ijk.ne.ij .and. t%length(ij).gt.0.d0) then
          kf=find_kf(ij,l,t,p)
          if (kf.lt.0.d0) then
          dh=t%h(ij)-p%sea_level
          t%h(ij)=p%sea_level
!          t%pmask(1,ij)=.true.
          else
          ss=kf*t%discharge(ij)**m*(t%h(ij)-t%h(ijk))/t%length(ij)
            if (ss.gt.tauc .and. t%discharge(ij).gt.l%area_min) then
            fact=kf*p%dt*t%discharge(ij)**m/t%length(ij)
            dh=(t%h(ij)+fact*t%h(ijk)+p%dt*tauc)/(1.d0+fact)-t%h(ij)
            t%h(ij)=t%h(ij)+dh
            t%pmask(1,ij)=.true.
            endif
          endif
          if (sc.gt.0.d0) then
          t%h(ij)=min(t%h(ij),t%h(ijk)+sc*(1.d0+dsc*scrd(ij))*t%length(ij))
          endif
          endif
        endif
      endif
    enddo
  enddo
!$OMP end do nowait
!$OMP end parallel

elseif (l%law.eq.-1) then
m=l%m
tauc=l%tauc
sc=tan(l%sc*3.141592654d0/180.d0)
t%pmask(1,:)=.false.
!$OMP parallel shared(g,t,m,p,tauc)
!$OMP do private(i,j,ij,ijk,ss,fact,dh,kf) schedule(dynamic,1)
  do j=1,p%num_threads
    do i=1,g%nstackp
    ij=g%stack((j-1)*g%nstackp+i)
      if (ij.ne.0) then
        if (.not.g%bc(ij).and.t%h(ij).gt.p%sea_level) then
        ijk=g%receiver(ij)
          if (ijk.ne.ij .and. t%length(ij).gt.0.d0) then
          kf=find_kf(ij,l,t,p)
          if (kf.lt.0.d0) then
          dh=t%h(ij)-p%sea_level
          t%h(ij)=p%sea_level
!          t%pmask(1,ij)=.true.
          else
          ss=kf*t%discharge(ij)**m*(t%h(ij)-t%h(ijk))/t%length(ij)
            if (ss.gt.tauc .and. t%discharge(ij).gt.l%area_min) then
            dh=-p%dt*(kf*t%discharge(ij)**m*(t%h(ij)-t%h(ijk))/t%length(ij)-tauc)
            t%h(ij)=t%h(ij)+dh
            t%pmask(1,ij)=.true.
            endif
          endif
          if (sc.gt.0.d0) then
          t%h(ij)=min(t%h(ij),t%h(ijk)+sc*t%length(ij))
          endif
          endif
        endif
      endif
    enddo
  enddo
!$OMP end do nowait
!$OMP end parallel

elseif (l%law.eq.2) then
m=l%m
n=l%n
nk=l%nk
tauc=l%tauc
sc=tan(l%sc*3.141592654d0/180.d0)
t%pmask(1,:)=.false.
if (nk.gt.0) then
!$OMP parallel shared(g,t,m,n,nk,p,tauc,sc)
!$OMP do private(i,j,ij,ijk,ss,fact,dh,hnn,hn,kf) schedule(dynamic,1)
  do j=1,p%num_threads
    do i=1,g%nstackp
    ij=g%stack((j-1)*g%nstackp+i)
      if (ij.ne.0) then
        if (.not.g%bc(ij).and.t%h(ij).gt.p%sea_level) then
        ijk=g%receiver(ij)
          if (ijk.ne.ij .and. t%length(ij).gt.0.d0) then
          kf=find_kf(ij,l,t,p)
          if (kf.lt.0.d0) then
          dh=t%h(ij)-p%sea_level
          t%h(ij)=p%sea_level
!          t%pmask(1,ij)=.true.
          else
          ss=kf*t%discharge(ij)**m*((t%h(ij)-t%h(ijk))/t%length(ij))**n
            if (ss.gt.tauc .and. t%discharge(ij).gt.l%area_min) then
            fact=p%dt*kf*t%discharge(ij)**m/t%length(ij)**n
            hnn=t%h(ij)
              do k=1,nk
              hn=max(hnn,t%h(ijk)+1.d-10)
              hnn=hn-(hn-t%h(ij)+fact*(hn-t%h(ijk))**n-tauc*p%dt)/(1.d0+fact*(hn-t%h(ijk))**(n-1.d0)*n)
              enddo
            dh=hn-t%h(ij)
            t%h(ij)=t%h(ij)+dh
!            t%pmask(1,ij)=.true.
            endif
          endif
          if (sc.gt.0.d0) then
          t%h(ij)=min(t%h(ij),t%h(ijk)+sc*t%length(ij))
          endif
          endif
        endif
      endif
    enddo
  enddo
!$OMP end do nowait
!$OMP end parallel
else
tol=l%tol
!$OMP parallel shared(g,t,m,n,p,tauc,tol)
!$OMP do private(i,j,ij,ijk,ss,fact,dh,hnn,hn,kf) schedule(dynamic,1)
  do j=1,p%num_threads
    do i=1,g%nstackp
    ij=g%stack((j-1)*g%nstackp+i)
      if (ij.ne.0) then
        if (.not.g%bc(ij).and.t%h(ij).gt.p%sea_level) then
        ijk=g%receiver(ij)
          if (ijk.ne.ij .and. t%length(ij).gt.0.d0) then
          kf=find_kf(ij,l,t,p)
          if (kf.lt.0.d0) then
          dh=t%h(ij)-p%sea_level
          t%h(ij)=p%sea_level
!          t%pmask(1,ij)=.true.
          else
          ss=kf*t%discharge(ij)**m*((t%h(ij)-t%h(ijk))/t%length(ij))**n
            if (ss.gt.tauc .and. t%discharge(ij).gt.l%area_min) then
            fact=p%dt*kf*t%discharge(ij)**m/t%length(ij)**n
            hnn=t%h(ij)
11111         continue
              hn=max(hnn,t%h(ijk)+1.d-10)
              hnn=hn-(hn-t%h(ij)+fact*(hn-t%h(ijk))**n-tauc*p%dt)/(1.d0+fact*(hn-t%h(ijk))**(n-1.d0)*n)
              if (abs(hn-hnn).gt.tol) goto 11111
            dh=hn-t%h(ij)
            t%h(ij)=t%h(ij)+dh
!            t%pmask(1,ij)=.true.
            endif
          endif
          if (sc.gt.0.d0) then
          t%h(ij)=min(t%h(ij),t%h(ijk)+sc*t%length(ij))
          endif
          endif
        endif
      endif
    enddo
  enddo
!$OMP end do nowait
!$OMP end parallel
endif

elseif (l%law.eq.-2) then
m=l%m
n=l%n
tauc=l%tauc
sc=tan(l%sc*3.141592654d0/180.d0)
t%pmask(1,:)=.false.
!$OMP parallel shared(g,t,m,n,tauc,p)
!$OMP do private(i,j,ij,ijk,ss,fact,dh,kf) schedule(dynamic,1)
  do j=1,p%num_threads
    do i=1,g%nstackp
    ij=g%stack((j-1)*g%nstackp+i)
      if (ij.ne.0) then
        if (.not.g%bc(ij).and.t%h(ij).gt.p%sea_level) then
        ijk=g%receiver(ij)
          if (ijk.ne.ij .and. t%length(ij).gt.0.d0) then
          kf=find_kf(ij,l,t,p)
          if (kf.lt.0.d0) then
          dh=t%h(ij)-p%sea_level
          t%h(ij)=p%sea_level
          t%pmask(1,ij)=.true.
          else
          ss=kf*t%discharge(ij)**m*((t%h(ij)-t%h(ijk))/t%length(ij))**n
            if (ss.gt.tauc .and. t%discharge(ij).gt.l%area_min) then
            dh=-p%dt*(kf*t%discharge(ij)**m*((t%h(ij)-t%h(ijk))/t%length(ij))**n-tauc)
            t%h(ij)=t%h(ij)+dh
            t%pmask(1,ij)=.true.
            endif
          endif
          if (sc.gt.0.d0) then
          t%h(ij)=min(t%h(ij),t%h(ijk)+sc*t%length(ij))
          endif
          endif
        endif
      endif
    enddo
  enddo
!$OMP end do nowait
!$OMP end parallel

elseif (l%law.eq.3) then
m=l%m
lf=l%lf
sc=tan(l%sc*3.141592654d0/180.d0)
t%pmask(1,:)=.false.
!$OMP parallel shared(g,t,m,n,p,tauc,tol)
!$OMP do private(i,j,ij,ijk,fact,dh,hnn,hn,kf) schedule(dynamic,1)
  do j=1,p%num_threads
    do i=1,g%nstackp
    ij=g%stack((j-1)*g%nstackp+i)
      if (ij.ne.0) then
        if (.not.g%bc(ij).and.t%h(ij).gt.p%sea_level) then
        ijk=g%receiver(ij)
          if (ijk.ne.ij) then
          kfp=find_kf(ij,l,t,p)
          kf=kfp
          dhs=t%h(ij)-t%hb(ij)
          if (dhs.gt.2.d-10) kf=l%kfs
          if (kf.lt.0.d0) then
          dh=t%h(ij)-p%sea_level
          t%h(ij)=p%sea_level
          t%pmask(1,ij)=.true.
          else
          fact=0.d0
          if (t%length(ij).gt.0.d0) fact=kf*p%dt*t%discharge(ij)**m/t%length(ij)
          factp=t%sediment(ij)*p%dt/(lf*t%discharge(ij))
          dh=(t%h(ij)+fact*t%h(ijk)+factp)/(1.d0+fact)-t%h(ij)
            if (dhs.gt.2.d-10 .and. -dh.gt.dhs) then
            t%h(ij)=t%h(ij)-dhs
            kf=kfp
            fact=0.d0
            if (t%length(ij).gt.0.d0) fact=kf*p%dt*t%discharge(ij)**m/t%length(ij)
            factp=t%sediment(ij)*p%dt/(lf*t%discharge(ij))
            dh=(t%h(ij)+fact*t%h(ijk)+factp)/(1.d0+fact)-t%h(ij)
            t%h(ij)=t%h(ij)+dh
            else
            t%h(ij)=t%h(ij)+dh
            endif
          t%pmask(1,ij)=.true.
          endif
          if (sc.gt.0.d0) then
          t%h(ij)=min(t%h(ij),t%h(ijk)+sc*t%length(ij))
          endif
          endif
        endif
      endif
    enddo
  enddo
!$OMP end do nowait
!$OMP end parallel

endif

if (tauc.lt.tiny(tauc) .and. l%area_min.lt.tiny(l%area_min)) t%pmask=.false.

return

end subroutine fluvial_erosion
