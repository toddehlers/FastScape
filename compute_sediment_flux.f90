subroutine compute_sediment_flux (t,g,p,l)

! to compute sediment flux for Beaumont's law (3)

use omp_lib
use definitions

implicit none

type (topography) t
type (topology) g
type (param) p
type (law) l

integer i,j,ij,ijk,k
double precision dhdt,kf,dhs,dh,x,y,dist
double precision, dimension(:), allocatable :: sediment

!!$OMP parallel shared(t)
!!$OMP workshare
!t%sediment=0.d0
!!$OMP end workshare nowait
!!$OMP end parallel

!$OMP parallel shared(t,p) private(ij)
!$OMP do
  do ij=1,p%nn
  t%sediment(ij)=0.d0
  enddo
!$OMP end do nowait
!$OMP end parallel

!$OMP parallel shared(t,g,p,l) private(i,j,ij,ijk,dhs,kf,dhdt,dh,x,y,dist)
!$OMP do
  do j=1,p%num_threads
    do i=g%nstackp,1,-1
    ij=g%stack((j-1)*g%nstackp+i)
      if (ij.ne.0) then
      ijk=g%receiver(ij)
        if (ijk.ne.ij) then
        dhs=t%hp(ij)-t%hb(ij)
        kf=l%kf
        if (dhs.gt.2.d-10) kf=l%kfs
          x=p%dx*(1+mod(ij-1,p%nx))
          y=p%dy*(1+(ij-1)/p%nx)
            do k=1,p%granite_n
            dist=sqrt((x-p%granite_x(k))**2/p%granite_rx(k)**2+(y-p%granite_y(k))**2/p%granite_ry(k)**2)
            if (dist.lt.1.d0.and. &
                t%href(ij)-t%h(ij).ge.p%granite_top(k).and.t%href(ij)-t%h(ij).le.p%granite_bottom(k)) kf=kf*p%granite_dk(k)
            enddo
          if (kf.lt.0.d0) then
          dhdt=(t%h(ij)-p%sea_level)/p%dt
          else
          if (t%length(ij).gt.0.d0.and.t%hp(ij).gt.t%hp(ijk)) then
          dhdt=-(kf*(t%hp(ij)-t%hp(ijk))/t%length(ij)*t%discharge(ij)**l%m-t%sediment(ij)/(l%lf*t%discharge(ij)))
          dh=dhdt*p%dt
            if (dhs.gt.2.d-10 .and. -dh.gt.dhs) then
            dhdt=-dhs/p%dt-(l%kf*(t%hp(ij)-dhs-t%hp(ijk))/t%length(ij)*t%discharge(ij)**l%m-t%sediment(ij)/(l%lf*t%discharge(ij)))
            endif
          else
          dhdt=t%sediment(ij)/(l%lf*t%discharge(ij))
          endif
          endif
        t%sediment(ijk)=t%sediment(ijk)+max(0.d0,t%sediment(ij)-dhdt*p%dx*p%dy)
        endif
      endif
    enddo
  enddo
!$OMP end do nowait
!$OMP end parallel

return

end subroutine compute_sediment_flux
