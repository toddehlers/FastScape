subroutine Chi (g,t,l,p)

use omp_lib
use definitions

implicit none

type (topography) t
type (topology) g
type (param) p
type (law) l

double precision discharge0,uplift0,kf0
double precision fact,kf,dhs,find_kf
integer i,j,ij,ijk
external find_kf

discharge0=sum(t%discharge(1:p%nn))/p%nn
uplift0=sum(t%u(1:p%nn))/p%nn
kf0=l%kf
fact=uplift0/kf0/discharge0**l%m

t%chi=0.d0

!$OMP parallel shared(g,t,p,fact)
!$OMP do private(i,j,ij,ijk,kf,dhs)
  do j=1,p%num_threads
    do i=1,g%nstackp
    ij=g%stack((j-1)*g%nstackp+i)
      if (ij.ne.0) then
        if (.not.g%bc(ij).and.t%h(ij).gt.p%sea_level) then
        ijk=g%receiver(ij)
          if (ijk.ne.ij) then
            if (p%plot_chi.eq.1) then
              if (t%discharge(ij).ge.p%discharge_min_chi) then
              t%chi(ij)=t%chi(ijk)+t%length(ij)*(discharge0/t%discharge(ij))**(l%m/l%n)
              endif
            elseif (p%plot_chi.eq.2) then
            kf=find_kf(ij,l,t,p)
            dhs=t%h(ij)-t%hb(ij)
            if (dhs.gt.2.d-10) kf=l%kfs
              if (kf.gt.0.d0.and.t%discharge(ij).ge.p%discharge_min_chi) then
              t%chi(ij)=t%chi(ijk)+t%length(ij)*(t%u(ij)/kf/t%discharge(ij)**l%m*fact)**(1.d0/l%n)
              endif
            else
            stop 'this value for chi is not implemented yet'
            endif
          endif
        endif
      endif
    enddo
  enddo
!$OMP end do nowait
!$OMP end parallel

t%chi=t%chi/maxval(t%chi)

return

end subroutine Chi
