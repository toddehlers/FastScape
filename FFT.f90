subroutine FFT (t,p)

use omp_lib
use definitions

implicit none

type (param) p
type (topography) t

real*8, dimension(:,:), allocatable :: flex,work
integer i,j,ii,jj,nflex

write (555,*) t%h

nflex=2
do i=1,15
if (nflex.ge.p%nx .and. nflex.ge.p%ny) goto 111
nflex=nflex*2
enddo
111 continue

print*,nflex

allocate (flex(nflex,nflex),work(nflex,nflex))

!$OMP parallel shared(flex,nflex,t,p) private(ii,jj,i,j)
!$OMP do schedule(dynamic)
do j=1,nflex
jj=1+((j-1)*(p%ny-1))/(nflex-1)
do i=1,nflex
ii=1+((i-1)*(p%nx-1))/(nflex-1)
flex(i,j)=t%h(ii+(jj-1)*p%nx)
enddo
enddo
!$OMP end do nowait
!$OMP end parallel

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

!print*,'work',minval(work),maxval(work)

!$OMP parallel shared(flex,nflex,t,p) private(ii,jj,i,j)
!$OMP do schedule(dynamic)
do j=1,p%ny
jj=1+((j-1)*(nflex-1))/(p%ny-1)
do i=1,p%nx
ii=1+((i-1)*(nflex-1))/(p%nx-1)
t%fft(i+(j-1)*p%nx)=work(ii,jj)*4.d0/p%xl/p%yl
enddo
enddo 
!$OMP end do nowait
!$OMP end parallel

!print*,'fft',minval(t%fft),maxval(t%fft)

end subroutine FFT
