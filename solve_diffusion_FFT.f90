!---------------------

subroutine solve_diffusion_FFT (x,nx,ny,dt,kappa,xl,yl,num_threads)

! subroutine to compute the solution to the diffusion equation using FFT
! beware that the method only works for uniform diffusivity

! in input:
!    x(n,n): double precision array of dimension (n,n) containing the initial distribution
!    n: integer, dimension of array x; must be a power of 2
!    dt: double precision time step length/increment
!    kappa: double precision diffusivity
!    l: dimension of the problem

! in output:
!    x(n,n): contains the final distribution

implicit none

double precision dt,kappa,xl,yl
integer nx,ny,nstep,i,j,num_threads
double precision x(nx,ny),rx,ry,r0x,r0y
double precision, dimension(:,:), allocatable :: w
double precision pi,rax,ray,rbx,rby,wi,wj,hk,ropt

nstep=100

j=1
do i=1,20
j=j*2
if (nx.eq.j) goto 999
enddo
print*,'solve_FFT only works for arrays dimensioned as powers of 2...'
print*,'Your nx value is ',nx
stop
999 continue

j=1
do i=1,20
j=j*2
if (ny.eq.j) goto 998
enddo
print*,'solve_FFT only works for arrays dimensioned as powers of 2...'
print*,'Your ny value is ',ny
stop
998 continue

if (dt.le.0.d0) then
print*,'dt must be positive number'
print*,'Your value is ',dt
endif

if (kappa.le.0.d0) then
print*,'kappa must be positive number'
print*,'Your value is ',kappa
endif

if (xl.le.0.d0) then
print*,'xl must be positive number'
print*,'your value is ',xl
endif

if (yl.le.0.d0) then
print*,'yl must be positive number'
print*,'your value is ',yl
endif

rx=dt*kappa*nx*nx/xl/xl

if (rx.gt.100.d0) then
print*,'your time step is too large'
print*,'it should not be greater than ',sqrt(5.d0)/2.d0/kappa/nx/nx/xl/xl*100.d0
stop
endif

ry=dt*kappa*ny*ny/yl/yl

if (ry.gt.100.d0) then
print*,'your time step is too large'
print*,'it should not be greater than ',sqrt(5.d0)/2.d0/kappa/ny/ny/yl/yl*100.d0
stop
endif

rx=rx/nstep
ry=ry/nstep

pi=atan(1.d0)*4.d0

rax=1.d0/12.d0+rx/2.d0
rbx=1.d0/12.d0-rx/2.d0
ray=1.d0/12.d0+ry/2.d0
rby=1.d0/12.d0-ry/2.d0

allocate (w(ny,nx))

!$OMP parallel shared(x,nx,ny) private(j)
!$OMP do schedule(dynamic,ny/num_threads)
do j=1,ny
call sinft (x(:,j),nx)
enddo
!$OMP end do nowait
!$OMP end parallel

!$OMP parallel shared(w,x)
!$OMP workshare
w=transpose(x)
!$OMP end workshare nowait
!$OMP end parallel

!$OMP parallel shared(w,nx,ny) private(i)
!$OMP do schedule(dynamic,nx/num_threads)
do i=1,nx
call sinft (w(:,i),ny)
enddo
!$OMP end do nowait
!$OMP end parallel

!$OMP parallel shared(w,nx,ny)
!$OMP workshare
w=w*4./nx/ny
!$OMP end workshare nowait
!$OMP end parallel

!$OMP parallel shared(w,nx,ny,pi,nstep,rax,rbx,ray,rby) private(i,j,wi,wj,hk)
!$OMP do schedule(dynamic,ny/num_threads)
do j=1,ny
wj=pi*j/ny
wj=pi*(j-1)/(ny)
  do i=1,nx
  wi=pi*i/nx
  wi=pi*(i-1)/(nx)
  hk=(1.d0-2.d0*rax*(1.d0-cos(wi)))/(1.d0-2.d0*rbx*(1.d0-cos(wi)))*(1.d0-2.d0*ray*(1.d0-cos(wj)))/(1.d0-2.d0*rby*(1.d0-cos(wj)))
  w(j,i)=w(j,i)*hk**nstep
  enddo
enddo
!$OMP end do nowait
!$OMP end parallel

!$OMP parallel shared(w,nx,ny) private(i)
!$OMP do schedule(dynamic,nx/num_threads)
do i=1,nx
call sinft (w(:,i),ny)
enddo
!$OMP end do nowait
!$OMP end parallel

!$OMP parallel shared(w,x)
!$OMP workshare
x=transpose(w)
!$OMP end workshare nowait
!$OMP end parallel

!$OMP parallel shared(w,nx,ny) private(j)
!$OMP do schedule(dynamic,ny/num_threads)
do j=1,ny
call sinft (x(:,j),nx)
enddo
!$OMP end do nowait
!$OMP end parallel

deallocate (w)

return
end subroutine solve_diffusion_FFT
