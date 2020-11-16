subroutine orography (t,p)

! not tested, not vectorized, in development

use omp_lib
use definitions

implicit none

type (param) p
type (topography) t

double precision, dimension(:,:), allocatable :: h,pr,w
complex*16, dimension(:,:), allocatable :: work
complex*16, dimension(:), allocatable :: d

integer nx,ny,i,j,ij
double precision dist,cw,u,v,Nm,tauf,taus,hw,pm,lx,ly,sigma,hmax
double precision k,l,mkl,Nmsigma
complex*16 im

! grid dimension
nx=p%noro
ny=p%noro
allocate (h(nx,ny),pr(nx,ny))
lx=p%xl
ly=p%yl

!cw is uplift sensitivity factor (0.001-0.02 kg/m3)
cw=0.02d0
cw=p%cw
!u is E-W wind velocity
u=1.d0
u=p%u
!v is S-N wind velocity
v=1.d0
v=p%v
!Nm is moist stability frequency (0-0.01 1/s)
Nm=0.005d0
Nm=p%nm
!taus is conversion time (200-2000 s)
taus=500.d0
taus=p%taus
!tauf is fallout time (200-2000 s)
tauf=500.d0
tauf=p%tauf
!pm is mean precipitation (m/s)
pm=3.d-3/3600.d0
pm=3.d-2/3600.d0
pm=p%pm
!hw is Water vapour scale height (1-5 km)
hw=2000.d0
hw=p%hw
!im is imaginary unit
im=cmplx(0.d0,1.d0)

allocate (w(p%nx,p%ny))
  do j=1,p%ny
    do i=1,p%nx
    ij=i+(j-1)*p%nx
    w(i,j)=t%h(ij)
    enddo
  enddo
call InterpolSmooth (w,p%nx,p%ny,h,nx,ny)
deallocate (w)

allocate (work(ny,nx))

! FFT

allocate (d(nx))
  do j=1,ny
    do i=1,nx
    d(i)=dcmplx(h(i,j),0.d0)
    enddo
  call four1 (d,nx,-1)
    do i=1,nx
    work(j,i)=d(i)
    enddo
  enddo
deallocate (d)

allocate (d(ny))
  do i=1,nx
    do j=1,ny
    d(j)=work(j,i)
    enddo
  call four1 (d,ny,-1)
    do j=1,ny
    work(j,i)=d(j)
    enddo
  enddo 
deallocate (d)

  do j=1,ny
  do i=1,nx
  work(j,i)=work(j,i)/nx/ny
  enddo
  enddo

!Orographic filter

  do j=1,ny
   do i=1,nx
   k=float(i-1)/lx
   if (i.gt.nx/2) k=-float(nx-i+1)/lx
   l=float(j-1)/ly
   if (j.gt.ny/2) l=-float(ny-j+1)/ly
   sigma=u*k+v*l
   mkl=0.d0
   Nmsigma=Nm**2-sigma**2
   if (abs(sigma).gt.tiny(sigma).and.Nmsigma.gt.tiny(Nmsigma)) mkl=sqrt((Nmsigma/sigma**2)*(k**2+l**2))
   work(j,i)=Cw*im*sigma*work(j,i)/(1.d0-im*mkl*Hw)/(1.d0+im*sigma*taus)/(1.d0+im*sigma*tauf)
   enddo
 enddo

!Inverse FFT

allocate (d(ny))
  do i=1,nx
    do j=1,ny
    d(j)=work(j,i)
    enddo
  call four1 (d,ny,1)
    do j=1,ny
    work(j,i)=d(j)
    enddo
  enddo
deallocate (d)

allocate (d(nx))
  do j=1,ny
    do i=1,nx
    d(i)=work(j,i)
    enddo
  call four1 (d,nx,1)
    do i=1,nx
    pr(i,j)=real(d(i))
    enddo
  enddo
deallocate (d)

deallocate (work)

pr=pr*3600.d0*24.d0*365.25d0/1.d3

pr=pr+pm
pr=max(pr,pm)

allocate (w(p%nx,p%ny))
call InterpolLinear (pr,nx,ny,w,p%nx,p%ny)
  do j=1,p%ny
    do i=1,p%nx
    ij=i+(j-1)*p%nx
    t%discharge(ij)=w(i,j)*p%dx*p%dy
    enddo
  enddo
deallocate (w)

deallocate (h,pr)

end subroutine orography
