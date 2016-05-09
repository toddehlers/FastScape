program run_FastScapeSolo

! example of how to call FastScapeSolo

implicit none

integer nx,ny,nstep,ibc
double precision, dimension(:,:), allocatable :: z,vx,vy,vz,kf,kd,p
double precision dx,dy,dt,m,n
integer i,j
double precision pi,x0,y0,x,y,rad
character example*10

nx=301
ny=301
dx=333.d0
dy=333.d0
dt=10000.d0
nstep=1000
m=0.4d0
n=1.d0

allocate (z(nx,ny),vx(nx,ny),vy(nx,ny),vz(nx,ny),kf(nx,ny),kd(nx,ny),p(nx,ny))

pi=atan(1.d0)*4.d0
call random_number (z)

example="NewZealand"

select case (example)

case ("advection ")
! advection towards a corner
vx=5.d-3
vy=5.d-3
vz=2.d-3
kf=1.d-4
ibc=1111

case ("swirl     ")
!swirl
  do j=1,ny
    do i=1,nx
    vx(i,j)=1.d-2*sin(dfloat(i-1)/(nx-1)*pi)*cos(dfloat(j-1)/(ny-1)*pi)
    vy(i,j)=-1.d-2*cos(dfloat(i-1)/(nx-1)*pi)*sin(dfloat(j-1)/(ny-1)*pi)
    enddo
  enddo
vz=2.d-3
kf=1.d-4
ibc=1111

case ("NewZealand")
! new zealand
  do j=1,ny
    do i=1,nx
      if (i.gt.nx/3) then
      vx(i,j)=-1.d-2
      vy(i,j)=-4.d-2*dfloat(i-nx/3)/(2*nx/3)
      vz(i,j)=1.d-2/2.d0*dfloat(nx-i)/(2*nx/3)
      else
      vx(i,j)=0.d0
      vy(i,j)=0.d0
      vz(i,j)=0.d0
      endif
    enddo
  enddo
kf=1.d-4
ibc=0101

case ("SpinVolcan")
! spinning volcano
x0=dx*nx/2.d0
y0=dy*ny/2.d0
rad=dx*ny/3.d0
  do j=1,ny
  y=dy*float(j-1)-y0
    do i=1,nx
    x=dx*float(i-1)-x0
      if (sqrt(x**2+y**2).lt.rad) then
      vx(i,j)=-y/rad*2.d-2
      vy(i,j)=x/rad*2.d-2
      else
      vx(i,j)=0.d0
      vy(i,j)=0.d0
      endif
    vz(i,j)=3.d-3
    enddo
  enddo
kf=3d-5
ibc=1111

case default
vx=0.d0
vy=0.d0
vz=2.d-3
kf=1.d-4
ibc=1111
end select

kd=2.d0
p=1.d0

do i=1,nstep/10
call FastScapeSolo (z,nx,ny,dx,dy,dt,10,vx,vy,vz,kf,kd,p,m,n,ibc)
write (*,'(i5,3f8.1)') i,minval(z),sum(z)/(nx*ny),maxval(z)
call vtk (z,nx,ny,dx,dy,i)
enddo

end
