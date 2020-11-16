program Extract_DEM

! simple program to read topography data from a world database in a given area
! uses GTOPO1 at approximately 1 km resolution

implicit none

double precision latmin,latmax,lonmin,lonmax
double precision con,dx,dy
integer nx,ny
double precision, dimension(:,:), allocatable :: h0

!lat-lon ranges for the model

latmin=20.d0
latmax=90.d0
lonmin=-130.d0
lonmax=-70.d0

con=atan(1.d0)/90.d0
dy=1.d0/60.d0
dx=dy

! readin geometry for topofile

call find_size (lonmin,lonmax,latmin,latmax,nx,ny)
print*,'Topo dataset size is ',nx,ny

! readin topo

allocate (h0(nx,ny))
call extract_topo (lonmin,lonmax,latmin,latmax,h0,nx,ny)
print*,'Topo range is ',minval(h0),sum(h0)/nx/ny,maxval(h0)

! create topo file

open (7,file='DEM.bin',status='unknown',form='unformatted',access='direct', &
      recl=2*(nx*ny))
write (7,rec=1) int2(h0)
close (7)

end program Extract_DEM

!----------------------

subroutine find_size (lonmin,lonmax,latmin,latmax,nx,ny)

implicit none

double precision lonmin,lonmax,latmin,latmax
integer nx,ny,minx,maxx,miny,maxy

minx=int((lonmin+180.d0)*60.d0)
maxx=int((lonmax+180.d0)*60.d0)
nx=maxx-minx+1

miny=int((90.d0-latmax)*60.d0)
maxy=int((90.d0-latmin)*60.d0)
ny=maxy-miny+1

return

end subroutine find_size

!---------------------

subroutine extract_topo (lonmin,lonmax,latmin,latmax,h,nx,ny)

implicit none

integer nx,ny,minx,maxx,miny,maxy
double precision h(nx,ny),lonmin,lonmax,latmin,latmax
integer nx0,ny0
integer*2, dimension(:,:), allocatable :: h2

nx0=21600
ny0=10800
allocate (h2(nx0,ny0))
open (77,file="etopo1_bed_c_i2.bin",status="old",form="unformatted", &
      access="direct",recl=nx0*ny0*2,convert='little_endian')
read (77,rec=1) h2
close (77)

minx=int((lonmin+180.d0)*60.d0)
maxx=int((lonmax+180.d0)*60.d0)
miny=int((90.d0-latmax)*60.d0)
maxy=int((90.d0-latmin)*60.d0)

h=h2(minx:maxx,maxy:miny:-1)
deallocate (h2)

return

end subroutine extract_topo
