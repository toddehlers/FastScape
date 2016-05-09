subroutine vtk (h,nx,ny,dx,dy,istep)

! sibrutine to create a simple VTK file for plotting

implicit none

double precision h(nx,ny),dx,dy
integer nx,ny,istep

integer iunit,i,j,ij
character*6 cstep

iunit=74

write (cstep,'(i6)') istep
if (istep.lt.10) cstep(1:5)='00000'
if (istep.lt.100) cstep(1:4)='0000'
if (istep.lt.1000) cstep(1:3)='000'
if (istep.lt.10000) cstep(1:2)='00'
if (istep.lt.100000) cstep(1:1)='0'


open(unit=iunit,file='VTK/Topo'//cstep//'.vtk')
write(iunit,'(a)')'# vtk DataFile Version 3.0'
write(iunit,'(a)')'Surface'
write(iunit,'(a)')'ASCII'
write(iunit,'(a)')'DATASET STRUCTURED_POINTS'
write(iunit,'(a10,3i10)')'DIMENSIONS',nx,ny,1
write(iunit,'(a10,3f15.5)')'ORIGIN',0.d0,0.d0,0.d0
write(iunit,'(a10,3f15.5)')'SPACING',dx,dy,1.d0

write(iunit,'(a11,i10)')'POINT_DATA ',nx*ny

write(iunit,'(a)')'SCALARS H float 1'
write(iunit,'(a)')'LOOKUP_TABLE default'
  do j=1,ny
    do i=1,nx
    write(iunit,'(f9.2)') h(i,j)
    enddo
  enddo

close (iunit)
return
end subroutine vtk
