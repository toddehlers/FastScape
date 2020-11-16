subroutine uplift_rate (t,p,g)

use omp_lib
use definitions

implicit none

type (topography) t
type (topology) g
type (param) p

integer i,j,ij
double precision x,y,time,dh,uplift,xl,yl,val_t,rat,dt
double precision w1,w2,w3,w4,w5,w6,w7,w8,w9,w0
integer i1,i2,i3,i4,i5,i6,i7,i8,i9,i0

! Adding advection parameters
double precision, dimension(:), allocatable :: diag,sup,inf,rhs,res
double precision vx_input(p%nn),vy_input(p%nn),vz_input(p%nn),h_advect(p%nx,p%ny)
double precision vx(p%nx,p%ny),vy(p%nx,p%ny),vz(p%nx,p%ny)
double precision hi_advect(p%nx,p%ny),hi_iso_advect(p%nx,p%ny),href_advect(p%nx,p%ny),hb_advect(p%nx,p%ny)
character cbc*4

! Adding velocity interface
integer :: VelocityInput, FileUnit = 392, LineCounter, ios, NumberOfVelocityFiles
character(255) VelocityFile
character(255), dimension(:), allocatable :: VelocityFiles
double precision, dimension(:), allocatable :: InputTimeStep, TimeStep
double precision TimeSpan, RegionalUplift

time=p%time
xl=p%xl
yl=p%yl
dt=p%dt

if (p%uplift_n.ne.-1) then

  if (p%uplift_nt.eq.0) then
  val_t=1.d0
  elseif (p%time.lt.p%uplift_t(1)) then
  val_t=p%uplift_f(1)
  elseif (p%time.gt.p%uplift_t(p%uplift_nt)) then
  val_t=p%uplift_f(p%uplift_nt)
  else
    do i=1,p%uplift_nt-1
      if ((p%time-p%uplift_t(i))*(p%time-p%uplift_t(i+1)).le.0.d0) then
        if (abs(p%uplift_t(i)-p%uplift_t(i+1)).lt.tiny(x)) then
        val_t=(p%uplift_f(i)+p%uplift_f(i+1))/2.d0
        else
        rat=(p%time-p%uplift_t(i))/(p%uplift_t(i)-p%uplift_t(i+1))
        val_t=p%uplift_f(i)+rat*(p%uplift_f(i+1)-p%uplift_f(i))
        endif
      endif
    enddo
  endif
!$OMP parallel num_threads(p%num_threads) shared(t,p,g,val_t) private (i,j,ij,x,y,uplift,dh)
!$OMP do schedule(dynamic,p%ny/p%num_threads)
  do j=1,p%ny
    do i=1,p%nx
    ij=i+(j-1)*p%nx
      if (.not.g%bc(ij)) then
      x=p%xl*float(i-1)/(p%nx-1)
      y=p%yl*float(j-1)/(p%ny-1)
      call interpolate (x/p%xl*2.d0-1.d0,y/p%yl*2.d0-1.d0, &
                        p%uplift_n,p%uplift_v,uplift)
      dh=uplift*p%dt*val_t
      t%u(ij)=uplift
      t%h(ij)=t%h(ij)+dh
      t%hb(ij)=t%hb(ij)+dh
      t%hi(ij)=t%hi(ij)+dh
      t%href(ij)=t%href(ij)+dh
      if (p%refflex) t%hi_iso(ij)=t%hi_iso(ij)+dh
      endif
    enddo
  enddo
!$OMP end do nowait
!$OMP end parallel

else

!$OMP parallel shared(t,p,g,time,xl,yl,dt)
!$OMP do private (i,j,ij,x,y,uplift,dh,w1,w2,w3,w4,w5,w6,w7,w8,w9,w0,i1,i2,i3,i4,i5,i6,i7,i8,i9,i0) &
!$OMP schedule(dynamic,p%ny/p%num_threads)
  do j=1,p%ny
    do i=1,p%nx
    ij=i+(j-1)*p%nx
      if (.not.g%bc(ij)) then
      x=p%xl*float(i-1)/(p%nx-1)
      y=p%yl*float(j-1)/(p%ny-1)
      !uplift=5.d-3
! user supplied
!user supplied
      t%u(ij)=uplift
      !dh=uplift*p%dt
      !t%h(ij)=t%h(ij)+dh
      !t%hb(ij)=t%hb(ij)+dh
      !t%hi(ij)=t%hi(ij)+dh
      !t%href(ij)=t%href(ij)+dh
      !if (p%refflex) t%hi_iso(ij)=t%hi_iso(ij)+dh
      endif
    enddo
  enddo
!$OMP end do nowait
!$OMP end parallel

! Adding velocity interface (PRE March2020) START 
  if (VelocityInput .eq. 1) then
    VelocityFile = p%run//trim(VelocityFile)
    call LengthOfFile(VelocityFile, LineCounter)
    open(FileUnit,file=VelocityFile,status='old',iostat=ios)
      do i=1,1
        read(FileUnit,*,iostat=ios)
      enddo
      NumberOfVelocityFiles = LineCounter - 1
      allocate (TimeStep(NumberOfVelocityFiles), VelocityFiles(NumberOfVelocityFiles), InputTimeStep(NumberOfVelocityFiles))
      do i = 1, NumberOfVelocityFiles
        read(FileUnit, *, iostat=ios) InputTimestep(i), VelocityFiles(i)
        VelocityFiles(i) = '/input/'//trim(VelocityFiles(i))
        TimeStep(i) = (InputTimeStep(1) - InputTimeStep(i)) * 1e6_8
        call ShortenString(VelocityFiles(i), 4)
      end do
    close(FileUnit)

    ! Identify time step and grids
    do i = 1,(NumberOfVelocityFiles+1)
      if (i == 1) then
        VelocityFile = p%run//trim(VelocityFiles(i))//'.dat'
      else
        VelocityFile = p%run//trim(VelocityFiles(i-1))//'.dat'
      end if
      if (TimeStep(i) >= time) then
        TimeSpan = TimeStep(i) - TimeStep(i-1)
        exit
      end if
    end do
    open(FileUnit, file = VelocityFile, status='old')
      do j=1,p%ny
        do i=1,p%nx
        ij=i+(j-1)*p%nx
          read(FileUnit, *, iostat=ios) vx_input(ij), vy_input(ij), vz_input(ij)
          vx(i,j) = vx_input(ij)
          vy(i,j) = vy_input(ij)
          vz(i,j) = vz_input(ij) + RegionalUplift
        enddo
      enddo
    close(FileUnit)
  end if

  ! Adding Velocity interface (PRE MArch2020) END

! Adding advection (PRE March2019) START
! x-advection using an implicit, second-order scheme to solve advection eq.

allocate (diag(p%nx),sup(p%nx),inf(p%nx),rhs(p%nx),res(p%nx))

write (cbc,'(i4.4)') p%ibc

! Create 2D array
  do j=1,p%ny
    do i=1,p%nx
      ij=i+(j-1)*p%nx
      h_advect(i,j) = t%h(ij)
      hb_advect(i,j) = t%hb(ij)
      hi_advect(i,j) = t%hi(ij)
      href_advect(i,j) = t%href(ij)
      if (p%refflex) hi_iso_advect(i,j) = t%hi_iso(ij)
    enddo
  enddo

  do j=1,p%ny

    diag=1.d0
    sup=0.d0
    inf=0.d0

    do i=1,p%nx
      if (vx(i,j).gt.0.d0) then
        diag(i)=1.d0+vx(i,j)*p%dt/p%dx
        inf(i)=-vx(i,j)*p%dt/p%dx
      elseif (vx(i,j).lt.0.d0) then
        diag(i)=1.d0-vx(i,j)*p%dt/p%dx
        sup(i)=vx(i,j)*p%dt/p%dx
      endif
    enddo
    sup(1)=0.d0
    diag(1)=1.d0
    diag(p%nx)=1.d0
    inf(p%nx)=0.d0

    rhs=h_advect(:,j)
    call tridag (inf,diag,sup,rhs,res,p%nx)
    h_advect(:,j)=res

    rhs=hb_advect(:,j)
    call tridag (inf,diag,sup,rhs,res,p%nx)
    hb_advect(:,j)=res

    rhs=hi_advect(:,j)
    call tridag (inf,diag,sup,rhs,res,p%nx)
    hi_advect(:,j)=res

    rhs=href_advect(:,j)
    call tridag (inf,diag,sup,rhs,res,p%nx)
    href_advect(:,j)=res    

    if (p%refflex) then
    rhs=hi_iso_advect(:,j)
    call tridag (inf,diag,sup,rhs,res,p%nx)
    hi_iso_advect(:,j)=res 
    end if

  enddo

  deallocate (diag,sup,inf,rhs,res)

! y-advection using an implicit, second-order scheme to solve advection eq.

 allocate (diag(p%ny),sup(p%ny),inf(p%ny),rhs(p%ny),res(p%ny))

  do i=1,p%nx

    diag=1.d0
    sup=0.d0
    inf=0.d0

    do j=1,p%ny
      if (vy(i,j).gt.0.d0) then
        diag(j)=1.d0+vy(i,j)*p%dt/p%dy
        inf(j)=-vy(i,j)*p%dt/p%dy
      elseif (vy(i,j).lt.0.d0) then
        diag(j)=1.d0-vy(i,j)*p%dt/p%dy
        sup(j)=vy(i,j)*p%dt/p%dy
      endif
    enddo
    sup(1)=0.d0
    diag(1)=1.d0
    diag(p%ny)=1.d0
    inf(p%ny)=0.d0

    rhs=h_advect(i,:)
    call tridag (inf,diag,sup,rhs,res,p%ny)
    h_advect(i,:)=res

    rhs=hb_advect(i,:)
    call tridag (inf,diag,sup,rhs,res,p%ny)
    hb_advect(i,:)=res

    rhs=hi_advect(i,:)
    call tridag (inf,diag,sup,rhs,res,p%ny)
    hi_advect(i,:)=res

    rhs=href_advect(i,:)
    call tridag (inf,diag,sup,rhs,res,p%ny)
    href_advect(i,:)=res    

    if (p%refflex) then
    rhs=hi_iso_advect(i,:)
    call tridag (inf,diag,sup,rhs,res,p%ny)
    hi_iso_advect(i,:)=res 
    end if

  enddo

  deallocate (diag,sup,inf,rhs,res)

  i1=1
  i2=p%nx
  if (cbc(4:4)=='1') i1=2
  if (cbc(2:2)=='1') i2=p%nx-1
  i3=1
  i4=p%ny
  if (cbc(1:1)=='1') i3=2
  if (cbc(3:3)=='1') i4=p%ny-1

  do j=i3,i4
    do i=i1,i2
    ij=i+(j-1)*p%nx
    h_advect(i,j)=h_advect(i,j)+vz(i,j)*p%dt
    hb_advect(i,j)=hb_advect(i,j)+vz(i,j)*p%dt
    hi_advect(i,j)=hi_advect(i,j)+vz(i,j)*p%dt
    href_advect(i,j)=href_advect(i,j)+vz(i,j)*p%dt
    if (p%refflex) hi_iso_advect(i,j)=hi_iso_advect(i,j)+vz(i,j)*p%dt
    enddo
  enddo

! Create 1D array
  do j=1,p%ny
    do i=1,p%nx
      ij=i+(j-1)*p%nx
      t%h(ij) = h_advect(i,j)
      t%hb(ij) = hb_advect(i,j)
      t%hb(ij)=min(t%h(ij),t%hb(ij))
      t%hi(ij) = hi_advect(i,j)
      t%href(ij) = href_advect(i,j)
	    if (p%refflex) t%hi_iso(ij) = hi_iso_advect(i,j)
      t%u(ij)=vz(i,j)
    enddo
  enddo

! Adding advection (PRE March2019) END

endif

return

end subroutine uplift_rate