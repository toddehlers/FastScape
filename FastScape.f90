program FastScape

! FastScape is a Landscape Evolution Model written by Jean BRAUN
! starting in early 2012
! It solves the basic fluvial incision equation (or stream power law)
! on a very high resolution regular mesh

! It is O(N) (number of operations strictly proportional to, N, the number of nodes
! used to discretize the landscape), implicit in time and parallel (using OpenMP)

! The author intends to make it evolve in the very near future so look for updates

use omp_lib
use definitions

implicit none

! for the definition of these derived types have a look at module_definitions.f90

type (topography) t
type (topology) g
type (param) p
type (law) l
type (box), dimension(:), allocatable :: b

integer,dimension(:), allocatable :: stackp,ierror
double precision, dimension(:), allocatable :: strati,dhdt,dudt,plot,al,sed_flux,tinit
integer nb_cycles_sec,nb_cycles_max,nb_cycles_initial,nb_cycles_final,nb_cycles
integer i,j,ij,ii,jj,ijk,k,icatch,nstack,numarg,istep,istep0,iii,jjj
integer chunkny,chunknn,num,nk,nfreq_screen,nsum,nsump
real temps_elapsed,tin,tdt
double precision slopemax,slopep,lengthp,x,y
double precision fact,dh,flux_in,flux_out,flux_true,sr,sumh,minh,maxh
double precision dist,mean,var,skew,kurt
character cstep*7,c5*5
logical first,second,inside
logical, dimension(:), allocatable :: coast
double precision, dimension(:), allocatable :: temperatureH,aft_age,zft_age,ahe_age,zhe_age
double precision, dimension(:), allocatable :: SurfaceGradientToday,timeH
double precision, dimension(:,:), allocatable :: eRate
integer, dimension(:), allocatable :: ij_age
character*1024, dimension(:), allocatable :: title
real, dimension(:,:), allocatable :: field
integer nf,nfmax

external inside

! uses RUN00 as default run name, unless a name is specified in the command line

numarg=command_argument_count()
if (numarg.eq.0) then
p%run='RUN00'
else
call getarg (1,p%run)
endif

! checks that the run directory exists, otherwise exits

open (7,file=p%run//'/FastScape.in',status='unknown',err=999)
goto 998
999 print*,'Cant open/create '//p%run//'/FastScape.in'
print*,'Check that the directory '//p%run//' exists'
stop
998 continue

p%infile=p%run//'/FastScape.in'

! reads the input file (or creates it if it does not exist)

call initialize (p,l)

! checks that number of requested threads is not larger than maximum allowed by system/hardware

if (p%num_threads.gt.omp_get_max_threads()) then
print*,'Maximum number of threads is ',omp_get_max_threads()
stop 'too many threads'
endif
call omp_set_num_threads (p%num_threads)
allocate (ierror(p%num_threads))
p%chunk=p%nn/p%num_threads

! resets the system clock

call system_clock (count_rate=nb_cycles_sec,count_max=nb_cycles_max)
call system_clock (count=nb_cycles_initial)
call system_clock (count=nb_cycles_final)
nb_cycles=nb_cycles_final - nb_cycles_initial
IF (nb_cycles_final.lt.nb_cycles_initial) &
        nb_cycles=nb_cycles+nb_cycles_max
temps_elapsed=real(nb_cycles)/nb_cycles_sec

! this is the length of the stack array that needs to be slightly large than N
! as the catchment topology might not allow for a perfect division of the nodes
! by the number of threads

chunknn=p%nn/p%num_threads/2
chunkny=chunkny/p%num_threads/2

! allocates (most of the) arrays

allocate (t%length(p%nn),g%receiver(p%nn),g%bc(p%nn),t%h(p%nn),t%hp(p%nn))
allocate (t%hiso(p%nn),t%hi(p%nn),t%pmask(1,p%nn),t%hi_iso(p%nn),t%href(p%nn))
allocate (g%ndon(p%nn+1),g%donors(p%nn),g%catchment(p%nn),t%discharge(p%nn))
allocate (t%width(p%nn),t%sediment(p%nn),t%hb(p%nn),t%u(p%nn),t%s(p%nn))
allocate (sed_flux(p%nn),coast(p%nn))
nfmax=20
if (p%vtk) allocate(field(p%nn,nfmax),title(nfmax+1))
g%nstackp=p%nn/p%num_threads+p%nn/p%num_threads/5
g%nstackp=p%nn

allocate (g%stack(g%nstackp*p%num_threads))

second=.false.

99999 continue
if (p%nage.gt.0) then
if (.not.second) then
allocate (timeH(p%nstep),eRate(p%nstep,p%nage),ij_age(p%nage))
  do i=1,p%nage
  ii=p%x_age(i)/p%xl*(p%nx-1)+1
  ii=min(ii,p%nx)
  jj=p%y_age(i)/p%yl*(p%ny-1)+1
  jj=min(jj,p%ny)
  ij_age(i)=(jj-1)*p%nx+ii
  ij_age(i)=min(ij_age(i),p%nn)
  ij_age(i)=max(ij_age(i),1)
  enddo
endif
endif

! set/reset/initialize variables where needed

p%time=0.d0
p%istep=0
nfreq_screen=19
t%s=0.d0

! prepares for Demoulin store if necessary

if (p%restart.ne.1.and.p%ndemoulin.ne.0) then
open (61,file=p%run//'/Demoulin.txt',status='unknown')
write (61,'(a7)') 'Time Sr'
close (61)
endif

! reads in or calculates the initial topography

if (p%debug) call timer (tin,0,nb_cycles_initial,'doing initial topography')
call initial_topography (t,p,g)
if (second) then
t%h=tinit
t%hi=t%h
t%href=t%h
t%hi_iso=t%h
t%hb=t%h
t%hp=t%h
t%hiso=0.d0
deallocate (tinit)
if (p%reference_surface.eq.0) t%href=0.
  do i=1,p%nbox
  b(i)%volume_predicted=0.d0
  enddo
else
allocate (tinit(p%nn))
tinit=t%h
endif
if (p%debug) call timer (tin,1,nb_cycles_initial,'initial topography done')

! imposes/defines the bonudary conditions

call set_bc (g,p,t)

! if fluxes are to be saved/stored the corresponding file is created

write (c5,'(i5)') p%imetric
if (c5(4:4).eq.'1') open (91,file=p%run//'/Slopes.txt',status='unknown')
if (c5(5:5).eq.'1') open (89,file=p%run//'/Fluxes.txt',status='unknown')

if (.not.second) p%nbox=0
if (c5(5:5).eq.'1' .and. .not.second) then
call system ('cd '//p%run//'; ls -1 Volume?* | wc > lsVolume')
open (7,file=p%run//'/lsVolume',status='unknown')
read (7,*,err=993) p%nbox
993 close (7)
call system ('rm '//p%run//'/lsVolume')
allocate (b(p%nbox))
  do i=1,p%nbox
  allocate (b(i)%mask(p%nn))
  b(i)%volume_predicted=0.d0
  enddo
call system ('cd '//p%run//'; ls -1 Volume?* > lsVolume')
open (7,file=p%run//'/lsVolume',status='old')
  do i=1,p%nbox
  read (7,'(a)') b(i)%fnme
  open (8,file=p%run//'/'//trim(b(i)%fnme),status='old')
  read (8,*) (b(i)%x(j),b(i)%y(j),j=1,4)
    do j=1,1000
    read (8,*,end=991) b(i)%timei(j),b(i)%timef(j),b(i)%volume(j),b(i)%dvolume(j)
    b(i)%ntime=j
    enddo
991 close (8)
  b(i)%mask=.false.
    do jj=1,p%ny
      do ii=1,p%nx
      ij=ii+(jj-1)*p%nx
      if (inside(p%xl*float(ii-1)/(p%nx-1),p%yl*float(jj-1)/(p%ny-1),b(i)%x,b(i)%y,4)) b(i)%mask(ij)=.true.
      enddo
    enddo
  enddo
close (7)
call system ('rm '//p%run//'/lsVolume')
endif

! beginning of time stepping

istep0=p%istep+1
flux_true=0.d0

do istep=istep0,p%nstep

call timer (tdt,0,nb_cycles_initial,'$')

p%istep=istep

p%time=p%time+p%dt
if (p%nage.gt.0) timeH(istep)=p%dt*p%nstep-p%time

write (cstep,'(i7)') p%istep
if (p%istep.lt.1000000) cstep(1:1)='0'
if (p%istep.lt.100000) cstep(1:2)='00'
if (p%istep.lt.10000) cstep(1:3)='000'
if (p%istep.lt.1000) cstep(1:4)='0000'
if (p%istep.lt.100) cstep(1:5)='00000'
if (p%istep.lt.10) cstep(1:6)='000000'

if (p%debug) print*,'Step ',p%istep
if (p%debug) print*,'------------'

if ((p%istep/p%nfreq)*p%nfreq.eq.p%istep .and. p%vtk) &
  call write_vtk (p%run//'/Step'//cstep//'.vtk',21,p,t%h,0,field,title,nf,nfmax)

!$OMP parallel shared(t,p) private(ij)
!$OMP do schedule(dynamic,p%chunk)
  do ij=1,p%nn
  t%hp(ij)=t%h(ij)
  enddo
!$OMP end do nowait
!$OMP end parallel

! initializes surface uplift rate
allocate (dudt(p%nn))
!$OMP parallel shared(dudt,t,p) private(ij)
!$OMP do schedule(dynamic,p%chunk)
  do ij=1,p%nn
  dudt(ij)=t%h(ij)
  enddo
!$OMP end do nowait
!$OMP end parallel

! calculates uplift rate and updates topographic height accordingly

if (p%debug) call timer (tin,0,nb_cycles_initial,'doing uplift')

flux_in=0.d0
!$OMP parallel shared(p,t)
!$OMP do private(ij) reduction(+:flux_in) schedule(dynamic,p%chunk)
do ij=1,p%nn
flux_in=flux_in+t%h(ij)
enddo
!$OMP end do nowait
!$OMP end parallel

call uplift_rate (t,p,g)

flux_out=0.d0
!$OMP parallel shared(p,t)
!$OMP do private(ij) reduction(+:flux_out) schedule(dynamic,p%chunk)
do ij=1,p%nn
flux_out=flux_out+t%h(ij)
enddo
!$OMP end do nowait
!$OMP end parallel

!$OMP workshare
flux_in=(flux_out-flux_in)/p%dt*p%dx*p%dy
!$OMP end workshare nowait

if (p%debug) call timer (tin,1,nb_cycles_initial,'uplift done')

if (p%Pecube_nfreq.ne.0.and.p%istep.eq.1) then
if (p%debug) call timer (tin,0,nb_cycles_initial,'writing Pecube output')
p%nPecube=0
call Interface_Pecube (t,p,p%nPecube)
if (p%debug) call timer (tin,1,nb_cycles_initial,'Pecube output written')
endif

!$OMP parallel shared(t,p) private(ij)
!$OMP do schedule(dynamic,p%chunk)
  do ij=1,p%nn
  t%length(ij)=0.d0
  enddo
!$OMP end do nowait
!$OMP end parallel

!$OMP parallel shared(g,p)
!$OMP do private (ij)  schedule(dynamic,p%chunk)
  do ij=1,p%nn
  g%receiver(ij)=ij
  enddo
!$OMP end do nowait
!$OMP end parallel

! finds receiver to each node and calculates the length between each node and its receiver

if (p%debug) call timer (tin,0,nb_cycles_initial,'finding receiver')
call find_receiver (g,p,t)
if (p%debug) call timer (tin,1,nb_cycles_initial,'receiver found')

! first is a flag that is used when the local minima algorithm is used and may require to iterate

first=.true.

111 continue

! from the receiver information builds the donor information

if (p%debug) call timer (tin,0,nb_cycles_initial,'finding ndon')
call find_donors (g,p)
if (p%debug) call timer (tin,1,nb_cycles_initial,'ndon found')

! builds the stack

if (p%debug) call timer (tin,0,nb_cycles_initial,'building stack')

1101 continue

nstack=0

!$OMP parallel shared(g,p) private(ij)
!$OMP do schedule(dynamic,g%nstackp)
  do ij=1,g%nstackp*p%num_threads
  g%stack(ij)=0
  enddo
!$OMP end do nowait
!$OMP end parallel

ierror=0
!$OMP parallel shared(p,g,ierror) firstprivate(nstack) private(stackp,num)
allocate (stackp(g%nstackp))
num=omp_get_thread_num()
!$OMP do private (ij) schedule(dynamic,g%nstackp/p%num_threads)
  do ij=1,g%nstackp
  stackp(ij)=0
  enddo
!$OMP end do nowait
!$OMP do private (ij,icatch) schedule(dynamic,p%chunk)
  do ij=1,p%nn
     if (g%receiver(ij).eq.ij) then
     icatch=ij
     if (.not.g%bc(ij)) icatch=-ij
     g%catchment(ij)=icatch
     call add_to_stack (ij,nstack,g%nstackp,stackp,g%ndon,g%donors,g%catchment,p%nn,icatch,ierror(num+1))
     endif
  enddo
!$OMP end do nowait
g%stack(num*g%nstackp+1:num*g%nstackp+nstack)=stackp(1:nstack)
deallocate (stackp)
!$OMP end parallel

if (sum(ierror).ne.0) then
deallocate (g%stack)
g%nstackp=g%nstackp+p%nn/p%num_threads/5
allocate (g%stack(g%nstackp*p%num_threads))
goto 1101
endif

if (p%debug) call timer (tin,1,nb_cycles_initial,'stack built')

! check if local minima have to be removed and if an iteration on finding the donors and
! receivers is needed

if (p%local_minima.and.first) then
if (p%debug) call timer (tin,0,nb_cycles_initial,'finding local_minima')
call local_minima (t,g,p)
if (p%debug) call timer (tin,1,nb_cycles_initial,'local minima found')
first=.false.
goto 111
endif

nstack=g%nstackp*p%num_threads

! calculates precipitation

if (p%debug) call timer (tin,0,nb_cycles_initial,'calculating precipitation')
  if (p%orography.eq.0) then
  call precipitation_function (t,p)
  else
  call orography (t,p)
  endif
if (p%debug) call timer (tin,1,nb_cycles_initial,'precipitation calculated')

! stores precipitation and plots it if requested

  if ((p%istep/p%nfreq)*p%nfreq.eq.p%istep .and. p%plot_precip.ne.0) then
  allocate (plot(p%nx*p%ny))
!$OMP parallel shared(t,plot,p) private(ij)
!$OMP do schedule(dynamic,p%chunk)
  do ij=1,p%nn
  plot(ij)=t%discharge(ij)/p%dx/p%dy
  enddo
!$OMP end do nowait
!$OMP end parallel
  if (p%plot_precip.eq.1) call bmp (t%h,plot,p%nx,p%ny, &
                     p%run//'/Precip'//cstep//'.bmp',23,256,0, &
                     p%precip_min,p%precip_max,1,p%vex)
  if (p%plot_precip.eq.-1) call bmp (t%h,plot,p%nx,p%ny, &
                     p%run//'/Precip'//cstep//'.bmp',23,256,0, &
                     p%precip_min,p%precip_max,0,p%vex)
  if (p%vtk) call write_vtk ('Discharge',9,p,plot,1,field,title,nf,nfmax)
  deallocate (plot)
  endif

! computes discharge

if (p%debug) call timer (tin,0,nb_cycles_initial,'calculating discharge')
!$OMP parallel num_threads(p%num_threads) shared(t,g,p) private(i,j,ij,ijk)
!$OMP do schedule(dynamic,1)
  do j=1,p%num_threads
    do i=g%nstackp,1,-1
    ij=g%stack((j-1)*g%nstackp+i)
      if (ij.ne.0) then
      ijk=g%receiver(ij)
        if (ijk.ne.ij) then
        t%discharge(ijk)=t%discharge(ijk)+t%discharge(ij)
        endif
      endif
    enddo
  enddo
!$OMP end do nowait
!$OMP end parallel
if (p%debug) call timer (tin,1,nb_cycles_initial,'discharge calculated')

! computes channel width (for some geomorphic laws only)

!!$OMP parallel num_threads(p%num_threads) shared(t,fact) private(ij)
!!$OMP do schedule(dynamic,p%chunk)
!  do ij=1,p%nn
!  t%width(ij)=10.d0*sqrt(t%discharge(ij)/fact+0.01d0)
!  enddo
!!$OMP end do nowait
!!$OMP end parallel
!if (p%debug) call timer (tin,1,nb_cycles_initial,'width calculated')

! and adjust discharge accordingly §if necessary)

if (l%law.eq.3) then
if (p%debug) call timer (tin,0,nb_cycles_initial,'adjusting discharge')
!!$OMP parallel shared(t,p) private(ij)
!!$OMP do schedule(dynamic,p%chunk)
!!  do ij=1,p%nn
!  t%discharge(ij)=t%discharge(ij)/t%width(ij)
!  enddo
!!$OMP end do nowait
!!$OMP end parallel
if (p%debug) call timer (tin,1,nb_cycles_initial,'discharge adjusted')
endif

if (p%debug) call timer (tin,0,nb_cycles_initial,'computing sediment flux')

if (l%law.eq.3) call compute_sediment_flux (t,g,p,l)

if (p%debug) call timer (tin,1,nb_cycles_initial,'sediment flux computed')

! initialize dhdt array (to store and plot/display the erosion rate)

allocate (dhdt(p%nn))
!$OMP parallel shared(dhdt,t,p) private(ij)
!$OMP do schedule(dynamic,p%chunk)
  do ij=1,p%nn
  dhdt(ij)=t%h(ij)
  enddo
!$OMP end do nowait
!$OMP end parallel

! calculates fluvial erosion

flux_true=0.d0
!$OMP parallel shared(p,t)
!$OMP do private(ij) reduction(+:flux_true) schedule(dynamic,p%chunk)
do ij=1,p%nn
flux_true=flux_true+t%h(ij)
enddo
!$OMP end do nowait
!$OMP end parallel

!$OMP parallel shared(sed_flux,t,p) private(ij)
!$OMP do schedule(dynamic,p%chunk)
  do ij=1,p%nn
  sed_flux(ij)=t%h(ij)
  enddo
!$OMP end do nowait
!$OMP end parallel

if (p%debug) call timer (tin,0,nb_cycles_initial,'computing fluvial erosion')
call fluvial_erosion (t,g,p,l)
if (p%debug) call timer (tin,1,nb_cycles_initial,'fluvial erosion computed')

sumh=0.d0
!$OMP parallel shared(p,t)
!$OMP do private(ij) reduction(+:sumh) schedule(dynamic,p%chunk)
do ij=1,p%nn
sumh=sumh+t%h(ij)
enddo
!$OMP end do nowait
!$OMP end parallel

flux_true=flux_true-sumh

!$OMP parallel shared(sed_flux,t,p) private(ij)
!$OMP do schedule(dynamic,p%chunk)
  do ij=1,p%nn
  sed_flux(ij)=sed_flux(ij)-t%h(ij)
  enddo
!$OMP end do nowait
!$OMP end parallel

if (l%tauc.lt.tiny(l%tauc)) then
!$OMP parallel shared(t,p) private(ij)
!$OMP do schedule(dynamic,p%chunk)
  do ij=1,p%nn
  t%pmask(1,ij)=.false.
  enddo
!$OMP end do nowait
!$OMP end parallel
endif

! calculate diffusion

flux_true=0.d0
!$OMP parallel shared(p,t)
!$OMP do private(ij) reduction(+:flux_true) schedule(dynamic,p%chunk)
do ij=1,p%nn
flux_true=flux_true+t%h(ij)
enddo
!$OMP end do nowait
!$OMP end parallel

!$OMP parallel shared(sed_flux,t,p) private(ij)
!$OMP do schedule(dynamic,p%chunk)
  do ij=1,p%nn
  sed_flux(ij)=sed_flux(ij)+t%h(ij)
  enddo
!$OMP end do nowait
!$OMP end parallel

if (p%debug) call timer (tin,0,nb_cycles_initial,'diffusing')
call diffusion (t,g,p,l)
if (p%debug) call timer (tin,1,nb_cycles_initial,'diffused')

sumh=0.d0
!$OMP parallel shared(p,t)
!$OMP do private(ij) reduction(+:sumh) schedule(dynamic,p%chunk)
do ij=1,p%nn
sumh=sumh+t%h(ij)
enddo
!$OMP end do nowait
!$OMP end parallel

flux_true=flux_true-sumh

!$OMP parallel shared(sed_flux,t,p) private(ij)
!$OMP do schedule(dynamic,p%chunk)
  do ij=1,p%nn
  sed_flux(ij)=sed_flux(ij)-t%h(ij)
  enddo
!$OMP end do nowait
!$OMP end parallel

! calculates accumulated flux

!$OMP parallel shared(sed_flux,t,p) private(ij)
!$OMP do schedule(dynamic,p%chunk)
  do ij=1,p%nn
  sed_flux(ij)=sed_flux(ij)*p%dx*p%dy/p%dt
  enddo
!$OMP end do nowait
!$OMP end parallel

! note that we make sure that negative values are not passed down the network
! these negative values come from the diffusion term

if (p%debug) print*,'calculating accumulated flux'
!$OMP parallel shared(g,p,sed_flux)
!$OMP do private(i,j,ij,ijk) schedule(dynamic,1)
  do j=1,p%num_threads
    do i=g%nstackp,1,-1
    ij=g%stack((j-1)*g%nstackp+i)
      if (ij.ne.0) then
      ijk=g%receiver(ij)
      if (ijk.ne.ij) sed_flux(ijk)=max(0.d0,sed_flux(ijk)+sed_flux(ij))
      endif
    enddo
  enddo
!$OMP end do nowait
!$OMP end parallel

! there might still be some negative values near the top of the stack where
! we may have started with a negative value; this is likely to happen for
! very small amplitude topography

!$OMP parallel shared(sed_flux,p) private(ij)
!$OMP do schedule(dynamic,p%chunk)
  do ij=1,p%nn
  sed_flux(ij)=max(0.d0,sed_flux(ij))
  enddo
!$OMP end do nowait
!$OMP end parallel

if (p%debug) print*,'accumulated flux calculated'

! computes erosion rate

!$OMP parallel shared(dhdt,t,p) private(ij)
!$OMP do schedule(dynamic,p%chunk)
do ij=1,p%nn
dhdt(ij)=(dhdt(ij)-t%h(ij))/p%dt
enddo
!$OMP end do nowait
!$OMP end parallel

if (p%nage.gt.0) then
!$OMP parallel shared(dhdt,p,eRate) private(ij)
!$OMP do schedule(dynamic,p%chunk)
do ij=1,p%nage
eRate(p%istep,ij)=-dhdt(ij_age(ij))
enddo
!$OMP end do nowait
!$OMP end parallel
endif

! calculates flexure

if (p%debug) call timer (tin,0,nb_cycles_initial,'calculating flexure')
if (p%flexure) call flexure (p,t,g)
if (p%debug) call timer (tin,1,nb_cycles_initial,'flexure calculated')

! computes surface uplift rate

!$OMP parallel shared(dudt,t,p) private(ij)
!$OMP do schedule(dynamic,p%chunk)
  do ij=1,p%nn
  dudt(ij)=(t%h(ij)-dudt(ij))/p%dt
  enddo
!$OMP end do nowait
!$OMP end parallel

! updates basement height to avoid confusion due to round off error

!$OMP parallel shared(t,p) private(ij)
!$OMP do schedule(dynamic,p%chunk)
  do ij=1,p%nn
  if (t%h(ij).lt.t%hb(ij)) t%hb(ij)=t%h(ij)+1.d-10
  enddo
!$OMP end do nowait
!$OMP end parallel

! computes coastline

!$OMP parallel shared(t,p,g,coast) private(ij)
!$OMP do schedule(dynamic,p%chunk)
  do ij=1,p%nn
  if (t%h(ij).gt.p%sea_level .and. t%h(g%receiver(ij)).lt.p%sea_level) coast(ij)=.true.
  enddo
!$OMP end do nowait
!$OMP end parallel

! computes slope statistics

if (c5(4:4).eq.'1') then
allocate (plot(p%nn))
call Slope (t,p,g,l,plot,mean,var,skew,kurt)
deallocate (plot)
write (91,'(128g14.6)') p%time,mean,var,skew,kurt
endif

! computes out flux

sumh=0.d0
!$OMP parallel shared(p,t)
!$OMP do private(ij) reduction(+:sumh) schedule(dynamic,p%chunk)
do ij=1,p%nn
sumh=sumh+t%h(ij)
enddo
!$OMP end do nowait
!$OMP end parallel
flux_out=-(sumh-flux_out)/p%dt*p%dx*p%dy

if (c5(5:5).eq.'1') write (89,'(128g14.6)') p%time,minval(t%h),sum(t%h)/p%nn,maxval(t%h),flux_in, &
sum(sed_flux(1:p%nx))+sum(sed_flux(p%nx:p%nn:p%nx))+sum(sed_flux(p%nx*(p%ny-1):p%nn))+sum(sed_flux(1:p%nn:p%nx)), &
sum(sed_flux(1:p%nx)),sum(sed_flux(p%nx:p%nn:p%nx)),sum(sed_flux(p%nx*(p%ny-1):p%nn)),sum(sed_flux(1:p%nn:p%nx))

do i=1,p%nbox
if (b(i)%fnme(7:7).eq.'$') then
  do j=1,b(i)%ntime
  if ((p%time-b(i)%timei(j))*(p%time-b(i)%timef(j)).le.0.d0) b(i)%volume_predicted(j)=b(i)%volume_predicted(j)+ &
                                                                  sum(dhdt,mask=b(i)%mask)*p%dt*p%dx*p%dy
  enddo
else
  do j=1,b(i)%ntime
  if ((p%time-b(i)%timei(j))*(p%time-b(i)%timef(j)).le.0.d0) b(i)%volume_predicted(j)=b(i)%volume_predicted(j)+ &
                                                                  sum(sed_flux,mask=b(i)%mask.and.coast)*p%dt
  enddo
endif
enddo
  
! saves:plots the results according to the value of various flags (see module_definitions.f90)

  if ((p%istep/p%nfreq)*p%nfreq.eq.p%istep) then

    if (p%plot_DEM.ne.0) call DEM (p%run//'/Step'//cstep,17,p,t%h,p%plot_DEM)

    if (p%plot_topo.ne.0) then
    if (p%plot_topo.eq.1) call bmp (t%h,t%h,p%nx,p%ny, &
          p%run//'/Topo'//cstep//'.bmp',21,256,0,p%topo_min,p%topo_max,1,p%vex)
    if (p%plot_topo.eq.-1) call bmp (t%h,t%h,p%nx,p%ny, &
          p%run//'/Topo'//cstep//'.bmp',21,256,0,p%topo_min,p%topo_max,0,p%vex)
    if (p%vtk) call write_vtk ('Topo',4,p,t%h,1,field,title,nf,nfmax)
    endif

    if (p%plot_strati.ne.0) then
    allocate (strati(p%nx*p%ny))
      do ij=1,p%nn
      strati(ij)=1.d0
        do i=1,l%strati_n
          if (t%href(ij)-t%h(ij).ge.l%strati_top(i).and.t%href(ij)-t%h(ij).le.l%strati_bottom(i)) &
              strati(ij)=l%strati_f(i)
        enddo
      enddo
    if (p%plot_strati.eq.1) call bmp (t%h,strati,p%nx,p%ny, &
          p%run//'/Strati'//cstep//'.bmp',23,256,0,p%strati_min,p%strati_max,1,p%vex)
    if (p%plot_strati.eq.-1) call bmp (t%h,strati,p%nx,p%ny, &
          p%run//'/Strati'//cstep//'.bmp',23,256,0,p%strati_min,p%strati_max,0,p%vex)
    if (p%vtk) call write_vtk ('Strati',6,p,strati,1,field,title,nf,nfmax)
    deallocate (strati)
    endif

    if (p%plot_hardness.ne.0) then
    allocate (strati(p%nx*p%ny))
      do ij=1,p%nn
      x=p%dx*(1+mod(ij-1,p%nx))
      y=p%dy*(1+(ij-1)/p%nx)
      strati(ij)=1.d0
        do k=1,p%granite_n
        dist=sqrt((x-p%granite_x(k))**2/p%granite_rx(k)**2+(y-p%granite_y(k))**2/p%granite_ry(k)**2)
        if (dist.lt.1.d0.and. &
            t%href(ij)-t%h(ij).ge.p%granite_top(k).and.t%href(ij)-t%h(ij).le.p%granite_bottom(k)) strati(ij)=p%granite_dk(k)
        enddo
      enddo
    if (p%plot_hardness.eq.1) call bmp (t%h,strati,p%nx,p%ny, &
          p%run//'/Hardness'//cstep//'.bmp',25,256,0,p%hardness_min,p%hardness_max,1,p%vex)
    if (p%plot_hardness.eq.-1) call bmp (t%h,strati,p%nx,p%ny, &
          p%run//'/Hardness'//cstep//'.bmp',25,256,0,p%hardness_min,p%hardness_max,0,p%vex)
    if (p%vtk) call write_vtk ('Hardness',8,p,strati,1,field,title,nf,nfmax)
    deallocate (strati)
    endif

    if (p%plot_density.ne.0) then
    allocate (strati(p%nx*p%ny))
      do ij=1,p%nn
      x=p%dx*(1+mod(ij-1,p%nx))
      y=p%dy*(1+(ij-1)/p%nx)
      strati(ij)=0.d0
        do k=1,p%granite_n
        dist=sqrt((x-p%granite_x(k))**2/p%granite_rx(k)**2+(y-p%granite_y(k))**2/p%granite_ry(k)**2)
!        if (dist.lt.1.d0) strati(ij)=t%h(ij)-t%href(ij)+p%granite_top(k)
        if (dist.lt.1.d0.and. &
            t%href(ij)-t%h(ij).ge.p%granite_top(k).and.t%href(ij)-t%h(ij).le.p%granite_bottom(k)) strati(ij)=p%granite_drho(k)
        enddo
      enddo
    if (p%plot_density.eq.1) call bmp (t%h,strati,p%nx,p%ny, &
          p%run//'/Density'//cstep//'.bmp',24,256,0,p%density_min,p%density_max,1,p%vex)
    if (p%plot_density.eq.-1) call bmp (t%h,strati,p%nx,p%ny, &
          p%run//'/Density'//cstep//'.bmp',24,256,0,p%density_min,p%density_max,0,p%vex)
    if (p%vtk) call write_vtk ('Density',7,p,strati,1,field,title,nf,nfmax)
    deallocate (strati)
    endif

    if (p%plot_erosion.ne.0) then
    allocate (plot(p%nx*p%ny))
    plot=max(0.d0,t%hi-t%h)
    if (p%plot_erosion.eq.1) call bmp (t%h,plot,p%nx,p%ny, &
          p%run//'/Erosion'//cstep//'.bmp',24,256,0,p%erosion_min,p%erosion_max,1,p%vex)
    if (p%plot_erosion.eq.-1) call bmp (t%h,plot,p%nx,p%ny, &
          p%run//'/Erosion'//cstep//'.bmp',24,256,0,p%erosion_min,p%erosion_max,0,p%vex)
    if (p%vtk) call write_vtk ('Erosion',7,p,plot,1,field,title,nf,nfmax)
    deallocate (plot)
    endif

    if (p%plot_rate.ne.0) then
    if (p%plot_rate.eq.1) call bmp (t%h,dhdt,p%nx,p%ny, &
          p%run//'/Rate'//cstep//'.bmp',21,256,0,p%rate_min,p%rate_max,1,p%vex)
    if (p%plot_rate.eq.-1) call bmp (t%h,dhdt,p%nx,p%ny, &
          p%run//'/Rate'//cstep//'.bmp',21,256,0,p%rate_min,p%rate_max,0,p%vex)
    if (p%vtk) call write_vtk ('Rate',4,p,dhdt,1,field,title,nf,nfmax)
    endif

     if (p%plot_uplift.ne.0) then
     if (p%plot_uplift.eq.1) call bmp (t%h,dudt,p%nx,p%ny, &
           p%run//'/Uplift'//cstep//'.bmp',23,256,0,p%uplift_min,p%uplift_max,1,p%vex)
     if (p%plot_uplift.eq.-1) call bmp (t%h,dudt,p%nx,p%ny, &
           p%run//'/Uplift'//cstep//'.bmp',23,256,0,p%uplift_min,p%uplift_max,0,p%vex)
     if (p%vtk) call write_vtk ('Uplift',6,p,dudt,1,field,title,nf,nfmax)
     endif

     if (p%plot_rock_uplift.ne.0) then
     if (p%plot_rock_uplift.eq.1) call bmp (t%h,t%u,p%nx,p%ny, &
           p%run//'/RockUplift'//cstep//'.bmp',27,256,0,p%rock_uplift_min,p%rock_uplift_max,1,p%vex)
     if (p%plot_rock_uplift.eq.-1) call bmp (t%h,t%u,p%nx,p%ny, &
           p%run//'/RockUplift'//cstep//'.bmp',27,256,0,p%rock_uplift_min,p%rock_uplift_max,0,p%vex)
     if (p%vtk) call write_vtk ('RockUplift',10,p,t%u,1,field,title,nf,nfmax)
     endif


    if (p%plot_sedim_flux.ne.0) then
!    call bmp (t%h,t%sediment,p%nx,p%ny, &
!          p%run//'/SedFlux'//cstep//'.bmp',23,256,0,0.d0,100.d0,1,p%vex)
!    if (p%vtk) call write_vtk ('Sedim_Flux',10,p,t%sediment,1,field,title,nf,nfmax)
    call bmp (t%h,log10(sed_flux),p%nx,p%ny, &
          p%run//'/SedFlux'//cstep//'.bmp',24,256,0,0.d0,100.d0,1,p%vex)
    if (p%vtk) call write_vtk ('Sedim_Flux',10,p,log10(sed_flux),1,field,title,nf,nfmax)
    endif

    if (p%plot_sedim.ne.0) then
    allocate (plot(p%nx*p%ny))
    plot=t%h-t%hb
    if (p%plot_sedim.eq.1) call bmp (t%h,plot,p%nx,p%ny, &
          p%run//'/Sediment'//cstep//'.bmp',25,256,0,p%sedim_min,p%sedim_max,1,p%vex)
    if (p%plot_sedim.eq.-1) call bmp (t%h,plot,p%nx,p%ny, &
          p%run//'/Sediment'//cstep//'.bmp',25,256,0,p%sedim_min,p%sedim_max,0,p%vex)
    if (p%vtk) call write_vtk ('Sedim',5,p,plot,1,field,title,nf,nfmax)
    deallocate (plot)
    endif

    if (p%plot_discharge.ne.0) then
    allocate (plot(p%nx*p%ny))
    plot=log10(t%discharge)
    if (p%plot_discharge.eq.1) call bmp (t%h,plot,p%nx,p%ny, &
          p%run//'/Discharge'//cstep//'.bmp',26,256,0,p%discharge_min,p%discharge_max,1,p%vex)
    if (p%plot_discharge.eq.-1) call bmp (t%h,plot,p%nx,p%ny, &
          p%run//'/Discharge'//cstep//'.bmp',26,256,0,p%discharge_min,p%discharge_max,0,p%vex)
    if (p%vtk.and.p%plot_discharge.ne.0) call write_vtk ('Discharge',9,p,plot,1,field,title,nf,nfmax)
    deallocate (plot)
    endif

    if (p%plot_catchment.eq.1) call bmp (t%h,dfloat(g%catchment),p%nx,p%ny, &
          p%run//'/Catchment'//cstep//'.bmp',26,256,1,0.d0,400.d0,0,p%vex)
    if (p%vtk.and.p%plot_catchment.ne.0) call write_vtk ('Catchment',9,p,dfloat(g%catchment),1,field,title,nf,nfmax)

    if (p%plot_steepnessindex.ne.0) then
    allocate (plot(p%nn))
    call SteepnessIndex (t,p,g,l,plot)
    if (p%plot_steepnessindex.eq.1) call bmp (t%h,plot,p%nx,p%ny, &
          p%run//'/SteepnessIndex'//cstep//'.bmp',31,256,0,p%steepnessindex_min,p%steepnessindex_max,1,p%vex)
    if (p%plot_steepnessindex.eq.-1) call bmp (t%h,plot,p%nx,p%ny, &
          p%run//'/SteepnessIndex'//cstep//'.bmp',31,256,0,p%steepnessindex_min,p%steepnessindex_max,0,p%vex)
    if (p%vtk.and.p%plot_steepnessindex.ne.0) call write_vtk ('SteepnessIndex',14,p,plot,1,field,title,nf,nfmax)
    deallocate (plot)
    endif

    if (p%plot_slope.ne.0) then
    allocate (plot(p%nn))
    call Slope (t,p,g,l,plot,mean,var,skew,kurt)
    if (p%plot_slope.eq.1) call bmp (t%h,plot,p%nx,p%ny, &
          p%run//'/Slope'//cstep//'.bmp',22,256,0,p%slope_min,p%slope_max,1,p%vex)
    if (p%plot_slope.eq.-1) call bmp (t%h,plot,p%nx,p%ny, &
          p%run//'/Slope'//cstep//'.bmp',22,256,0,p%slope_min,p%slope_max,0,p%vex)
    if (p%vtk.and.p%plot_slope.ne.0) call write_vtk ('Slope',5,p,plot,1,field,title,nf,nfmax)
    deallocate (plot)
    endif

    if (p%plot_curvature.ne.0) then
    allocate (plot(p%nn))
    call Curvature (t,p,g,l,plot)
    if (p%plot_curvature.eq.1) call bmp (t%h,plot,p%nx,p%ny, &
          p%run//'/Curvature'//cstep//'.bmp',26,256,0,p%curvature_min,p%curvature_max,1,p%vex)
    if (p%plot_curvature.eq.-1) call bmp (t%h,plot,p%nx,p%ny, &
          p%run//'/Curvature'//cstep//'.bmp',26,256,0,p%curvature_min,p%curvature_max,0,p%vex)
    if (p%vtk.and.p%plot_curvature.ne.0) call write_vtk ('Curvature',9,p,plot,1,field,title,nf,nfmax)
    deallocate (plot)
    endif

    if (p%plot_concavity.ne.0) then
    allocate (plot(p%nn))
    call Concavity (t,p,g,l,plot)
    if (p%plot_concavity.eq.1) call bmp (t%h,plot,p%nx,p%ny, &
          p%run//'/Concavity'//cstep//'.bmp',26,256,0,p%concavity_min,p%concavity_max,1,p%vex)
    if (p%plot_concavity.eq.-1) call bmp (t%h,plot,p%nx,p%ny, &
          p%run//'/Concavity'//cstep//'.bmp',26,256,0,p%concavity_min,p%concavity_max,0,p%vex)
    if (p%vtk.and.p%plot_concavity.ne.0) call write_vtk ('Concavity',9,p,plot,1,field,title,nf,nfmax)
    deallocate (plot)
    endif

    if (p%plot_chi.ne.0) then
    allocate (t%chi(p%nn))
    call Chi (g,t,l,p)
    call bmp (t%h,t%chi,p%nx,p%ny, &
          p%run//'/Chi'//cstep//'.bmp',20,256,0,0.d0,1.d0,1,p%vex)
    if (p%vtk.and.p%plot_chi.ne.0) call write_vtk ('Chi',3,p,t%chi,1,field,title,nf,nfmax)
    deallocate (t%chi)
    endif

    if (p%plot_fft.ne.0) then
    allocate (t%fft(p%nn))
    call FFT (t,p)
    call bmp (t%h,t%fft,p%nx,p%ny, &
          p%run//'/FFT'//cstep//'.bmp',20,256,0,0.d0,1.d0,0,p%vex)
    if (p%vtk.and.p%plot_fft.ne.0) call write_vtk ('FFT',3,p,t%fft,1,field,title,nf,nfmax)
    deallocate (t%fft)
    endif

    if (p%plot_volume_mask.eq.1 .and. p%nbox.gt.0) then
    allocate (plot(p%nn))
    plot=0.d0
      do i=1,p%nbox
      where (b(i)%mask) plot=i
      enddo
    where (coast) plot=p%nbox+1
    call bmp (t%h,plot,p%nx,p%ny, &
          p%run//'/MaskVolume'//cstep//'.bmp',27,256,0,0.d0,dfloat(p%nbox+1),1,p%vex)
    if (p%vtk.and.p%plot_volume_mask.ne.0) call write_vtk ('MaskVolume',10,p,plot,1,field,title,nf,nfmax)
    deallocate (plot)
    endif

    if (p%vtk) call write_vtk ('Nil',3,p,t%h,-1,field,title,nf,nfmax)

    open (7,file=p%run//'/RESTART',status='unknown',access='direct',recl=48*p%nx*p%ny+8+8)
    write (7,rec=1) p%time,p%istep,p%nPecube,t%h,t%hi,t%href,t%hiso,t%hb,t%hp
    close (7)

  endif

  if (p%ndemoulin.ne.0) then
  if ((p%istep/p%ndemoulin)*p%ndemoulin.eq.p%istep) then
  call demoulin (t,g,p,l,sr)
  open (61,file=p%run//'/Demoulin.txt',status='unknown',access='append')
  write (61,*) p%time,sr
  close (61)
  endif
  endif

  if ((p%istep/p%nmetric)*p%nmetric.eq.p%istep) then
    if (p%plot_topo.eq.1) then
    call metric (t,p,g,dhdt,p%run//'/Metric'//cstep,19,1)
    else
    call metric (t,p,g,dhdt,p%run//'/Metric'//cstep,19,-1)
    endif
  endif

call timer (tdt,1,nb_cycles_initial,'$')

nfreq_screen=nfreq_screen+1
if (nfreq_screen.eq.20) then
if (.not.p%quiet) write (*,'(a)') '  Step     hmin     hmean      hmax     CPU time'
nfreq_screen=0
endif

sumh=0.d0
!$OMP parallel shared(p,t)
!$OMP do private(ij) reduction(+:sumh) schedule(dynamic,p%chunk)
do ij=1,p%nn
sumh=sumh+t%h(ij)
enddo
!$OMP end do nowait
!$OMP end parallel

!$OMP workshare
minh=minval(t%h)
maxh=maxval(t%h)
!$OMP end workshare

if (.not.p%quiet) write (*,'(i7,3f10.3,a3,f9.4,a1)') p%istep,minh,sumh/p%nn,maxh,'  (',tdt,')'

if (p%Pecube_nfreq.ne.0) then
if (p%debug) call timer (tin,0,nb_cycles_initial,'writing Pecube output')
if ((p%istep/p%Pecube_nfreq)*p%Pecube_nfreq.eq.p%istep) call Interface_Pecube (t,p,p%nPecube)
if (p%debug) call timer (tin,1,nb_cycles_initial,'Pecube output writen')
endif

deallocate (dhdt,dudt)

!end of time stepping

enddo

  do i=1,p%nbox
  open (1000,file=p%run//'/V'//trim(b(i)%fnme),status='unknown')
    do j=1,b(i)%ntime
    write (1000,*) b(i)%timei(j),b(i)%timef(j),b(i)%volume(j),b(i)%dvolume(j),b(i)%volume_predicted(j)
    enddo
  close (1000)
  enddo

if (c5(4:4).eq.'1') then
close (91)
endif

if (c5(5:5).eq.'1') then
close (89)
open (3,file=p%run//'/Metric.R',status='unknown',access='append')
write (3,'(a)') 'pdf("'//p%run//'/Fluxes.pdf")'
write (3,'(a)') 'x=read.table("'//p%run//'/Fluxes.txt")'
write (3,'(a)') 'plot(x[,1],x[,6],type="l",main="Fluxes",'// &
                'xlab="Time (yr)",'// &
                'ylab="Flux (m^3/yr)")'
write (3,'(a)') 'lines(x[,1],x[,5])'
write (3,'(a)') 'dev.off()'
write (3,'(a)')
close (3)
open (3,file=p%run//'/Metric.R',status='unknown',access='append')
write (3,'(a)') 'pdf("'//p%run//'/TopoEvolution.pdf")'
write (3,'(a)') 'x=read.table("'//p%run//'/Fluxes.txt")'
write (3,'(a)') 'plot(x[,1],x[,4],type="l",main="Maximum Topography",'// &
                'xlab="Time (yr)",'// &
                'ylab="Height (m)")'
write (3,'(a)') 'dev.off()'
write (3,'(a)')
close (3)
endif

if (p%nage.gt.0) then
  if ((p%dx_age.gt.0.d0 .or. p%dy_age.gt.0.d0 .or. maxval(p%h_age(1:p%nage)).gt.-99998.d0) &
      .and. .not.second) then
  ii=int(p%dx_age*p%nx/p%xl)
  jj=int(p%dy_age*p%ny/p%yl)
    do i=1,p%nage
    ijk=ij_age(i)
      do jjj=-jj,jj
        do iii=-ii,ii
        ij=ijk+iii+jjj*p%nx
        ij=1+modulo(ij-1,p%nn)
        if (abs(t%h(ij)-p%h_age(i)).lt.abs(t%h(ijk)-p%h_age(i))) ijk=ij
        enddo
      enddo
    ij_age(i)=ijk
    enddo
  second=.true.
  goto 99999
  else
  allocate (temperatureH(p%nstep),aft_age(p%nage),zft_age(p%nage),ahe_age(p%nage),zhe_age(p%nage))
  allocate (SurfaceGradientToday(p%nage))
  open (71,file=p%run//'/Ages.txt',status='unknown')
  write (71,*) p%nage,p%limit_age
    do i=1,p%nage
    call compute_temperature (timeH,eRate(1,i),temperatureH,p%nstep,p%l_age, &
                              p%tmin_age,p%tmax_age,p%heat_age,p%zheat_age,p%cond_age,p%zcond_age, &
                              SurfaceGradientToday(i))
    call compute_age (timeH,temperatureH,p%nstep,'AFT',p%size_age(i),p%limit_age,aft_age(i))
    call compute_age (timeH,temperatureH,p%nstep,'ZFT',p%size_age(i),p%limit_age,zft_age(i))
    call compute_age (timeH,temperatureH,p%nstep,'AHE',p%size_age(i),p%limit_age,ahe_age(i))
    call compute_age (timeH,temperatureH,p%nstep,'ZHE',p%size_age(i),p%limit_age,zhe_age(i))
    write (71,*) p%x_age(i),p%y_age(i),p%h_age(i),t%h(ij_age(i)), &
                 aft_age(i),p%aft_age(i),p%aft_dage(i),zft_age(i),p%zft_age(i),p%zft_dage(i), &
                 ahe_age(i),p%ahe_age(i),p%ahe_dage(i),zhe_age(i),p%zhe_age(i),p%zhe_dage(i)
    enddo
  close (71)
  deallocate (temperatureH,aft_age,zft_age,ahe_age,zhe_age,SurfaceGradientToday)
  deallocate (timeH,eRate,ij_age)
  endif
endif

call system_clock (count=nb_cycles_final)
nb_cycles=nb_cycles_final-nb_cycles_initial
IF (nb_cycles_final.lt.nb_cycles_initial) &
        nb_cycles=nb_cycles+nb_cycles_max
temps_elapsed=real(nb_cycles)/nb_cycles_sec
if (.not.p%quiet) print*,'time elapsed',temps_elapsed

deallocate (t%length,g%receiver,g%bc,t%h,g%ndon,g%donors,g%catchment,t%discharge,g%stack)
deallocate (t%width,t%sediment,t%hb,t%hp,t%pmask,ierror)
if (p%vtk) deallocate (title,field)

end program FastScape

! ---

recursive subroutine add_to_stack (i,nstack,nstackp,stack,ndon,donors,catchment,n,icatch,ierror)

! Recursive routine to go through the node upstream and keeping track of confluences

implicit none

integer n,nstack,nstackp,stack(n),ndon(n+1),donors(n),catchment(n)
integer i,j,icatch,ierror

if (ierror.ne.0) return

nstack=nstack+1
  if (nstack.gt.nstackp) then
  ierror=1
  return
  endif
stack(nstack)=i
catchment(i)=icatch
  do j=ndon(i),ndon(i+1)-1
  call add_to_stack (donors(j),nstack,nstackp,stack,ndon,donors,catchment,n,icatch,ierror)
  enddo

return

end subroutine add_to_stack
