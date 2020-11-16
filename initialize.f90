subroutine initialize (p,l)

! routine to initialize all variables and/or read their value from the input file
! this routine uses the scanfile routine to do so

! note that if the input file does not exist it is created and filled with default
! values

use omp_lib
use definitions

implicit none

type (param) p
type (law) l

integer ires,iflexure,imeanflex,local_minima,refflex
integer plot_all,debug,vtk,plot_dem,quiet
character c1*1,c2*2,c3*3,c4*4
integer i,vocal

vocal=1

debug=0
call scanfile (p%infile,'debug',debug,ires,vocal)
p%debug=.TRUE.
if (debug.eq.0) p%debug=.FALSE.

quiet=0
call scanfile (p%infile,'quiet',quiet,ires,vocal)
p%quiet=.TRUE.
if (quiet.eq.0) p%quiet=.FALSE.

p%reset_random=0
call scanfile (p%infile,'reset_random',p%reset_random,ires,vocal)

p%nx=1000
call scanfile (p%infile,'nx',p%nx,ires,vocal)
p%nx0=p%nx
  if ((p%nx/4)*4.ne.p%nx) then
  p%nx=(p%nx/4)*4
  print*,'BEWARE NX HAS BEEN REDUCED TO ',p%nx
  endif
p%ny=1000
call scanfile (p%infile,'ny',p%ny,ires,vocal)
p%ny0=p%ny
  if ((p%ny/4)*4.ne.p%ny) then
  p%ny=(p%ny/4)*4
  print*,'BEWARE NY HAS BEEN REDUCED TO ',p%ny
  endif
p%nn=p%nx*p%ny
p%xl=100.d3
call scanfile (p%infile,'xl',p%xl,ires,vocal)
p%yl=100.d3
call scanfile (p%infile,'yl',p%yl,ires,vocal)
p%dx=p%xl/p%nx
p%dy=p%yl/p%ny
p%dt=100000.d0
call scanfile (p%infile,'dt',p%dt,ires,vocal)
p%nstep=100
call scanfile (p%infile,'nstep',p%nstep,ires,vocal)
p%nfreq=p%nstep
call scanfile (p%infile,'nfreq',p%nfreq,ires,vocal)
p%nmetric=p%nfreq
call scanfile (p%infile,'nmetric',p%nmetric,ires,vocal)

p%ndemoulin=0
call scanfile (p%infile,'ndemoulin',p%ndemoulin,ires,vocal)
p%stream_threshold_area=0.5d6
call scanfile (p%infile,'stream_threshold_area',p%stream_threshold_area,ires,vocal)
p%minimum_area=10.d6
call scanfile (p%infile,'minimum_area',p%minimum_area,ires,vocal)

plot_all=1
call scanfile (p%infile,'plot_all',plot_all,ires,vocal)
p%plot_all=.TRUE.
if (plot_all.eq.0) p%plot_all=.FALSE.
vtk=0
call scanfile (p%infile,'vtk',vtk,ires,vocal)
p%vtk=.FALSE.
if (vtk.eq.1) p%vtk=.TRUE.
p%plot_DEM=0
call scanfile (p%infile,'plot_DEM',p%plot_DEM,ires,vocal)
p%vex=2.
call scanfile (p%infile,'vex',p%vex,ires,vocal)

  if (p%plot_all) then
  p%plot_topo=1
  p%plot_strati=1
  p%plot_hardness=1
  p%plot_density=1
  p%plot_erosion=1
  p%plot_rate = 1
  p%plot_uplift = 1
  p%plot_rock_uplift = 1
  p%plot_sedim=1
  p%plot_sedim_flux = 1
  p%plot_discharge=1
  p%plot_catchment=1
  p%plot_precip=1
  p%plot_slope=1
  p%plot_curvature=1
  p%plot_steepnessindex=1
  p%plot_concavity=1
  p%plot_chi=1
  p%plot_fft=1
  p%plot_volume_mask=1
  p%imetric=11111
  else
  p%plot_topo=0
  p%plot_strati=0
  p%plot_hardness=0
  p%plot_density=0
  p%plot_erosion=0
  p%plot_rate = 0
  p%plot_uplift = 0
  p%plot_rock_uplift = 0
  p%plot_sedim=0
  p%plot_sedim_flux = 0
  p%plot_discharge=0
  p%plot_catchment=0
  p%plot_precip=0
  p%plot_slope=0
  p%plot_curvature=0
  p%plot_steepnessindex=0
  p%plot_concavity=0
  p%plot_chi=0
  p%plot_fft=0
  p%plot_volume_mask=0
  p%imetric=00000
  endif
call scanfile (p%infile,'plot_topo',p%plot_topo,ires,vocal)
call scanfile (p%infile,'plot_strati',p%plot_strati,ires,vocal)
call scanfile (p%infile,'plot_hardness',p%plot_hardness,ires,vocal)
call scanfile (p%infile,'plot_density',p%plot_density,ires,vocal)
call scanfile (p%infile,'plot_erosion',p%plot_erosion,ires,vocal)
call scanfile (p%infile,'plot_rate',p%plot_rate,ires,vocal)
call scanfile (p%infile,'plot_uplift',p%plot_uplift,ires,vocal)
call scanfile (p%infile,'plot_rock_uplift',p%plot_rock_uplift,ires,vocal)
call scanfile (p%infile,'plot_sedim',p%plot_sedim,ires,vocal)
call scanfile (p%infile,'plot_sedim_flux',p%plot_sedim_flux,ires,vocal)
call scanfile (p%infile,'plot_discharge',p%plot_discharge,ires,vocal)
call scanfile (p%infile,'plot_catchment',p%plot_catchment,ires,vocal)
call scanfile (p%infile,'plot_precip',p%plot_precip,ires,vocal)
call scanfile (p%infile,'plot_slope',p%plot_slope,ires,vocal)
call scanfile (p%infile,'plot_curvature',p%plot_curvature,ires,vocal)
call scanfile (p%infile,'plot_steepnessindex',p%plot_steepnessindex,ires,vocal)
call scanfile (p%infile,'plot_concavity',p%plot_concavity,ires,vocal)
call scanfile (p%infile,'plot_chi',p%plot_chi,ires,vocal)
call scanfile (p%infile,'plot_fft',p%plot_fft,ires,vocal)
call scanfile (p%infile,'plot_volume_mask',p%plot_volume_mask,ires,vocal)
call scanfile (p%infile,'metric',p%imetric,ires,vocal)
p%topo_min=0.d0
call scanfile (p%infile,'topo_min',p%topo_min,ires,vocal)
p%topo_max=1000.d0
call scanfile (p%infile,'topo_max',p%topo_max,ires,vocal)
p%strati_min=0.d0
call scanfile (p%infile,'strati_min',p%strati_min,ires,vocal)
p%strati_max=1.d0
call scanfile (p%infile,'strati_max',p%strati_max,ires,vocal)
p%hardness_min=0.d0
call scanfile (p%infile,'hardness_min',p%hardness_min,ires,vocal)
p%hardness_max=1.d0
call scanfile (p%infile,'hardness_max',p%hardness_max,ires,vocal)
p%density_min=0.d0
call scanfile (p%infile,'density_min',p%density_min,ires,vocal)
p%density_max=400.d0
call scanfile (p%infile,'density_max',p%density_max,ires,vocal)
p%erosion_min=0.d0
call scanfile (p%infile,'erosion_min',p%erosion_min,ires,vocal)
p%erosion_max=p%topo_max
call scanfile (p%infile,'erosion_max',p%erosion_max,ires,vocal)
p%rate_min=-1.d-4
call scanfile (p%infile,'rate_min',p%rate_min,ires,vocal)
p%rate_max=1.d-4
call scanfile (p%infile,'rate_max',p%rate_max,ires,vocal)
p%uplift_min=-1.d-4
call scanfile (p%infile,'uplift_min',p%uplift_min,ires,vocal)
p%uplift_max=1.d-4
call scanfile (p%infile,'uplift_max',p%uplift_max,ires,vocal)
p%rock_uplift_min=-1.d-4
call scanfile (p%infile,'rock_uplift_min',p%rock_uplift_min,ires,vocal)
p%rock_uplift_max=1.d-4
call scanfile (p%infile,'rock_uplift_max',p%rock_uplift_max,ires,vocal)
p%sedim_min=0.d0
call scanfile (p%infile,'sedim_min',p%sedim_min,ires,vocal)
p%sedim_max=p%topo_max
call scanfile (p%infile,'sedim_max',p%sedim_max,ires,vocal)
p%precip_min=0.d0
call scanfile (p%infile,'precip_min',p%precip_min,ires,vocal)
p%precip_max=1.d0
call scanfile (p%infile,'precip_max',p%precip_max,ires,vocal)
p%discharge_min=log10(p%dx*p%dy*p%precip_max)
call scanfile (p%infile,'discharge_min',p%discharge_min,ires,vocal)
p%discharge_max=log10(p%xl*p%yl*p%precip_max)
call scanfile (p%infile,'discharge_max',p%discharge_max,ires,vocal)
p%slope_min=0.d0
call scanfile (p%infile,'slope_min',p%slope_min,ires,vocal)
p%slope_max=60.d0
call scanfile (p%infile,'slope_max',p%slope_max,ires,vocal)
p%curvature_min=-5.d-3
call scanfile (p%infile,'curvature_min',p%curvature_min,ires,vocal)
p%curvature_max=5.d-3
call scanfile (p%infile,'curvature_max',p%curvature_max,ires,vocal)
p%steepnessindex_min=0.d0
call scanfile (p%infile,'steepnessindex_min',p%steepnessindex_min,ires,vocal)
p%steepnessindex_max=1.d4
call scanfile (p%infile,'steepnessindex_max',p%steepnessindex_max,ires,vocal)
p%concavity_min=0.d0
call scanfile (p%infile,'concavity_min',p%concavity_min,ires,vocal)
p%concavity_max=1.d0
call scanfile (p%infile,'concavity_max',p%concavity_max,ires,vocal)
p%discharge_min_chi=0.5d6
call scanfile (p%infile,'discharge_min_chi',p%discharge_min_chi,ires,vocal)
p%HeightDis_min=-999.d0
call scanfile (p%infile,'heightdis_min',p%HeightDis_min,ires,vocal)
p%HeightDis_max=-999.d0
call scanfile (p%infile,'heightdis_max',p%HeightDis_max,ires,vocal)
p%SlopeDis_min=-999.d0
call scanfile (p%infile,'slopedis_min',p%SlopeDis_min,ires,vocal)
p%SlopeDis_max=-999.d0
call scanfile (p%infile,'slopedis_max',p%SlopeDis_max,ires,vocal)
p%CurveDis_min=-999.d0
call scanfile (p%infile,'curvedis_min',p%CurveDis_min,ires,vocal)
p%CurveDis_max=-999.d0
call scanfile (p%infile,'curvedis_max',p%CurveDis_max,ires,vocal)
p%DischargeDis_min=-999.d0
call scanfile (p%infile,'dischargedis_min',p%DischargeDis_min,ires,vocal)
p%DischargeDis_max=-999.d0
call scanfile (p%infile,'dischargedis_max',p%DischargeDis_max,ires,vocal)
iflexure=0
call scanfile (p%infile,'flexure',iflexure,ires,vocal)
p%flexure=.FALSE.
if (iflexure.eq.1) p%flexure=.TRUE.
refflex=1
call scanfile (p%infile,'refflex',refflex,ires,vocal)
p%refflex=.FALSE.
if (refflex.eq.1) p%refflex=.TRUE.
p%thickflex=20.d3
call scanfile (p%infile,'thickflex',p%thickflex,ires,vocal)
p%ym=1.d11
call scanfile (p%infile,'ym',p%ym,ires,vocal)
p%pratio=0.25d0
call scanfile (p%infile,'pratio',p%pratio,ires,vocal)
p%rhocflex=2750.d0
call scanfile (p%infile,'rhocflex',p%rhocflex,ires,vocal)
p%rhoaflex=3300.d0
call scanfile (p%infile,'rhoaflex',p%rhoaflex,ires,vocal)
imeanflex=0
p%meanflex=.FALSE.
call scanfile (p%infile,'meanflex',imeanflex,ires,vocal)
if (imeanflex.eq.1) p%meanflex=.TRUE.

p%num_threads=1
call scanfile (p%infile,'num_threads',p%num_threads,ires,vocal)

p%sea_level=0.d0
call scanfile (p%infile,'sea_level',p%sea_level,ires,vocal)

local_minima=1
call scanfile (p%infile,'local_minima',local_minima,ires,vocal)
p%local_minima=.FALSE.
if (local_minima.eq.1) p%local_minima=.TRUE.

l%law=1
call scanfile (p%infile,'law',l%law,ires,vocal)
  if (iabs(l%law).eq.1) then
  l%m=0.4d0
  l%n=1.d0
  l%kf=1.d-5
  l%tauc=0.d0
  l%area_min=0.d0
  elseif (iabs(l%law).eq.2) then
  l%m=1.d0/3.d0
  l%n=2.d0/3.d0
  l%kf=1.d-3
  l%nk=0
  l%tol=1.d-6
  l%tauc=0.d0
  l%area_min=0.d0
  elseif (iabs(l%law).eq.3) then
  l%m=1.d0
  l%n=1.d0
  l%tauc=0.d0
  l%area_min=0.d0
  l%kf=1.d-5
  l%kfs=l%kf
  l%lf=1.d0
  else
  print*,'fluvial erosion law not implemented'
  stop
  endif
call scanfile (p%infile,'m',l%m,ires,vocal)
call scanfile (p%infile,'n',l%n,ires,vocal)
call scanfile (p%infile,'kf',l%kf,ires,vocal)
l%kfs=l%kf
call scanfile (p%infile,'kfs',l%kfs,ires,vocal)
call scanfile (p%infile,'alpha',l%lf,ires,vocal)
call scanfile (p%infile,'nk',l%nk,ires,vocal)
call scanfile (p%infile,'tol',l%tol,ires,vocal)
call scanfile (p%infile,'tauc',l%tauc,ires,vocal)
call scanfile (p%infile,'area_min',l%area_min,ires,vocal)
l%sc=-1.d0
call scanfile (p%infile,'sc',l%sc,ires,vocal)
l%dsc=-1.d0
call scanfile (p%infile,'dsc',l%dsc,ires,vocal)
l%kd=0.d0
call scanfile (p%infile,'kd',l%kd,ires,vocal)

l%p0=0.d0
call scanfile (p%infile,'p0',l%p0,ires,vocal)
l%s0=0.5d0
call scanfile (p%infile,'s0',l%s0,ires,vocal)

l%slope=0.d0
call scanfile (p%infile,'slope',l%slope,ires,vocal)

l%strati_n=0
do i=1,9
l%strati_f(i)=1.d0
l%strati_fd(i)=1.d0
l%strati_top(i)=0.d0
l%strati_bottom(i)=0.d0
write (c1,'(i1)') i
call scanfile (p%infile,'strati_f'//c1,l%strati_f(i),ires,vocal)
if (ires.eq.1) l%strati_n=i
if (ires.lt.0) l%strati_f(i)=l%strati_f(-ires)
call scanfile (p%infile,'strati_fd'//c1,l%strati_fd(i),ires,vocal)
if (ires.eq.1) l%strati_n=i
if (ires.lt.0) l%strati_fd(i)=l%strati_fd(-ires)
call scanfile (p%infile,'strati_top'//c1,l%strati_top(i),ires,vocal)
if (ires.eq.1) l%strati_n=i
if (ires.lt.0) l%strati_top(i)=l%strati_top(-ires)
call scanfile (p%infile,'strati_bottom'//c1,l%strati_bottom(i),ires,vocal)
if (ires.eq.1) l%strati_n=i
if (ires.lt.0) l%strati_bottom(i)=l%strati_bottom(-ires)
enddo

p%granite_n=0
do i=1,100
p%granite_x(i)=0.d0
p%granite_y(i)=0.d0
p%granite_top(i)=0.d0
p%granite_bottom(i)=1.d6
p%granite_rx(i)=0.d0
p%granite_ry(i)=0.d0
p%granite_drho(i)=0.d0
p%granite_dk(i)=1.d0
p%granite_dkd(i)=1.d0
if (i.lt.10) then
write (c1,'(i1)') i
call scanfile (p%infile,'granite_x'//c1,p%granite_x(i),ires,vocal)
if (ires.eq.1) p%granite_n=i
if (ires.lt.0) p%granite_x(i)=p%granite_x(-ires)
call scanfile (p%infile,'granite_y'//c1,p%granite_y(i),ires,vocal)
if (ires.eq.1) p%granite_n=i
if (ires.lt.0) p%granite_y(i)=p%granite_y(-ires)
call scanfile (p%infile,'granite_top'//c1,p%granite_top(i),ires,vocal)
if (ires.eq.1) p%granite_n=i
if (ires.lt.0) p%granite_top(i)=p%granite_top(-ires)
call scanfile (p%infile,'granite_bottom'//c1,p%granite_bottom(i),ires,vocal)
if (ires.eq.1) p%granite_n=i
if (ires.lt.0) p%granite_bottom(i)=p%granite_bottom(-ires)
call scanfile (p%infile,'granite_rx'//c1,p%granite_rx(i),ires,vocal)
if (ires.eq.1) p%granite_n=i
if (ires.lt.0) p%granite_rx(i)=p%granite_rx(-ires)
p%granite_ry(i)=p%granite_rx(i)
call scanfile (p%infile,'granite_ry'//c1,p%granite_ry(i),ires,vocal)
if (ires.eq.1) p%granite_n=i
if (ires.lt.0) p%granite_ry(i)=p%granite_ry(-ires)
call scanfile (p%infile,'granite_drho'//c1,p%granite_drho(i),ires,vocal)
if (ires.eq.1) p%granite_n=i
if (ires.lt.0) p%granite_drho(i)=p%granite_drho(-ires)
call scanfile (p%infile,'granite_dk'//c1,p%granite_dk(i),ires,vocal)
if (ires.eq.1) p%granite_n=i
if (ires.lt.0) p%granite_dk(i)=p%granite_dk(-ires)
call scanfile (p%infile,'granite_dkd'//c1,p%granite_dkd(i),ires,vocal)
if (ires.eq.1) p%granite_n=i
if (ires.lt.0) p%granite_dkd(i)=p%granite_dkd(-ires)
elseif (i.lt.100) then
write (c2,'(i2)') i
call scanfile (p%infile,'granite_x'//c2,p%granite_x(i),ires,vocal)
if (ires.eq.1) p%granite_n=i
if (ires.lt.0) p%granite_x(i)=p%granite_x(-ires)
call scanfile (p%infile,'granite_y'//c2,p%granite_y(i),ires,vocal)
if (ires.eq.1) p%granite_n=i
if (ires.lt.0) p%granite_y(i)=p%granite_y(-ires)
call scanfile (p%infile,'granite_top'//c2,p%granite_top(i),ires,vocal)
if (ires.eq.1) p%granite_n=i
if (ires.lt.0) p%granite_top(i)=p%granite_top(-ires)
call scanfile (p%infile,'granite_bottom'//c2,p%granite_bottom(i),ires,vocal)
if (ires.eq.1) p%granite_n=i
if (ires.lt.0) p%granite_bottom(i)=p%granite_bottom(-ires)
call scanfile (p%infile,'granite_rx'//c2,p%granite_rx(i),ires,vocal)
if (ires.eq.1) p%granite_n=i
if (ires.lt.0) p%granite_rx(i)=p%granite_rx(-ires)
p%granite_ry(i)=p%granite_rx(i)
call scanfile (p%infile,'granite_ry'//c2,p%granite_ry(i),ires,vocal)
if (ires.eq.1) p%granite_n=i
if (ires.lt.0) p%granite_ry(i)=p%granite_ry(-ires)
call scanfile (p%infile,'granite_drho'//c2,p%granite_drho(i),ires,vocal)
if (ires.eq.1) p%granite_n=i
if (ires.lt.0) p%granite_drho(i)=p%granite_drho(-ires)
call scanfile (p%infile,'granite_dk'//c2,p%granite_dk(i),ires,vocal)
if (ires.eq.1) p%granite_n=i
if (ires.lt.0) p%granite_dk(i)=p%granite_dk(-ires)
call scanfile (p%infile,'granite_dkd'//c2,p%granite_dkd(i),ires,vocal)
if (ires.eq.1) p%granite_n=i
if (ires.lt.0) p%granite_dkd(i)=p%granite_dkd(-ires)
else
write (c3,'(i3)') i
call scanfile (p%infile,'granite_x'//c3,p%granite_x(i),ires,vocal)
if (ires.eq.1) p%granite_n=i
if (ires.lt.0) p%granite_x(i)=p%granite_x(-ires)
call scanfile (p%infile,'granite_y'//c3,p%granite_y(i),ires,vocal)
if (ires.eq.1) p%granite_n=i
if (ires.lt.0) p%granite_y(i)=p%granite_y(-ires)
call scanfile (p%infile,'granite_top'//c3,p%granite_top(i),ires,vocal)
if (ires.eq.1) p%granite_n=i
if (ires.lt.0) p%granite_top(i)=p%granite_top(-ires)
call scanfile (p%infile,'granite_bottom'//c3,p%granite_bottom(i),ires,vocal)
if (ires.eq.1) p%granite_n=i
if (ires.lt.0) p%granite_bottom(i)=p%granite_bottom(-ires)
call scanfile (p%infile,'granite_rx'//c3,p%granite_rx(i),ires,vocal)
if (ires.eq.1) p%granite_n=i
if (ires.lt.0) p%granite_rx(i)=p%granite_rx(-ires)
p%granite_ry(i)=p%granite_rx(i)
call scanfile (p%infile,'granite_ry'//c3,p%granite_ry(i),ires,vocal)
if (ires.eq.1) p%granite_n=i
if (ires.lt.0) p%granite_ry(i)=p%granite_ry(-ires)
call scanfile (p%infile,'granite_drho'//c3,p%granite_drho(i),ires,vocal)
if (ires.eq.1) p%granite_n=i
if (ires.lt.0)p%granite_drho(i)=p%granite_drho(-ires)
call scanfile (p%infile,'granite_dk'//c3,p%granite_dk(i),ires,vocal)
if (ires.eq.1) p%granite_n=i
if (ires.lt.0)p%granite_dk(i)=p%granite_dk(-ires)
call scanfile (p%infile,'granite_dkd'//c3,p%granite_dkd(i),ires,vocal)
if (ires.eq.1) p%granite_n=i
if (ires.lt.0)p%granite_dkd(i)=p%granite_dkd(-ires)
endif
enddo

p%ibc=0001
call scanfile (p%infile,'boundary_condition',p%ibc,ires,vocal)
p%restart=0
call scanfile (p%infile,'restart',p%restart,ires,vocal)
p%convert='native'
call scanfile (p%infile,'convert',p%convert,ires,vocal)
p%DEM='DEM'
call scanfile (p%infile,'DEM',p%DEM,ires,vocal)
p%base_level=1.d-6
call scanfile (p%infile,'base_level',p%base_level,ires,vocal)

p%uplift_n=-1
call scanfile (p%infile,'uplift_n',p%uplift_n,ires,vocal)
  if (p%uplift_n.eq.0) then
  p%uplift_v(1)=0.d0
  call scanfile (p%infile,'uplift_v1',p%uplift_v(1),ires,vocal)
  elseif (p%uplift_n.eq.1) then
  p%uplift_v(1:4)=0.d0
    do i=1,4
    write (c1,'(i1)') i
    call scanfile (p%infile,'uplift_v'//c1,p%uplift_v(i),ires,vocal)
    if (ires.lt.0) p%uplift_v(i)=p%uplift_v(-ires)
    enddo
  elseif (p%uplift_n.eq.2) then
  p%uplift_v(1:16)=0.d0
    do i=1,9
    write (c1,'(i1)') i
    call scanfile (p%infile,'uplift_v'//c1,p%uplift_v(i),ires,vocal)
    if (ires.lt.0) p%uplift_v(i)=p%uplift_v(-ires)
    enddo
    do i=10,16
    write (c2,'(i2)') i
    call scanfile (p%infile,'uplift_v'//c2,p%uplift_v(i),ires,vocal)
    if (ires.lt.0) p%uplift_v(i)=p%uplift_v(-ires)
    enddo
  endif
p%uplift_nt=0
call scanfile (p%infile,'uplift_nt',p%uplift_nt,ires,vocal)
  if (p%uplift_nt.gt.0) then
  if (p%uplift_nt.lt.2) stop 'uplift_nt needs to be gt 1'
  if (p%uplift_nt.gt.9) stop 'uplift_nt needs to be lt 10'
    do i=1,p%uplift_nt
    write (c1,'(i1)') i
    p%uplift_t(i)=p%nstep*p%dt*float(i)/p%uplift_nt
    call scanfile (p%infile,'uplift_t'//c1,p%uplift_t(i),ires,vocal)
    if (ires.lt.0) p%uplift_t(i)=p%uplift_t(-ires)
    p%uplift_f(i)=1.d0
    call scanfile (p%infile,'uplift_f'//c1,p%uplift_f(i),ires,vocal)
    if (ires.lt.0) p%uplift_f(i)=p%uplift_f(-ires)
    enddo
  endif

p%precipitation_n=-1
call scanfile (p%infile,'precipitation_n',p%precipitation_n,ires,vocal)
  if (p%precipitation_n.eq.0) then
  p%precipitation_v(1)=0.d0
  call scanfile (p%infile,'precipitation_v1',p%precipitation_v(1),ires,vocal)
  elseif (p%precipitation_n.eq.1) then
  p%precipitation_v(1:4)=0.d0
    do i=1,4
    write (c1,'(i1)') i
    call scanfile (p%infile,'precipitation_v'//c1,p%precipitation_v(i),ires,vocal)
    if (ires.lt.0) p%precipitation_v(i)=p%precipitation_v(-ires)
    enddo
  elseif (p%precipitation_n.eq.2) then
  p%precipitation_v(1:16)=0.d0
    do i=1,9
    write (c1,'(i1)') i
    call scanfile (p%infile,'precipitation_v'//c1,p%precipitation_v(i),ires,vocal)
    if (ires.lt.0) p%precipitation_v(i)=p%precipitation_v(-ires)
    enddo
    do i=10,16
    write (c2,'(i2)') i
    call scanfile (p%infile,'precipitation_v'//c2,p%precipitation_v(i),ires,vocal)
    if (ires.lt.0) p%precipitation_v(i)=p%precipitation_v(-ires)
    enddo
  endif
p%precipitation_nt=0
call scanfile (p%infile,'precipitation_nt',p%precipitation_nt,ires,vocal)
  if (p%precipitation_nt.gt.0) then
  if (p%precipitation_nt.lt.2) stop 'precipitation_nt needs to be gt 1'
  if (p%precipitation_nt.gt.9) stop 'precipitation_nt needs to be lt 10'
    do i=1,p%precipitation_nt
    write (c1,'(i1)') i
    p%precipitation_t(i)=p%nstep*p%dt*float(i)/p%precipitation_nt
    call scanfile (p%infile,'precipitation_t'//c1,p%precipitation_t(i),ires,vocal)
    if (ires.lt.0) p%precipitation_t(i)=p%precipitation_t(-ires)
    p%precipitation_f(i)=1.d0
    call scanfile (p%infile,'precipitation_f'//c1,p%precipitation_f(i),ires,vocal)
    if (ires.lt.0) p%precipitation_f(i)=p%precipitation_f(-ires)
    enddo
  endif

p%initial_topography_n=-1
call scanfile (p%infile,'initial_topography_n',p%initial_topography_n,ires,vocal)
  if (p%initial_topography_n.eq.0) then
  p%initial_topography_v(1)=0.d0
  call scanfile (p%infile,'initial_topography_v1',p%initial_topography_v(1),ires,vocal)
  elseif (p%initial_topography_n.eq.1) then
  p%initial_topography_v(1:4)=0.d0
    do i=1,4
    write (c1,'(i1)') i
    call scanfile (p%infile,'initial_topography_v'//c1,p%initial_topography_v(i),ires,vocal)
    if (ires.lt.0) p%initial_topography_v(i)=p%initial_topography_v(-ires)
    enddo
  elseif (p%initial_topography_n.eq.2) then
  p%initial_topography_v(1:16)=0.d0
    do i=1,9
    write (c1,'(i1)') i
    call scanfile (p%infile,'initial_topography_v'//c1,p%initial_topography_v(i),ires,vocal)
    if (ires.lt.0) p%initial_topography_v(i)=p%initial_topography_v(-ires)
    enddo
    do i=10,16
    write (c2,'(i2)') i
    call scanfile (p%infile,'initial_topography_v'//c2,p%initial_topography_v(i),ires,vocal)
    if (ires.lt.0) p%initial_topography_v(i)=p%initial_topography_v(-ires)
    enddo
  endif

p%orography = 0
call scanfile (p%infile,'orography',p%orography,ires,vocal)
p%cw=0.02d0
call scanfile (p%infile,'cw',p%cw,ires,vocal)
p%u=1.d0
call scanfile (p%infile,'u',p%u,ires,vocal)
p%v=1.d0
call scanfile (p%infile,'v',p%v,ires,vocal)
p%Nm=0.005d0
call scanfile (p%infile,'Nm',p%Nm,ires,vocal)
p%taus=500.d0
call scanfile (p%infile,'taus',p%taus,ires,vocal)
p%tauf=500.d0
call scanfile (p%infile,'tauf',p%tauf,ires,vocal)
p%pm=0.25d0
call scanfile (p%infile,'pm',p%pm,ires,vocal)
p%hw=2000.d0
call scanfile (p%infile,'hw',p%hw,ires,vocal)
p%noro=256
call scanfile (p%infile,'noro',p%noro,ires,vocal)

p%Pecube_nfreq=0
call scanfile (p%infile,'Pecube_nfreq',p%Pecube_nfreq,ires,vocal)
p%Pecube_nx=p%nx
call scanfile (p%infile,'Pecube_nx',p%Pecube_nx,ires,vocal)
p%Pecube_ny=p%ny
call scanfile (p%infile,'Pecube_ny',p%Pecube_ny,ires,vocal)
p%Pecube_xmin=0.d0
call scanfile (p%infile,'Pecube_xmin',p%Pecube_xmin,ires,vocal)
p%Pecube_xmax=p%xl
call scanfile (p%infile,'Pecube_xmax',p%Pecube_xmax,ires,vocal)
p%Pecube_ymin=0.d0
call scanfile (p%infile,'Pecube_ymin',p%Pecube_ymin,ires,vocal)
p%Pecube_ymax=p%yl
call scanfile (p%infile,'Pecube_ymax',p%Pecube_ymax,ires,vocal)

p%nage=0
call scanfile (p%infile,'nage',p%nage,ires,vocal)
p%L_age=120.d3
call scanfile (p%infile,'L_age',p%L_age,ires,vocal)
p%tmin_age=0.d0
call scanfile (p%infile,'tmin_age',p%tmin_age,ires,vocal)
p%tmax_age=1350.d0
call scanfile (p%infile,'tmax_age',p%tmax_age,ires,vocal)
p%heat_age=0.d0
call scanfile (p%infile,'heat_age',p%heat_age,ires,vocal)
p%zheat_age=35.d3
call scanfile (p%infile,'zheat_age',p%zheat_age,ires,vocal)
p%cond_age=1.d0
call scanfile (p%infile,'cond_age',p%cond_age,ires,vocal)
p%zcond_age=35.d3
call scanfile (p%infile,'zcond_age',p%zcond_age,ires,vocal)
p%limit_age=p%dt*p%nstep
call scanfile (p%infile,'limit_age',p%limit_age,ires,vocal)
p%dx_age=0.d0
call scanfile (p%infile,'dx_age',p%dx_age,ires,vocal)
p%dy_age=0.d0
call scanfile (p%infile,'dy_age',p%dy_age,ires,vocal)
  do i=1,p%nage
  p%x_age(i)=p%xl/2.d0
  p%y_age(i)=p%yl/2.d0
  p%h_age(i)=-99999.d0
  p%size_age(i)=100.d-6
  p%aft_age(i)=-999.d0
  p%aft_dage(i)=-999.d0
  p%zft_age(i)=-999.d0
  p%zft_dage(i)=-999.d0
  p%ahe_age(i)=-999.d0
  p%ahe_dage(i)=-999.d0
  p%zhe_age(i)=-999.d0
  p%zhe_dage(i)=-999.d0
  if (i.lt.10) then
  write (c1,'(i1)') i
  call scanfile (p%infile,'age_x'//c1,p%x_age(i),ires,vocal)
  call scanfile (p%infile,'age_y'//c1,p%y_age(i),ires,vocal)
  call scanfile (p%infile,'age_h'//c1,p%h_age(i),ires,vocal)
  call scanfile (p%infile,'age_size'//c1,p%size_age(i),ires,vocal)
  call scanfile (p%infile,'age_aft'//c1,p%aft_age(i),ires,vocal)
  call scanfile (p%infile,'age_daft'//c1,p%aft_dage(i),ires,vocal)
  call scanfile (p%infile,'age_zft'//c1,p%zft_age(i),ires,vocal)
  call scanfile (p%infile,'age_dzft'//c1,p%zft_dage(i),ires,vocal)
  call scanfile (p%infile,'age_ahe'//c1,p%ahe_age(i),ires,vocal)
  call scanfile (p%infile,'age_dahe'//c1,p%ahe_dage(i),ires,vocal)
  call scanfile (p%infile,'age_zhe'//c1,p%zhe_age(i),ires,vocal)
  call scanfile (p%infile,'age_dzhe'//c1,p%zhe_dage(i),ires,vocal)
  elseif (i.lt.100) then
  write (c2,'(i2)') i
  call scanfile (p%infile,'age_x'//c2,p%x_age(i),ires,vocal)
  call scanfile (p%infile,'age_y'//c2,p%y_age(i),ires,vocal)
  call scanfile (p%infile,'age_h'//c2,p%h_age(i),ires,vocal)
  call scanfile (p%infile,'age_size'//c2,p%size_age(i),ires,vocal)
  call scanfile (p%infile,'age_aft'//c2,p%aft_age(i),ires,vocal)
  call scanfile (p%infile,'age_daft'//c2,p%aft_dage(i),ires,vocal)
  call scanfile (p%infile,'age_zft'//c2,p%zft_age(i),ires,vocal)
  call scanfile (p%infile,'age_dzft'//c2,p%zft_dage(i),ires,vocal)
  call scanfile (p%infile,'age_ahe'//c2,p%ahe_age(i),ires,vocal)
  call scanfile (p%infile,'age_dahe'//c2,p%ahe_dage(i),ires,vocal)
  call scanfile (p%infile,'age_zhe'//c2,p%zhe_age(i),ires,vocal)
  call scanfile (p%infile,'age_dzhe'//c2,p%zhe_dage(i),ires,vocal)
  elseif (i.lt.1000) then
  write (c3,'(i3)') i
  call scanfile (p%infile,'age_x'//c3,p%x_age(i),ires,vocal)
  call scanfile (p%infile,'age_y'//c3,p%y_age(i),ires,vocal)
  call scanfile (p%infile,'age_h'//c3,p%h_age(i),ires,vocal)
  call scanfile (p%infile,'age_size'//c3,p%size_age(i),ires,vocal)
  call scanfile (p%infile,'age_aft'//c3,p%aft_age(i),ires,vocal)
  call scanfile (p%infile,'age_daft'//c3,p%aft_dage(i),ires,vocal)
  call scanfile (p%infile,'age_zft'//c3,p%zft_age(i),ires,vocal)
  call scanfile (p%infile,'age_dzft'//c3,p%zft_dage(i),ires,vocal)
  call scanfile (p%infile,'age_ahe'//c3,p%ahe_age(i),ires,vocal)
  call scanfile (p%infile,'age_dahe'//c3,p%ahe_dage(i),ires,vocal)
  call scanfile (p%infile,'age_zhe'//c3,p%zhe_age(i),ires,vocal)
  call scanfile (p%infile,'age_dzhe'//c3,p%zhe_dage(i),ires,vocal)
  else
  write (c4,'(i4)') i
  call scanfile (p%infile,'age_x'//c4,p%x_age(i),ires,vocal)
  call scanfile (p%infile,'age_y'//c4,p%y_age(i),ires,vocal)
  call scanfile (p%infile,'age_h'//c4,p%h_age(i),ires,vocal)
  call scanfile (p%infile,'age_size'//c4,p%size_age(i),ires,vocal)
  call scanfile (p%infile,'age_aft'//c4,p%aft_age(i),ires,vocal)
  call scanfile (p%infile,'age_daft'//c4,p%aft_dage(i),ires,vocal)
  call scanfile (p%infile,'age_zft'//c4,p%zft_age(i),ires,vocal)
  call scanfile (p%infile,'age_dzft'//c4,p%zft_dage(i),ires,vocal)
  call scanfile (p%infile,'age_ahe'//c4,p%ahe_age(i),ires,vocal)
  call scanfile (p%infile,'age_dahe'//c4,p%ahe_dage(i),ires,vocal)
  call scanfile (p%infile,'age_zhe'//c4,p%zhe_age(i),ires,vocal)
  call scanfile (p%infile,'age_dzhe'//c4,p%zhe_dage(i),ires,vocal)
  endif
  enddo

p%reference_surface=0
call scanfile (p%infile,'reference_surface',p%reference_surface,ires,vocal)

open (7,file=p%infile,status='unknown')
read (7,*,end=999)
goto 998
999 print*,'Empty or no input file'
rewind (7)
print*,'It has been filled with default parameter values'
write (7,'(a42)') 'DEFAULT INPUT FILE GENERATED BY FastScape'
write (7,'(a)') ' '
write (7,'(a)') 'Basic geometry'
write (7,'(a)') ' '
write (7,'(a5,i6)') 'nx = ',p%nx
write (7,'(a5,i6)') 'ny = ',p%ny
write (7,'(a5,g14.6)') 'xl = ',p%xl
write (7,'(a5,g14.6)') 'yl = ',p%yl
write (7,'(a)') ' '
write (7,'(a)') 'Time stepping'
write (7,'(a)') ' '
write (7,'(a5,g14.6)') 'dt = ',p%dt
write (7,'(a8,i6)') 'nstep = ',p%nstep
write (7,'(a8,i6)') 'nfreq = ',p%nstep/10
write (7,'(a)') ' '
write (7,'(a)') 'Number of threads/cores to use'
write (7,'(a)') ' '
write (7,'(a14,i6)') 'num_threads = ',p%num_threads
write (7,'(a)') ' '
write (7,'(a)') 'Incision law'
write (7,'(a)') ' '
write (7,'(a6,i1)') 'law = ',l%law
write (7,'(a4,g14.6)') 'm = ',l%m
write (7,'(a5,g14.6)') 'kf = ',l%kf
write (7,'(a)') ' '
write (7,'(a)') 'Boundary conditions'
write (7,'(a)') ' '
write (c4,'(i4)') p%ibc
if (p%ibc.lt.10) c4(1:3)='000'
if (p%ibc.lt.100) c4(1:2)='00'
if (p%ibc.lt.1000) c4(1:1)='0'
write (7,'(a,i4)') 'boundary_condition = '//c4
write (7,'(a)') ' '
write (7,'(a)') 'Precipitations'
write (7,'(a)') ' '
write (7,'(a,i1)') 'precipitation_n = ',0
write (7,'(a,g14.6)') 'precipitation_v1 = ',1.d0
write (7,'(a)') ' '
write (7,'(a)') 'Initial topography'
write (7,'(a)') ' '
write (7,'(a,i1)') 'initial_topography_n = ',1
write (7,'(a,g14.6)') 'initial_topography_v1 = ',0.d0
write (7,'(a,g14.6)') 'initial_topography_v2 = ',100.d0
write (7,'(a,g14.6)') 'initial_topography_v3 = ',100.d0
write (7,'(a,g14.6)') 'initial_topography_v4 = ',0.d0
write (7,'(a)') ' '
write (7,'(a)') 'Uplift rate'
write (7,'(a)') ' '
write (7,'(a,i1)') 'uplift_n = ',1
write (7,'(a,g14.6)') 'uplift_v1 = ',0.d0
write (7,'(a,g14.6)') 'uplift_v2 = ',1.d-4
write (7,'(a,g14.6)') 'uplift_v3 = ',1.d-4
write (7,'(a,g14.6)') 'uplift_v4 = ',0.d0
write (7,'(a)') ' '
write (7,'(a)') 'Plotting options'
write (7,'(a)') ' '
write (7,'(a,i1)') 'plot_all = ',0
write (7,'(a,i1)') 'plot_topo = ',1

write (*,'(a)') p%infile//' created - Ready to run FastScape now'
stop

998 close (7)

end subroutine initialize
