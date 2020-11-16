module definitions

!=====[topography]====================================================================
! this type carries the topography with its related arrays
! h                = topographic height (m)
! hb               = basement height (m)
! hiso             = cumulated isostatic response (m)
! hi               = initial topography (m)
! hi_iso           = initial topography for isostasy (m)
! hp               = topography at the beginning of previous time step
! discharge        = volumetric discharge (m^3/yr)
! length           = length to receiver node (m)
! width            = channel width (m)
! sediment         = sediment volumetric discharge (m^3/yr)
! pmask            = active process mask
! s                = soil thickness (m)

      type topography
      double precision, dimension(:), allocatable :: h,hb,hiso,hi,hp,s,u,hi_iso,href
      double precision, dimension(:), allocatable :: discharge,length,width,sediment
      double precision, dimension(:), allocatable :: chi,fft
      logical , dimension(:,:), allocatable :: pmask
      end type topography

!=====[topology]====================================================================
! this type carries the topology or connectivity of the landscape
! receiver         = receiver node (where the water flows to)
! ndon             = number of donor nodes (from which the water is coming from)
! donors           = list of donor nodes
! catchment        = catchment number
! stack            = stack or order in which the nodes are processed for fluvial processes
! nstackp          = length of stack
! bc               = boundary condition (=0 means free node; =1 means base level node)

      type topology
      integer, dimension(:), allocatable :: receiver,ndon,donors,catchment,stack
      logical, dimension(:), allocatable :: bc
      integer nstackp
      end type topology

!=====[law]====================================================================
! this type carries the erosion law parameters
! law              = law to be used (0 is linear stream power law)
!                                   (1 is non linear stream power law)
!                                   (2 is Davy's Q-Ksi law)
! nk               = number of iterations to use in the implicit scheme
!                    when the law is non linear in slope; note that a
!                    convergence criterion is used when nk=<0
! m                = discharge exponent
! n                = slope exponent
! kf               = fluvial K
! sc               = maximum slope angle
! dsc              = variability in sc as a fraction of sc
! kfs              = fluvial kf for sediments (only for law=3)
! tol              = tolerance used in the Taylor algortithm to find root (m)
! tauc             = critical erosion rate (shear stress) used in Davy's algorithm
! lf               = reaction constant in Beaumont's algorithm
! kd               = diffusion coefficient (m^2/yr)
! p0               = soil production at zero soil thickness (m/yr)
! s0               = exponential decay depth for soil production function (m)
! slope            = threshold slope (Josh's algorithm to mimic landsliding) (in deg)
! strati_n         = number of stratigraphic layers
! strati_f         = multiplying factors for fluvial erosion constant (kf)
! strati_fd        = multiplying factors for diffusivity (kd, transport coefficient)
! strati_top       = top and bottom define the stratigraphic intervals
! strati_bottom    = over which the multiplying factor should be applied to the
!                    fluvial K constant
! area_min         = minimum drainage area for fluvial erosion

      type law
      integer law,nk
      double precision m,n,kf,kfs,tol,tauc,lf,kd,sc,p0,s0,slope,dsc,area_min
      double precision strati_top(9),strati_bottom(9),strati_f(9),strati_fd(9)
      integer strati_n
      end type law

!=====[param]====================================================================
! this type carries all "general" parameters defining the model run
! nx              = grid dimension in x-direction
! ny              = grid dimension in y-direction
! nx0             = size given in the input file (different from nx if it is odd)
! ny0             = size given in the input file (different from ny if it is odd)
! nn              = total number of nodes (=nx*ny)
! nstep           = number of time steps
! nfreq           = output frequency
! nmetric         = metric calculation frequency
! ndemoulin       = Demoulin's index calculation frequency
! dx              = mesh step length in x-direction (m)
! dy              = mesh step length in y-direction (m)
! xl              = mesh length in x-direction (=(nx-1)*dx) (m)
! yl              = mesh length in y-direction (=(ny-1)*dy) (m)
! dt              = time step length (yr)
! time            = current time (yr)
! num_threads     = number of OPENMP threads
! chunk           = number of loop instances to be performed per each trhead
! debug           = debugging flag
! quiet           = set true to kill all output to window (useful in inversion mode)
! infile          = name of input file
! run             = name of run directory
! ibc             = integer construction to decide on the boundary conditions
! imetric         = integer construction to decide on which metrics to output
! flexure         = logical flag to turn on flexure calculations
! thickflex       = elastic plate thicness (m)
! ym              = Young Modulus (Pa)
! pratio          = Poisson's ratio
! rhocflex        = crustal density (kg/m^3)
! rhoaflex        = asthenospheric density (kg/m^3)
! meanflex        = flag to control averageing in either x- or y- direction depending
!                   on boundary conditions (if two opposite sides are set to base level
!                   in the bc, the load is automatically averaged in the corresponding
!                   direction if meanflex = .TRUE. only
! refflex         = flag (1 for thicknening rate - 0 for uplift rate)
! sea_level       = sea-level (something different happens below sea-level
!                              in terms of erosion/deposition) (m)
! local_minima    = logical flag to decide whether to resolve the local minima
! plot_topo       = integer flag to decide whether to build topo image files
! plot_strati     = integer flag to decide whether to build strati image files
! plot_hardness   = integer flag to decide whether to build hardness image files (granite hardness)
! plot_density    = integer flag to decide whether to build density image files (granite density)
! plot_erosion    = integer flag to decide whether to build total erosion image files
! plot_rate       = integer flag to decide whether to build erosion rate image files
! plot_uplift     = integer flag to decide whether to build surface uplift rate image files
! plot_rock_uplift= integer flag to decide whether to build rock uplift rate image files
! plot_sedim_flux = sedimentary flux (law=3 only)
! plot_sedim      = integer flag to decide whether to build total sediment thickness image files
! plot_discharge  = integer flag to decide whether to build discharge image files
! plot_catchment  = integer flag to decide whether to build catchment image files
! plot_precip     = integer flag to decide whether to build precipitation image files
! plot_slope      = integer flag to decide whether to build slope image files
! plot_curvature  = integer flag to decide whether to build curvature image files
! plot_steepnessindex = integer flag to decide whether to build steepnessindex image files
! plot_concavity  = integer flag to decide whether to build concavity image files
!                   concerning all the plot_xxxx flag they can be either -1 0 or 1
!                   if 0, no plot is produced, if 1 or -1 the corresponding plot is produced
!                   if negative, the range used to color the field is derived from the
!                   min-max values of field at each time step
! plot_volume_mask = integer flag to decide whether to build volume mask image files (only applies
!                    when a box is specified in a Volume* file
! plot_all        = logical flag to make all plots unless specified otherwise or to make only the
!                   plot specified in the input file
! vtk             = flag to decide whether vtk output is generated (=1) or not (=0)
! plot_DEM        = flag to decide whether topographic dem is generated (=1) or not (=0)
!                   (if =2 then the DEM is stored as int*4; =3 as real*4; =4 as real*8)
! vex             = vertical exaggeration when creating the 3D perspective vtk files
!                   or for the shading of the image files (.bmp)
! topo_min        = minimum value for topography in bmp image file
! topo_max        = maximum value for topography in bmp image file
! strati_min      = minimum value for strati in bmp image file
! strati_max      = maximum value for strati in bmp image file
! hardness_min    = minimum value for hardness in bmp image file
! hardness_max    = maximum value for hardness in bmp image file
! density_min     = minimum value for density in bmp image file
! density_min     = maximum value for density in bmp image file
! erosion_min     = minimum value for total erosion in bmp image file
! erosion_max     = maximum value for total erosion in bmp image file
! rate_min        = minimum value for the erosion rate in bmp image file
! rate_max        = maximum value for the erosion rate in bmp image file
! uplift_min      = minimum value for the surface uplift rate in bmp image file
! uplift_max      = maximum value for the surface uplift rate in bmp image file
! rock_uplift_min = minimum value for the rock uplift rate in bmp image file
! rock_uplift_max = maximum value for the rock uplift rate in bmp image file
! sedim_min       = minimum value for sediment thickness in bmp image file
! sedim_max       = maximum value for sediment thickness in bmp image file
! discharge_min   = minimum value for discharge in bmp image file
! discharge_max   = maximum value for discharge in bmp image file
! precip_max      = minimum value for precipitation in bmp image file
! precip_max      = maximum value for precipitation in bmp image file
! slope_min       = minimum slope value for metric (slope vs area and flux vs slope)
! slope_max       = maximum slope value for metric (slope vs area and flux vs slope)
! curvature_min   = minimum curvature (in per meter)
! curvature_max   = maximum curvature (in per meter)
! steepnessindex_min = minimum steepness index
! steepnessindex_max = maximum steepness index
! concavity_min   = minimum concavity
! concavity_max   = maximum concavity
! HeightDis_min   = minimum value for Height CDF plot
! HeightDis_max   = maximum value for Height CDF plot
! SlopeDis_min    = minimum value for Slope CDF plot
! SlopeDis_max    = maximum value for Slope CDF plot
! CurveDis_min    = minimum value for Curvature CDF plot
! CurveDis_max    = maximum value for Curvature CDF plot
! DischargeDis_min = minimum value for Discharge CDF plot
! DischargeDis_max = maximum value for Discharge CDF plot
! stream_threshold_area = minimum drainage area for a pixel to be considered a stream
! minimum_area    = minimum drainage area needed to compute Demoulin's R1 index
! restart         = flag to decide whether it is a restart job or not
!                   =0 means not a restart (new job)
!                   =1 means restart job and FastScape will read the initial topography from the
!                      RESTART file stored in the working directory
!                   It is also used to read in external initial topographies from DEM files
!                   =-1 reads in initial topo in integer*2 format
!                   =-2 reads in initial topo in integer*4 format
!                   =-3 reads in initial topo in real*4 format
!                   =-4 reads in initial topo in double precision format
! DEM             nale of DEM file containing initial topography
! convert         = conversion option for reading in direct access binary files:
!                   = native (default): no conversion
!                   = other choices are (big_endian, cray, fdx, fgx, ibm, little_endian, vaxd, vaxg)
!                   see ifort man page for further details
! reset_random    = to enable reseting the random number generator creating the initial topography
!                   (=0, default, same random series is used every time; =1 the random set is different
!                    for each run)
! istep           = current time step
! orography       = flag to turn on orographic control of precipitation
!                   (=0 no orography - =1 orography)
! cw              = uplift sensitivity factor (kg/m3) (orographic model)
! u               = E-W wind velocity (m/s) (orographic model)
! v               = S-N wind velocity (m/s)
! Nm              = moist stability frequency (1/s)
! taus            = conversion time (s)
! tauf            = fallout time (s)
! pm              = mean precipitation (m/yr)
! hw              = water vapour scale height (m)
! noro            = resolution of FFT
! granite_xi      = x-position of center of granite i (in m)
! granite_yi      = y-position of center of granite i (in m)
! granite_topi    = top of granite (in positive depth) (in m)
! granite_bottomi = bottom of granite (in positive depth) (in m)
! granite_rxi     = radius of granite i in x-direction (a rectangle if negative) (in m)
! granite_ryi     = radius of granite i in y-direction (a rectangle if negative) (in m) (if omitted ry=rx is assumed)
! granite_drhoi   = density anomaly (added to reference density for the crust) for granite i (in kg/m3)
! granite_dki     = erodibility factor (multiplying the reference erodibility kf) for granite i (no dim)
! granite_dkdi    = diffusivity (transport coefficient) (multiplying ref diffusivity) for granite i (no dim)
! granite_n       = number of granites (max is 100)
! n_Pecube        = number of Pecube output so far (internal variable)
! Pecube_nfreq    = flag to indicate the frequency at which output is saved for use as input to Pecube
!                   if =0 (default) no output is saved
! Pecube_nx       = x-discretization used for output to Pecube (default is nx)
! Pecube_ny       = y-discretization used for output to Pecube (default is ny)
! Pecube_xmin     = x-coordinate of bottom left corner of window extracted for output to Pecube (default is 0)
! Pecube_xmax     = x-coordinate of top right corner of window extracted for output to Pecube (default is xl)
! Pecube_ymin     = y-coordinate of bottom left corner of window extracted for output to Pecube (default is 0)
! Pecube_ymax     = y-coordinate of top right corner of window extracted for output to Pecube (default is xl)
! plot_chi        = flag to compute and plot chi (see Willett et al, Science (2014) for definition)
!                 = 0 chi is not calculated
!                 = 1 simplest form of chi is computed every nfreq time step
!                 = 2 full definition of chi is computed (ie taking variations in both U and kf into account)
!                     every nfreq time step
! plot_fft        = flag to compute and plot 2D fft of landscape
! discharge_min = minimum discharge (or drainage area if precipitation rate is unifor andset to 1)
!                     for the calculation of chi
! base_level      = elevation below which all nodes are set at base level (only works when restart<0
!                   ie when a DEM is used as initial topography
! nage            = number of locations where low-T ages are to be calculated
! x_age(i)        = x-position of location where age(i) is to be calculated
! y_age(i)        = y-position of location where age(i) is to be calculated
! size_age(i)     = grain size of sample i in m
! aft_age(i)      = observed/measured apatite fission track age at location (age_xi,age_yi) in yr
! aft_dage(i)     = uncertainty/error on observed apatite fission track age at location (age_xi,age_yi) in yr
! zft_age(i)      = observed/measured zircon fission track age at location (age_xi,age_yi) in yr
! zft_dage(i)     = uncertainty/error on observed zircon fission track age at location (age_xi,age_yi) in yr
! ahe_age(i)      = observed/measured apatite helium age at location (age_xi,age_yi) in yr
! ahe_dage(i)     = uncertainty/error on observed apatite helium age at location (age_xi,age_yi) in yr
! zhe_age(i)      = observed/measured zircon helium age at location (age_xi,age_yi) in yr
! zhe_dage(i)     = uncertainty/error on observed zircon helium age at location (age_xi,age_yi) in yr
! L_age           = thickness of thermal model in m (where the basal boundary condition is imposed)
! tmin_age        = surface temperature in deg C
! tmax_age        = basal temperature in deg C
! heat_age        = heat production in deg C/Myr
! zheat_age       = initial thickness of heat producing layer in m
! cond_age        = anomalous diffusivity/conductivity (multiplying factor)
! zcond_age       = initial thickness of anomalous conductivity layer
! limit_age       = maximum age (in yr)
! reference_surface = if 1 the reference surface is the initial topography
!                     if 0 the reference surface is a flat Earth, regardless of initial topography
! nbox            = number of boxes where sed flux is computed

      type param
      integer nage
      double precision x_age(1000),y_age(1000),size_age(1000),h_age(1000)
      double precision aft_age(1000),zft_age(1000),ahe_age(1000),zhe_age(1000)
      double precision aft_dage(1000),zft_dage(1000),ahe_dage(1000),zhe_dage(1000)
      double precision L_age,tmin_age,tmax_age,heat_age,zheat_age,cond_age,zcond_age,limit_age,dx_age,dy_age
      integer nx,ny,nx0,ny0,nn,nstep,istep,nfreq,num_threads,ibc,nmetric,imetric,restart,ndemoulin,chunk
      double precision dx,dy,xl,yl,dt,time,stream_threshold_area,minimum_area
      double precision thickflex,ym,pratio,rhocflex,rhoaflex,sea_level
      logical debug,flexure,meanflex,local_minima,quiet,refflex
      integer plot_topo,plot_erosion,plot_sedim,plot_discharge,plot_catchment,plot_precip
      integer plot_strati,plot_rate,plot_uplift,plot_rock_uplift,plot_sedim_flux,plot_slope,plot_curvature
      integer plot_hardness,plot_density,plot_chi,plot_fft
      integer plot_concavity,plot_steepnessindex,plot_DEM,plot_volume_mask
      double precision HeightDis_min,HeightDis_max,SlopeDis_min,SlopeDis_max
      double precision CurveDis_min,CurveDis_max,DischargeDis_min,DischargeDis_max
      logical plot_all,vtk
      double precision topo_min,topo_max,erosion_min,erosion_max,sedim_min,sedim_max
      double precision discharge_min,discharge_max,precip_min,precip_max
      double precision strati_min,strati_max,rate_min,rate_max,uplift_min,uplift_max,rock_uplift_min,rock_uplift_max
      double precision hardness_min,hardness_max,density_min,density_max
      double precision slope_min,slope_max,curvature_min,curvature_max,steepnessindex_min,steepnessindex_max
      double precision concavity_min,concavity_max
      character infile*40,run*5,convert*13,DEM*40
      integer uplift_n,uplift_nt,precipitation_n,precipitation_nt
      integer initial_topography_n
      double precision uplift_v(16),uplift_t(9),uplift_f(9)
      double precision precipitation_v(16),precipitation_t(9),precipitation_f(9)
      double precision initial_topography_v(16)
      double precision vex
      integer orography,noro
      double precision cw,u,v,Nm,taus,tauf,pm,hw
      integer granite_n,reset_random
      double precision granite_x(100),granite_y(100),granite_rx(100),granite_ry(100),granite_drho(100),granite_dk(100)
      double precision granite_top(100),granite_bottom(100),granite_dkd(100)
      integer Pecube_nfreq,Pecube_nx,Pecube_ny,nPecube
      double precision Pecube_xmin,Pecube_xmax,Pecube_ymin,Pecube_ymax
      double precision discharge_min_chi,base_level
      integer reference_surface,nbox
      end type param

!=====[box]====================================================================
! this type carries masks for sed flux computation
! mask            = mask for sed flux computation

      type box
      logical , dimension(:), allocatable :: mask
      character*128 fnme
      double precision x(4),y(4),timei(1000),timef(1000),volume(1000),dvolume(1000)
      double precision volume_predicted(1000)
      integer ntime
      end type box

!=====[INTERFACE SCANFILE]======================================================
! following is a general interface to read stuff from a file

      interface scanfile
        subroutine iscanfile (fnme,text,res,ires,vocal)
        character*(*) fnme,text
        integer,intent(out)::res
        integer,intent(out)::ires
        integer,intent(in)::vocal
        end subroutine iscanfile
        subroutine dscanfile (fnme,text,res,ires,vocal)
        character*(*) fnme,text
        double precision,intent(out)::res
        integer,intent(out)::ires
        integer,intent(in)::vocal
        end subroutine dscanfile
        subroutine cscanfile (fnme,text,res,ires,vocal)
        character*(*) fnme,text
        character*(*),intent(out)::res
        integer,intent(out)::ires
        integer,intent(in)::vocal
        end subroutine cscanfile
      end interface

end module definitions
