subroutine compute_age (time,temperature,nstep,sage,size,age_limit,age)

! subroutine to compute ages for four low-T systems (AFT,ZFT,AHE,ZHE)
! time and temperature contain the temperature history
! made of nstep starting by the largest time (oldest age)
! sage is the age system to consider ("AFT","ZFT","AHE","ZHE")
! size is the grain size
! age_limit is maximum age
! age is resulting age

! time, age_limit and age are in years
! temperature is in deg C
! size is in m

implicit none

double precision time(nstep),temperature(nstep),age_limit,age,size
character sage*3
integer nstep
integer i,ntime
double precision tc,r
real*4, dimension(:), allocatable :: time_i,temp_i
real*4 age_i,size_i
real*4  ftld(17),ftldmean,ftldsd

allocate (time_i(nstep+1),temp_i(nstep+1))

if (time(1).lt.age_limit) then
time_i(1)=age_limit/1.d6
temp_i(1)=temperature(1)
  do i=1,nstep
  time_i(1+i)=time(i)/1.d6
  temp_i(1+i)=temperature (i)
  enddo
ntime=nstep+1
else
  do i=1,nstep
  time_i(i)=time(i)/1.d6
  temp_i(i)=temperature(i)
  enddo
ntime=nstep
endif

size_i=size

select case (sage)
  case("AFT")
  call Mad_Trax (time_i,temp_i,ntime,0,2, &
                 age_i,ftld,ftldmean,ftldsd)
  case("ZFT")
  call Mad_Zirc (time_i,temp_i,ntime,0,2, &
                 age_i,ftld,ftldmean,ftldsd)
  case("AHE")
!  do i=1,ntime
!  print*,time_i(i),temp_i(i)
!  enddo
  call Mad_He (time_i,temp_i,ntime,age_i,1,size_i)
!  print*,age_i,size_i
  case("ZHE")
  call Mad_He (time_i,temp_i,ntime,age_i,2,size_i)
  case default
  stop 'this thermochronometer is not implemented yet'
end select

age=age_i*1.d6

deallocate (time_i,temp_i)

end subroutine compute_age
