  subroutine InterpolStepwise(a,nxa,nya,b,nx,ny)

! Stepwise interpolation (fast but rough)

  use omp_lib

  integer nxa,nya,nx,ny
  double precision, dimension(nxa,nya), intent(in) :: a
  double precision, dimension(nx,ny), intent(out) :: b

  integer i,j,ia,ja

!$OMP parallel shared(nx,ny,nxa,nya,a,b)
!$OMP do private(i,j,ia,ja)
  do j=1,ny
  ja=1+((nya-1)*(j-1))/(ny-1)
    do i=1,nx
    ia=1+((nxa-1)*(i-1))/(nx-1)
    b(i,j)=a(ia,ja)
    enddo
  enddo
!$OMP end do nowait
!$OMP end parallel

  end subroutine InterpolStepwise

!!!!!!!!!!!!!!!!!!

  subroutine InterpolLinear(a,nxa,nya,b,nx,ny)

! Linear interpolation (less fast but smoother)

  use omp_lib

  integer nxa,nya,nx,ny
  double precision, dimension(nxa,nya), intent(in) :: a
  double precision, dimension(nx,ny), intent(out) :: b

  integer i,j,ia,ja
  double precision r,s

!$OMP parallel shared(nx,ny,nxa,nya,a,b)
!$OMP do private(i,j,ia,ja,r,s)
  do j=1,ny
  ja=1+((nya-1)*(j-1))/(ny-1)
  s=(float(j-1)/(ny-1)-float(ja-1)/(nya-1))/(nya-1)*2.d0-1.d0
    if (ja.eq.nya) then
    ja=nya-1
    s=1.d0
    endif
    do i=1,nx
    ia=1+((nxa-1)*(i-1))/(nx-1)
    r=(float(i-1)/(nx-1)-float(ia-1)/(nxa-1))/(nxa-1)*2.d0-1.d0
      if (ia.eq.nxa) then
      ia=nxa-1
      r=1.d0
      endif
    b(i,j)=(a(ia,ja)*(1.d0-r)*(1.d0-s) &
          +a(ia+1,ja)*(1.d0+r)*(1.d0-s) &
          +a(ia+1,ja+1)*(1.d0+r)*(1.d0+s) &
          +a(ia,ja+1)*(1.d0-r)*(1.d0+s))/4.d0
    enddo
  enddo
!$OMP end do nowait
!$OMP end parallel

  end subroutine InterpolLinear

!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine InterpolLinearBox(a,nxa,nya,b,nx,ny,x1,x2,y1,y2)

! Linear interpolation (less fast but smoother)

  use omp_lib

  integer nxa,nya,nx,ny
  double precision, dimension(nxa,nya), intent(in) :: a
  double precision, dimension(nx,ny), intent(out) :: b
  double precision, intent(in) :: x1,x2,y1,y2
  integer i,j,ia,ja
  double precision r,s,x,y

!$OMP parallel shared(nx,ny,nxa,nya,a,b)
!$OMP do private(i,j,ia,ja,r,s)
  do j=1,ny
  y=y1+(y2-y1)*float(j-1)/(ny-1)
  ja=1+int((nya-1)*y)
  s=(y-float(ja-1)/(nya-1))/(nya-1)*2.d0-1.d0
    if (ja.eq.nya) then
    ja=nya-1
    s=1.d0
    endif
    do i=1,nx
    x=x1+(x2-x1)*float(i-1)/(nx-1)
    ia=1+int((nxa-1)*x)
    r=(x-float(ia-1)/(nxa-1))/(nxa-1)*2.d0-1.d0
      if (ia.eq.nxa) then
      ia=nxa-1
      r=1.d0
      endif
    b(i,j)=(a(ia,ja)*(1.d0-r)*(1.d0-s) &
          +a(ia+1,ja)*(1.d0+r)*(1.d0-s) &
          +a(ia+1,ja+1)*(1.d0+r)*(1.d0+s) &
          +a(ia,ja+1)*(1.d0-r)*(1.d0+s))/4.d0
    enddo
  enddo
!$OMP end do nowait
!$OMP end parallel

  end subroutine InterpolLinearBox

!!!!!!!!!!!!!!!!!!!!!

  subroutine InterpolSmooth(a,nxa,nya,b,nx,ny)

! Smoothing interpolation (less fast but smoother)

  use omp_lib

  integer nxa,nya,nx,ny
  double precision, dimension(nxa,nya), intent(in) :: a
  double precision, dimension(nx,ny), intent(out) :: b

  integer i,j,ia,ja,ib,jb
  double precision r,s

!$OMP parallel shared(nx,ny,nxa,nya,a,b)
!$OMP do private(i,j,ia,ja,r,s)
  do j=1,ny
  ja=1+((nya-1)*(j-1))/(ny-1)
  if (ja.eq.nya) ja=nya-1
  jb=1+((nya-1)*j)/(ny-1)
  if (jb.gt.nya) jb=nya
    do i=1,nx
    ia=1+((nxa-1)*(i-1))/(nx-1)
    if (ia.eq.nxa) ia=nxa-1
    ib=1+((nxa-1)*i)/(nx-1)
    if (ib.gt.nxa) ib=nxa
    b(i,j)=sum(a(ia:ib,ja:jb))/size(a(ia:ib,ja:jb))
    enddo
  enddo
!$OMP end do nowait
!$OMP end parallel

  end subroutine InterpolSmooth
