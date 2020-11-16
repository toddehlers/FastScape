subroutine bmp (h,f,nx,ny,fnme,nfnme,ncol,cmap,hmi,hma,flag,vex)

! subroutine to create a bmp image file
! suffers from the need for nx and ny to be even numbers (I need to improve this)
! an option has been added to create a vtk file is requested

character fnme*(*)
integer nfnme,nx,ny,i,j,k,icol,jp,jm,ip,im,flag,cmap
double precision h(nx,ny),f(nx,ny),hmi,hma,fmin,fmax
character*1,dimension(:,:,:),allocatable::ima
double precision hmin,hmax,grey,cose,cosi,dhx,dhy,den,sqrt2,dx,dy,fact
integer,dimension(:,:),allocatable::rgb_col
double precision rann(3),vex

allocate (ima(3,nx,ny),rgb_col(ncol,3))

sqrt2=sqrt(2.d0)/2.d0
  if (flag.eq.1) then
  fmin=hmi
  fmax=hma
  else
  fmin=minval(f)
  fmax=maxval(f)
  endif
hmin=minval(h)
hmax=maxval(h)
dx=(hmax-hmin)/vex
dy=dx

!call rgb_spectrum (rgb_col,ncol)
if (cmap.eq.0) call topo_cmap (rgb_col,ncol)
if (cmap.eq.1) then
  do i=1,ncol
  call random_number (rann)
  rgb_col(i,:)=int(rann(:)*255.)
  enddo
endif

  do j=1,ny
  jp=min(j+1,ny)
  jm=max(j-1,1)
  do i=1,nx
  ip=min(i+1,nx)
  im=max(i-1,1)
    if (fmin.lt.0. .and. fmax.gt.0.) then
    icol=1+(ncol-1)*f(i,j)/fmax
    else
    icol=1+(ncol-1)*(f(i,j)-fmin)/(fmax-fmin)
    endif
  icol=min(icol,ncol)
  icol=max(icol,1)
  dhx=h(ip,j)-h(im,j)
  dhy=h(i,jp)-h(i,jm)
  den=sqrt(dhx**2*dy**2+dhy**2*dx**2+dx**2*dy**2)
  cose=dx*dy/den
  cosi=(sqrt2*dx*dy-dhx*dy/2.d0-dhy*dx/2.d0)/den
  grey=(1.d0+1.d0/(1.d0+cose/cosi))/2.d0
  grey=max(0.d0,grey)
  grey=min(1.d0,grey)
  do k=1,3
  ima(k,i,j)=char(int(rgb_col(icol,k)*grey))
  enddo
    if (fmin.lt.0. .and. fmax.gt.0. .and. f(i,j).lt.0.) then
    fact=(f(i,j)-fmin)/(0.-fmin)
    ima(1,i,j)=char(int(50.*fact*grey))
    ima(2,i,j)=char(int(70.*fact*grey))
    ima(3,i,j)=char(int(255.*fact*grey))
    endif
  enddo
  enddo

open (3,file=fnme(1:nfnme),status='unknown')
call pixout (ima,3,nx,ny)
close (3)

end subroutine bmp

!-----------
subroutine rgb_spectrum (rgb,M)

integer rgb (M,3)

       MAX=255
       GAMMA=.80

        DO I=1,M

            WL = 380. + REAL(I * 400. / M)

            IF ((WL.GE.380.).AND.(WL.LE.440.)) THEN 
              R = -1.*(WL-440.)/(440.-380.)
              G = 0.
              B = 1.
            ENDIF
            IF ((WL.GE.440.).AND.(WL.LE.490.)) THEN
              R = 0.
              G = (WL-440.)/(490.-440.)
              B = 1.
            ENDIF
            IF ((WL.GE.490.).AND.(WL.LE.510.)) THEN 
              R = 0.
              G = 1.
              B = -1.*(WL-510.)/(510.-490.)
            ENDIF
            IF ((WL.GE.510.).AND.(WL.LE.580.)) THEN 
              R = (WL-510.)/(580.-510.)
              G = 1.
              B = 0.
            ENDIF
            IF ((WL.GE.580.).AND.(WL.LE.645.)) THEN
              R = 1.
              G = -1.*(WL-645.)/(645.-580.)
              B = 0.
            ENDIF
            IF ((WL.GE.645.).AND.(WL.LE.780.)) THEN
              R = 1.
              G = 0.
              B = 0.
            ENDIF

         IF (WL.GT.700.) THEN
            SSS=.3+.7* (780.-WL)/(780.-700.)
         ELSE IF (WL.LT.420.) THEN
            SSS=.3+.7*(WL-380.)/(420.-380.)
         ELSE
            SSS=1.
         ENDIF

         rgb(I,1)=((SSS*R)**GAMMA)*MAX
         rgb(I,2)=((SSS*G)**GAMMA)*MAX
         rgb(I,3)=((SSS*B)**GAMMA)*MAX
        ENDDO

return
end

!---

       subroutine pixout (rgb,iunit,ihpixf,jvpixf)

       implicit none
       integer ihpixf, jvpixf
       character*1 rgb(3,ihpixf,jvpixf)      ! RGB data array
       integer iunit
       integer i, j, k
       integer itmp, icnt
       character*15 frmtstr
       character*54 headmsw
       character*4  byt4
       character*2  byt2

         itmp = mod(ihpixf, 4)
         if (itmp .NE. 0) then
           write(*,*) 'width must be multiple of 4'
           return
         endif

         headmsw( 1: 2) = 'BM'             ! declaring this is BMP file
         itmp = 54 + ihpixf * jvpixf * 3 ! total file size = header + data
         call num2bit4(itmp,byt4)
         headmsw( 3: 6) = byt4(1:4)
         itmp = 0                        ! may be 0
         call num2bit2(itmp,byt2)
         headmsw( 7: 8) = byt2(1:2)
         itmp = 0                        ! may be 0
         call num2bit2(itmp,byt2)
         headmsw( 9:10) = byt2(1:2)
         itmp = 54                       ! must be 54 : total length of header
         call num2bit4(itmp,byt4)
         headmsw(11:14) = byt4(1:4)
         itmp = 40                       ! must be 40 : length of bit-map header
         call num2bit4(itmp,byt4)
         headmsw(15:18) = byt4(1:4)
         itmp = ihpixf                   ! width
         call num2bit4(itmp,byt4)
         headmsw(19:22) = byt4(1:4)
         itmp = jvpixf                   ! height
         call num2bit4(itmp,byt4)
         headmsw(23:26) = byt4(1:4)
         itmp = 1                        ! must be 1
         call num2bit2(itmp,byt2)
         headmsw(27:28) = byt2(1:2)
         itmp = 24                       ! must be 24 : color depth in bit.
         call num2bit2(itmp,byt2)
         headmsw(29:30) = byt2(1:2)
         itmp = 0                        ! may be 0 : compression method index
         call num2bit4(itmp,byt4)
         headmsw(31:34) = byt4(1:4)
         itmp = 0                        ! may be 0 : file size if compressed
         call num2bit4(itmp,byt4)
         headmsw(35:38) = byt4(1:4)
         itmp = 0                        ! arbit. : pixel per meter, horizontal
         call num2bit4(itmp,byt4)
         headmsw(39:42) = byt4(1:4)
         itmp = 0                        ! arbit. : pixel per meter, vertical
         call num2bit4(itmp,byt4)
         headmsw(43:46) = byt4(1:4)
         itmp = 0                        ! may be 0 here : num. of color used
         call num2bit4(itmp,byt4)
         headmsw(47:50) = byt4(1:4)
         itmp = 0                        ! may be 0 here : num. of important color
         call num2bit4(itmp,byt4)
         headmsw(51:54) = byt4(1:4)
         write(iunit,'(a54,$)') headmsw(1:54)
         itmp = ihpixf * jvpixf * 3
         write(frmtstr,'(''('',i9.9,''A,$)'')') itmp
         write(iunit,fmt=frmtstr) &
           (((rgb(k,i,j),k=3,1,-1),i=1,ihpixf),j=1,jvpixf) ! writing in BGR order, not RGB.

       return
       end subroutine pixout

       subroutine num2bit4(inum,byt4)
       implicit none
       integer inum
       character*4 byt4
       integer itmp1, itmp2
       itmp1 = inum
       itmp2 = itmp1 / 256**3
       byt4(4:4) = char(itmp2)
       itmp1 =-itmp2 * 256**3 +itmp1
       itmp2 = itmp1 / 256**2
       byt4(3:3) = char(itmp2)
       itmp1 =-itmp2 * 256**2 +itmp1
       itmp2 = itmp1 / 256
       byt4(2:2) = char(itmp2)
       itmp1 =-itmp2 * 256    +itmp1
       byt4(1:1) = char(itmp1)
       return
       end subroutine num2bit4

       subroutine num2bit2(inum,byt2)
       implicit none
       integer inum
       character*2 byt2
       integer itmp1, itmp2
       itmp1 = inum
       itmp2 = itmp1 / 256
       byt2(2:2) = char(itmp2)
       itmp1 =-itmp2 * 256 + itmp1
       byt2(1:1) = char(itmp1)
       return
       end subroutine num2bit2

!---
subroutine topo_cmap (rgb,n)

implicit none

integer n,i,k,j
double precision h,rat
double precision vrgb(4,14)
integer rgb(n,3)

data vrgb / &
  0,0,97,71,50,16,122,47,&
  50,16,122,47,500,232,215,125,&
  500,232,215,125,1200,161,67,0,&
  1200,161,67,0,1700,130,30,30,&
  1700,130,30,30,2800,110,110,110,&
  2800,110,110,110,4000,255,255,255,&
  4000,255,255,255,6000,255,255,255 /

do i=1,n
h=6000.*float(i-1)/(n-1)
j=1
  do k=1,13
  if ((h-vrgb(1,k))*(h-vrgb(1,k+1)).le.0) then
  j=k
  rat=(h-vrgb(1,k))/(vrgb(1,k+1)-vrgb(1,k))
  goto 1
  endif
  enddo
1 continue
rgb(i,:)=int(vrgb(2:4,j)+rat*(vrgb(2:4,j+1)-vrgb(2:4,j)))
enddo

return
end

! -----------

subroutine write_vtk (fnme,nfnme,p,x,flag,field,title,nf,nfmax)

use definitions

implicit none

type (param) p

character fnme*(*),header*1024,footer*1024,part1*1024,part2*1024,nxc*6,nyc*6,nnc*12
integer nfnme,flag,nheader,nfooter,npart1,npart2,nf,ntitle,nfmax
double precision x(p%nn)
integer ij,i,j
real field(p%nn,nfmax)
character*1024 title(nfmax+1)

if (flag.eq.0) then
nf=0
title(1)=''
title(1)=fnme(1:nfnme)
elseif (flag.eq.1) then
nf=nf+1
if (nf.gt.nfmax) stop 'VTK called too many times'
field(:,nf)=sngl(x)
title(nf+1)=''
title(nf+1)=fnme(1:nfnme)
else
write (nxc,'(i6)') p%nx
write (nyc,'(i6)') p%ny
write (nnc,'(i12)') p%nn
header(1:1024)=''
header='# vtk DataFile Version 3.0'//char(10)//'FastScape'//char(10) &
       //'BINARY'//char(10)//'DATASET STRUCTURED_GRID'//char(10) &
       //'DIMENSIONS '//nxc//' '//nyc//' 1'//char(10)//'POINTS' &
       //nnc//' float'//char(10)
nheader=len_trim(header)
footer(1:1024)=''
footer='POINT_DATA'//nnc//char(10)
nfooter=len_trim(footer)
part1(1:1024)=''
part1='SCALARS '
npart1=len_trim(part1)+1
part2(1:1024)=''
part2=' float 1'//char(10)//'LOOKUP_TABLE default'//char(10)
npart2=len_trim(part2)
ntitle=0
  do i=2,nf+1
  ntitle=max(len_trim(title(i)),ntitle)
  enddo
  do i=2,nf+1
  title(i)(len_trim(title(i))+1:ntitle)=''
  enddo
call system ('rm -f '//trim(title(1)))
open (77,file=trim(title(1)),status='unknown',form='unformatted',access='direct', &
      recl=nheader+3*4*p%nn+nfooter+nf*(npart1+ntitle+npart2+4*p%nn),convert='big_endian')
write (77,rec=1) header(1:nheader), &
((sngl(p%dx*(i-1)/1000.),sngl(p%dy*(j-1)/1000.),sngl(x(i+(j-1)*p%nx)/1000./p%vex),i=1,p%nx),j=1,p%ny), &
footer(1:nfooter),(part1(1:npart1)//title(i+1)(1:ntitle)//part2(1:npart2),(field(ij,i),ij=1,p%nn),i=1,nf)
close(77)
endif

return
end subroutine write_vtk
