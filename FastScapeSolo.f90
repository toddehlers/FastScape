subroutine FastScapeSolo (z,nx,ny,dx,dy,dt,nstep,vx,vy,vz,kf,kd,p,m,n,ibc)

! stand alone version of FastScape for hooking up to DOUAR
! user is free to use his/her own units but it is recommended that
! all distances be in m and times in yr
! in input:
! z(nx,ny) double precision array of initial topography
! nx,ny are dimension of z
! dx and dy are x- and y-spacing
! dt is time step
! nstep is number of time steps
! vx(nx,ny),vy(nx,ny),vz(nx,ny) are the velocity field defining
!         the tectonic motion of the landscape
! kf(nx,ny) is the spatially variable fluvial erodibility constant
! kd(nx,ny) is the spatially variable hillslope diffusion constant
! p(nx,ny) is the spatially variable precipitation rate
! m and n are the stream power law exponent
! ibc defines the boundary conditions as a number made of four digits that are
! either 0 or 1 (1111 or 0101 or 1110, for example). Each number corresponds to
! a side of the model, the first to the bottom (y=0), the second to the right 
! (x=xl), the third to the top (y=yl), the fourth to the left (x=0) boundary.
! A value of 1 means the corresponding boundary is at base level
! A value of 0 means the corresponding boundary is reflective
! unless both opposing bondaries are set to 0 in which case they are both cyclic
! in output:
! z(nx,ny) eroded and advected topography
! kf(nx,ny) advected erodibility constant
! kd(nx,ny) advected hillslope diffusion constant

implicit none

double precision z(nx,ny),vx(nx,ny),vy(nx,ny),vz(nx,ny)
double precision kf(nx,ny),kd(nx,ny),p(nx,ny)
double precision dt,dx,dy,m,n
integer nstep,ibc,nx,ny

! declaring arrays

double precision, dimension(:), allocatable :: h,a,length,k,kp,precip
integer, dimension(:), allocatable :: rec,ndon,stack
integer, dimension(:,:), allocatable :: donor

integer nn,nstack
integer i,j,ij,ii,jj,iii,jjj,ijk,ijr,istep,iout,ip,ipmax,i1,i2,j1,j2
double precision l,slope,smax
double precision diff,fact,h0,hp,tol
character cbc*4

! defining size of the problem

nn=nx*ny

! allocating memory

allocate (h(nn),a(nn),length(nn),rec(nn),ndon(nn),stack(nn))
allocate (donor(8,nn),k(nn),kp(nn),precip(nn))

! generating initial topography

  do j=1,ny
    do i=1,nx
    ij=i+(j-1)*nx
    precip(ij)=p(i,j)
    enddo
  enddo

! initializing other parameters

tol=1.e-3
write (cbc,'(i4)') ibc
if (ibc.lt.1000) cbc(1:1)='0'
if (ibc.lt.100) cbc(2:2)='0'
if (ibc.lt.10) cbc(3:3)='0'

! begining of time stepping

  do istep=1,nstep

! advects landscape and material properties

  call advect (z,kf,kd,vx,vy,vz,nx,ny,dt,dx,dy,cbc)

! puts topography in a linear array

  do j=1,ny
    do i=1,nx
    ij=i+(j-1)*nx
    h(ij)=z(i,j)
    enddo
  enddo

! initializing rec and length

    do ij=1,nn
    rec(ij)=ij
    length(ij)=0.
    enddo

! computing receiver array

  i1=1
  i2=nx
  if (cbc(4:4)=='1') i1=2
  if (cbc(2:2)=='1') i2=nx-1
  j1=1
  j2=ny
  if (cbc(1:1)=='1') j1=2
  if (cbc(3:3)=='1') j2=ny-1
    do j=j1,j2
      do i=i1,i2
      ij=i+(j-1)*nx
      smax=tiny(smax)
        do jj=-1,1
          do ii=-1,1
          iii=i+ii
          if (cbc(2:2)=='0'.and.cbc(4:4)=='0') iii=modulo(iii-1,nx)+1
          if (cbc(4:4)=='0') iii=max(iii,1)
          if (cbc(2:2)=='0') iii=min(iii,nx)
          jjj=j+jj
          if (cbc(1:1)=='0'.and.cbc(3:3)=='0') jjj=modulo(jjj-1,ny)+1
          if (cbc(1:1)=='0') jjj=max(jjj,1)
          if (cbc(3:3)=='0') jjj=min(jjj,ny)
          ijk=iii+(jjj-1)*nx
            if (ijk.ne.ij) then
            l=sqrt((dx*ii)**2+(dy*jj)**2)
            slope=(h(ij)-h(ijk))/l
              if (slope.gt.smax) then
              smax=slope
              rec(ij)=ijk
              length(ij)=l
              endif
            endif
          enddo
        enddo
      enddo
    enddo

! initialising number of donors per node to 0

  ndon=0

! computing donor arrays

    do ij=1,nn
      if (rec(ij).ne.ij) then
      ijk=rec(ij)
      ndon(ijk)=ndon(ijk)+1
      donor(ndon(ijk),ijk)=ij
      endif
    enddo

! computing stack

  nstack=0
    do ij=1,nn
      if (rec(ij).eq.ij) then
      nstack=nstack+1
      stack(nstack)=ij
      call find_stack (ij,donor,ndon,nn,stack,nstack)
      endif
    enddo

! computing drainage area

  a=dx*dy*precip
    do ij=nn,1,-1
    ijk=stack(ij)
      if (rec(ijk).ne.ijk) then
      a(rec(ijk))=a(rec(ijk))+a(ijk)
      endif
    enddo

! transpose rock properties onto a linear array

  do j=1,ny
    do i=1,nx
    ij=i+(j-1)*nx
    k(ij)=kf(i,j)
    kp(ij)=kd(i,j)
    enddo
  enddo

! computing erosion

    do ij=1,nn
    ijk=stack(ij)
    ijr=rec(ijk)
      if (ijr.ne.ijk) then
      fact=k(ij)*dt*a(ijk)**m/length(ijk)**n
      h0=h(ijk)
      hp=h0
      diff=tol*2.
        do while (abs(diff).gt.tol)
        h(ijk)=h(ijk)-(h(ijk)-h0+ &
        fact*(h(ijk)-h(ijr))**n)/(1.+fact*n*(h(ijk)-h(ijr))**(n-1))
        diff=h(ijk)-hp
        hp=h(ijk)
        enddo
      endif
    enddo

! computing diffusion

  call diffusion (h,nx,ny,kp,dt,dx,dy,ibc)

! transpose height into a square matrix

  do j=1,ny
    do i=1,nx
    ij=i+(j-1)*nx
    z(i,j)=h(ij)
    enddo
  enddo

  enddo

deallocate (h,a,length,rec,donor,ndon,stack,k,kp)

return

end

!----------

! recursive routine to compute the stack

recursive subroutine find_stack	(ij,donor,ndon,nn,stack,nstack)

implicit none

integer ij,k,ijk,nn,nstack
integer donor(8,nn),ndon(nn),stack(nstack)

  do k=1,ndon(ij)
  ijk=donor(k,ij)
  nstack=nstack+1
  stack(nstack)=ijk
  call find_stack (ijk,donor,ndon,nn,stack,nstack)
  enddo

return
end

!-----------------------------------

subroutine diffusion (h,nx,ny,kd,dt,dx,dy,ibc)

! to compute hillslope diffusion

! here we solve the diffusion equation using an ADI (Alternating Direction Implicit) scheme
! that is implicit and O(n)
! The solution has been checked against an analytical solution for a point source problem and a sinusoidal pertubation

implicit none

double precision h(nx*ny),kd(nx*ny)
double precision dx,dy,dt
integer nx,ny,ibc

logical xcycle,ycycle
character c4*4
double precision, dimension(:), allocatable :: f,diag,sup,inf,res,tint
integer i,j,i1,i2,j1,j2,ij,ipj,imj,ijp,ijm,ijk,iter,niter,nn
double precision fact,factxp,factxm,factyp,factym,s0,p0,dh,slope
double precision f1,s1,f2,s2,f3,s3,f4,s4

xcycle=.false.
ycycle=.false.
write (c4,'(i4)') ibc
if (c4(1:1).eq.'0' .and. c4(3:3).eq.'0') ycycle=.true.
if (c4(2:2).eq.'0' .and. c4(4:4).eq.'0') xcycle=.true.

nn=nx*ny

allocate (tint(nn))

tint=h

! half a time step in the x-direction

j1=1
if (c4(1:1).eq.'1') j1=2
j2=ny
if (c4(3:3).eq.'1') j2=ny-1
f4=1.d0
s4=0.d0
  if (c4(4:4).ne.'1') then
  f4=0.d0
  s4=-1.d0
  endif
f2=1.d0
s2=0.d0
  if (c4(2:2).ne.'1') then
  f2=0.d0
  s2=-1.d0
  endif
allocate (f(nx),diag(nx),sup(nx),inf(nx),res(nx))
do j=j1,j2
  do i=2,nx-1
  ij=(j-1)*nx+i
  ipj=(j-1)*nx+i+1
  if (i.eq.nx) ipj=(j-1)*nx+1
  imj=(j-1)*nx+i-1
  if (i.eq.1) imj=(j-1)*nx+nx
  ijp=(j)*nx+i
  if (j.eq.ny) ijp=i
  ijm=(j-2)*nx+i
  if (j.eq.1) ijm=(ny-1)*nx+i
  factxp=(kd(ipj)+kd(ij))/2.d0*(dt/2.d0)/dx**2
  factxm=(kd(imj)+kd(ij))/2.d0*(dt/2.d0)/dx**2
  factyp=(kd(ijp)+kd(ij))/2.d0*(dt/2.d0)/dy**2
  factym=(kd(ijm)+kd(ij))/2.d0*(dt/2.d0)/dy**2
  diag(i)=1.d0+factxp+factxm
  sup(i)=-factxp
  inf(i)=-factxm
  f(i)=tint(ij)+factyp*tint(ijp)-(factyp+factym)*tint(ij)+factym*tint(ijm)
  enddo
diag(1)=1.d0
sup(1)=s4
f(1)=tint((j-1)*nx+1)*f4
diag(nx)=1.d0
inf(nx)=s2
f(nx)=tint((j-1)*nx+nx)*f2
call tridiagonal (inf,diag,sup,f,res,nx)
  do i=1,nx
  ij=(j-1)*nx+i
  tint(ij)=res(i)
  enddo
enddo
deallocate (f,diag,sup,inf,res)

! half a time step in the y-direction

i1=1
if (c4(4:4).eq.'1') i1=2
i2=nx
if (c4(2:2).eq.'1') i2=nx-1
f1=1.d0
s1=0.d0
  if (c4(1:1).ne.'1') then
  f1=0.d0
  s1=-1.d0
  endif
f3=1.d0
s3=0.d0
  if (c4(3:3).ne.'1') then
  f3=0.d0
  s3=-1.d0
  endif
allocate (f(ny),diag(ny),sup(ny),inf(ny),res(ny))
do i=i1,i2
  do j=2,ny-1
  ij=(j-1)*nx+i
  ipj=(j-1)*nx+i+1
  if (i.eq.nx) ipj=(j-1)*nx+1
  imj=(j-1)*nx+i-1
  if (i.eq.1) imj=(j-1)*nx+nx
  ijp=(j)*nx+i
  if (j.eq.ny) ijp=i
  ijm=(j-2)*nx+i
  if (j.eq.1) ijm=(ny-1)*nx+i
  factxp=(kd(ipj)+kd(ij))/2.d0*(dt/2.d0)/dx**2
  factxm=(kd(imj)+kd(ij))/2.d0*(dt/2.d0)/dx**2
  factyp=(kd(ijp)+kd(ij))/2.d0*(dt/2.d0)/dy**2
  factym=(kd(ijm)+kd(ij))/2.d0*(dt/2.d0)/dy**2
  diag(j)=1.+factyp+factym
  sup(j)=-factyp
  inf(j)=-factym
  f(j)=tint(ij)+factxp*tint(ipj)-(factxp+factxm)*tint(ij)+factxm*tint(imj)
  enddo
diag(1)=1.d0
sup(1)=s1
f(1)=tint(i)*f1
diag(ny)=1.d0
inf(ny)=s3
f(ny)=tint((ny-1)*nx+i)*f3
call tridiagonal (inf,diag,sup,f,res,ny)
  do j=1,ny
  ij=(j-1)*nx+i
  tint(ij)=res(j)
  enddo
enddo
deallocate (f,diag,sup,inf,res)

h=tint

deallocate (tint)

return

end subroutine diffusion

!---------------------------------

      SUBROUTINE tridiagonal(a,b,c,r,u,n)
! subroutine to solve a tri-diagonal system of linear equations
! from numerical recipes
      INTEGER n
      double precision a(n),b(n),c(n),r(n),u(n)
      INTEGER j
      double precision bet
      double precision,dimension(:),allocatable::gam
      allocate (gam(n))
      if(b(1).eq.0.d0) stop 'in tridiagonal'
      bet=b(1)
      u(1)=r(1)/bet
      do 11 j=2,n
        gam(j)=c(j-1)/bet
        bet=b(j)-a(j)*gam(j)
        if(bet.eq.0.) then
        print*,'tridiagonal failed'
        stop
        endif
        u(j)=(r(j)-a(j)*u(j-1))/bet
11    continue
      do 12 j=n-1,1,-1
        u(j)=u(j)-gam(j+1)*u(j+1)
12    continue
      deallocate (gam)
      return
      END

!-----------------------------------------

subroutine advect (z,kf,kd,vx,vy,vz,nx,ny,dt,dx,dy,cbc)

! routine to advect the topography (z), the two erosional constant (kf,kd)
! by a velocity field vx,vy,vz known at the same locations on a rectangular
! grid of nx by ny points separated by dx,dy over a time step dt
! cbc is a character holding the boundary conditions

implicit none

double precision z(nx,ny),kf(nx,ny),kd(nx,ny),vx(nx,ny),vy(nx,ny),vz(nx,ny)
double precision dt,dx,dy
integer nx,ny,nn
character cbc*4

double precision, dimension(:), allocatable :: diag,sup,inf,rhs,res

integer i,j,ij,i1,i2,j1,j2

nn=nx*ny

! x-advection using an implicit, second-order scheme to solve advection eq.

allocate (diag(nx),sup(nx),inf(nx),rhs(nx),res(nx))

  do j=1,ny

  diag=1.d0
  sup=vx(:,j)*dt/dx/2.d0
  inf=-vx(:,j)*dt/dx/2.d0

  diag(1)=1.d0
  sup(1)=0.d0
    if (cbc(4:4)=='0') then
    sup(1)=-1.d0
    endif
  diag(nx)=1.d0
  inf(nx)=0.d0
    if (cbc(2:2)=='0') then
    inf(nx)=-1.d0
    endif

  rhs=z(:,j)
    if (cbc(4:4)=='0') then
    rhs(1)=0.d0
    endif
    if (cbc(2:2)=='0') then
    rhs(nx)=0.d0
    endif
  call tridiagonal (inf,diag,sup,rhs,res,nx)
  z(:,j)=res

  rhs=kf(:,j)
    if (cbc(4:4)=='0') then
    rhs(1)=0.d0
    endif
    if (cbc(2:2)=='0') then
    rhs(nx)=0.d0
    endif
  call tridiagonal (inf,diag,sup,rhs,res,nx)
  kf(:,j)=res

  rhs=kd(:,j)
    if (cbc(4:4)=='0') then
    rhs(1)=0.d0
    endif
    if (cbc(2:2)=='0') then
    rhs(nx)=0.d0
    endif
  call tridiagonal (inf,diag,sup,rhs,res,nx)
  kd(:,j)=res

  enddo

deallocate (diag,sup,inf,rhs,res)

! y-advection using an implicit, second-order scheme to solve advection eq.

allocate (diag(ny),sup(ny),inf(ny),rhs(ny),res(ny))

  do i=1,nx
  
  diag=1.d0
  sup=vy(i,:)*dt/dy/2.d0
  inf=-vy(i,:)*dt/dy/2.d0

  diag(1)=1.d0
  sup(1)=0.d0
    if (cbc(1:1)=='0') then
    sup(1)=-1.d0
    endif
  diag(ny)=1.d0
  inf(ny)=0.d0
    if (cbc(3:3)=='0') then
    inf(ny)=-1.d0
    endif

  rhs=z(i,:)
    if (cbc(1:1)=='0') then
    rhs(1)=0.d0
    endif
    if (cbc(3:3)=='0') then
    rhs(ny)=0.d0
    endif
  call tridiagonal (inf,diag,sup,rhs,res,ny)
  z(i,:)=res

  rhs=kf(i,:)
    if (cbc(1:1)=='0') then
    rhs(1)=0.d0
    endif
    if (cbc(3:3)=='0') then
    rhs(ny)=0.d0
    endif
  call tridiagonal (inf,diag,sup,rhs,res,ny)
  kf(i,:)=res

  rhs=kd(i,:)
    if (cbc(1:1)=='0') then
    rhs(1)=0.d0
    endif
    if (cbc(3:3)=='0') then
    rhs(ny)=0.d0
    endif
  call tridiagonal (inf,diag,sup,rhs,res,ny)
  kd(i,:)=res
  
  enddo

deallocate (diag,sup,inf,rhs,res)

! z-component by simply moving material upwards

i1=1
i2=nx
if (cbc(4:4)=='1') i1=2
if (cbc(2:2)=='1') i2=nx-1
j1=1
j2=ny
if (cbc(1:1)=='1') j1=2
if (cbc(3:3)=='1') j2=ny-1

  do j=j1,j2
    do i=i1,i2
    ij=i+(j-1)*nx
    z(i,j)=z(i,j)+vz(i,j)*dt
    enddo
  enddo

return
end
