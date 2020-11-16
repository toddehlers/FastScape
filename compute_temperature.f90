      subroutine compute_temperature (timeH,erosionRate,temperatureH,nstep,xl,tmin,tmax,heat,xheat,cond,xcond, &
                                      SurfaceGradientToday)

! subroutine to compute the temperature history of a rock particle exhumed according to
! an erosion history stored in erosionRate(nstep)/timeH(nstep). The result is stored in
! temperatureH(nstep). The routine solves the 1D heat conduction, advection and production
! equation. It assumes that the temperature is fixed at tmin at the surface and tmax at a depth
! xl. Heat production is heat applied over the top xheat layer. Note that the heat producing
! layer is affected by erosion too. Heat diffusivity is assumed to be 25 km^2/Myr but one
! can specify the thickness of a top layer xcond that is cond times more diffusive. In output
! the routine also provides the present-day (final) surface temperature gradient 
! (SurfaceGradientToday)

! time is in years
! erosionRate in m/yr
! temperature in deg C
! xl in m
! tmin, tmax in deg C
! heat in deg C/yr
! xheat in m
! cond is a multiplication factor
! xcond in m
! SurfaceGradietnToday is in deg C/m

      implicit none

      integer nnode,nstep,istep,i
      double precision, dimension(:), allocatable :: x,t,error
      double precision, dimension(:,:), allocatable :: tint
      double precision timeH(nstep),temperatureH(nstep),erosionRate(nstep)
      double precision tmin,tmax,depth,temp,r,f,h,z,Pe,diffusivity,heat,xheat,xheatp,SurfaceGradientToday,dt,xl
      double precision cond,xcond,xcondp

      nnode=101

      allocate (x(nnode),t(nnode),tint(nnode,nstep),error(nstep))

      diffusivity=25.d0
      Pe=-erosionRate(1)*xl/diffusivity
      f=heat*xl*xl/diffusivity/tmax

        do i=1,nnode
        x(i)=xl*float(i-1)/(nnode-1)
        enddo

      Pe=-erosionRate(1)*xl/diffusivity

      t=(tmax-tmin)*x/xl+tmin

      call Heat1D (nnode,-1.d9,x,t,erosionRate(1),diffusivity,heat,xheat,cond,xcond)

      tint(:,1)=t

      xheatp=xheat
      xcond=xcond
        do istep=2,nstep
        dt=timeH(istep-1)-timeH(istep)
        xheatp=max(0.d0,xheatp+erosionRate(istep)*dt)
        xcondp=max(0.d0,xcondp+erosionRate(istep)*dt)
        call Heat1D (nnode,dt,x,t,erosionRate(istep),diffusivity,heat,xheatp,cond,xcondp)
        tint(:,istep)=t
        enddo

      depth=0.d0
      temperatureH(nstep)=tmin
        do istep=nstep-1,1,-1
        dt=timeH(istep)-timeH(istep+1)
        depth=depth-erosionRate(istep+1)*dt
        i=1+int(((nnode-1)*depth)/xl)
        temp=tmax
          if (i.lt.nnode) then
          r=(depth-x(i))/(x(i+1)-x(i))
          temp=tint(i,istep)+r*(tint(i+1,istep)-tint(i,istep))
          endif
        temperatureH(istep)=temp
        enddo

      SurfaceGradientToday=(tint(2,nstep)-tint(1,nstep))/(x(2)-x(1))*1.d3
      deallocate (x,t,tint,error)

      return
      end subroutine compute_temperature

!---------------------------------------------

      subroutine Heat1D (nnode,tfinal,x,cinit,er,diffusivity,heat,xheat,cond,xcond)

! Finite element program to solve:
!  rc dc/dt = -v dc/dz + d/dz(cond dc/dz) + heat


      implicit none

      double precision, dimension(:), allocatable :: c,f,tstep,xint,a4p
      double precision, dimension(:,:), allocatable :: a
      integer, dimension(:), allocatable :: ipvt

      double precision cinit(nnode),x(nnode)
      double precision bel(3),ael1(3,3),ael2(3,3)
      double precision r(2),h(2,3),b(2,3)
      double precision tfinal,er,heat,xlength,hh,diffusivity,xheat,cond,xcond
      double precision bc1,bc2,bc1p,bc2p,alpha,timen,timep,det,dt,a1,a2,a3,a4
      double precision acond,amass,dynamic
      integer nstep,ibc1,ibc2,i,j,nelem,nnode,iint,ielem,jint
      integer istep,i1,i2,i3,ic,info

      data r/-.57735,.57735/

      xlength=x(nnode)-x(1)

      if (tfinal.lt.0.d0) then
      nstep=1
      dynamic=0.d0
      else
      nstep=max(int(er*tfinal/xlength*10.d0),1)
      dynamic=1.d0
      endif

      allocate (c(nnode),f(nnode),a(7,nnode))
      allocate (tstep(nstep),xint(nnode))
      allocate (a4p(nnode))
      allocate (ipvt(nnode))

      ibc1=1
      ibc2=1

      bc1=cinit(1)
      bc2=cinit(nnode)

      alpha=1.d0
      hh=0.d0

        do i=1,nnode
        c(i)=cinit(i)
        enddo

      if (tfinal.gt.0.d0) then
        do i=1,nstep
        tstep(i)=float(i)/float(nstep)*abs(tfinal)
        enddo
      else
        tstep(1)=1.d0
      endif

      nelem=(nnode-1)/2

        do iint=1,2
        h(iint,1)=-r(iint)/2.d0*(1.d0-r(iint))
        h(iint,2)=1.d0-r(iint)*r(iint)
        h(iint,3)=r(iint)/2.d0*(1.d0+r(iint))
        b(iint,1)=(-.5d0+r(iint))
        b(iint,2)=(-2.d0*r(iint))
        b(iint,3)=(.5d0+r(iint))
        enddo

      timep=0.d0
        do ielem=1,nelem 
        i1=(ielem-1)*2+1
          do iint=1,2
          jint=(ielem-1)*2+iint
          xint(jint)=x(i1+1)+r(iint)*(x(i1+1)-x(i1))
          a4p(jint)=0.d0
          if (xint(jint).lt.xheat) a4p(jint)=heat
          enddo
        enddo

! calculates the position through time of a particle that will end up
! at the surface of the Earth (x=0) at the end of the numerical experiment
! (time=tfinal)

        do istep=1,nstep

        timen=tstep(istep)
        dt=timen-timep

          do i=1,nnode
          f(i)=0.d0
            do j=1,7
            a(j,i)=0.d0
            enddo
          enddo

          do ielem=1,nelem

            do i=1,3
            bel(i)=0.d0
               do j=1,3
               ael1(i,j)=0.d0
               ael2(i,j)=0.d0
               enddo
            enddo

          i1=(ielem-1)*2+1
          i3=i1+2
          det=(x(i3)-x(i1))/2.d0

            do iint=1,2                          
            jint=(ielem-1)*2+iint
                           
            a1=1.d0
            a2=er*a1
            a3=diffusivity
            if (xint(jint).lt.xcond) a3=diffusivity*cond
            a4=0.d0
            if (xint(jint).lt.xheat) a4=heat

               do i=1,3
               bel(i)=bel(i)+h(iint,i)*dt*det* &
                            ((1.d0-alpha)*a4p(jint)+alpha*a4)
                  do j=1,3
                  acond=b(iint,i)*a3*b(iint,j)/det/det &
                      +h(iint,i)*a2*b(iint,j)/det
                  amass=h(iint,i)*a1*h(iint,j)*dynamic
                  ael1(i,j)=ael1(i,j) &
                          +(amass+dt*alpha*acond)*det
                  ael2(i,j)=ael2(i,j) &
                          +(amass-dt*(1.d0-alpha)*acond)*det
                  enddo
               enddo            

            a4p(jint)=a4

            enddo

            do i=1,3
               do j=1,3
               ic=(ielem-1)*2+j
               bel(i)=bel(i)+ael2(i,j)*c(ic)
               enddo
            enddo
!
          if (ielem.eq.1) then

            if (ibc1.eq.1) then
! fixed temperature at top = bc1 (node 1)
            ael1(1,1)=1.d0
            bel(1)=bc1
              do j=2,3
              bel(j)=bel(j)-ael1(j,1)*bc1
              ael1(1,j)=0.d0
              ael1(j,1)=0.d0
              enddo
            
            elseif (ibc1.eq.2) then
! fixed flux at top = bc1 (node 1)
            bel(1)=bel(1)+dt*((1.d0-alpha)*bc1p+alpha*bc1)

            elseif (ibc1.eq.3) then
! radiative b.c. at top (node 1)
            ael1(1,1)=ael1(1,1)+dt*hh
            bel(1)=bel(1)+dt*hh*((1.d0-alpha)*bc1p+alpha*bc1)
            endif

          endif

          if (ielem.eq.nelem) then

            if (ibc2.eq.1) then
! fixed temperature at bottom = bc2 (node nnode)
            ael1(3,3)=1.d0
            bel(3)=bc2
              do j=1,2
              bel(j)=bel(j)-ael1(j,3)*bc2
              ael1(3,j)=0.d0
              ael1(j,3)=0.d0
              enddo
            
            elseif (ibc2.eq.2) then
! fixed flux at bottom = bc2 (node nnode)
            bel(3)=bel(3)+dt*((1.d0-alpha)*bc2p+alpha*bc2)

            elseif (ibc2.eq.3) then
            ael1(3,3)=ael1(3,3)+dt*hh
            bel(3)=bel(3)+dt*hh*((1.d0-alpha)*bc2p+alpha*bc2)
            endif

          endif

            do i=1,3
            i1=(ielem-1)*2+i
            f(i1)=f(i1)+bel(i)
               do j=1,3
               i2=(ielem-1)*2+j
               a(5+(i2-i1),i1)=a(5+(i2-i1),i1)+ael1(j,i)
               enddo
            enddo

         enddo

        call sgbfa (a,7,nnode,2,2,ipvt,info)

          if (info.ne.0) then
          print*,'Negative pivot in position: ',info
          stop 'end of processing...'
          endif

        call sgbsl (a,7,nnode,2,2,ipvt,f,0)

          do i=1,nnode
          c(i)=f(i)
          enddo

        timep=timen
        bc1p=bc1
        bc2p=bc2

        enddo

      cinit=c

      deallocate (c,f,a,tstep,xint,a4p,ipvt)

      return

      end subroutine Heat1D

!-----------------------------------------------------

      SUBROUTINE SGBFA (ABD,LDA,N,ML,MU,IPVT,INFO)

      IMPLICIT NONE

      INTEGER LDA,N,ML,MU,IPVT(*),INFO
      DOUBLE PRECISION ABD(LDA,*)

      DOUBLE PRECISION T
      INTEGER I,ISAMAX,I0,J,JU,JZ,J0,J1,K,KP1,L,LM,M,MM,NM1

      M = ML + MU + 1
      INFO = 0

      J0 = MU + 2
      J1 = MIN0(N,M) - 1
      IF (J1 .LT. J0) GO TO 30
      DO 20 JZ = J0, J1
         I0 = M + 1 - JZ
         DO 10 I = I0, ML
            ABD(I,JZ) = 0.0E0
   10    CONTINUE
   20 CONTINUE
   30 CONTINUE
      JZ = J1
      JU = 0

      NM1 = N - 1
      IF (NM1 .LT. 1) GO TO 130
      DO 120 K = 1, NM1
         KP1 = K + 1

         JZ = JZ + 1
         IF (JZ .GT. N) GO TO 50
         IF (ML .LT. 1) GO TO 50
            DO 40 I = 1, ML
               ABD(I,JZ) = 0.0E0
   40       CONTINUE
   50    CONTINUE

         LM = MIN0(ML,N-K)
         L = ISAMAX(LM+1,ABD(M,K),1) + M - 1
         IPVT(K) = L + K - M

         IF (ABD(L,K) .EQ. 0.) GO TO 100

            IF (L .EQ. M) GO TO 60
               T = ABD(L,K)
               ABD(L,K) = ABD(M,K)
               ABD(M,K) = T
   60       CONTINUE

            T = -1.0E0/ABD(M,K)
            CALL SSCAL(LM,T,ABD(M+1,K),1)

            JU = MIN0(MAX0(JU,MU+IPVT(K)),N)
            MM = M
            IF (JU .LT. KP1) GO TO 90
            DO 80 J = KP1, JU
               L = L - 1
               MM = MM - 1
               T = ABD(L,J)
               IF (L .EQ. MM) GO TO 70
                  ABD(L,J) = ABD(MM,J)
                  ABD(MM,J) = T
   70          CONTINUE
               CALL SAXPY(LM,T,ABD(M+1,K),1,ABD(MM+1,J),1)
   80       CONTINUE
   90       CONTINUE
         GO TO 110
  100    CONTINUE
            INFO = K
  110    CONTINUE
  120 CONTINUE
  130 CONTINUE
      IPVT(N) = N
      IF (ABD(M,N) .EQ. 0.) INFO = N
      RETURN
      END

!-----------------------------------------------------------------------
      SUBROUTINE SGBSL (ABD,LDA,N,ML,MU,IPVT,B,JOB)

      IMPLICIT NONE

      INTEGER LDA,N,ML,MU,IPVT(*),JOB
      DOUBLE PRECISION ABD(LDA,*),B(*)
 
      DOUBLE PRECISION SDOT,T
      INTEGER K,KB,L,LA,LB,LM,M,NM1

      M = MU + ML + 1
      NM1 = N - 1
      IF (JOB .NE. 0) GO TO 50

         IF (ML .EQ. 0) GO TO 30
         IF (NM1 .LT. 1) GO TO 30
            DO 20 K = 1, NM1
               LM = MIN0(ML,N-K)
               L = IPVT(K)
               T = B(L)
               IF (L .EQ. K) GO TO 10
                  B(L) = B(K)
                  B(K) = T
   10          CONTINUE
               CALL SAXPY(LM,T,ABD(M+1,K),1,B(K+1),1)
   20       CONTINUE
   30    CONTINUE

         DO 40 KB = 1, N
            K = N + 1 - KB
            B(K) = B(K)/ABD(M,K)
            LM = MIN0(K,M) - 1
            LA = M - LM
            LB = K - LM
            T = -B(K)
            CALL SAXPY(LM,T,ABD(LA,K),1,B(LB),1)
   40    CONTINUE
      GO TO 100
   50 CONTINUE

         DO 60 K = 1, N
            LM = MIN0(K,M) - 1
            LA = M - LM
            LB = K - LM
            T = SDOT(LM,ABD(LA,K),1,B(LB),1)
            B(K) = (B(K) - T)/ABD(M,K)
   60    CONTINUE

         IF (ML .EQ. 0) GO TO 90
         IF (NM1 .LT. 1) GO TO 90
            DO 80 KB = 1, NM1
               K = N - KB
               LM = MIN0(ML,N-K)
               B(K) = B(K) + SDOT(LM,ABD(M+1,K),1,B(K+1),1)
               L = IPVT(K)
               IF (L .EQ. K) GO TO 70
                  T = B(L)
                  B(L) = B(K)
                  B(K) = T
   70          CONTINUE
   80       CONTINUE
   90    CONTINUE
  100 CONTINUE
      RETURN
      END

!-----------------------------------------------------------------------

      DOUBLE PRECISION FUNCTION SDOT(N,SX,INCX,SY,INCY)

      IMPLICIT NONE

      DOUBLE PRECISION SX(*),SY(*),STEMP
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N

      STEMP = 0.0E0
      SDOT = 0.0E0
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20

      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        STEMP = STEMP + SX(IX)*SY(IY)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      SDOT = STEMP
      RETURN

   20 M = MOD(N,5)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        STEMP = STEMP + SX(I)*SY(I)
   30 CONTINUE
      IF( N .LT. 5 ) GO TO 60
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
        STEMP = STEMP + SX(I)*SY(I) + SX(I + 1)*SY(I + 1) + &
         SX(I + 2)*SY(I + 2) + SX(I + 3)*SY(I + 3) + SX(I + 4)*SY(I + 4)
   50 CONTINUE
   60 SDOT = STEMP
      RETURN
      END

!-------------------------------------------------------------------------

      SUBROUTINE  SSCAL(N,SA,SX,INCX)

      IMPLICIT NONE

      DOUBLE PRECISION SA,SX(*)
      INTEGER I,INCX,M,MP1,N,NINCX

      IF(N.LE.0)RETURN
      IF(INCX.EQ.1)GO TO 20

      NINCX = N*INCX
      DO 10 I = 1,NINCX,INCX
        SX(I) = SA*SX(I)
   10 CONTINUE
      RETURN

   20 M = MOD(N,5)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        SX(I) = SA*SX(I)
   30 CONTINUE
      IF( N .LT. 5 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
        SX(I) = SA*SX(I)
        SX(I + 1) = SA*SX(I + 1)
        SX(I + 2) = SA*SX(I + 2)
        SX(I + 3) = SA*SX(I + 3)
        SX(I + 4) = SA*SX(I + 4)
   50 CONTINUE
      RETURN
      END

!-------------------------------------------------------------------------

      INTEGER FUNCTION ISAMAX(N,SX,INCX)

      IMPLICIT NONE

      DOUBLE PRECISION SX(*),SMAX
      INTEGER I,INCX,IX,N

      ISAMAX = 0
      IF( N .LT. 1 ) RETURN
      ISAMAX = 1
      IF(N.EQ.1)RETURN
      IF(INCX.EQ.1)GO TO 20

      IX = 1
      SMAX = ABS(SX(1))
      IX = IX + INCX
      DO 10 I = 2,N
         IF(ABS(SX(IX)).LE.SMAX) GO TO 5
         ISAMAX = I
         SMAX = ABS(SX(IX))
    5    IX = IX + INCX
   10 CONTINUE
      RETURN

   20 SMAX = ABS(SX(1))
      DO 30 I = 2,N
         IF(ABS(SX(I)).LE.SMAX) GO TO 30
         ISAMAX = I
         SMAX = ABS(SX(I))
   30 CONTINUE
      RETURN
      END

!-----------------------------------------------------------------------

      SUBROUTINE SAXPY(N,SA,SX,INCX,SY,INCY)

      IMPLICIT NONE

      DOUBLE PRECISION SX(*),SY(*),SA
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N

      IF(N.LE.0)RETURN
      IF (SA .EQ. 0.) RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20

      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        SY(IY) = SY(IY) + SA*SX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN

   20 M = MOD(N,4)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        SY(I) = SY(I) + SA*SX(I)
   30 CONTINUE
      IF( N .LT. 4 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,4
        SY(I) = SY(I) + SA*SX(I)
        SY(I + 1) = SY(I + 1) + SA*SX(I + 1)
        SY(I + 2) = SY(I + 2) + SA*SX(I + 2)
        SY(I + 3) = SY(I + 3) + SA*SX(I + 3)
   50 CONTINUE
      RETURN
      END
