!---
      logical function inside (x0,y0,x,y,n)
 
! subroutine to check whether a point (x0,y0) is inside a polygon
! defined by n vertices x(n),y(n)

! note that the polygone does not have to be convex but the result
! is umpredictable if x0,y0 is on the polygon

      implicit none

      integer n,i,ip
      double precision  x0,y0,angle,x1,x2,y1,y2,beta,cbeta,sbeta,xp2,yp2
      double precision x(n),y(n)
      double precision atanm2
      external atanm2

      angle=0.
      inside=.true.

        do i=1,n
        ip=mod(i,n)+1
        x1=x(i)
        y1=y(i)
        x2=x(ip)
        y2=y(ip)
        beta=atanm2(y1-y0,x1-x0)
        cbeta=cos(beta)
        sbeta=sin(beta)
        xp2=(x2-x0)*cbeta+(y2-y0)*sbeta
        yp2=-(x2-x0)*sbeta+(y2-y0)*cbeta
        angle=angle+atanm2(yp2,xp2)
        enddo

      if (abs(angle).lt.1.) inside=.false.

      return
      end

!---
      double precision function atanm2 (x,y)

      implicit none

      double precision x,y

      atanm2=0.
      if (x.ne.0. .or. y.ne.0.) atanm2=atan2(x,y)

      return
      end
