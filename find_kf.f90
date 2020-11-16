double precision function find_kf (ij,l,t,p)

use omp_lib
use definitions

implicit none

type (topography) t
type (topology) g
type (param) p
type (law) l
integer ij

double precision kf,kfp,x,y,dist
integer is,k

kfp=l%kf
kf=kfp    
  do is=1,l%strati_n
  if (t%href(ij)-t%h(ij).ge.l%strati_top(is).and.t%href(ij)-t%h(ij).le.l%strati_bottom(is)) &
      kf=kfp*l%strati_f(is)
  enddo     
x=p%dx*(1+mod(ij-1,p%nx))
y=p%dy*(1+(ij-1)/p%nx)
  do k=1,p%granite_n
  dist=sqrt((x-p%granite_x(k))**2/p%granite_rx(k)**2+(y-p%granite_y(k))**2/p%granite_ry(k)**2)
  if (dist.lt.1.d0.and. &
      t%href(ij)-t%h(ij).ge.p%granite_top(k).and.t%href(ij)-t%h(ij).le.p%granite_bottom(k)) kf=kf*p%granite_dk(k)
  enddo     

find_kf=kf

return

end function find_kf
