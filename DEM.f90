subroutine DEM (fnme,nfnme,p,h,flag)

use definitions

implicit none

type (param) p
double precision h(p%nx,p%ny)
integer nfnme,flag
character fnme*(*)

open (7,file=fnme(1:nfnme)//'.hdr',status='unknown')
write (7,'(A14,a1)') 'BYTEORDER     ','I'
write (7,'(A14,a3)') 'LAYOUT        ','BIL'
write (7,'(A14,i5)') 'NROWS         ',p%ny
write (7,'(A14,i5)') 'NCOLS         ',p%nx
write (7,'(A14,i1)') 'NBANDS        ',1
if (flag.eq.1) then
write (7,'(A14,i2)') 'NBITS         ',16
elseif (flag.eq.2) then
write (7,'(A14,i2)') 'NBITS         ',32
elseif (flag.eq.3) then
write (7,'(A14,i2)') 'NBITS         ',32
elseif (flag.eq.4) then
write (7,'(A14,i2)') 'NBITS         ',64
else
stop 'value not acceptable for plot_DEM'
endif
write (7,'(A14,i5)') 'BANDROWBYTES  ',p%nx*2
write (7,'(A14,i5)') 'TOTALROWBYTES ',p%nx*2
write (7,'(A14,i1)') 'BANDGAPBYTES  ',0
write (7,'(A14,i5)') 'NODATA        ',-9999
write (7,'(A14,f3.0)') 'ULXMAP        ',0.
write (7,'(A14,f3.0)') 'ULYMAP        ',0.
write (7,'(A14,f16.13)') 'XDIM          ',p%dx/111111.111d0
write (7,'(A14,f16.13)') 'YDIM          ',p%dy/111111.111d0
close (7)

if (flag.eq.1) then
open (7,file=fnme(1:nfnme)//'.dem',status='unknown',form='unformatted', &
      access='direct',recl=p%nx*p%ny*2)
write (7,rec=1) int2(h)
close (7)

elseif (flag.eq.2) then
open (7,file=fnme(1:nfnme)//'.dem',status='unknown',form='unformatted', &
      access='direct',recl=p%nx*p%ny*4)
write (7,rec=1) int(h)
close (7)

elseif (flag.eq.3) then
open (7,file=fnme(1:nfnme)//'.dem',status='unknown',form='unformatted', &
      access='direct',recl=p%nx*p%ny*4)
write (7,rec=1) sngl(h)
close (7)

elseif (flag.eq.4) then
open (7,file=fnme(1:nfnme)//'.dem',status='unknown',form='unformatted', &
      access='direct',recl=p%nx*p%ny*8)
write (7,rec=1) h
close (7)

else
stop 'value not acceptable for plot_DEM'
endif

return

end subroutine DEM
