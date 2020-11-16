!program test_extract

!implicit none

!integer n
!character*256 line
!double precision v1(10),v2(10)

!print*,'enter line'
!read (*,'(a)') line
!call extract (line,':',n,v1,v2)
!call replace (line,n,(v1+v2)/2.d0)
!print*,line

!end

subroutine extract (line,symbol,n,v1,v2)

implicit none

character*(*) line
character*1 symbol
integer n,nline,np
double precision v1(*),v2(*)
integer ipos(10),i,j,i1,i2

nline=len_trim(line)

np=n

do i=1,10
ipos(i)=index(line,symbol)
if (ipos(i).eq.0) goto 999
n=n+1
  do j=ipos(i)-1,1,-1
  if (line(j:j).eq.'[') then
  read (line(j+1:ipos(i)-1),*) v1(n)
  goto 998
  endif
  enddo
stop 'error 1'
998 continue
  do j=ipos(i)+1,nline
  if (line(j:j).eq.']') then
  read (line(ipos(i)+1:j-1),*) v2(n)
  goto 997
  endif
  enddo
stop 'error 2'
997 continue
line(ipos(i):ipos(i))='#'
enddo

999 continue

do i=1,n-np
line(ipos(i):ipos(i))=symbol
enddo

return
end

!-----

subroutine replace (fnme,v)

implicit none

character*(*) fnme
character line*256,linep*256
real v
integer nfnme,i1,i2,nline,nlinep,iflag

nfnme=len_trim(fnme)
open (71,file=fnme(1:nfnme),status='old')
open (72,status='scratch')
iflag=0
111 continue
read (71,'(a)',end=999) line
nline=len_trim(line)
i1=index(line,'[')
i2=index(line,']')
  if (i1.ne.0.and.iflag.eq.0) then
  linep(1:i1-1)=line(1:i1-1)
  write (linep(i1:i1+15-1),'(g15.10)') v
  linep(i1+15:nline-(i2-i1+1)+15)=line(i2+1:nline)
  nlinep=nline-(i2-i1+1)+15
  iflag=1
  else
  linep=line
  nlinep=nline
  endif
write (72,'(a)') linep(1:nlinep)
goto 111

999 rewind (72)
rewind (71)

112 continue
read (72,'(a)',end=998) line
write (71,'(a)') trim(line)
goto 112

998 close (72)
close (71)

return
end
