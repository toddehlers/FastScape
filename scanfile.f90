!------------------------------------------------------------------------------|

subroutine iscanfile (fnme,text,res,ires,vocal)

!------------------------------------------------------------------------------|
! subroutine to read the value of an integer parameter whose name is stored in 
! text  from file fnme the result is stored in res and the flag ires is set 
! to 1 if all went well. if vocal is set yo 1 the routine echoes its find to
! standard output

character*(*) fnme,text
integer,intent(out)::res,ires
integer,intent(in)::vocal

character line*256,record*256
integer ios

ires=0

open (7,file=fnme(1:len_trim(fnme)),status='old',iostat=ios) 
if (ios/=0) write(*,*) 'pb opening input.txt in iscanfile'

111 continue

read (7,'(a)',end=112) line

ieq=scan(line,'=')
if (ieq/=0) then
   if (trim(adjustl(line(1:ieq-1))).eq.text(1:len_trim(text)))  then
      nrec=len_trim(adjustl(line(ieq+1:len_trim(line))))
      record(1:nrec)=trim(adjustl(line(ieq+1:len_trim(line))))
        if (record(1:1).eq.text(1:1)) then
        res=0
        read (record(nrec-1:nrec),'(i2)',err=888) ires
        ires=-ires
        goto 889
888     read (record(nrec:nrec),'(i1)',err=999) ires
        ires=-ires
889     continue
        else
        read (record(1:nrec),*,err=999) res
        ires=1
        endif
      if (vocal.eq.1) write (*,*) 'Value read for '//text//': '//record(1:nrec)
   end if
end if

goto 111

112 close (7)

return

999 print*,''
    print*,'Error in input file'
    print*,'You have entered the following value : (',record(1:nrec),')'
    print*,'while parameter ('//text//') needs to be specified as an integer'
    stop ' Terminating job'
end subroutine iscanfile



!---------------------------------------------------------------------------

subroutine dscanfile (fnme,text,res,ires,vocal)

!------------------------------------------------------------------------------|
! subroutine to read the value of a double precision  parameter whose name 
! is stored in  text  from file fnme the result is stored in res and the flag 
! ires is set  to 1 if all went well. if vocal is set yo 1 the routine echoes its find to
! standard output

character*(*) fnme,text
double precision res
integer ires,vocal

character line*256,record*256
integer ios

ires=0

open (7,file=fnme(1:len_trim(fnme)),status='old',iostat=ios) 
if (ios/=0) write(*,*) 'pb opening input.txt in iscanfile'

111     continue

read (7,'(a)',end=112) line

ieq=scan(line,'=')
if (ieq/=0) then
   if (trim(adjustl(line(1:ieq-1))).eq.text(1:len_trim(text))) then
      nrec=len_trim(adjustl(line(ieq+1:len_trim(line))))
      record(1:nrec)=trim(adjustl(line(ieq+1:len_trim(line))))
        if (record(1:1).eq.text(1:1)) then
        res=0.d0
        read (record(nrec-1:nrec),'(i2)',err=888) ires
        ires=-ires
        goto 889
888     read (record(nrec:nrec),'(i1)',err=999) ires
        ires=-ires
889     continue
        else
        read (record(1:nrec),*,err=999) res
        ires=1
        endif
      if (vocal.eq.1) write (*,*) 'Value read for '//text//': '//record(1:nrec)
   endif
end if

goto 111

112   close (7)

return

999 print*,''
    print*,'Error in input file'
    print*,'You have entered the following value : (',record(1:nrec),')'
    print*,'while parameter ('//text//') needs to be specified as a floating point number'
    stop ' Terminating job'
end subroutine dscanfile


!------------------------------------------------------------------------------|

subroutine cscanfile (fnme,text,res,ires,vocal)

!------------------------------------------------------------------------------|
! subroutine to read the value of a character string parameter whose name 
! is stored in  text  from file fnme the result is stored in res and the flag 
! ires is set  to 1 if all went well. if vocal is set yo 1 the routine echoes its find to
! standard output
!------------------------------------------------------------------------------|
      
character*(*) fnme,text,res
integer ires,vocal

character line*256,record*256
integer ios

ires=0

open (7,file=fnme(1:len_trim(fnme)),status='old',iostat=ios) 
if (ios/=0) write(*,*) 'pb opening input.txt in iscanfile'

111     continue

read (7,'(a)',end=112) line

ieq=scan(line,'=')
if (ieq/=0) then
   if (trim(adjustl(line(1:ieq-1))).eq.text(1:len_trim(text))) then
      nrec=len_trim(adjustl(line(ieq+1:len_trim(line))))
      record(1:nrec)=trim(adjustl(line(ieq+1:len_trim(line))))
      res=record(1:nrec)
      if (vocal.eq.1) write (*,*) 'Value read for '//text//': '//record(1:nrec)
      ires=1
   endif
end if

goto 111

112   close (7)

return
end subroutine cscanfile


!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
