!------

subroutine timer (tin,flag,nb_cycles_initial,message)

implicit none

real tin,temps_elapsed
integer flag
integer nb_cycles_sec,nb_cycles_max,nb_cycles_initial,nb_cycles_final,nb_cycles
character message*(*)

call system_clock (count_rate=nb_cycles_sec,count_max=nb_cycles_max)
call system_clock (count=nb_cycles_final)
nb_cycles=nb_cycles_final - nb_cycles_initial
IF (nb_cycles_final.lt.nb_cycles_initial) &
        nb_cycles=nb_cycles+nb_cycles_max
temps_elapsed=real(nb_cycles)/nb_cycles_sec
  if (flag.eq.0) then
  tin=temps_elapsed
  if (message.ne.'$') print*,message
  else
  tin=temps_elapsed-tin
  if (message.ne.'$') print*,message,tin
  endif


return
end

