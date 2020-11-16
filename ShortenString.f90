subroutine ShortenString(FileName, ReduceBy)
	implicit none

	character(255), intent(out) :: FileName
	integer(4), intent(in) :: ReduceBy
	integer(4) :: FileName_StringLength

	FileName_StringLength = len_trim(FileName)

	if (FileName_StringLength < ReduceBy) then
		print*, 'Subroutine ShortenString: String is too short to be shortened by ', ReduceBy, ' characters.'
	else
		FileName = FileName(1:(FileName_StringLength - ReduceBy))
	end if
end subroutine ShortenString