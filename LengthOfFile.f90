! Subroutine that determines the length of a file
subroutine LengthOfFile(FileName, LineCounter)
    implicit none

    character(255), intent(in) :: FileName
	integer(4) :: FileUnit = 11, ios
	integer(4), intent(out) :: LineCounter

	open(FileUnit, file=FileName, status='old')

		LineCounter = 0

		do
			read(FileUnit, *, iostat=ios)
			if (ios/=0) then
				exit
			end if
			LineCounter = LineCounter + 1
		end do
	close(FileUnit)
end subroutine LengthOfFile