subroutine count_xanes(c,x,p,itrain,conv,idx)
! How many spectra with corresponding xyz files are there?

  implicit none
  integer :: i,c,x,p
  integer, dimension(x) :: idx
  logical :: xyz_file_exists,xanes_file_exists
  character(len=28) :: convfile
  character(len=13) :: xyzfile
  character(len=2) :: itrain 
  logical :: conv

  104 format(a5,a2,a1,I0.5,a)
  105 format(a4,I0.5,a)
  ! Outer loop: over inputs
  c = 0
  do i=1,x
    if (conv) then
      write(convfile,104) "train",itrain,"/",i,"_tddft_conv.txt"
    else
      write(convfile,104) "train",itrain,"/",i,"_tddft.txt"
    endif
    write(xyzfile,105) "xyz/",i,".xyz"
    INQUIRE(FILE=convfile, EXIST=xanes_file_exists)
    INQUIRE(FILE=xyzfile, EXIST=xyz_file_exists)
    if (xyz_file_exists .and. xanes_file_exists) then
      c = c + 1
      idx(c)=i
    endif
  enddo
  p=c  ! total number of XANES files named xxxxx_tddft.txt in results directory

  return
end
