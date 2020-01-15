subroutine read_xanes(filename,ev,xanes,emax)
  implicit none
   
  character(len=28) :: filename         
  integer :: n,emax
  real*8, dimension (emax) :: ev,xanes
    
  !     READ IN XANES
  open(2,file=filename,form='formatted',status='old')
  read(2,*)
  do n=1,emax
    read(2,110) ev(n),xanes(n)
  enddo
  close(2)
  110 format(3x,f7.2,2x,f13.8)
  
  return
end
