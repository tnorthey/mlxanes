subroutine reducekgrid(ev,xanes,kev,kxanes,kmax,n1,n2)

  ! Reduce size of energy grid
  
  implicit none
  integer :: kmax,n1,n2,deltak,k,c
  real*8 :: remax,rkmax
  real*8 :: ev(kmax),xanes(kmax)
  real*8 :: kev(n2-n1+1),kxanes(n2-n1+1)
  
  c=0
  do k=n1,n2
    c=c+1
    kev(c) = ev(k)
    kxanes(c) = xanes(k)
  enddo
  
  return
end
