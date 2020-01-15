module functions
  contains
   
  function hypothesis_function(tik,xij) result(hkj)
 
    implicit none
    integer :: i,j,k,m,n,p
    real*8, dimension (:,:), intent(in) :: tik,xij
    real*8, dimension (size(tik,2),size(xij,2)) :: hkj
    
    m=size(xij,1)
    p=size(xij,2)
    n=size(tik,2)
    do i=1,m+1
    do k=1,n
    do j=1,p
      hkj(k,j) = hkj(k,j) + tik(i,k)*xij(i,j)
    enddo
    enddo
    enddo

  end function hypothesis_function


  function cost_function(hkj,yjk,percent) result(Jth)
 
    implicit none
    integer :: n,p
    real*8, dimension (:,:), intent(in) :: hkj,yjk
    real*8 :: Jth
    logical :: percent
    
    n=size(hkj,1)
    p=size(hkj,2)
    if (percent) then
      Jth = 0.5d0*sum( ((hkj-yjk)**2)/yjk )/dble(p*n)
    else
      Jth = 0.5d0*sum((hkj-yjk)**2)/dble(p*n)
    endif

  end function cost_function


  function gradient_theta(hkj,yjk,xij,percent) result(gradJ)
    
    implicit none
    
    real*8, dimension (:,:), intent(in) :: yjk,hkj,xij
    integer :: i,j,k,m,n,p
    real*8, dimension (size(xij,1),size(yjk,1)) :: gradJ
    real*8, dimension (size(hkj,1),size(hkj,2)) :: tmp
    logical :: percent
    
    m=size(xij,1)
    p=size(xij,2)
    n=size(yjk,1)
    ! PARTIAL DERIVATIVES OF COST FUNCTION
    if (percent) then
      tmp = (hkj - yjk) / yjk
    else
      tmp = (hkj - yjk) 
    endif
    do j=1,p
    do k=1,n
    do i=1,m+1
      gradJ(i,k) = gradJ(i,k) + tmp(k,j)*xij(i,j)
    enddo
    enddo
    enddo
    gradJ = gradJ / dble(p*n)
    ! END PARTIAL DERIVATIVES OF COST FUNCTION
  end function gradient_theta


  function fdmnes_gamma(E,xanes) result(gam)
    ! Apply convolution to XANES
    implicit none
    
    REAL, PARAMETER :: Pi = 3.1415927
    real*8, dimension(:), intent(in) :: E,xanes
    real*8, dimension(:), allocatable :: gam,Eshift
    real*8 :: gamma_hole,gamma_max,E_larg,Efermi,E_cent,E1,h
    integer :: j,nx

    ! fixed values (eV)
    gamma_hole = 2.0  ! the fwhm at the core-hole
    gamma_max = 15  ! the maximum value of the fwhm at higher energy (15 is FDMNES default)
    E_larg = 30  ! Changes the curvature of the arctangent (30 is FDMNES default)
    Efermi = 0   ! The Fermi energy (this isn't default, it is set by FDMNES and is element dependent; I don't know what it is for Fe)
    E_cent = 30  ! The centre of the arctan function (30 is FDMNES default)
    E1 = 7112    ! The edge energy (for Fe ~= 7112)

    ! fdmnes broadening equation
    nx = size(xanes)
    allocate(gam(nx),Eshift(nx))
    Eshift = E - E1   ! shift origin to edge
    do j=1,nx
        h = (Eshift(j) - Efermi) / E_cent
        gam(j) = gamma_max*(0.5 + (1/pi)*atan(((pi/3)*(gamma_max/E_larg))*(h-h**(-2))))
    enddo
    gam = gam + gamma_hole

  end function fdmnes_gamma
  

  function lorentzian_broaden(E,xanes,fwhm) result(bxanes)

    implicit none

    real, parameter :: pi = 3.1415927
    real*8, dimension(:), intent(in) :: E,xanes,fwhm
    real*8, dimension(:), allocatable :: bxanes,l
    real*8 :: dE,tmp
    integer :: i,j,k,nx

    nx = size(xanes)
    allocate (bxanes(nx),l(nx))
      
    do i=1,nx
      bxanes(i) = 0.0d0
      do j=1,nx
            
        dE = E(i)-E(j)
        if (dE==0) then
          l(i) = 1
        else
          tmp = 0.5*fwhm(i)
          l(i) = tmp / ( dE**2 + tmp**2 )  ! Lorentzian
        endif
        bxanes(i) = bxanes(i) + xanes(j) * l(i)
      enddo
    enddo
    bxanes = sum(xanes) * bxanes / sum(bxanes)  ! normalise to original xanes
    
  end function lorentzian_broaden


end module functions
