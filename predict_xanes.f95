subroutine predict_xanes(wxas,kev,tik,rcm,idx,m,n,itrain,convolute)

  use functions
  implicit none
  
  integer :: i,k,m,n,idx
  character(len=36) :: predicted_xanes_file
  character(len=2) :: itrain
  real*8, dimension (m+1,n) :: tik
  real*8, dimension (m+1,1) :: xi
  real*8, dimension (m) :: rcm
  real*8, dimension (n) :: kev,fwhm,wxas
  real*8, dimension (n,1) :: tmp
  logical :: convolute
 
  ! The predicted spectrum is the hypothesis function with the optimised coefficients
  ! hypothesis_function expects rank 2 arrays, fix: make 2nd dim have 1 column
  xi(1,1) = 1.0
  xi(2:m+1,1)=rcm
  tmp = hypothesis_function(tik,xi)
  wxas = tmp(:,1)

  ! write predicted XANES to file
  !print*,'FINISHED PREDICTING.'
  !print*,'WRITING TO PREDICTED XANES FILE.'
  117 format(a19,I0.5,a6,a2,a4)
  write(predicted_xanes_file,117) "predicted/predicted",idx,"_train",itrain,".dat"
  OPEN(66,FILE=predicted_xanes_file,FORM="FORMATTED",STATUS="replace")
  110 format(3x,f7.2,2x,f13.8)
  write(66,*) 'Energy(eV)   Absorbance(a.u.)'
  do k=1,n
    write(66,110) kev(k),wxas(k)
  enddo
  close(66)

  ! Convolution (if enabled)
  if (convolute) then
    fwhm = fdmnes_gamma(kev,wxas) 
    wxas = lorentzian_broaden(kev,wxas,fwhm) 
    write(predicted_xanes_file,117) "predicted/pred_conv",idx,"_train",itrain,".dat"
    OPEN(66,FILE=predicted_xanes_file,FORM="FORMATTED",STATUS="replace")
    write(66,*) 'Energy(eV)   Absorbance(a.u.)'
    do k=1,n
      write(66,110) kev(k),wxas(k)
    enddo
    close(66)
  endif

  return
end
