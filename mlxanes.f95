program mlxanes
  ! Title: MLXANES Machine Learning XANES Program with OpenMP parallelisation 
  ! Author: T. Northey
  ! Date: 2018-2019
  use omp_lib  ! Use OPENMP libraries 
  implicit none
  integer :: nmax,mmax,p,n,n1,n2,m,j,k,natm,maxiter,x,c,stat,mx
  integer, allocatable, dimension(:) :: idx,tmp
  character(len=17) :: xyzfile
  character(len=18) :: coeffs_file
  character(len=10) :: info_file
  character(len=28) :: xasfile
  character(len=2), dimension(1000) :: atmtyp
  character(len=2) :: itrain
  real*8, dimension(:), allocatable :: xanes,kxanes,ev,kev,wxas,rcm
  real*8, dimension(:), allocatable :: Jth,totgrad,dJ,dtotgrad
  real*8, dimension(1000) :: xcoo,ycoo,zcoo
  real*8, dimension(:,:), allocatable :: cm,xij,yjk,tik
  real*8 :: tol,step,t1,t2,Jtarget
  logical :: conv,jprint,rcm_on,learn_on,file_exists,vary_stepsize,percent
  logical, dimension(:), allocatable :: converged
  ! begin timer
  call cpu_time(t1)
  ! Read from setup.txt file
  print*,'Beginning mlxanes.'
  print*,'Reading setup.txt'
  x = 10000  ! The number of files the program will try to read in order of 00001,00002,...,x
  call read_setup(learn_on,tol,Jtarget,maxiter,step,mmax,m,nmax,n1,n2,itrain,conv,jprint,rcm_on,vary_stepsize,percent)
  print*,'Finished reading setup.txt'
  call itrain_exists(itrain)
  n=n2-n1+1  ! number of energy grid points
  ! Allocate arrays
  allocate (ev(nmax),xanes(nmax))
  allocate (kev(n),kxanes(n),wxas(n))
  allocate (cm(mmax,mmax))
  ! Count the number of spectrum files in 'itrain' with corresponding xyz files in xyz/
  x = 10000  ! The number of files the program will try to read in order of 00001,00002,...,x
  allocate (tmp(x))
  call count_xanes(c,x,p,itrain,conv,tmp)
  print*,'Number of spectra files in "training data": ',p
  allocate (yjk(p,n))
  allocate (idx(p))
  idx = tmp(1:p)
  if (conv) then
    print*,'TRAINING FROM CONVOLUTED SPECTRA.'
  else
    print*,'TRAINING FROM NOT CONVOLUTED SPECTRA.'
  endif
  ! Allocate rcm or full cm
  if (rcm_on) then
    print*,'Using reduced Coulomb matrix of size ',m
  else
    mx = m  ! store original m
    m=int((m*m-m)/2)
    print*,'Using half the non-diagonals of Coulomb matrix; total elements: ',m
  endif
  allocate (rcm(m),xij(m+1,p),tik(m+1,n))
  allocate (Jth(n),totgrad(n),dJ(n),dtotgrad(n),converged(n))
  ! print array information
  print*,'Size of RCM = ',size(rcm)  
  print*,'Size of xij = ',size(xij,1),size(xij,2)  
  print*,'Size of tik = ',size(tik,1),size(tik,2)  
  ! Initiate files
  open(unit=1234, iostat=stat, file='rcm.dat', status='old')
  if (stat == 0) close(1234, status='delete')
  ! Formatting
  103 format(a,I0.5,a)
  104 format(a5,a2,a1,I0.5,a)
  109 format(a12,a2,a4)
  110 format(a5,a12,a14,a14,a14,a14)
  111 format(I0.5,l12,e14.6,e14.6,e14.6,e14.6)
  write(coeffs_file,109) "coeffs_train",itrain,".dat"
  info_file = "final.info"
  !============================================!
  ! LEARN SECTION                              !
  !============================================!
  if (learn_on) then
    print*,'Enter learn gradient descent section.'
    open(unit=66,file=info_file,form='formatted')
    write(66,*) 'Final information:'
    write(66,*) 'Training set = train',itrain
    write(66,*) 'Maxiter = ',maxiter
    write(66,*) 'Step-size = ',step
    write(66,*) 'Tolerance = ',tol
    write(66,*) ' '
    print*,'Maxiter = ',maxiter
    print*,'Step-size = ',step
    print*,'Tolerance = ',tol
    do j=1,p  ! Loop over training set of size p, j=1,p
      write(xyzfile,103) "xyz/",idx(j),".xyz"
      ! READ XYZ FILE
      call read_xyz(xyzfile,natm,atmtyp,xcoo,ycoo,zcoo)
      ! CALCULATE THE COULOMB MATRIX (CM)
      call cmcalc(natm,atmtyp,xcoo,ycoo,zcoo,cm,rcm,mx,m,rcm_on)
      ! DEFINE FEATURE VECTORS AS ELEMENTS OF CM
      xij(1,j)=1.0
      xij(2:m+1,j)=rcm
      ! READ XANES from x_tddft_conv.txt or x_tddft.txt
      if (conv) then
        write(xasfile,104) "train",itrain,"/",idx(j),"_tddft_conv.txt"
      else
        write(xasfile,104) "train",itrain,"/",idx(j),"_tddft.txt"
      endif
      call read_xanes(xasfile,ev,xanes,nmax)
      call reducekgrid(ev,xanes,kev,kxanes,nmax,n1,n2)
      yjk(j,:)=kxanes
    enddo  ! End loop over training set of size p, j=1,p
    ! If coeffs files exists replace it with blank file, if not new file
    !call file_creator(coeffs_file)
    ! If info files exists replace it with blank file, if not new file
    !call file_creator(info_file)
    ! Gradient descent to optimise "theta_ik" coefficients
    print*,'Writing coefficients to file ...'
    ! 0pen file and append coeffs
    open(unit=55,file=coeffs_file,form='formatted')
    write(66,110) "k","converged?","J","dJ","gradJ","dgradJ"
    ! ======= Begin OpenMP parallel section ==================
    !$omp parallel
    !$omp do
    do k=1,n  ! Loop over energy grid
      print*,'Optimising energy grid point: ',k,'/',n
      call x2_gradient_descent(xij,yjk(:,k),tik(:,k),p,1,m,tol,Jtarget,step,maxiter,jprint,itrain, & 
                               converged(k),Jth(k),totgrad(k),dJ(k),dtotgrad(k),vary_stepsize,percent)
    enddo
    !$omp end do
    !$omp end parallel
    ! ======= End OpenMP parallel section ====================
    ! Write coeffs to file and convergence info to info.out
    do k=1,n
      write(66,111) k+n1,converged(k),Jth(k),dJ(k),totgrad(k),dtotgrad(k)
      write(55,*) tik(:,k)
    enddo
    close(55)  ! close file
    close(66)  ! close file
  else
    print*,'Skipping learn...'
  endif  ! end if learn==1
  !============================================!
  ! PREDICTION SECTION                         !
  !============================================!
  inquire(file=coeffs_file, exist=file_exists)
  if (file_exists) then
    open(44,file=coeffs_file,form="FORMATTED",status="old")
    print*,'PREDICTING from coefficient file ',coeffs_file
  else
    print*,'The file: ',coeffs_file,' does not exist to read.'
    print*,'Exiting program.'
    stop
  endif
  do k=1,n
    read(44,*) tik(:,k)
  enddo
  close(44)
  print*,'All xyz files in "xyz" directory will be predicted...'
  ! read energy range from first existing xanes file
  write(xasfile,104) "train",itrain,"/",idx(1),"_tddft.txt"
  call read_xanes(xasfile,ev,xanes,nmax)
  call reducekgrid(ev,xanes,kev,kxanes,nmax,n1,n2)
  c=0
  do j=1,x  ! loop over all the xyz files and predict all of them
    write(xyzfile,103) "xyz/",j,".xyz"
    inquire(file=xyzfile, exist=file_exists)
    if (file_exists) then
      c=c+1
      ! READ XYZ FILE
      print*,'Predicting: ',xyzfile
      call read_xyz(xyzfile,natm,atmtyp,xcoo,ycoo,zcoo)
      ! CALCULATE THE COULOMB MATRIX (CM)
      call cmcalc(natm,atmtyp,xcoo,ycoo,zcoo,cm,rcm,mx,m,rcm_on)
      call predict_xanes(wxas,kev,tik,rcm,j,m,n,itrain,.true.)
    endif
  enddo  ! end loop over files
  print*,c,' spectra were predicted. Files written to predicted/ directory.'
  call cpu_time(t2)
  write ( *, * ) 'Elapsed CPU time = ', t2 - t1 ,' seconds.'
  print*,'END OF PROGRAM.'
end


! SUBROUTINES
subroutine itrain_exists(itrain)

  implicit none
  character(len=2) :: itrain
  character(len=10) :: traindir
  logical :: file_exists

  113 format(a7,a2,a1)
  write(traindir,113) './train',itrain,'/'
  inquire(file=traindir, exist=file_exists)
  ! Check if the itrain directory exists
  if (file_exists) then
    print*,'Reading from training set train',itrain
  else
    print*,'WARNING: ',traindir,' does not exist! Exiting.'
    stop
  endif

return
end subroutine itrain_exists

