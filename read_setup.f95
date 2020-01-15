subroutine read_setup(learn,tol,Jtarget,maxiter,step,mmax,m,nmax,n1,n2,itrain,conv,jprint,rcm,vary_stepsize,percentin)
  
!     SUBROUTINE: read_setup
!     DESCRIPTION: Reads values from setup.txt
!     INPUTS: learn (logical), if true learning is on
!             tol (real), tolerance value for the gradient descent
!             Jtarget (real), target value for cost function
!             maxiter (integer), maximum number of gradient descent iterations
!             step (real), gradient descent step-size
!             mmax (integer), number of atoms
!             m (integer), size of reduced Coulomb matrix
!             nmax (integer), total number of energy grid points in output files
!             n1 (integer), use energy grid starting at n1
!             n2 (integer), use energy grid ending at n2
!             itrain (character), defines which training set to use (01 means read from train01 directory)
!             conv (logical), use convoluted spectra (true) or not (false)
!             jprint (logical), print gradient descent converence information every iteration if true
!             rcm (logical), use reduced Coulomb matrix (true) or full Coulomb matrix (false)
!             vary_stepsize (logical), whether or not to use a variable step-size
  
  implicit none
  integer :: learnin,maxiter,mmax,m,nmax,n1,n2,convin,jprintin,rcmin,varyin,percentin
  character(len=2) :: itrain
  character :: del
  real*8 :: tol,step,Jtarget
  logical :: conv,jprint,rcm,learn,vary_stepsize,percent
  
  ! Read from setup file to define variables
  open(3,file='setup.txt',form='formatted',status='old')
  read(3,*) del,learnin
  read(3,*) del,tol
  read(3,*) del,Jtarget
  read(3,*) del,maxiter
  read(3,*) del,step
  read(3,*) del,mmax
  read(3,*) del,m
  read(3,*) del,nmax
  read(3,*) del,n1
  read(3,*) del,n2
  read(3,*) del,itrain
  read(3,*) del,convin
  read(3,*) del,jprintin
  read(3,*) del,rcmin
  read(3,*) del,varyin
  read(3,*) del,percentin
  close(3)

  if (learnin==1) then
    learn=.true.
  else
    learn=.false.
  endif
  
  if (convin==1) then
    conv=.true.
  else
    conv=.false.
  endif
   
  if (jprintin==1) then
    jprint=.true.
  else
    jprint=.false.
  endif

  if (rcmin==1) then
    rcm=.true.
  else
    rcm=.false.
  endif

  if (varyin==1) then
    vary_stepsize=.true.
  else
    vary_stepsize=.false.
  endif

  if (percentin==1) then
    percent=.true.
  else
    percent=.false.
  endif

  return
end
