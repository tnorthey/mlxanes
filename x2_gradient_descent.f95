subroutine x2_gradient_descent(xij,ykj,tik,p,n,m,tol,Jtarget,stepin,maxiter,jprint,itrain, & 
                               converged,Jth,totgrad,dJ,dgradJ,vary_stepsize,percent)
! SUBROUTINE : x2_gradient_descent
! DESCRIPTION: Solves iteratively a gradient descent for
!              rank 2 hypothesis function defined as: h(th) = th0*x0 + th1*x1 + ... + thp*xp,
!              the rank 2 array "xij" containing the feature vectors,
!              The cost function is defined as: J(th) = (1/2np) * sum_l(sum_k(|hkj - xij|)**2),
!              The cost function is minimised with a gradient descent of constant stepsize.
! INPUTS     : xij     (rank 2 array), contains the p (j=1,p) feature vectors, each with length m+1 (i=1,m+1)
!            : p       (integer), the number of feature vectors
!            : m       (integer), the length of each feature vector
!            : n       (integer), the number of grid points in energy
!            : ykj     (rank 2 array), the supervised learning dataset (i.e. XANES spectra)
!            : hkj     (rank 2 array), the hypothesis function
!            : tol     (real), the tolerance value that defines the convergence of the gradient descent
!            : Jtarget (real), the target value of J, once reached stops the loop
!            : stepin  (real), the stepsize of each step in the gradient descent
!            : maxiter (integer), the maximum number of iterations of the gradient descent
!            : jprint  (boolean), .true. or .false. if true will print convergence information every iteration
!            : vary_stepsize (boolean), whether or not to vary the step-size 
!            : percent (boolean), whether or not to use a percentage based cost function
!=============================================
  use functions ! Use functions module containing the hypothesis and cost funcitons, and gradient
  implicit none
  integer :: m,n,p,c,i,j,k,l,maxiter,nstep
  real*8, dimension (m+1,p) :: xij
  real*8, dimension (n,p) :: ykj,hkj
  real*8, dimension (m+1,n) :: gradJ,tik
  real*8 :: tol,stepin,step,Jth,Jth_prev,totgrad,grad_prev,sumth,dJ,dgradJ,beta,Jtarget
  character(len=20) :: gradfile
  character(len=2) :: itrain
  logical :: jprint,converged,vary_stepsize,percent
  !=============================================
  !print*,'Enter gradient descent...'
  ! Initial values
  step = stepin ! save step-size to another memory slot (makes parallelisation work somehow)
  c=0          ! intialise counter
  Jth = 0.01d0   ! arbitrary inital value for cost function Jth
  dJ = 0.01d0  ! arbitrary inital value for deltaJ
  tik=0.0d0
  totgrad = 0.0d0
  beta = 2.0d0  ! factor to increase or decrease step-size by
  nstep = 5    ! step-size increases by factor beta every nstep steps
  !=============================================
  ! Write gradient descent iterations to file
  109 format(a14,a2,a4)
  !write(gradfile,109) "gradconv_train",itrain,".dat"
  !open(33,file=gradfile,form='formatted',status='replace')
  do while ( ( dJ > tol .or. Jth > Jtarget ) .and. c < maxiter )  ! begin gradient descent loop
    c=c+1  ! iteration count
    ! previous iteration values
    Jth_prev = Jth
    grad_prev = totgrad
    ! INITIALISE VALUES
    hkj = 0.0d0
    Jth = 0.0d0
    gradJ = 0.0d0
    totgrad = 0.0d0
    ! HYPOTHESIS ARRAY
    hkj = hypothesis_function(tik,xij)
    ! COST FUNCTION
    Jth = cost_function(hkj,ykj,percent)
    ! PARTIAL DERIVATIVES OF COST FUNCTION
    gradJ = gradient_theta(hkj,ykj,xij,percent)
    ! Armijo rule to adjust step-size (probably wrong currently)
    !if (Jth > Jth_prev + beta*step*totgrad) then
    !  step = beta*step
    !endif
    ! GRADIENT
    totgrad = sum(gradJ)
    ! ITERATE COEFFICIENTS
    tik = tik - step*gradJ
    ! Deltas
    dJ = abs(Jth-Jth_prev)
    dgradJ = abs(totgrad-grad_prev)
    ! Check for divergence (to very large value or infinity)
    if (dJ-1 == dJ) then
      print('(I14,4x,E12.4,4x,E12.4,4x,E12.4,4x,E12.4)'), c,Jth,dJ,totgrad,dgradJ
      print*,'WARNING: J has diverged !'
      print*,'Suggestion: Decrease step-size.'
      print*,'Exiting program...'
      stop
    endif
    ! Vary step-size
    ! Every n steps increase step-size
    if ( vary_stepsize ) then
      if ( modulo(c,nstep)==0 ) then
        step = step*beta
      endif
      ! if J increases reduce step-size
      if ( Jth > Jth_prev ) then
        step = step/beta
      endif
    endif
    ! Print information
    if (jprint) then
      if (c==1) then
        print('(a14,a14,a14,a14,a14)'),'ITERATION','Jth','|J(th)-prev|','Grad(J)','|Grad-prev|'
        !write(33,'(a14,a14,a14,a14,a14)'),'ITERATION','Jth','|J(th)-prev|','Grad(J)','|Grad-prev|'
      endif
      if ( (modulo(c,10000)==0) .or. (c==maxiter) ) then
        print('(I14,4x,E12.4,4x,E12.4,4x,E12.4,4x,E12.4)'), c,Jth,dJ,totgrad,dgradJ
        !write(33,'(I14,4x,E12.4,4x,E12.4,4x,E12.4,4x,E12.4)'), c,Jth,dJ,totgrad,dgradJ
      endif
    endif
  enddo  ! end gradient descent while loop
  !=============================================
  ! Print final information
  !print*,'FINAL INFORMATION:'
  if (c == maxiter) then
    converged = .false.
    !print*,'MAXIMUM ITERATIONS REACHED: STOPPING...'
    !print*,'(NOT CONVERGED) |Jth - prev| = ', dJ
    !print*,'(NOT CONVERGED) Grad(Jth) ', totgrad
  else
    converged = .true.
    !print*,'CONVERGENCE ACHIEVED!'
    !print*,'Converged Jth = ', Jth
    !print*,'Converged |Jth - prev| = ', dJ
    !print*,'Converged Grad(Jth) ', totgrad
  endif
  print('(I14,4x,E12.4,4x,E12.4,4x,E12.4,4x,E12.4)'), c,Jth,dJ,totgrad,dgradJ
  !write(33,'(I14,4x,E12.4,4x,E12.4,4x,E12.4,4x,E12.4)'), c,Jth,dJ,totgrad,dgradJ
  !close(33)  ! close gradient descent output file
  !print*,'Exiting gradient descent...'
  !=============================================
  return
end
