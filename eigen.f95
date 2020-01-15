!  =============================================================================
!
!  DSYEV Example.
!  ==============
!
!  Program computes all eigenvalues and eigenvectors of a real symmetric
!  matrix A:
!
!    1.96  -6.49  -0.47  -7.20  -0.65
!   -6.49   3.80  -6.39   1.50  -6.34
!   -0.47  -6.39   4.17  -1.51   2.67
!   -7.20   1.50  -1.51   5.70   1.80
!   -0.65  -6.34   2.67   1.80  -7.10
!
!  Description.
!  ============
!
!  The routine computes all eigenvalues and, optionally, eigenvectors of an
!  n-by-n real symmetric matrix A. The eigenvector v(j) of A satisfies
!
!  A*v(j) = lambda(j)*v(j)
!
!  where lambda(j) is its eigenvalue. The computed eigenvectors are
!  orthonormal.
!
!  Example Program Results.
!  ========================
!
! DSYEV Example Program Results
!
! Eigenvalues
! -11.07  -6.23   0.86   8.87  16.09
!
! Eigenvectors (stored columnwise)
!  -0.30  -0.61   0.40  -0.37   0.49
!  -0.51  -0.29  -0.41  -0.36  -0.61
!  -0.08  -0.38  -0.66   0.50   0.40
!   0.00  -0.45   0.46   0.62  -0.46
!  -0.80   0.45   0.17   0.31   0.16
!  =============================================================================
!

!  =============================================================================
      subroutine eigen(cm,cme,n)

      implicit none

!     .. Parameters ..
      INTEGER          N,M,P
!      INTEGER          LDA
!      PARAMETER        (LDA=N)
      INTEGER          LWMAX
      PARAMETER        (LWMAX=10000)
      double precision cm(n,n),cme(n)
!
!     .. Local Scalars ..
      INTEGER          INFO,LWORK
!
!     .. Local Arrays ..
!      DOUBLE PRECISION A(LDA,N),W(N),WORK(LWMAX)
      DOUBLE PRECISION A(N,N),W(N),WORK(LWMAX)
!
!     .. External Subroutines ..
      EXTERNAL         DSYEV
      EXTERNAL         PRINT_MATRIX
!
!     .. Intrinsic Functions ..
      INTRINSIC        INT, MIN
!
!     .. Executable Statements ..
!     WRITE(*,*)'DSYEV Example Program Results'
!
!     Query the optimal workspace.
!      DO M=1,LDA
      DO M=1,N
        DO P=1,N
           A(M,P)=CM(M,P)
        ENDDO
      ENDDO
!
!
      LWORK = -1
!      CALL DSYEV( 'Vectors', 'Upper', N, A, LDA, W, WORK, LWORK, INFO )
      CALL DSYEV( 'Vectors', 'Upper', N, A, N, W, WORK, LWORK, INFO )
      LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
!
!     Solve eigenproblem.
!
!      CALL DSYEV( 'Vectors', 'Upper', N, A, LDA, W, WORK, LWORK, INFO )
      CALL DSYEV( 'Vectors', 'Upper', N, A, N, W, WORK, LWORK, INFO )
!
!     Check for convergence.
!
      IF( INFO.GT.0 ) THEN
         WRITE(*,*)'The algorithm failed to compute eigenvalues.'
         STOP
      END IF
!
!     Print eigenvalues.
!
!     CALL PRINT_MATRIX( 'Eigenvalues', 1, N, W, 1 )
      P=0
       DO M=N,1,-1
        P=P+1
        cme(P)=W(M)
      ENDDO
!
!     Print eigenvectors.
!
!     CALL PRINT_MATRIX( 'Eigenvectors (stored columnwise)', N, N, A,
!    $                   LDA )
      
      return
      END
!
!     End of DSYEV Example.
!
!  =============================================================================
!
!     Auxiliary routine: printing a matrix.
!
      SUBROUTINE PRINT_MATRIX( DESC, M, N, A, LDA )
      CHARACTER*(*)    DESC
      INTEGER          M, N, LDA
      DOUBLE PRECISION A( LDA, * )
!
      INTEGER          I, J
!
      WRITE(*,*)
      WRITE(*,*) DESC
      DO I = 1, M
         WRITE(*,9998) ( A( I, J ), J = 1, N )
      END DO
!
 9998 FORMAT( 11(:,1X,F6.2) )
      RETURN
      END
