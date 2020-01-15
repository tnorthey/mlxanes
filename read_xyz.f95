subroutine read_xyz(xyzfile,natm,atmtyp,xcoo,ycoo,zcoo)
  
  ! Reads xyz file if it contains less than 1000 atoms...
  ! Input: xyzfile (name of xyz file)
  ! Outputs: natm (number of atoms)
  !        : atmtyp (list of atom types: H, He, Li, ...) 
  !        : xcoo,ycoo,zcoo (list of coordinates)
  implicit none
  
  integer :: j,natm
  character(len=17) :: xyzfile
  character(len=2) :: atmtyp(1000)
  real*8 :: xcoo(1000),ycoo(1000),zcoo(1000)
  
  !print*,'reading ',xyzfile
  open(22,file=xyzfile,form='formatted',status='old')
  read(22,'(i4)') natm
  read(22,*)
  do j=1,natm
    read(22,*) atmtyp(j),xcoo(j),ycoo(j),zcoo(j) 
  enddo
  close(22)
  
  return
end
