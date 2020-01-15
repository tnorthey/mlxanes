subroutine cmcalc(natm,atmtyp,xc,yc,zc,cm,rcm,mx,m,rcm_on)
  
  implicit none
  
  integer :: i,j,m,natm,mx,c
  character(len=2) :: atmtyp(natm)
  real*8 :: xc(natm),yc(natm),zc(natm),z(natm),dist
  real*8, dimension (natm,natm) :: cm
  real*8, dimension (m) :: rcm
  logical :: rcm_on

  ! convert atom type list to atomic number Z
  call atmtyp2z(natm,atmtyp,z)
 
  ! sort by atomic number (z) then loop over sorted indices to create CM
  ! Implement this here <----------

  ! diagonal elements
  do i=1,natm
    cm(i,i) = 0.5d0*z(i)**2.4
  enddo
  
  ! OFF DIAGONALS
  do i=1,natm
    do j=i+1,natm    
      dist = sqrt( (xc(i)-xc(j))**2 + (yc(i)-yc(j))**2 + (zc(i)-zc(j))**2 )
      if (dist .gt. 1.0e-8) then
        cm(i,j) = z(i)*z(j)/dist
      else
        cm(i,j) = 0.0d0
      endif
      ! other side of diagonal is equivalent
      cm(j,i) = cm(i,j) 
    enddo
  enddo
  
  ! CALCULATE REDUCED CM OR MODIFIED FULL CM
  106 format(f12.8)
  !OPEN(7,FILE='rcm.dat',FORM="FORMATTED",STATUS="replace")
  !write(7,*) 'CM: ',idx
  if (rcm_on) then
    rcm = cm(1:m,1)
  else
    c=0
    do i=1,mx-1
      do j=i+1,mx
        c=c+1
        rcm(c) = cm(i,j)
      enddo
    enddo
  endif
  rcm(1) = 0.0d0  ! first element is absorbing atom Z^2, set to 0 (for reduced or full CM)

  !do i=1,m
    !write(7,106) rcm(i)
  !enddo
  !close(7)

  return
end


! SUBROUTINES
subroutine atmtyp2z(natm,atmtyp,z)
  
  ! convert atom type (atmtyp) i.e. H, He, Li, etc. to atomic number Z
  ! inputs: number of atoms (natm)
  !       : list of atom types of length natm
  ! output: list of atomic numbers of length natm
  implicit none

  integer :: natm,i
  character(len=2) :: atmtyp(natm)
  real*8 :: z(natm)
  
  do i=1,natm
    if      (atmtyp(i) .eq. ' H' .or. atmtyp(i) .eq. 'H ') then
       z(i) = 1.0d0
    else if (atmtyp(i) .eq. ' B' .or. atmtyp(i) .eq. 'B ') then
       z(i) = 5.0d0
    else if (atmtyp(i) .eq. ' C' .or. atmtyp(i) .eq. 'C ') then
       z(i) = 6.0d0
    else if (atmtyp(i) .eq. ' N' .or. atmtyp(i) .eq. 'N ') then
       z(i) = 7.0d0
    else if (atmtyp(i) .eq. ' O' .or. atmtyp(i) .eq. 'O ') then
       z(i) = 8.0d0
    else if (atmtyp(i) .eq. ' P' .or. atmtyp(i) .eq. 'P ') then
       z(i) = 15.0d0
    else if (atmtyp(i) .eq. ' F' .or. atmtyp(i) .eq. 'F ') then
       z(i) = 9.0d0
    else if (atmtyp(i) .eq. ' I' .or. atmtyp(i) .eq. 'I ') then
       z(i) = 53.0d0
    else if (atmtyp(i) .eq. ' K' .or. atmtyp(i) .eq. 'K ') then
       z(i) = 19.0d0
    else if (atmtyp(i) .eq. ' S' .or. atmtyp(i) .eq. 'S ') then
       z(i) = 16.0d0
    else if (atmtyp(i) .eq. ' Y' .or. atmtyp(i) .eq. 'Y ') then
       z(i) = 39.0d0
    else if (atmtyp(i) .eq. ' W' .or. atmtyp(i) .eq. 'W ') then
       z(i) = 74.0d0
    else if (atmtyp(i) .eq. ' U' .or. atmtyp(i) .eq. 'U ') then
       z(i) = 92.0d0
    else if (atmtyp(i) .eq. ' V' .or. atmtyp(i) .eq. 'V ') then
       z(i) = 23.0d0
    else if (atmtyp(i) .eq. 'Si') then
       z(i) = 14.0d0
    else if (atmtyp(i) .eq. 'Zn') then
       z(i) =30.0d0
    else if (atmtyp(i) .eq. 'Nb') then
       z(i) =41.0d0
    else if (atmtyp(i) .eq. 'Mo') then
       z(i) =42.0d0
    else if (atmtyp(i) .eq. 'Ru') then
       z(i) =44.0d0
    else if (atmtyp(i) .eq. 'Rh') then
       z(i) =45.0d0
    else if (atmtyp(i) .eq. 'Pd') then
       z(i) =46.0d0
    else if (atmtyp(i) .eq. 'Ag') then
       z(i) =47.0d0
    else if (atmtyp(i) .eq. 'Cd') then
       z(i) =48.0d0
    else if (atmtyp(i) .eq. 'In') then
       z(i) =49.0d0
    else if (atmtyp(i) .eq. 'Hf') then
       z(i) =72.0d0
    else if (atmtyp(i) .eq. 'Ta') then
       z(i) =73.0d0
    else if (atmtyp(i) .eq. 'Re') then
       z(i) =75.0d0
    else if (atmtyp(i) .eq. 'Os') then
       z(i) =76.0d0
    else if (atmtyp(i) .eq. 'Ir') then
       z(i) =77.0d0
    else if (atmtyp(i) .eq. 'Pt') then
       z(i) =78.0d0
    else if (atmtyp(i) .eq. 'Au') then
       z(i) =79.0d0
    else if (atmtyp(i) .eq. 'Pb') then
       z(i) =82.0d0
    else if (atmtyp(i) .eq. 'Bi') then
       z(i) =83.0d0
    else if (atmtyp(i) .eq. 'Po') then
       z(i) =84.0d0
    else if (atmtyp(i) .eq. 'At') then
       z(i) =85.0d0
    else if (atmtyp(i) .eq. 'Rb') then
       z(i) =37.0d0
    else if (atmtyp(i) .eq. 'Cs') then
       z(i) =55.0d0
    else if (atmtyp(i) .eq. 'Fr') then
       z(i) =87.0d0
    else if (atmtyp(i) .eq. 'Ce') then
       z(i) =58.0d0
    else if (atmtyp(i) .eq. 'Pr') then
       z(i) =59.0d0
    else if (atmtyp(i) .eq. 'Nd') then
       z(i) =60.0d0
    else if (atmtyp(i) .eq. 'Sm') then
       z(i) =62.0d0
    else if (atmtyp(i) .eq. 'Eu') then
       z(i) =63.0d0
    else if (atmtyp(i) .eq. 'Gd') then
       z(i) =64.0d0
    else if (atmtyp(i) .eq. 'Tb') then
       z(i) =65.0d0
    else if (atmtyp(i) .eq. 'Dy') then
       z(i) =66.0d0
    else if (atmtyp(i) .eq. 'Ho') then
       z(i) =67.0d0
    else if (atmtyp(i) .eq. 'Er') then
       z(i) =68.0d0
    else if (atmtyp(i) .eq. 'Tm') then
       z(i) =69.0d0
    else if (atmtyp(i) .eq. 'Yb') then
       z(i) =70.0d0
    else if (atmtyp(i) .eq. 'Lu') then
       z(i) =71.0d0
    else if (atmtyp(i) .eq. 'Th') then
       z(i) =90.0d0
    else if (atmtyp(i) .eq. 'Pa') then
       z(i) =91.0d0
    else if (atmtyp(i) .eq. 'Np') then
       z(i) =93.0d0
    else if (atmtyp(i) .eq. 'Pu') then
       z(i) =94.0d0
    else if (atmtyp(i) .eq. 'Am') then
       z(i) =95.0d0
    else if (atmtyp(i) .eq. 'Cm') then
       z(i) =96.0d0
    else if (atmtyp(i) .eq. 'Zr') then
       z(i) =40.0d0
    else if (atmtyp(i) .eq. 'Cr') then
       z(i) = 24.0d0
    else if (atmtyp(i) .eq. 'Mn') then
       z(i) = 25.0d0
    else if (atmtyp(i) .eq. 'Ni') then
       z(i) = 28.0d0 
    else if (atmtyp(i) .eq. 'Cu') then
       z(i) = 29.0d0
    else if (atmtyp(i) .eq. 'As') then
       z(i) = 33.0d0
    else if (atmtyp(i) .eq. 'Na') then
       z(i) = 11.0d0
    else if (atmtyp(i) .eq. 'Be') then
       z(i) = 4.0d0
    else if (atmtyp(i) .eq. 'Li') then
       z(i) = 3.0d0
    else if (atmtyp(i) .eq. 'Mg') then
       z(i) = 21.0d0
    else if (atmtyp(i) .eq. 'Ca') then
       z(i) = 20.0d0
    else if (atmtyp(i) .eq. 'Co') then
       z(i) = 27.0d0
    else if (atmtyp(i) .eq. 'Ba') then
       z(i) = 56.0d0
    else if (atmtyp(i) .eq. 'Sr') then
       z(i) = 38.0d0
    else if (atmtyp(i) .eq. 'Sc') then
       z(i) = 21.0d0
    else if (atmtyp(i) .eq. 'Ra') then
       z(i) = 88.0d0
    else if (atmtyp(i) .eq. 'Ge') then
       z(i) = 32.0d0
    else if (atmtyp(i) .eq. 'Sn') then
       z(i) = 50.0d0
    else if (atmtyp(i) .eq. 'Ti') then
       z(i) = 81.0d0
    else if (atmtyp(i) .eq. 'Ga') then
       z(i) = 31.0d0
    else if (atmtyp(i) .eq. 'Cl') then
       z(i) = 17.0d0
    else if (atmtyp(i) .eq. 'Br') then
       z(i) = 35.0d0
    else if (atmtyp(i) .eq. 'Se') then
       z(i) = 34.0d0
    else if (atmtyp(i) .eq. 'Al') then
        z(i) = 13.0d0
    else if (atmtyp(i) .eq. 'Sb') then        
        z(i) = 51.0d0
    else if (atmtyp(i) .eq. 'Te') then        
        z(i) = 52.0d0
    else if (atmtyp(i) .eq. 'Fe') then
        z(i) = 26.0d0
    else 
        print*,'Element undefined. Manually add to cmcalc.'
        print*,atmtyp(i)
        print*,i
        stop
    endif
  enddo

end subroutine atmtyp2z
