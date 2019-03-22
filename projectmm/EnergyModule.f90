!!!____________________________________________________________________________
!!! Author: Aaron Eisses
!!! Studentnumber: 10593527
!!! University: University of Amsterdam
!!! Date changed: 11-02-2019
!!!____________________________________________________________________________
!!!
!!! PROJECT SCIENTIFIC COMPUTING AND PROGRAMMING: MOLECULAR MECHANICS
!!!
!!! EnergyModule.f90
!!!
!!! This file is the module energy. Here the functions/subroutines are given
!!! which are used to finally calculate the total energy of a certain 
!!! configuration. The total energy function is split into stretch, bending,
!!! torsional and non-bonding (electronic and van der Waals) energy.
!!!____________________________________________________________________________

module EnergyModule
use DataModule

   implicit none
   save
   private
   public               :: StretchEnergy, BendEnergy, NonBondedEnergy, TotalEnergy 

contains

   subroutine StretchEnergy(EnStretch, Molecule, Variables, Bonds, NumberofAtoms)
      real, intent(inout)                       :: EnStretch 
      type (Atom), intent(inout), allocatable   :: Molecule(:)
      type (Parameters), intent(in)             :: Variables
      type (Binding), intent(inout)             :: Bonds(:)
      integer                                   :: i
      integer, intent(in)                       :: NumberofAtoms

      EnStretch = 0   
     
      !!! Calculate the energy of bonds with E=K(r-r0)^2. It loops over all the bonds and the if statement
      !!! is there to seperate the bonds that are C-C and C-H to make sure that the proper constants are
      !!! used here.

      do i=1,size(Bonds)
         if (Bonds(i)%elements == 'CH') then
            EnStretch = EnStretch + Variables%ForceConstantCH*((Bonds(i)%length - Variables%EquiBondCH)**2)
         else
            EnStretch = EnStretch + Variables%ForceConstantCC*((Bonds(i)%length - Variables%EquiBondCC)**2)
         endif
      enddo     

   endsubroutine

   subroutine BendEnergy(EnBend, Molecule, Variables, Bonds, NumberofAtoms)
      real, intent(inout)                       :: EnBend
      type (Atom), intent(inout), allocatable   :: Molecule(:)
      type (Parameters), intent(in)             :: Variables
      type (Binding), intent(inout)             :: Bonds(:) 
      real, allocatable                         :: Angles(:)
      integer                                   :: i, j, k, l
      integer, intent(in)                       :: NumberofAtoms
      real                                      :: Angle

      EnBend = 0
      l = 0

      !!! The angles are calculated between the atoms of interest. The structure beneath is quite
      !!! complicated so here is the explanation. First, it loops over all of the atoms in the 
      !!! molecule. In the first if statement, only the C atoms are selected in the molecule. Since 
      !!! H-atoms can only have one bond here, they can't be the connecting atom between two bonds.
      !!! Then, there is a double loop with an if statement in it. In here, two bonds are taken and
      !!! there is being checked if they are connected by the same atom but with another atom at the
      !!! other side of the bond. If so, an angle(l) is created with the corresponding angle between
      !!! the two bonds. 
        
      



      do i = 1,NumberofAtoms
         if (Molecule(i)%element == 'C') then
            do j = 1,size(Bonds)
               do k = j,size(Bonds)
                  if (j /= k .and. Bonds(j)%FirstAtom == i .and. Bonds(k)%FirstAtom == i) then
                !     print *, '------------First-----------------------'
                !     print *, 'i =', i
                !     print *, 'j =', j, 'Bond(j) first atom =', Bonds(j)%FirstAtom, 'Bond(j) second atom =', Bonds(j)%SecondAtom
                !     print *, 'k =', k, 'Bond(k) first atom =', Bonds(k)%FirstAtom, 'Bond(k) second atom =', Bonds(k)%SecondAtom
                     
                     
                     
                     Angle = -((BondLength(Bonds(j)%SecondAtom, Bonds(k)%SecondAtom, Molecule))**2 - &
                             (Bonds(j)%length)**2 - (Bonds(k)%length)**2) / (2*(Bonds(j)%length)*   &
                             (Bonds(k)%length))
                     call AddToList_Angle(Angles, Angle)
                  endif
                  if (j /= k .and. Bonds(j)%SecondAtom == i .and. Bonds(k)%FirstAtom == i) then
                !     print *, '-----------------Second--------------------'
                !     print *, 'i =', i
                !     print *, 'j =', j, 'Bond(j) first atom =', Bonds(j)%FirstAtom, 'Bond(j) second atom =', Bonds(j)%SecondAtom
                !     print *, 'k =', k, 'Bond(k) first atom =', Bonds(k)%FirstAtom, 'Bond(k) second atom =', Bonds(k)%SecondAtom
                     
                          
                          
                     Angle = -((BondLength(Bonds(j)%SecondAtom, Bonds(k)%SecondAtom, Molecule))**2 - &
                             ((Bonds(j)%length)**2 - (Bonds(k)%length)**2) / (2*(Bonds(j)%length)*  &
                             (Bonds(k)%length)))
                     call AddToList_Angle(Angles, Angle)
                  endif
               enddo
            enddo
         endif
      enddo
      !!! The bending energy is calculated by looping over all the angles 
      
      do l = 1, size(Angles)
         if (abs(acos(Angles(l)) - Variables%EquiAngle) < 1) then
            EnBend = EnBend + Variables%ForceConstantAngle*((acos(Angles(l)) - Variables%EquiAngle)**2)
         endif
      enddo
   endsubroutine

   subroutine TorsionalEnergy(EnTors, Molecule, Variables, Bonds, NumberofAtoms)
      real, intent(inout)                       :: EnTors
      type (Atom), intent(inout), allocatable   :: Molecule(:)
      type (Parameters), intent(in)             :: Variables
      type (Binding), intent(inout)             :: Bonds(:)
      integer                                   :: i, j, k, l, m, n
      integer, intent(in)                       :: NumberofAtoms
      real                                      :: VectorX1, VectorY1, VectorZ1, VectorX2, VectorY2, VectorZ2, phi
      type (Planes), allocatable                :: PlanesList(:)
      type (Planes)                             :: Plane

      EnTors = 0
      l = 0
        
      !!! This structure is quite similar to the one in the subroutine 'BendEnergy'. In the same way,
      !!! the different planes in the molecule are figured out, since the bonds desribe the
      !!! planes of the molecule. In the inner if statements, the planes are calculated by using
      !!! the coordinateds of three points on that plane.  
      
      do i = 1,NumberofAtoms
         if (Molecule(i)%element == 'C') then
            do j = 1,size(Bonds)
               do k = j,size(Bonds)
                  if (j /= k .and. Bonds(j)%FirstAtom == i .and. Bonds(k)%FirstAtom == i ) then
                     VectorX1 = Molecule(Bonds(j)%SecondAtom)%x - Molecule(i)%x
                     VectorY1 = Molecule(Bonds(j)%SecondAtom)%y - Molecule(i)%y
                     VectorZ1 = Molecule(Bonds(j)%SecondAtom)%z - Molecule(i)%z
                     VectorX2 = Molecule(Bonds(k)%SecondAtom)%x - Molecule(i)%x
                     VectorY2 = Molecule(Bonds(k)%SecondAtom)%y - Molecule(i)%y
                     VectorZ2 = Molecule(Bonds(k)%SecondAtom)%z - Molecule(i)%z
                     !!! Determine the values of a, b and c in the standard formula
                     !!! a(x-x0) + b(y-y0) + c(z-z0) = 0
                     
                     Plane%a = (VectorY1*VectorZ2)-(VectorZ1*VectorY2)
                     Plane%b = (VectorX1*VectorZ2)-(VectorZ1*VectorX2)
                     Plane%c = (VectorX1*VectorY2)-(VectorY1*VectorX2)

                     call AddToList_Plane(Planeslist, Plane)
                     
             
                  elseif (j /= k .and. Bonds(j)%SecondAtom == i .and. Bonds(k)%FirstAtom == i) then
                     VectorX1 = Molecule(Bonds(j)%SecondAtom)%x - Molecule(i)%x
                     VectorY1 = Molecule(Bonds(j)%SecondAtom)%x - Molecule(i)%y
                     VectorZ1 = Molecule(Bonds(j)%SecondAtom)%x - Molecule(i)%z
                     VectorX2 = Molecule(Bonds(k)%SecondAtom)%x - Molecule(i)%x
                     VectorY2 = Molecule(Bonds(k)%SecondAtom)%x - Molecule(i)%y
                     VectorZ2 = Molecule(Bonds(k)%SecondAtom)%x - Molecule(i)%z

                     !!! Determine the values of a, b and c in the standard formula
                     !!! a(x-x0) + b(y-y0) + c(z-z0) = 0
                     
                     Plane%a = (VectorY1*VectorZ2)-(VectorZ1*VectorY2)
                     Plane%b = (VectorX1*VectorZ2)-(VectorZ1*VectorX2)
                     Plane%c = (VectorX1*VectorY2)-(VectorY1*VectorX2)

                     call AddToList_Plane(PlanesList, Plane)
                  endif
               enddo
            enddo
         endif
      enddo
      
      !!! Here, there is a double do loop with an if statement inside. There is being looped over two
      !!! different planes and the angle between the planes is calculated. With this angle, the 
      !!! torsional energy can be calculated.
      
      do m = 1, size(PlanesList)
         do n = m, size(PlanesList)
            if (n /= m) then
               phi = (abs((PlanesList(m)%a*PlanesList(n)%a) + (PlanesList(m)%b*PlanesList(n)%b) + &
                       (PlanesList(m)%c*PlanesList(n)%c)) /     &
                     (sqrt(PlanesList(m)%a**2 + PlanesList(m)%b**2 + PlanesList(m)%c**2) * sqrt(PlanesList(m)%a**2 +       &
                     PlanesList(m)%b**2 + PlanesList(m)%c**2)))
               if (abs(phi) <= 1) then
                  phi = acos(phi)
                  EnTors = EnTors + 0.5*Variables%V1*(1 + cos(Variables%n*phi - Variables%gama))
               endif
            endif 
         enddo
      enddo
   
   endsubroutine

   subroutine NonBondedEnergy(EnNonBond, Molecule, Variables, NumberofAtoms)
      real, intent(inout)                       :: EnNonBond
      type (Atom), intent(inout), allocatable   :: Molecule(:)
      type (Parameters), intent(in)             :: Variables
      integer                                   :: i, j
      real(8)                                   :: Kelec, Rij, Avogadro
      integer, intent(in)                       :: NumberofAtoms
      
      EnNonBond = 0
      Kelec = 8.987551787*(10**9) ! Constant in Law of Coulomb
      Avogadro = 6.022*(10**23) ! Constant of Avogadro
      
      !!! In this double loop, the non bonded energy is calculated between all the atoms
      !!! with both the electronic effect and the Van Der Waals effect. The double loop
      !!! is there since all the atoms have these interactions with each other. Then there
      !!! are three if statements to determine what kind of atoms interact with each other
      !!! (C-C, C-H, or H-H).  

      do i = 1, NumberofAtoms
         do j = 1, NumberofAtoms
            if (i /= j) then
               Rij = BondLength(i, j, Molecule)
               if (Molecule(i)%element == 'C' .and. Molecule(j)%element == 'C') then
                  EnNonBond = EnNonBond + Avogadro*((Variables%ChargeC**2) * Kelec / (Rij/(10**10))) 
                  EnNonBond = EnNonBond + Variables%WellDepthC*(((Variables%RstarC/Rij)**12)- 2*((Variables%RstarC/Rij)**6))
               elseif (Molecule(i)%element == 'H' .and. Molecule(j)%element == 'H') then
                  EnNonBond = EnNonBond + Avogadro*((Variables%ChargeH**2) * Kelec / (Rij/(10**10))) 
                  EnNonBond = EnNonBond + Variables%WellDepthH*(((Variables%RstarH/Rij)**12)- 2*((Variables%RstarH/Rij)**6)) 
               else    
                  EnNonBond = EnNonBond + Avogadro*((Variables%ChargeH*Variables%ChargeC) * Kelec / (Rij/(10**10)))
                  EnNonBond = EnNonBond + sqrt(Variables%WellDepthH*Variables%WellDepthC) *             &
                  (((((Variables%RstarC/2)+(Variables%RstarH/2))/Rij)**12) - 2*((((Variables%RstarC/2)+ &
                  (Variables%RstarH/2))/Rij)**6))
               endif
            endif
         enddo
      enddo
   
   endsubroutine 
   
   subroutine TotalEnergy(TotEner, Molecule, Variables, Bonds, NumberofAtoms)
      real, intent(inout)                               :: TotEner
      type (Atom), intent(inout), allocatable           :: Molecule(:)               
      type (Parameters), intent(in)                     :: Variables
      type (Binding), intent(inout), allocatable        :: Bonds(:)
      real                                              :: EnStretch, EnBend, EnNonBond, EnTors
      integer, intent(in)                               :: NumberofAtoms

       
      call AssigningBonds(Bonds, NumberofAtoms, Molecule, Variables)
      call StretchEnergy(EnStretch, Molecule, Variables, Bonds, NumberofAtoms)
      call BendEnergy(EnBend, Molecule, Variables, Bonds, NumberofAtoms)
      call TorsionalEnergy(EnTors, Molecule, Variables, Bonds, NumberofAtoms)
      call NonBondedEnergy(EnNonBond, Molecule, Variables, NumberofAtoms)
      
      TotEner = EnStretch + EnBend + EnNonBond

      deallocate(Bonds)
!      print *, 'TotEner =', TotEner, 'EnStretch =', EnStretch, 'EnBend =', EnBend, 'EnTors =', EnTors, 'EnNonBond =', EnNonBond
   endsubroutine


endmodule
