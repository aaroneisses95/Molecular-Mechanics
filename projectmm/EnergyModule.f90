!!!____________________________________________________________________________
!!! Author: Aaron Eisses
!!! Studentnumber: 10593527
!!! University: University of Amsterdam
!!! Date changed: 31-03-2019
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
use CalculationModule

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
     
      !!! Calculate the energy of bonds with E=K(r-r0)^2. It loops over all the bonds and the if-statement
      !!! is there to separate the C-C and C-H bonds to make sure that the proper constants are used.

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

      !!! The angles are calculated between the atoms of interest. In the first if-statement, only 
      !!! the C atoms are selected in the molecule because H-atoms can only have one bond here, and 
      !!! therefore they can't be the connecting atom between two bonds. Then, two bonds are taken 
      !!! and there is being checked if they are connected by the same atom but with another atom at 
      !!! the other side of the bond. If so, an angle is created with the corresponding angle between 
      !!! the two bonds and added to the list. 

      do i = 1,NumberofAtoms
         if (Molecule(i)%element == 'C') then
            do j = 1,size(Bonds)
               do k = j,size(Bonds)
                  if (j /= k .and. Bonds(j)%FirstAtom == i .and. Bonds(k)%FirstAtom == i) then   
                     Angle = -((BondLength(Bonds(j)%SecondAtom, Bonds(k)%SecondAtom, Molecule))**2 - &
                             (Bonds(j)%length)**2 - (Bonds(k)%length)**2) / (2*(Bonds(j)%length)*   &
                             (Bonds(k)%length))
                     call AddToList_Angle(Angles, Angle)
                  endif
                  if (j /= k .and. Bonds(j)%SecondAtom == i .and. Bonds(k)%FirstAtom == i) then
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
      integer                                   :: i, j, k, m, n
      integer, intent(in)                       :: NumberofAtoms
      real                                      :: VectorX1, VectorY1, VectorZ1, VectorX2, VectorY2, &
                                                   VectorZ2, phi
      type (Planes), allocatable                :: PlanesList(:)
      type (Planes)                             :: Plane

      EnTors = 0
        
      !!! The different planes in the molecule are figured out, since every bond describes a plane 
      !!! of the molecule. The formulas of the planes are calculated by using the coordinates of 
      !!! three points on that plane.  
      
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
 
                     Plane%a = (VectorY1*VectorZ2)-(VectorZ1*VectorY2)
                     Plane%b = (VectorX1*VectorZ2)-(VectorZ1*VectorX2)
                     Plane%c = (VectorX1*VectorY2)-(VectorY1*VectorX2)

                     call AddToList_Plane(PlanesList, Plane)
                  endif
               enddo
            enddo
         endif
      enddo
      
      !!! There is being looped over two different planes and the angle between the planes is calculated. 
      !!! With this angle, the torsional energy can be calculated.
      
      do m = 1, size(PlanesList)
         do n = m, size(PlanesList)
            if (n /= m) then
               phi = (abs((PlanesList(m)%a*PlanesList(n)%a) + (PlanesList(m)%b*PlanesList(n)%b) + &
                     (PlanesList(m)%c*PlanesList(n)%c)) / (sqrt(PlanesList(m)%a**2 + PlanesList(m)%b**2 &
                     + PlanesList(m)%c**2) * sqrt(PlanesList(m)%a**2 + PlanesList(m)%b**2 + PlanesList(m)%c**2)))
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
      real                                      :: Rij
      integer, intent(in)                       :: NumberofAtoms
      
      EnNonBond = 0
      
      !!! The non-bonded energy is calculated between all the atoms, considering the electronic effect and the Van Der 
      !!! Waals effect. There are three if-statements to determine what kind of atoms interact with each other (C-C, C-H, 
      !!! or H-H).
 
      do i = 1, NumberofAtoms
         do j = i, NumberofAtoms
            if (i /= j) then
               Rij = BondLength(i, j, Molecule)
               if (Molecule(i)%element == 'C' .and. Molecule(j)%element == 'C') then

                  !!! This if-statement makes sure that the electronic effect between two atoms is only counted once.
                  !!! However, the Van der Waals interactions should be counted double since they can be different in an
                  !!! atom pair. 

                  if (i < j) then
                     EnNonBond = EnNonBond + Variables%Avogadro*((Variables%ChargeC**2)*Variables%Kelec/(Rij/(10**10))) 
                  endif                
                  EnNonBond = EnNonBond + Variables%WellDepthC*(((Variables%RstarC/Rij)**12)-2*((Variables%RstarC/Rij)**6))
               
               elseif (Molecule(i)%element == 'H' .and. Molecule(j)%element == 'H') then
                  if (i < j) then      
                     EnNonBond = EnNonBond + Variables%Avogadro*((Variables%ChargeH**2)*Variables%Kelec/(Rij/(10**10))) 
                  endif
                  EnNonBond = EnNonBond + Variables%WellDepthH*(((Variables%RstarH/Rij)**12)-2*((Variables%RstarH/Rij)**6)) 
               
               else
                  if (i < j) then    
                     EnNonBond = EnNonBond + Variables%Avogadro*((Variables%ChargeH*Variables%ChargeC)*Variables%Kelec/ &
                                 (Rij/(10**10)))
                  endif
                  EnNonBond = EnNonBond + sqrt(Variables%WellDepthH*Variables%WellDepthC) * &
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

   endsubroutine

endmodule
