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

   subroutine StretchEnergy(EnStretch, Molecule, Variables, Bond)
      real, intent(inout)               :: EnStretch 
      type (Atom), intent(inout)        :: Molecule(8)
      type (Parameters), intent(in)     :: Variables
      real, intent(inout)               :: Bond(7)
      integer                           :: i

      EnStretch = 0   

      !!! It must loop over all the bonds in the system and save whether C-C or
      !!! C-H and save their Cartesian coordinates and calculate the distance
      !!! between them. Bond(1) up until Bond(6) is c-h and bond(7) is c-c 

      Bond(1) = BondLength(2,3, Molecule)
      Bond(2) = BondLength(2,4, Molecule)
      Bond(3) = BondLength(2,5, Molecule)   
      Bond(4) = BondLength(1,6, Molecule)
      Bond(5) = BondLength(1,7, Molecule)
      Bond(6) = BondLength(1,8, Molecule)
      Bond(7) = BondLength(1,2, Molecule)

      !!! Calculate the energy of bonds with E=K(r-r0)^2. The first loop is for
      !!! the C-H bonds. The other line is for the only C-C bond in ethane.
    
      do i=1,6
         EnStretch = EnStretch + Variables%ForceConstantCH*((Bond(i) -         &
         Variables%EquiBondCH)**2)
      enddo

      EnStretch = EnStretch + Variables%ForceConstantCC*((Bond(7) -            &
      Variables%EquiBondCC)**2) 


   endsubroutine

   subroutine BendEnergy(EnBend, Molecule, Variables, Bond)
      real, intent(inout)               :: EnBend
      type (Atom), intent(inout)        :: Molecule(8)
      type (Parameters), intent(in)     :: Variables
      real, intent(in)                  :: Bond(7) 
      real                              :: Angle(12)
      integer                           :: i

      EnBend = 0
      
      !!! The angles are calculated between the atoms of interest. 

      Angle(1) = (BondLength(6,7, Molecule)**2 - Bond(4)**2 - Bond(5)**2) / (2*Bond(4)*Bond(5))
      
      Angle(2) = (BondLength(6,8, Molecule)**2 - Bond(4)**2 - Bond(6)**2) / (2*Bond(4)*Bond(6))
      
      Angle(3) = (BondLength(7,8, Molecule)**2 - Bond(5)**2 - Bond(6)**2) / (2*Bond(5)*Bond(6))
      
      Angle(4) = (BondLength(3,4, Molecule)**2 - Bond(1)**2 - Bond(2)**2) / (2*Bond(1)*Bond(2))
      
      Angle(5) = (BondLength(3,5, Molecule)**2 - Bond(1)**2 - Bond(3)**2) / (2*Bond(1)*Bond(3))      
      
      Angle(6) = (BondLength(4,5, Molecule)**2 - Bond(2)**2 - Bond(3)**2) / (2*Bond(2)*Bond(3))
      
      Angle(7) = (BondLength(2,6, Molecule)**2 - Bond(4)**2 - Bond(7)**2) / (2*Bond(4)*Bond(7))
      
      Angle(8) = (BondLength(2,7, Molecule)**2 - Bond(5)**2 - Bond(7)**2) / (2*Bond(5)*Bond(7))
      
      Angle(9) = (BondLength(2,8, Molecule)**2 - Bond(6)**2 - Bond(7)**2) / (2*Bond(6)*Bond(7))
      
      Angle(10) = (BondLength(1,3, Molecule)**2 - Bond(1)**2 - Bond(7)**2) / (2*Bond(1)*Bond(7))
      
      Angle(11) = (BondLength(1,4, Molecule)**2 - Bond(2)**2 - Bond(7)**2) / (2*Bond(2)*Bond(7))
      
      Angle(12) = (BondLength(1,5, Molecule)**2 - Bond(3)**2 - Bond(7)**2) / (2*Bond(3)*Bond(7))

      !!! The bending energy is calculated 

      do i = 1, 12
         EnBend = EnBend + Variables%ForceConstantAngle*((Angle(i) - Variables%EquiAngle)**2)
      enddo
   
   endsubroutine

 !  subroutine TorsionalEnergy
 
 !  endsubroutine
  
   subroutine NonBondedEnergy(EnNonBond, Molecule, Variables)
      real, intent(inout)               :: EnNonBond
      type (Atom), intent(inout)        :: Molecule(8)
      type (Parameters), intent(in)     :: Variables
      integer                           :: i, j
      real(8)                           :: Kelec, Rij
      
      EnNonBond = 0
      Kelec = 8.987551787*(10**9) ! Constant in Law of Coulomb
      
      !!! In this double loop, the non bonded energy is calculated between all the atoms
      !!! with both the electronic effect and the Van Der Waals effect. 

      do i = 1, 8
         do j = 1, 8
            if (i /= j) then
               Rij = BondLength(i,j, Molecule)
               EnNonBond = EnNonBond + ((Variables%Aij / (Rij**12)) - (Variables%Bij / (Rij**6)))
               if (Molecule(i)%element == 'C' .and. Molecule(j)%element == 'C') then
                  EnNonBond = EnNonBond + ((Variables%ChargeC**2) * Kelec / Rij)
               elseif (Molecule(i)%element == 'H' .and. Molecule(j)%element == 'H') then
                  EnNonBond = EnNonBond + ((Variables%ChargeH**2) * Kelec / Rij)
               else
                  EnNonBond = EnNonBond + ((Variables%ChargeH*Variables%ChargeC) * Kelec / Rij)
               endif
            endif
         enddo
      enddo

   endsubroutine 
   
   subroutine TotalEnergy(TotEner, Molecule, Variables)
      real, intent(inout)               :: TotEner
      type (Atom), intent(inout)        :: Molecule(8)               
      type (Parameters), intent(in)     :: Variables
      real                              :: EnStretch, EnBend, Bond(7), EnNonBond
      
      call StretchEnergy(EnStretch, Molecule, Variables, Bond)
      call BendEnergy(EnBend, Molecule, Variables, Bond)
!      call TorsionalEnergy(EnTors, Molecule, Variables)
      call NonBondedEnergy(EnNonBond, Molecule, Variables)

      TotEner = EnStretch + EnBend + EnNonBond
   endsubroutine

endmodule
