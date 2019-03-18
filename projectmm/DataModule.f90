!!!____________________________________________________________________________
!!! Author: Aaron Eisses
!!! Studentnumber: 10593527
!!! University: University of Amsterdam
!!! Date changed: 11-02-2019
!!!____________________________________________________________________________
!!!
!!! PROJECT SCIENTIFIC COMPUTING AND PROGRAMMING: MOLECULAR MECHANICS
!!!
!!! DataModule.f90
!!!
!!! This file is the module energy. Here, a type is created called Atom and a
!!! subroutine called ReadData. The type is there to store the element of each
!!! particle and its coordinates (x,y,x). The subroutine is used to store all 
!!! of the information about the system and parameters from an input file.
!!!____________________________________________________________________________

module DataModule
   implicit none
   save
   private
   public               :: ReadData, ReadParameters, Atom, BondLength, Parameters

   type Atom     
      character(1)      :: element
      real              :: x, y, z
   endtype
   
   type Parameters
      real              :: TEMP, PRESSURE, ForceConstantCC, ForceConstantCH,   &
                           EquiBondCC, EquiBondCH, ForceConstantAngle,         &
                           EquiAngle, gama, n, V1, ChargeH, ChargeC, Aij, Bij, &
                           Boltzmann
   endtype


contains

   subroutine ReadData(Molecule)
      type (Atom), intent(inout)   :: Molecule(8)
      integer                      :: i, waste
      
      !!! First, all of the atoms are put into the array molecule(:) with
      !!! their coordinates
      

      open(32,file='ethan_distorted.xyz')
      read(32,*) 
      do i = 1,8 
         read(32,*) waste, Molecule(i)%element, Molecule(i)%x, Molecule(i)%y, Molecule(i)%z
      enddo
      close(32) 
     
   endsubroutine



   subroutine ReadParameters(Variables)
      type (Parameters), intent(inout) :: Variables 

      open(33,file='parameter.txt')
      read(33,*) Variables%TEMP
      read(33,*) Variables%PRESSURE
      read(33,*) Variables%ForceConstantCC
      read(33,*) Variables%ForceConstantCH
      read(33,*) Variables%EquiBondCC
      read(33,*) Variables%EquiBondCH
      read(33,*) Variables%ForceConstantAngle
      read(33,*) Variables%EquiAngle
      read(33,*) Variables%gama
      read(33,*) Variables%n
      read(33,*) Variables%V1
      read(33,*) Variables%ChargeH
      read(33,*) Variables%ChargeC
      read(33,*) Variables%Aij
      read(33,*) Variables%Bij
      Variables%Boltzmann = 1.38064852
      Variables%Boltzmann = Variables%Boltzmann/(10**23)
      Variables%ChargeH = Variables%ChargeH*(1.60217653*(10**(-19)))
      Variables%ChargeC = Variables%ChargeC*(1.60217653*(10**(-19)))
      close(33)


   endsubroutine
 
   real function BondLength(atom1, atom2, Molecule)
      integer, intent(in)       :: atom1, atom2
      type (Atom), intent(inout)   :: Molecule(8)

      BondLength = sqrt((Molecule(atom1)%x - Molecule(atom2)%x)**2 + (Molecule(atom1)%y - Molecule(atom2)%y)**2 &
              + (Molecule(atom1)%z - Molecule(atom2)%z)**2) 
   endfunction

   real function RandomNumber()
       
      call random_number(RandomNumber)
  
   endfunction



endmodule
