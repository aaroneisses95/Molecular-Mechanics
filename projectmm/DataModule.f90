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
   public               :: ReadData, ReadParameters, Atom, BondLength, Parameters, &
                           AssigningBonds, Binding, Planes

   type Atom     
      character(1)      :: element
      real              :: x, y, z
   endtype
   
   type Parameters
      real              :: TEMP, PRESSURE, ForceConstantCC, ForceConstantCH,    &
                           EquiBondCC, EquiBondCH, ForceConstantAngle,          &
                           EquiAngle, gama, n, V1, ChargeH, ChargeC, WellDepthC,&
                           WellDepthH, RstarC, RstarH, Boltzmann
   endtype

   type Binding
      character(2)      :: elements
      real              :: length
      integer           :: FirstAtom, SecondAtom
   endtype

   type Planes
      real              :: a, b, c
   endtype
      
contains

   !!! In the subroutine ReadData, the coordinates and the type of all the atoms
   !!! are stored in an array of type 'Atom' from the input file. The variable 
   !!! 'waste' is there to get rid of the first word in the input file in each line.    
   
   subroutine ReadData(Molecule, NumberofAtoms)
      type (Atom), intent(inout), allocatable           :: Molecule(:)
      integer, intent(inout)                            :: NumberofAtoms    
      integer                                           :: i, waste
      
      open(32,file='ethan_distorted.xyz')
      read(32,*) NumberofAtoms 
      allocate(Molecule(NumberofAtoms))
      do i = 1,8 
         read(32,*) waste, Molecule(i)%element, Molecule(i)%x, Molecule(i)%y, Molecule(i)%z
      enddo
      close(32) 
     
   endsubroutine

   !!! In the subroutine ReadParameters, all the parameters of the system are stored
   !!! into 'Variables' with type 'Parameters', by reading them all from an input file

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
      read(33,*) Variables%WellDepthC
      read(33,*) Variables%WellDepthH
      read(33,*) Variables%RstarC
      read(33,*) Variables%RstarH
      Variables%Boltzmann = 1.38064852
      Variables%Boltzmann = Variables%Boltzmann/(10**23)
      Variables%ChargeH = Variables%ChargeH*1.60217653/(10**19)
      Variables%ChargeC = Variables%ChargeC*1.60217653/(10**19)
      close(33)

   endsubroutine
   
   !!! In the function BondLength, the distance is calculated between two atoms, whether
   !!! they are bonded or not. 

   real function BondLength(atom1, atom2, Molecule)
      integer, intent(in)       :: atom1, atom2
      type (Atom), intent(inout), allocatable   :: Molecule(:)

      BondLength = sqrt((Molecule(atom1)%x - Molecule(atom2)%x)**2 + (Molecule(atom1)%y - Molecule(atom2)%y)**2 &
              + (Molecule(atom1)%z - Molecule(atom2)%z)**2) 
   
   endfunction
   
   subroutine AssigningBonds(Bonds, NumberofAtoms, Molecule, Variables)
      type (Binding), intent(inout)             :: Bonds(:)
      integer, intent(in)                       :: NumberofAtoms
      type (Atom), intent(inout), allocatable   :: Molecule(:)
      type (Parameters), intent(in)             :: Variables
      integer                                   :: atom1, atom2, k   
      integer, allocatable                      :: Check(:,:)

      allocate(Check(NumberofAtoms,NumberofAtoms))
    !  allocate(Bonds(NumberofAtoms-1))

      k = 0
      Check = 0

      do atom1 = 1,NumberofAtoms
         do atom2 = 1,NumberofAtoms
            if (atom1 /= atom2) then     
               if (Molecule(atom1)%element == 'C' .and. Molecule(atom2)%element == 'C' .and. abs(BondLength(atom1, atom2,&
               Molecule)-Variables%EquiBondCC) < 0.2 .and. Check(atom1,atom2) == 0) then
                  k = k + 1
                  Bonds(k)%length = BondLength(atom1, atom2, Molecule)
                  Bonds(k)%elements = 'CC'
                  Check(atom1,atom2) = 1
                  Check(atom2,atom1) = 1
                  Bonds(k)%FirstAtom = atom1
                  Bonds(k)%SecondAtom = atom2
               elseif (Molecule(atom1)%element == 'C' .and. Molecule(atom2)%element == 'H' .and. abs(BondLength(atom1, atom2,&
               Molecule)-Variables%EquiBondCH) < 0.2) then
                  k = k + 1
                  Bonds(k)%length = BondLength(atom1, atom2, Molecule)
                  Bonds(k)%elements = 'CH'
                  Bonds(k)%FirstAtom = atom1
                  Bonds(k)%SecondAtom = atom2
               endif
            endif
         enddo
      enddo
       
       
   endsubroutine
   

endmodule
