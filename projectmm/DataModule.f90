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
                           AssigningBonds, Binding, Planes, AddToList_Angle,       &
                           AddToList_Plane

   type Atom     
      character(1)      :: element
      real              :: x, y, z
   endtype
   
   type Parameters
      real              :: TEMP, PRESSURE, ForceConstantCC, ForceConstantCH,            &
                           EquiBondCC, EquiBondCH, ForceConstantAngle,                  & 
                           EquiAngle, gama, n, V1, ChargeH, ChargeC, WellDepthC,        &
                           WellDepthH, RstarC, RstarH, Boltzmann, ForceConstantAngleCCC
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
       
      open(32,file='c4_cyclic.xyz')
      read(32,*) NumberofAtoms 
      allocate(Molecule(NumberofAtoms))
      do i = 1,NumberofAtoms 
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
      read(33,*) Variables%ForceConstantCC
      read(33,*) Variables%ForceConstantCH
      read(33,*) Variables%EquiBondCC
      read(33,*) Variables%EquiBondCH
      read(33,*) Variables%ForceConstantAngle
      read(33,*) Variables%ForceConstantAngleCCC
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
      type (Binding), intent(inout), allocatable        :: Bonds(:)
      integer, intent(in)                               :: NumberofAtoms
      type (Atom), intent(inout), allocatable           :: Molecule(:)
      type (Parameters), intent(in)                     :: Variables
      integer                                           :: atom1, atom2, signal   
      integer, allocatable                              :: Check(:,:)
      type (Binding)                                    :: Bond
      
      allocate(Check(NumberofAtoms,NumberofAtoms))
      Check = 0

      do atom1 = 1,NumberofAtoms
         do atom2 = 1,NumberofAtoms
            if (atom1 /= atom2) then     
               if (Molecule(atom1)%element == 'C' .and. Molecule(atom2)%element == 'C' .and. abs(BondLength(atom1, atom2,&
               Molecule)-Variables%EquiBondCC) < 0.2 .and. Check(atom1,atom2) == 0) then
                  
                  Bond%length = BondLength(atom1, atom2, Molecule)
                  Bond%elements = 'CC'
                  Check(atom1,atom2) = 1
                  Check(atom2,atom1) = 1
                  Bond%FirstAtom = atom1
                  Bond%SecondAtom = atom2
                  call AddToList_Bonds(Bonds, Bond)
               elseif (Molecule(atom1)%element == 'C' .and. Molecule(atom2)%element == 'H' .and. abs(BondLength(atom1, atom2,&
               Molecule)-Variables%EquiBondCH) < 0.2) then
                  
                  Bond%length = BondLength(atom1, atom2, Molecule)
                  Bond%elements = 'CH'
                  Bond%FirstAtom = atom1
                  Bond%SecondAtom = atom2
                  call AddToList_Bonds(Bonds, Bond)
               endif
            endif
         enddo
      enddo
   
   endsubroutine
   
   subroutine AddToList_Bonds(List, Element)
      implicit none
      type (Binding), intent(inout), allocatable        :: List(:)
      type (Binding), intent(in)                        :: Element
      type (Binding), allocatable                       :: ListDummy(:)
      integer                                           :: i, isize


      
      if (allocated(List)) then
         isize = size(List)
         allocate(ListDummy(isize+1))
         do i=1,isize
            ListDummy(i) = List(i)
         enddo
         ListDummy(isize+1) = Element

         deallocate(List)
         call move_alloc(ListDummy, List)
      else
         allocate(List(1))
         List(1) = Element
      endif

   endsubroutine

   subroutine AddToList_Angle(List, Element)
      real, intent(inout), allocatable  :: List(:) 
      real, intent(in)                  :: Element
      real, allocatable                 :: ListDummy(:)
      integer                           :: i, isize

      if(allocated(List)) then
         isize = size(List)
         allocate(ListDummy(isize+1))
         do i=1,isize  
             ListDummy(i) = List(i)
         end do
         ListDummy(isize+1) = Element

         deallocate(List)
         call move_alloc(ListDummy, List)
       else
          allocate(List(1))            ! When you list is not allocated yet allocate it with size 1.
          List(1) = Element
       end if

   endsubroutine


   subroutine AddToList_Plane(List, Element)
      type (Planes), intent(inout), allocatable :: List(:)
      type (Planes), intent(in)                 :: Element
      type (Planes), allocatable                :: ListDummy(:)
      integer                                   :: i, isize

      if(allocated(List)) then
         isize = size(List)
         allocate(ListDummy(isize+1))
         do i=1,isize
            ListDummy(i) = List(i)
         end do
         ListDummy(isize+1) = Element
 
         deallocate(List)
         call move_alloc(ListDummy, List)
      else
         allocate(List(1))
         List(1) = Element
      endif
   
   endsubroutine

endmodule
