!!!_____________________________________________________________________________________________________
!!! Author: Aaron Eisses
!!! Studentnumber: 10593527
!!! University: University of Amsterdam
!!! Date changed: 31-03-2019
!!!_____________________________________________________________________________________________________
!!!
!!! PROJECT SCIENTIFIC COMPUTING AND PROGRAMMING: MOLECULAR MECHANICS
!!!
!!! Calculations.f90
!!!
!!! This file is the module Calculation. Here, multiple types and subroutines are created that are being
!!! used throughout the whole program. The types are 'Binding' and 'Planes' and the subroutines are 
!!! 'BondLength', 'AssigningBonds','AddToList_Bonds', 'AddToList_Angle' and 'AddToList_Plane'. 
!!! The type 'Binding' stores the elements, length and the number of the atoms between a bond. The 
!!! type plane stores the variables 'a', 'b' and 'c', which will be used to calculate the angle 
!!! between two planes in the molecule. The subroutine 'BondLength' calculates the distance between 
!!! two atoms. The subroutine 'AssigningBonds' iterates over all the atoms and determines which atoms
!!! are bonded and which ones are not. If the calculated distance between two atoms is similar to the 
!!! equilibrium bondlength, they are considered to be a bond. The three subroutines 'AddToList_Bonds',
!!! 'AddToList_Angle' and 'AddToList_Plane' do the same thing but for different data types. It increases 
!!! a linked list of that type with one element.
!!! At the end of the module, there are also two functions called 'Units_kJ' and 'Units_kcal'. These are
!!! used to transfer units from J to kJ and kcal
!!!_____________________________________________________________________________________________________

module CalculationModule
use DataModule

   implicit none
   save
   private
   public               :: Binding, Planes, BondLength, AssigningBonds, AddToList_Bonds, AddToList_Angle, &
                           AddToList_Plane, Units_kJ, Units_kcal

   type Binding
      character(2)      :: elements
      real              :: length
      integer           :: FirstAtom, SecondAtom
   endtype

   type Planes
      real              :: a, b, c
   endtype
      
contains

   real function BondLength(atom1, atom2, Molecule)
      integer, intent(in)                       :: atom1, atom2
      type (Atom), intent(inout), allocatable   :: Molecule(:)
      real                                      :: Deltax, Deltay, Deltaz
      
      Deltax = Molecule(atom1)%x - Molecule(atom2)%x
      Deltay = Molecule(atom1)%y - Molecule(atom2)%y
      Deltaz = Molecule(atom1)%z - Molecule(atom2)%z

      BondLength = sqrt(Deltax**2 + Deltay**2 + Deltaz**2)
 
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

      !!! It iterates over all atoms. In the if statements is being check if it is a C-C bond or an 
      !!! C-H bond and whether it can be considered as a real bond. The Check-statement is in there 
      !!! to make sure that if there is already a bond between two atoms, it does not make another 
      !!! bond when they are checked the other way around. In the end, the new bond is added to the 
      !!! array called 'Bonds'. 

      do atom1 = 1,NumberofAtoms
         do atom2 = 1,NumberofAtoms
            if (atom1 /= atom2) then     
               if (Molecule(atom1)%element == 'C' .and. Molecule(atom2)%element == 'C' .and. &
                   abs(BondLength(atom1, atom2, Molecule)-Variables%EquiBondCC) < 0.2  .and. &
                   Check(atom1,atom2) == 0) then
                  
                  Bond%length = BondLength(atom1, atom2, Molecule)
                  Bond%elements = 'CC'
                  Check(atom1,atom2) = 1
                  Check(atom2,atom1) = 1
                  Bond%FirstAtom = atom1
                  Bond%SecondAtom = atom2
                  call AddToList_Bonds(Bonds, Bond)

               elseif (Molecule(atom1)%element == 'C' .and. Molecule(atom2)%element == 'H' .and. &
                       abs(BondLength(atom1, atom2, Molecule)-Variables%EquiBondCH) < 0.2) then
                  
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

      !!! It checks if the input array is already allocated. If so, the list is copied in a dummy 
      !!! variable which has an extra memory space which is then filled up with the new element. 
      !!! If not, the list is initiated and allocated. 
      
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
          allocate(List(1))
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

   real function Units_kJ(Energy) result(NewEnergy)
      real, intent(inout)       :: Energy

      NewEnergy = Energy/1000

   endfunction

   real function Units_kcal(Energy) result(NewEnergy)
      real, intent(inout)       :: Energy

      NewEnergy = Energy/4184

   endfunction

endmodule
