!!!________________________________________________________________________________________________
!!! Author: Aaron Eisses
!!! Studentnumber: 10593527
!!! University: University of Amsterdam
!!! Date changed: 31-03-2019
!!!________________________________________________________________________________________________
!!!
!!! PROJECT SCIENTIFIC COMPUTING AND PROGRAMMING: MOLECULAR MECHANICS
!!!
!!! DataModule.f90
!!!
!!! This file is the module Data. Here, some types and subroutines are created that are being
!!! used throughout the whole program. The types are 'Atom' and 'Parameters' and the subroutines 
!!! are 'ReadData' and 'ReadParameters'. The type 'Atom' is there to store the element of each 
!!! particle and its coordinates (x,y,x). The type 'Parameters' is there to store all the parameters 
!!! that are being used. The subroutine 'ReadData' is used to store all of the information about the 
!!! atoms of the system in an array of type 'Atom', where the length of the array is determined by 
!!! the amount of atoms. The subroutine 'ReadParameters' reads the parameters from an input file and 
!!! stores the values in 'Variables', which is of type 'Parameters'.  
!!!________________________________________________________________________________________________

module DataModule
   implicit none
   save
   private
   public               :: Atom, Parameters, ReadData, ReadParameters                  

   type Atom     
      character(1)      :: element
      real              :: x, y, z
   endtype
   
   type Parameters
      real              :: TEMP, ForceConstantCC, ForceConstantCH, EquiBondCC, EquiBondCH,            & 
                           ForceConstantAngle, EquiAngle, gama, n, V1, ChargeH, ChargeC, WellDepthC,  &
                           WellDepthH, RstarC, RstarH, Boltzmann, ForceConstantAngleCCC, maxcycle, r, &
                           Kelec, Avogadro
   endtype

contains
 
   subroutine ReadData(Molecule, NumberofAtoms)
      type (Atom), intent(inout), allocatable           :: Molecule(:)
      integer, intent(inout)                            :: NumberofAtoms    
      integer                                           :: i, waste
       
      open(32,file='ethan_distorted.xyz')
      read(32,*) NumberofAtoms 
      allocate(Molecule(NumberofAtoms))
      do i = 1,NumberofAtoms 
         read(32,*) waste, Molecule(i)%element, Molecule(i)%x, Molecule(i)%y, Molecule(i)%z
      enddo
      close(32) 

   endsubroutine

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
      read(33,*) Variables%maxcycle
      read(33,*) Variables%r
      read(33,*) Variables%Boltzmann
      read(33,*) Variables%Kelec
      read(33,*) Variables%Avogadro
      Variables%Avogadro = Variables%Avogadro*(10**23)
      Variables%Boltzmann = Variables%Boltzmann/(10**23)
      Variables%ChargeH = Variables%ChargeH*1.60217653/(10**19)
      Variables%ChargeC = Variables%ChargeC*1.60217653/(10**19)
      close(33)

   endsubroutine
   
endmodule
