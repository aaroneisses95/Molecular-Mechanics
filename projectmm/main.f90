!!!____________________________________________________________________________
!!! Author: Aaron Eisses
!!! Studentnumber: 10593527
!!! University: University of Amsterdam
!!! Date changed: 31-03-2019
!!!____________________________________________________________________________
!!!
!!! PROJECT SCIENTIFIC COMPUTING AND PROGRAMMING: MOLECULAR MECHANICS
!!!
!!! main.f90
!!!
!!! In this project, code will be written for the parametrization of an energy
!!! function from a given molecular structure, consisting of only hydrogen and
!!! carbon atoms, and is minimized. This is the main program called main.f90
!!! , and is combined with three module files called atom.f90, minimization.f90
!!! and energy.f90.
!!!____________________________________________________________________________

program main
use DataModule
use MinimizationModule
use EnergyModule
use CalculationModule
   implicit none

   type (Atom), allocatable     :: MoleculeInit(:), MoleculeRef(:)
   type (Parameters)            :: Variables
   type (Binding), allocatable  :: Bonds(:)
   integer                      :: i, NumberofAtoms
   real                         :: CurrentEnergy, EnergyRef
   logical                      :: Check
   
   print *, 
   print *, '-------------------Molecular Mechanics--------------------'

   !!! Initialize system by reading all the date from input
   
   Call ReadParameters(Variables)
   Call ReadData(MoleculeInit, NumberofAtoms)

   allocate(MoleculeRef(NumberofAtoms))
   
   Call TotalEnergy(CurrentEnergy, MoleculeInit, Variables, Bonds, NumberofAtoms)
 
   !!! Start the Metropolis algorithm
   
   EnergyRef = CurrentEnergy
   MoleculeRef = MoleculeInit

   call Minimization(CurrentEnergy, MoleculeInit, Variables, Check, NumberofAtoms, Bonds)

   !!! Print everything that is useful to know (energy, amount of cycles etc.)

   print *, 'Initial coordinates of atoms'
   print *, 'Number of Atom ', 'Type of Element ', 'X-Coordinate ', 'Y-Coordinate ', 'Z-Coordinate'
   
   do i = 1,NumberofAtoms
      print *, i, MoleculeRef(i)%element, MoleculeRef(i)%x, MoleculeRef(i)%y, MoleculeRef(i)%z
   enddo
   
   print *, 'Initial Energy (kJ/mol) = ', Units_kJ(EnergyRef)
   print *, 'Initial Energy (kcal/mol) =', Units_kcal(EnergyRef)   
   print *, '-----------------------------------------------------------'
   print *, 'End coordinates of atoms' 
   print *, 'Number of Atom ', 'Type of Element ', 'X-Coordinate ', 'Y-Coordinate ', 'Z-Coordinate'
   
   do i = 1,NumberofAtoms
      print *, i, MoleculeInit(i)%element, MoleculeInit(i)%x, MoleculeInit(i)%y, MoleculeInit(i)%z
   enddo
   
   print *, 'End Energy (kJ/mol) =', Units_kJ(CurrentEnergy)
   print *, 'End Energy (kcal/mol) =', Units_kcal(CurrentEnergy)
   print *,

   deallocate(MoleculeRef)

endprogram 
