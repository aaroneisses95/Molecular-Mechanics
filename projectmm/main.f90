!!!____________________________________________________________________________
!!! Author: Aaron Eisses
!!! Studentnumber: 10593527
!!! University: University of Amsterdam
!!! Date changed: 11-02-2019
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
!!!
!!!____________________________________________________________________________

program main
use DataModule
use MinimizationModule
use EnergyModule
   implicit none

   type (Atom), allocatable     :: MoleculeInit(:), MoleculeRef(:)
   type (Parameters)            :: Variables
   type (Binding), allocatable  :: Bonds(:)
   integer                      :: ncycle, i, NumberofAtoms
   real                         :: CurrentEnergy, r, EnergyRef
   logical                      :: Check
   
   print *, '-------------------Molecular Mechanics--------------------'

   !!! Initialize system by reading all the date from input
   
   Call ReadParameters(Variables)
   Call ReadData(MoleculeInit, NumberofAtoms)

   allocate(MoleculeRef(NumberofAtoms))
   
   Call TotalEnergy(CurrentEnergy, MoleculeInit, Variables, Bonds, NumberofAtoms)

   
   !!! Start the Metropolis algorithm
   EnergyRef = CurrentEnergy
   MoleculeRef = MoleculeInit
   ncycle = 0
   r = 0.0001

   !!! Here, the move subroutine 'Move' is called a number of times until the exit
   !!! strategy is reached. If ncycle reaches 5000, it means that the subroutine
   !!! has been called 5000 consecutive times with as a result unsuccesful changing
   !!! of the configuration of the molecule 

   do while (ncycle < 5000)
      Call Move(CurrentEnergy, MoleculeInit, Variables, Check, r, NumberofAtoms, Bonds)
      if (Check .eqv. .true.) then
         ncycle = 0
      elseif (Check .eqv. .false.) then
         ncycle = ncycle + 1
      endif
   enddo

   !!! Print everything that is useful to know (energy, amount of cycles etc.)

   print *, 'Initial coordinates of atoms'
   print *, 'Number of Atom ', 'Type of Element ', 'X-Coordinate ', 'Y-Coordinate ', 'Z-Coordinate'
   do i = 1,NumberofAtoms
      print *, i, MoleculeRef(i)%element, MoleculeRef(i)%x, MoleculeRef(i)%y, MoleculeRef(i)%z
   enddo
   print *, 'Initial Energy (kJ/mol) = ', EnergyRef/1000
   print *, 'Initial Energy (kcal/mol) =', EnergyRef/4184   
   print *, '-----------------------------------------------------------'
   print *, 'End coordinates of atoms' 
   print *, 'Number of Atom ', 'Type of Element ', 'X-Coordinate ', 'Y-Coordinate ', 'Z-Coordinate'
   do i = 1,NumberofAtoms
      print *, i, MoleculeInit(i)%element, MoleculeInit(i)%x, MoleculeInit(i)%y, MoleculeInit(i)%z
   enddo
   print *, 'End Energy (kJ/mol) =', CurrentEnergy/1000
   print *, 'End Energy (kcal/mol) =', CurrentEnergy/4184
   deallocate(MoleculeRef)

endprogram 
