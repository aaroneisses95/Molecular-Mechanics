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

   type (Atom)                   :: MoleculeInit(8)
   type (Parameters)             :: Variables
   integer                       :: ncycle, i
   real                          :: CurrentEnergy, r
   logical                       :: Check
   
   print *, '-------------------Molecular Mechanics--------------------'

   !!! Initialize system by reading all the date from input
   
   Call ReadParameters(Variables)
   Call ReadData(MoleculeInit)
   Call TotalEnergy(CurrentEnergy, MoleculeInit, Variables)
   
   !!! Start the Metropolis algorithm
   
   ncycle = 0
   r = 0.0001

   do while (ncycle < 5000)
      Call Move(CurrentEnergy, MoleculeInit, Variables, Check, r)
      if (Check .eqv. .true.) then
         ncycle = 0
      elseif (Check .eqv. .false.) then
         ncycle = ncycle + 1
      endif
   enddo

   !!! Print everything that is useful to know (energy, amount of cycles etc.)
   
   print *, 'The lowest Energy =', CurrentEnergy

endprogram 
