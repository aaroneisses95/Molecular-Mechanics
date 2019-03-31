!!!_____________________________________________________________________________________________________
!!! Author: Aaron Eisses
!!! Studentnumber: 10593527
!!! University: University of Amsterdam
!!! Date changed: 31-03-2019
!!!_____________________________________________________________________________________________________
!!!
!!! PROJECT SCIENTIFIC COMPUTING AND PROGRAMMING: MOLECULAR MECHANICS
!!!
!!! MinimizationModule.f90
!!!
!!! This file is the module Minimization. Here, there are two subroutines. The 'Move' subroutine does 
!!! an energy calculation of the current configuration, moves an atom with a random displacement and  
!!! calculates the energy of the new configuration. With a certain probability function, the new 
!!! configuration gets either accepted or rejected. The second subroutine, 'Minimization', calls the
!!! subroutine 'Move' an amount of times until a minimum energy value is reached. 
!!!_____________________________________________________________________________________________________

module MinimizationModule
use DataModule
use EnergyModule
use CalculationModule
   
   implicit none
   save 
   private
   public               :: Move, Minimization

contains

   subroutine Move(CurrentEnergy, MoleculeInit, Variables, Check, NumberofAtoms, Bonds)
      real, intent(inout)                               :: CurrentEnergy
      type (Atom), intent(inout), allocatable           :: MoleculeInit(:)
      type (Atom), allocatable                          :: MoleculeDummy(:)
      type (Parameters), intent(inout)                  :: Variables
      real                                              :: Denergy, ChangeX, ChangeY, ChangeZ,   &
                                                           VectorX, VectorY, VectorZ, NewEnergy, &   
                                                           Chance, q, random    
      logical, intent(inout)                            :: Check
      integer                                           :: k
      integer, intent(in)                               :: NumberofAtoms
      type (Binding), intent(inout), allocatable        :: Bonds(:)
      
      !!! Select a particle randomly (k is the index of the atoms)
      
      k = 10 
      do while (k >= 9 .or. k == 0)
         call random_number(random)
         k = int(10*random)
      enddo
      
      !!! Save the configuration of the molecule in a dummy variable
      
      allocate(MoleculeDummy(NumberofAtoms))
      MoleculeDummy = MoleculeInit

      !!! Give that particle a random displacement
      
      call random_number(ChangeX) 
      call random_number(ChangeY)
      call random_number(ChangeZ)
      call random_number(VectorX)
      call random_number(VectorY)
      call random_number(VectorZ)

      if (ChangeX >= 0.5) then
         VectorX = -VectorX
      endif
      
      if (ChangeY >= 0.5) then
         VectorY = -VectorY
      endif
      
      if (ChangeZ >= 0.5) then
         VectorZ = -VectorZ
      endif

      MoleculeDummy(k)%x = MoleculeDummy(k)%x + Variables%r*VectorX
      MoleculeDummy(k)%y = MoleculeDummy(k)%y + Variables%r*VectorY
      MoleculeDummy(k)%z = MoleculeDummy(k)%z + Variables%r*VectorZ

      !!! Calculate the energy of new configuration with dummy variable and variable NewEnergy

      call TotalEnergy(NewEnergy, MoleculeDummy, Variables, Bonds, NumberofAtoms)

      !!! Check if the new configuration is accepted. If the outcome is false, it returns the logical 
      !!! 'Check' as false and the new configuration is not accepted. Otherwise, 'Check' is returned 
      !!! as true and the current configuration is replaced with the new one. 

      Denergy = NewEnergy - CurrentEnergy

      if (Denergy < 0) then
         check = .TRUE.
         MoleculeInit = MoleculeDummy
         CurrentEnergy = NewEnergy
      else
         Chance = exp((-Denergy)/(Variables%Boltzmann*Variables%TEMP))
         call random_number(q)
         if (q < Chance) then
            Check = .TRUE.
            MoleculeInit = MoleculeDummy
            CurrentEnergy = NewEnergy   
            if (abs(NewEnergy - CurrentEnergy) <= 0.0000001) then
               Check = .FALSE.
            else
               Check = .TRUE.
            endif
         else
            Check = .FALSE.
         endif
      endif

   endsubroutine

   subroutine Minimization(CurrentEnergy, MoleculeInit, Variables, Check, NumberofAtoms, Bonds)
      real, intent(inout)                               :: CurrentEnergy
      type (Atom), intent(inout), allocatable           :: MoleculeInit(:)
      type (Parameters), intent(inout)                  :: Variables
      logical, intent(inout)                            :: Check
      integer, intent(in)                               :: NumberofAtoms
      type (Binding), intent(inout), allocatable        :: Bonds(:)
      integer                                           :: ncycle

      ncycle = 0
      
      !!! The subroutine 'Move' is called a certain number of times. If 'ncycle' reaches the same 
      !!! number as 'maxcycle', it means that the subroutine has been called so many consecutive 
      !!! times with as a result an unsuccesful attempt of changing the configuration of the molecule. 

      do while (ncycle < Variables%maxcycle)
         Call Move(CurrentEnergy, MoleculeInit, Variables, Check, NumberofAtoms, Bonds)
         if (Check .eqv. .true.) then
            ncycle = 0
         elseif (Check .eqv. .false.) then
            ncycle = ncycle + 1
         endif
      enddo

   endsubroutine

endmodule
