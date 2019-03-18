!!!____________________________________________________________________________
!!! Author: Aaron Eisses
!!! Studentnumber: 10593527
!!! University: University of Amsterdam
!!! Date changed: 11-02-2019
!!!____________________________________________________________________________
!!!
!!! PROJECT SCIENTIFIC COMPUTING AND PROGRAMMING: MOLECULAR MECHANICS
!!!
!!! MinimizationModule.f90
!!!
!!! This file is the module minimization. Here, there is only one subroutine
!!! called Move. This subroutine does an energy calculation of the current
!!! configuration, moves an atom with a random displacement and the calculates
!!! the energy of the new configuration. With a certain probility function, 
!!! the new configuration gets accepted or rejected. 
!!!____________________________________________________________________________

module MinimizationModule
use DataModule
use EnergyModule
   implicit none
   save 
   private
   public               :: Move

contains

   subroutine Move(CurrentEnergy, MoleculeInit, Variables, Check, r)
      real, intent(inout)               :: CurrentEnergy
      type (Atom), intent(inout)        :: MoleculeInit(8)
      type (Atom)                       :: MoleculeDummy(8)
      type (Parameters), intent(inout)  :: Variables
      real                              :: Denergy, ChangeX, ChangeY, ChangeZ, &
                                           VectorX, VectorY, VectorZ, NewEnergy&
                                           , Chance, q, random    
      logical, intent(inout)            :: Check
      integer                           :: k
      real, intent(in)                  :: r
      
      

      !!! Select a particle randomly (k is the index of the atoms)
      
      k = 10 
      do while (k >= 9 .or. k == 0)
         call random_number(random)
         k = int(10*random)
      enddo
      !!! Save the configuration of the molecule in a dummy variable
      
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

      MoleculeDummy(k)%x = MoleculeDummy(k)%x + r*VectorX
      MoleculeDummy(k)%y = MoleculeDummy(k)%y + r*VectorY
      MoleculeDummy(k)%z = MoleculeDummy(k)%z + r*VectorZ

      !!! Calculate the energy of new configuration with dummy variable
      !!! and variable NewEnergy

      call TotalEnergy(NewEnergy, MoleculeDummy, Variables)

      !!! Check whether the new configuration is accepted


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
         else
            if (abs(NewEnergy - CurrentEnergy) >= 0.0000001) then
               Check = .FALSE.
            else
               Check = .TRUE.
            endif
         endif
      endif
      
      !!! update new configuration if check = .TRUE.
      !!! for the exit strategy, it could return the value 1 and
      !!! add it to a dummy variable in the mainprogram. In that way
      !!! the amount of consecutive cycles can be checked. 

      

   endsubroutine


endmodule