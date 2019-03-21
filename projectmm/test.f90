real function RandomNumber()
    
      call random_number(RandomNumber)

endfunction





program test
        real, allocatable :: r1(:)
        allocate(r1(3))

        r1(1) = RandomNumber()

        print *, r1
        deallocate(r1)

        print *, r1









endprogram

