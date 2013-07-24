program main
    use eigen
    implicit none
    call runMain()
contains
    subroutine runMain()
        integer error
        error = EigenMnk12()
    end subroutine runMain
end program main
