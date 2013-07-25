program main
    use eigen, only: EigenMnk12
    implicit none
    call runMain()
contains
    subroutine runMain()
        double precision, dimension(1:4) :: eigenvals
        complex*16, dimension(1:4, 1:4) :: eigenvectors
        integer :: error
        error = EigenMnk12(1, eigenvals, eigenvectors)
        print *, error
        print *, eigenvals(1), eigenvals(2), eigenvals(3), eigenvals(4)
        print *, eigenvectors
    end subroutine runMain
end program main
