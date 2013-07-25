module eigen
    implicit none
    external ZHEEV
contains
    ! Call LAPACK routine ZHEEV to get eigenvalues and eigenvectors of the
    ! Hermitian, upper triangular matrix H of order N.
    ! Returns 0 on success.
    function EigenDecomp(H, N, eigenvals)
        complex*16, dimension(1:N, 1:N), intent(inout) :: H
        integer, intent(in) :: N
        double precision, dimension(1:N), intent(out) :: eigenvals
        complex*16, dimension(1:1) :: first_work
        double precision, dimension(1:3*N-2) :: rwork
        complex*16, dimension(:), allocatable :: work
        integer :: EigenDecomp, opt_lwork
        if (N <= 0) then
            ! TODO - handle invalid argument
            return
        end if
        ! get the optimal value for lwork
        call ZHEEV('V', 'U', N, H, N, eigenvals, first_work, -1, rwork, EigenDecomp)
        if (EigenDecomp /= 0) then
            return
        end if
        opt_lwork = first_work(1)
        ! Make optimal sized work array.
        ! allocate() crashes on error; TODO handle?
        ! Can get error code with allocate(work(1:opt_lwork), STAT=code)
        ! error if code = 0
        allocate(work(1:opt_lwork)) 
        ! make the real call to ZHEEV
        call ZHEEV('V', 'U', N, H, N, eigenvals, work, opt_lwork, rwork, EigenDecomp)
    end function EigenDecomp
    
    function EigenMnk12(numLayers, eigenvals, eigenvectors)
        integer, intent(in) :: numLayers
        double precision, dimension(1:4*numLayers), intent(out) :: eigenvals
        complex*16, dimension(1:4*numLayers, 1:4*numLayers), intent(out) :: eigenvectors
        integer :: EigenMnk12
        ! initialize Hamiltonian in eigenvectors, since call to EigenDecomp
        ! will overwrite H with eigenvectors
        ! TODO - replace temp eigenvectors - assumes numLayers = 1
        eigenvectors(:, :) = (0.0D0, 0.0D0)
        eigenvectors(1, 1) = (1.0D0, 0.0D0)
        eigenvectors(2, 2) = (2.0D0, 0.0D0)
        eigenvectors(3, 3) = (3.0D0, 0.0D0)
        eigenvectors(4, 4) = (4.0D0, 0.0D0)
        print *, eigenvectors
        EigenMnk12 = EigenDecomp(eigenvectors, 4*numLayers, eigenvals)
        if (EigenMnk12 /= 0) then
            ! TODO handle error
        end if
    end function EigenMnk12
end module eigen
