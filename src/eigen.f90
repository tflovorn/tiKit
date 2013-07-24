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
        integer :: EigenDecomp, lwork
    end function EigenDecomp

    ! Helper function to make the call to ZHEEV - required since we may need
    ! to adjust the length of work.
    recursive function EigenDecompHelper(H, N, eigenvals, lwork) result (info)
        complex*16, dimension(1:N, 1:N), intent(inout) :: H
        integer, intent(in) :: N
        double precision, dimension(1:N), intent(out) :: eigenvals
        integer, intent(in) :: lwork
        integer :: info
        double precision, dimension(1:lwork) :: work
        ! 'V': return eigenvalues and eigenvectors
        ! 'U': input H should be upper triangular
        !call ZHEEV('V', 'U', N, H, N, eigenvals, work, lwork, info)
        ! work(1) contains optimal value of lwork, which is the length of work
        info = 0
    end function EigenDecompHelper

    function EigenMnk12()
        integer :: EigenMnk12, info
        info = 0
        if (info /= 0) then
            ! failure - TODO handle error
        end if
    end function EigenMnk12
end module eigen
