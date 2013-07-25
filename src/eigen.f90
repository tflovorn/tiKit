module eigen
    implicit none
    external ZHEEV

    complex*16, dimension(1:4, 1:4) :: Gamma1, Gamma2, Gamma3, Gamma4, Gamma5, Ident

    type :: Mnk12Params
        double precision :: gamma_fm, Delta_n, M, A, B, E_F, C, gamma_c
    end type Mnk12Params
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
    
    subroutine InitGammas()
        ! Liu (2010) notation for Gamma matrices
        Gamma1(:, 1) = (/(0.0D0, 0.0D0), (0.0D0, 0.0D0), (0.0D0, 0.0D0), (1.0D0, 0.0D0)/)
        Gamma1(:, 2) = (/(0.0D0, 0.0D0), (0.0D0, 0.0D0), (1.0D0, 0.0D0), (0.0D0, 0.0D0)/)
        Gamma1(:, 3) = (/(0.0D0, 0.0D0), (1.0D0, 0.0D0), (0.0D0, 0.0D0), (0.0D0, 0.0D0)/)
        Gamma1(:, 4) = (/(1.0D0, 0.0D0), (0.0D0, 0.0D0), (0.0D0, 0.0D0), (0.0D0, 0.0D0)/)

        Gamma2(:, 1) = (/(0.0D0, 0.0D0), (0.0D0, 0.0D0), (0.0D0, 0.0D0), (0.0D0, -1.0D0)/)
        Gamma2(:, 2) = (/(0.0D0, 0.0D0), (0.0D0, 0.0D0), (0.0D0, 1.0D0), (0.0D0, 0.0D0)/)
        Gamma2(:, 3) = (/(0.0D0, 0.0D0), (0.0D0, -1.0D0), (0.0D0, 0.0D0), (0.0D0, 0.0D0)/)
        Gamma2(:, 4) = (/(0.0D0, 1.0D0), (0.0D0, 0.0D0), (0.0D0, 0.0D0), (0.0D0, 0.0D0)/)

        Gamma3(:, 1) = (/(0.0D0, 0.0D0), (1.0D0, 0.0D0), (0.0D0, 0.0D0), (0.0D0, 0.0D0)/)
        Gamma3(:, 2) = (/(1.0D0, 0.0D0), (0.0D0, 0.0D0), (0.0D0, 0.0D0), (0.0D0, 0.0D0)/)
        Gamma3(:, 3) = (/(0.0D0, 0.0D0), (0.0D0, 0.0D0), (0.0D0, 0.0D0), (-1.0D0, 0.0D0)/)
        Gamma3(:, 4) = (/(0.0D0, 0.0D0), (0.0D0, 0.0D0), (-1.0D0, 0.0D0), (0.0D0, 0.0D0)/)

        Gamma4(:, 1) = (/(0.0D0, 0.0D0), (0.0D0, -1.0D0), (0.0D0, 0.0D0), (0.0D0, 0.0D0)/)
        Gamma4(:, 2) = (/(0.0D0, 1.0D0), (0.0D0, 0.0D0), (0.0D0, 0.0D0), (0.0D0, 0.0D0)/)
        Gamma4(:, 3) = (/(0.0D0, 0.0D0), (0.0D0, 0.0D0), (0.0D0, 0.0D0), (0.0D0, -1.0D0)/)
        Gamma4(:, 4) = (/(0.0D0, 0.0D0), (0.0D0, 0.0D0), (0.0D0, 1.0D0), (0.0D0, 0.0D0)/)

        Gamma5(:, 1) = (/(1.0D0, 0.0D0), (0.0D0, 0.0D0), (0.0D0, 0.0D0), (0.0D0, 0.0D0)/)
        Gamma5(:, 2) = (/(0.0D0, 0.0D0), (-1.0D0, 0.0D0), (0.0D0, 0.0D0), (0.0D0, 0.0D0)/)
        Gamma5(:, 3) = (/(0.0D0, 0.0D0), (0.0D0, 0.0D0), (1.0D0, 0.0D0), (0.0D0, 0.0D0)/)
        Gamma5(:, 4) = (/(0.0D0, 0.0D0), (0.0D0, 0.0D0), (0.0D0, 0.0D0), (-1.0D0, 0.0D0)/)

        Ident(:, 1) = (/(1.0D0, 0.0D0), (0.0D0, 0.0D0), (0.0D0, 0.0D0), (0.0D0, 0.0D0)/)
        Ident(:, 2) = (/(0.0D0, 0.0D0), (1.0D0, 0.0D0), (0.0D0, 0.0D0), (0.0D0, 0.0D0)/)
        Ident(:, 3) = (/(0.0D0, 0.0D0), (0.0D0, 0.0D0), (1.0D0, 0.0D0), (0.0D0, 0.0D0)/)
        Ident(:, 4) = (/(0.0D0, 0.0D0), (0.0D0, 0.0D0), (0.0D0, 0.0D0), (1.0D0, 0.0D0)/)
    end subroutine

    function DefaultMnk12Params()
        type(Mnk12Params) :: DefaultMnk12Params
        DefaultMnk12Params%gamma_fm = 1.0D0
        DefaultMnk12Params%Delta_n = 1.0D0
        DefaultMnk12Params%M = 0.3D0
        DefaultMnk12Params%A = 0.5D0
        DefaultMnk12Params%B = 0.25D0
        DefaultMnk12Params%E_F = 3.1D0
        DefaultMnk12Params%C = 3.0D0
        DefaultMnk12Params%gamma_c = 0.25D0
    end function DefaultMnk12Params

    function EigenMnk12(p, numLayers, k, eigenvals, eigenvectors)
        type(Mnk12Params), intent(in) :: p
        integer, intent(in) :: numLayers
        double precision, dimension(1:2) :: k
        double precision, dimension(1:4*numLayers), intent(out) :: eigenvals
        complex*16, dimension(1:4*numLayers, 1:4*numLayers), intent(out) :: eigenvectors
        integer :: EigenMnk12, i
        complex*16, dimension(1:4, 1:4) :: diagonal, cross, cross_conj
        call InitGammas()

        ! block-diagonal piece
        diagonal = p%C*Ident + (p%M + 2.0D0*p%B*(cos(k(1)) + cos(k(2)) - 3.0D0))*Gamma5 &
                     + p%A*(sin(k(1))*Gamma2 + sin(k(2))*Gamma1)
        ! hop layer n+1 to n
        cross = p%B*Gamma5 - (0.0D0, 0.5D0)*p%A*Gamma4
        ! hop layer n to n+1
        cross_conj = conjg(transpose(cross))
        
        ! initialize Hamiltonian
        eigenvectors(:, :) = (0.0D0, 0.0D0)
        do i=1, 4*numLayers, 4  ! i = 1, 5, 9, ..., 4*numLayers - 3
            eigenvectors(i:i+3, i:i+3) = diagonal
            if (i < 4*numLayers - 3) then
                print *, "cross"
                eigenvectors(i+4:i+7, i:i+3) = cross
            end if
            if (i > 1) then
                print *, "cross_conj"
                eigenvectors(i-4:i-1, i:i+3) = cross_conj
            end if
        end do
        print *, eigenvectors(:, 1)
        print *, eigenvectors(:, 2)
        print *, eigenvectors(:, 3)
        print *, eigenvectors(:, 4)
        EigenMnk12 = EigenDecomp(eigenvectors, 4*numLayers, eigenvals)
        if (EigenMnk12 /= 0) then
            ! TODO handle error
        end if
    end function EigenMnk12
end module eigen
