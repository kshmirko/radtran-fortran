module atmos
use mathutils
use miev0mod

implicit none

real(kind=dp), allocatable :: Evans_mol(:,:) 
contains

subroutine init_evansmol(momdim)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Создаем матрицу эванса для молекулярного рассеяния
! для дальнейшего использования
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    integer, intent(in) ::  momdim
    integer ::  npts, ierr
    real(kind=dp)   ::  wl
    real(kind=dp), allocatable  :: xs(:), ys(:), pmom(:,:)


    wl = 0.355
    npts = 101
    allocate(xs(npts))
    call linspace(0.0001_dp, 0.0005_dp, xs)
    ! здесь происходит автоматическое выделение памяти
    ys = xs**(-7.0)
    if (allocated(Evans_mol)) deallocate(Evans_mol)

    allocate(Evans_mol(0:momdim,6))
    allocate(pmom(0:momdim, 4))

    call miesdist(xs, ys, dcmplx(1.5,0.0), wl, pmom, ierr=ierr)
    call wiscombe2evans(pmom, Evans_mol)
    deallocate(pmom, xs, ys)

end subroutine init_evansmol

subroutine atmos_evams(evansa, taua, taum, evanstot, ierr)
    real(kind=dp), intent(in)       :: evansa(:,:)
    real(kind=dp), intent(in)       :: taua, taum
    real(kind=dp), intent(out)      :: evanstot(:,:)
    integer, intent(out), optional  :: ierr
    integer                         :: N1, M1, N2, M2

    N1 = size(evansa,1)
    M1 = size(evansa,2)

    N2 = size(evanstot,1)
    M2 = size(evanstot,2)

    if ((N1/=N2).or.(M1/=M2)) then
        if(present(ierr))ierr=1
        return
    end if
    evanstot = 0.0_dp
    evanstot = (taua*evansa+taum*Evans_mol) / (taua+taum)
    if(present(ierr))ierr=0
    return
    
end subroutine atmos_evams

subroutine write_sca_file(fname, ext, sca, evans)
    real(kind=dp), intent(in)       :: ext, sca
    real(kind=dp), intent(in)       :: evans(:,:)
    character(len=*), intent(in)        :: fname
    integer I
    
    open(100, FILE=fname, status='replace')
    write(100, '(E12.4)') ext
    write(100, '(E12.4)') sca
    write(100, '(E12.4)') sca/ext
    write(100, '(I3)') size(evans,1)-1

    do i=1, size(evans,1)
        write(100, '(I3,6F10.6)') i-1, evans(I,:)
    end do
    close(100)

end subroutine write_sca_file


subroutine exta_at_h(taua, hpbl, h, ext)
    real(kind=dp), intent(in)   ::  taua, hpbl
    real(kind=dp), intent(in)   ::  h(:)
    real(kind=dp), intent(out)  ::  ext(:)
    real(kind=dp)               ::  exta0

    exta0 = taua/hpbl

    ext = exta0*dexp(-h/hpbl)
    return
end subroutine exta_at_h

subroutine extm_at_h(wl, h, ext)
    use rayleigh, only : ext_m, h_mol
    real(kind=dp), intent(in)   ::  wl
    real(kind=dp), intent(in)   ::  h(:)
    real(kind=dp), intent(out)  ::  ext(:)
    real(kind=dp)               ::  extm0

    extm0 = ext_m(wl)

    ext = extm0*dexp(-h/h_mol)
    return
end subroutine extm_at_h
    
end module atmos