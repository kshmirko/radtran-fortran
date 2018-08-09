module miev0mod
  
  use mathutils
  implicit none
  
  type AerosolDistrParams
    integer         ::  npts
    real(kind=dp)   ::  r0, r1, gamma, dens
  end type AerosolDistrParams

  type, extends(AerosolDistrParams) :: SDistParams
    real(kind=dp)       ::  wl
    complex(kind=dp)    ::  midx
    integer             ::  momdim
  end type SDistParams 
contains
  
  subroutine miev0easydriver(xx, midx, qext, qsca, gqsc, pmom, mxmdm)
    ! dummy parameters
    real(kind=dp)             ::  xx
    complex(kind=dp)          ::  midx
    integer, optional         ::  mxmdm
    real(kind=dp)             ::  qext, qsca, gqsc
    real(kind=dp)             ::  pmom(0:,:)
    
    ! actually used parametrs
    LOGICAL  ANYANG, PERFCT, PRNT(2)
    INTEGER  IPOLZN, MOMDIM, NUMANG, NMOM
    REAL     GQSC_S, MIMCUT, QEXT_S, QSCA_S, SPIKE, XMU
    COMPLEX  CREFIN, SFORW, SBACK, S1(1), S2(1), TFORW(2), TBACK(2)
    REAL, ALLOCATABLE ::  PMOM_TMP(:,:)
    
    EXTERNAL MIEV0
    INTEGER I, J
    
    if(.not. present(mxmdm)) then
      MOMDIM = size(pmom, 1)-1
    else
      MOMDIM = mxmdm
    end if
    
    ANYANG = .TRUE.
    PERFCT = .FALSE.
    PRNT   = .FALSE.
    NUMANG = 1
    NMOM   = MOMDIM!3*INT(XX)
    CREFIN = CMPLX(midx)
    MIMCUT = 1.0E-8
    IPOLZN = +1234
    ALLOCATE(PMOM_TMP(0:MOMDIM, 4))
    
    
    XMU = 0.0
    
    call MIEV0( XX, CREFIN, PERFCT, MIMCUT, ANYANG, NUMANG, XMU,&
                NMOM, IPOLZN, MOMDIM, PRNT, QEXT_S, QSCA_S, GQSC_S,&
                PMOM_TMP, SFORW, SBACK, S1, S2, TFORW, TBACK,&
                SPIKE )
    
    
    ! COPY VALUES
    qext = DBLE(QEXT_S)
    qsca = DBLE(QSCA_S)
    gqsc = DBLE(GQSC_S)

    ! COPY MOMENTS
    DO I=1, 4
      DO J=0, MOMDIM
        pmom(J,I) = DBLE(PMOM_TMP(j,I))
      END DO
    END DO
    
   
    DEALLOCATE(PMOM_TMP)
    
  end subroutine miev0easydriver

  subroutine miesdist1(params, pmom, ext, sca, asy, vol, ierr)
    type(SDistParams), intent(in)   ::  params
    real(kind=dp), intent(out), optional  ::  ext, sca, asy, vol
    integer, intent(out), optional        ::  ierr
    real(kind=dp), intent(inout)          ::  pmom(0:,:)

    real(kind=dp), allocatable            ::  xs(:), ys(:)

    allocate(xs(params%npts))
    call linspace(params%r0, params%r1, xs)
    ys = xs**(params%gamma)

    call miesdist(xs, ys, params%midx, params%wl, pmom, ext, sca, asy, vol, ierr)

    deallocate(xs, ys)
  end subroutine miesdist1
  
  subroutine miesdist(xs, ys, midx, wl, pmom, ext, sca, asy, vol, ierr)
    real(kind=dp), intent(in)   ::  xs(:), ys(:), wl
    complex(kind=dp), intent(in)::  midx
    real(kind=dp), intent(inout)::  pmom(0:, :)
    real(kind=dp), intent(out), optional  ::  ext, sca, asy, vol
    integer, intent(out), optional        ::  ierr
    
    integer     ::  nsize, I,J,L
    real(kind=dp) ::  k, xx, xsquared, xcubed, qext, qsca, gqsc
    real(kind=dp), allocatable  ::  pmom_tmp(:,:,:), extinction(:),&
                                    scattering(:), asymetry(:), volume(:),&
                                    pmom_i(:,:)
    real(kind=dp) :: norm
    integer       :: momdim    
    
    momdim = size(pmom, 1)-1
    nsize = size(xs,1)
    if (nsize/=size(ys,1)) then
      if(present(ierr)) then
        ierr=-1
      end if
      return
    end if
    
    k = 2.0_dp*pi/wl
    
    allocate(pmom_tmp(nsize, 0:momdim, 4), pmom_i(0:momdim, 4))
    allocate(extinction(nsize), scattering(nsize), asymetry(nsize),&
              volume(nsize))
              
    do I=1, nsize
      xx = k*xs(I)
      call miev0easydriver(xx, midx, qext, qsca, gqsc, pmom_i)
      xsquared = xs(I)*xs(I)*1.0D-12
      xcubed = xsquared*xs(I)*1.0D-6
      
      extinction(I) = qext*xsquared*ys(I)
      scattering(I) = qsca*xsquared*ys(I)
      asymetry(I) = gqsc*scattering(I)
      volume(I) = xcubed*ys(I)
      
      do J=1, 4
        do L=0, momdim
          pmom_tmp(I, L, J) = pmom_i(L,J)*ys(I)
        end do
      end do
    
    end do
    
    norm = trapz(xs, ys)

    if (present(ext)) ext = trapz(xs, extinction) / norm
    
    if (present(sca)) sca = trapz(xs, scattering) / norm
    
    if (present(asy)) asy = trapz(xs, asymetry) / norm
    
    if (present(vol)) vol = trapz(xs, volume) / norm
    
    DO I=1, 4
      DO J=0, momdim
        pmom(J,I) = trapz(xs, pmom_tmp(:, J, I)) / norm
      END DO
    END DO
    
    deallocate(pmom_tmp, pmom_i, extinction, scattering, asymetry, volume)
    ierr = 0
    return
    
  end subroutine miesdist
  
  subroutine wiscombe2evans(wisc, evans)
    real(kind=dp), intent(in) ::  wisc(0:, :)
    real(kind=dp), intent(out)::  evans(0:, :)
    real(kind=dp)             ::  norm, factor
    integer ::  I
    norm = wisc(0,1)+wisc(0,2)
    
    DO I=lbound(wisc,1), ubound(wisc,1)
      factor = dble(2*I+1)
      evans(I,1) = (wisc(I,1)+wisc(I,2))/norm*factor
      evans(I,2) = (wisc(I,2)-wisc(I,1))/norm*factor
      evans(I,3) = 2.0_dp*wisc(I,3)/norm*factor
      evans(I,4) = 2.0_dp*wisc(I,4)/norm*factor
      evans(I,5) = evans(I,1)
      evans(I,6) = evans(I,3)
    END DO 
    
  end subroutine wiscombe2evans
end module miev0mod
