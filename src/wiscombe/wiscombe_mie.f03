module wiscombe
implicit none
contains

subroutine calculate_mie(WL, R0, R1, NR, GAMMA, CREFIN, NUMANG, PMOMS, omega)
    use legendre
    use trapez
    use distribution
    
    implicit none

    real, intent(in)    ::  R0, R1, GAMMA, WL
    complex,intent(in) ::  CREFIN
    integer,intent(in) ::  NR, NUMANG
    real,intent(out), allocatable ::  PMOMS(:, :)
    real, intent(out)   ::  omega

    real(8) ::  xx

    logical  ANYANG, PERFCT, PRNT(2)
    integer  IPOLZN, NMOM, MAXMOM
    real     GQSC, MIMCUT
    real, allocatable   ::  PMOM( :, : ), QEXT(:), QSCA(:), XMU(:), TMOM( :, :, : ), fd(:), ri(:)
    complex, allocatable::  S1(:), S2(:)
    complex  SFORW, SBACK, TFORW(2), TBACK(2)
    real dang, norm, SPIKE

    external  miev0
    integer   i, j
    real k, dr

    k = real(2.0*PI/WL)
    mimcut = 1.e-6
    dang=2.0/real(numang-1)

    xx = k*R1

    ! определяем количество коэффициентов разложения
    MAXMOM = int(3.0d0*xx)
    if (MAXMOM<2) MAXMOM=2

    allocate(PMOM(0:MAXMOM, 4))
    allocate(xmu(NUMANG))
    allocate(S1(NUMANG))
    allocate(S2(NUMANG))
    allocate(QEXT(NR))
    allocate(QSCA(NR))
    allocate(ri(NR), fd(NR))
    allocate(PMOMS(0:MAXMOM,6))
    allocate(TMOM(NR,0:MAXMOM,4))



    xmu = [(-1.0+I*dang, I=0, NUMANG-1)]
    
    PMOM = 0.0
    qext = 0.0
    qsca = 0.0
    S1 = 0.0
    S2 = 0.0
    xmu = 0.0

    perfct = .false.
    anyang = .true.
    ipolzn = +1234
  

    ! Отключаем возможность печати служебной информации
    prnt( 1 ) = .false.
    prnt( 2 ) = .false.

    ! инициализируем массив накоплений
    TMOM = 0.0 

    dr = (r1-r0)/(NR-1)
    ri = [(r0+I*dr, I=0, NR-1)]
    norm=0.0

    do i=1, NR

        xx = dble(k*ri(i))

        nmom = int(3.d0*xx)
        if(nmom<2) nmom=2
          

        SFORW=0.0
        SBACK=0.0
        TFORW=0.0
        TBACK=0.0

        call MIEV0(XX, CREFIN, PERFCT, MIMCUT, ANYANG, NUMANG, XMU, &
          & NMOM, IPOLZN, MAXMOM, PRNT, QEXT(I), QSCA(I), GQSC,&
          & PMOM, SFORW, SBACK, S1, S2, TFORW, TBACK,&
          & SPIKE )
        
        !print *, 'Omega_i', QSCA(I)/QEXT(I)
        
        fd(i) = power_law(real(ri(i)), gamma)
        
        QEXT(I)=QEXT(I)*real(fd(I)*PI*ri(I)**2)
        QSCA(I)=QSCA(I)*real(fd(I)*PI*ri(I)**2)
        
        TMOM(I,:,:)=PMOM*fd(i)

    end do

    norm = trapz1d(ri, fd, NR)
    omega = trapz1d(ri, QSCA, NR)/trapz1d(ri, QEXT, NR)

    PMOM = 0.0

    do j=1, 4
        do i=0, NMOM
            PMOM(i,j) = trapz1d(ri, TMOM(:,i,j), NR) / norm*real(2*i+1)
        end do
    end do
    
    
    ! wiscombe to evans conversion

    norm=PMOM(0,1)+PMOM(0,2)
    
    PMOMS(:,1) = (PMOM(:,1)+PMOM(:,2))/norm
    PMOMS(:,2) = (PMOM(:,2)-PMOM(:,1))/norm
    PMOMS(:,3) = 2.0*PMOM(:,3)/norm
    PMOMS(:,4) = 2.0*PMOM(:,4)/norm
    PMOMS(:,5) = PMOMS(:,1)
    PMOMS(:,6) = PMOMS(:,3)

    !PMOMS(:,2) = PMOMS(:,2)/PMOMS(0,1)
    !PMOMS(:,3) = PMOMS(:,3)/PMOMS(0,1)
    !PMOMS(:,4) = PMOMS(:,4)/PMOMS(0,1)
    !PMOMS(:,5) = PMOMS(:,5)/PMOMS(0,1)
    !PMOMS(:,6) = PMOMS(:,6)/PMOMS(0,1)
    !PMOMS(:,1) = PMOMS(:,1)/PMOMS(0,1)

    deallocate(PMOM, xmu, S1, S2, QEXT, QSCA, ri, fd, TMOM)

end subroutine calculate_mie

subroutine make_sca_file_aer(WL, R0, R1, NR, GAMMA, CREFIN, FNAME, NUMANG, EXT)
    implicit none
    real, intent(in)    ::  R0, R1, GAMMA, WL, EXT
    complex,intent(in) ::  CREFIN
    integer,intent(in) ::  NR, NUMANG
    real, allocatable ::  PMOMS(:, :)
    real    ::  omega, sca
    CHARACTER(*)    ::  FNAME
    integer I, NMOM(2)

    call calculate_mie(WL, R0, R1, NR, GAMMA, CREFIN, NUMANG, PMOMS, omega)
    SCA=EXT*omega

    open (101, file=FNAME, status='replace')
    write(101, '(E12.6)') EXT
    write(101, '(E12.6)') SCA
    write(101, '(E12.6)') omega
    NMOM = shape(PMOMS)

    write(101, '(I5)') NMOM(1)-1
    DO I=0, NMOM(1)-1
        write(101, '(I5, 6E15.4)') I, PMOMS(I,:)
    END DO
    close(101)

end subroutine make_sca_file_aer

subroutine make_sca_file_tot(WL, R0, R1, NR, GAMMA, CREFIN, FNAME, NUMANG, EXTA, EXTM)
    implicit none
    real, intent(in)    ::  R0, R1, GAMMA, WL, EXTA, EXTM
    complex,intent(in) ::  CREFIN
    integer,intent(in) ::  NR, NUMANG
    real, allocatable ::  PMOMS(:, :), PMOMSM(:,:)
    real    ::  omega, sca, EXT, omegam
    CHARACTER(*)    ::  FNAME
    integer I, NMOM, shp1(2), shp2(2), funit1

    call calculate_mie(WL, R0, R1, NR, GAMMA, CREFIN, NUMANG, PMOMS, omega)
    ! молекулярные коэффициенты лежандра мы можем получить задав маленький размер частиц
    ! и нулемую мнимую часть показателя преломления, показатель ангстрема положим равным 0 
    ! (нет спектральной зависимости). В этом случае PMOMSM будет содержать всего 3 коэффициента 
    ! разложения в ряд по полиномам Лежандра
    call calculate_mie(WL, 0.001, 0.0015, 10, 0.0, (1.4, 0.0), NUMANG, PMOMSM, omegam)

    shp1=shape(PMOMS)
    shp2=shape(PMOMSM)

    ! Вычисляем суммарную экстинкцию и рассеяние 
    EXT=EXTA+EXTM
    SCA=(EXTA*omega+EXTM)
    ! новое альбедо
    omega = SCA/EXT

    ! вычислем коэфициенты разложения для средневзвешенной фазовой функции (аэрозоль и молекулы) 
    DO I=0, shp2(1)-1
        PMOMS(I,:)=(PMOMS(I,:)*EXTA+PMOMSM(I,:)*EXTM)/EXT
    END DO

    DO I=shp2(1), shp1(1)-1
        PMOMS(I,:)=(PMOMS(I,:)*EXTA)/EXT
    END DO

    open (newunit=funit1, file=FNAME, status='replace')
    write(funit1, '(E12.6)') EXT
    write(funit1, '(E12.6)') SCA
    write(funit1, '(E12.6)') omega
    NMOM = shp1(1)

    write(funit1, '(I5)') NMOM-1
    DO I=0, NMOM-1
        write(funit1, '(I5, 6E15.4)') I, PMOMS(I,:)
    END DO
    close(funit1)

end subroutine make_sca_file_tot

end module wiscombe