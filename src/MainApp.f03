! program MainApp
!   use miev0mod
!   use rayleigh  
!   use aerosol
!   use atmos
!   implicit none
  
!   real(kind=dp) ::  hpbl, taua, galbedo
!   integer numazim, nmu
  
!   type(SDistParams)           :: params
!   CHARACTER*64 OUT_FILE

!   params%r0 = 0.1_dp
!   params%r1 = 1.0_dp
!   params%npts = 101
!   params%wl = 0.750_dp
!   params%midx = DCMPLX(1.4_dp, 0.01_dp)
!   params%gamma = -3.5_dp
!   params%dens = 300
!   hpbl = 3000.0_dp
!   numazim = 2!10
!   taua=0.1
!   galbedo = 0.3_dp
!   nmu = 32

!   write(OUT_FILE,'(A11,F3.1,A4)') 'rt3_rho_eq_', galbedo, '.dat'
  
!   call DO_CALC(params, hpbl, taua, numazim, galbedo, nmu, OUT_FILE)
!   galbedo = 0.0_dp

!   write(OUT_FILE,'(A11,F3.1,A4)') 'rt3_rho_eq_', galbedo, '.dat'
!   call DO_CALC(params, hpbl, taua, numazim, galbedo, nmu, OUT_FILE)


! end program


subroutine DO_CALC(params, hpbl, taua, numazim, galbedo, nmu, OUT_FILE) 
  use miev0mod
  use rayleigh  
  use aerosol
  use atmos
  implicit none

  real(kind=dp) :: taua, hpbl, galbedo
  integer                     :: numazim
  type(SDistParams)           :: params

  real(kind=dp) ::  ext, sca, asy, extm, exta, extt,&
                    omegat, scat
  integer mxmdm, ierr

  
  real(kind=dp)               :: omegaa
  

  !character(len=13)           :: fname
  integer I

  INTEGER   MAXV, MAXA, MAXLAY, NLAYS, nmu
  PARAMETER (MAXV=64, MAXA=32, NLAYS=40, mxmdm=60)
  PARAMETER (MAXLAY=200)

  real(kind=dp)  :: pmom(0:mxmdm,4), evans(0:mxmdm,6), evanstot(0:mxmdm,6)

  INTEGER   NSTOKES, NUMMU, AZIORDER
  PARAMETER(NSTOKES=2)
  INTEGER   NUM_LAYERS, SRC_CODE
  INTEGER   NOUTLEVELS, OUTLEVELS(MAXLAY), NUMAZIMUTHS
  REAL*8    GROUND_TEMP, GROUND_ALBEDO
  COMPLEX*16  GROUND_INDEX
  REAL*8    SKY_TEMP, WAVELENGTH, MAX_DELTA_TAU, THETA
  REAL*8    DIRECT_FLUX, DIRECT_MU
  REAL*8    MU_VALUES(MAXV)
  REAL*8    HEIGHT(MAXLAY), TEMPERATURES(MAXLAY)
  REAL*8    GAS_EXTINCT(MAXLAY)
  REAL*8    UP_FLUX(4*MAXLAY), DOWN_FLUX(4*MAXLAY)
  REAL*8    UP_RAD(MAXLAY*MAXA*MAXV), DOWN_RAD(MAXLAY*MAXA*MAXV)
  CHARACTER QUAD_TYPE*1, DELTAM*1, UNITS*1, OUTPOL*2, GROUND_TYPE*1
  CHARACTER*64 LAYER_FILE, OUT_FILE
  CHARACTER*64 SCAT_FILES(MAXLAY)
  REAL*8    PHIVAL(numazim)
  REAL*8    UP(NSTOKES*MAXA*MAXV*MAXLAY), DOWN(NSTOKES*MAXA*MAXV*MAXLAY)
  external RADTRAN, OUTPUT_FILE

  
  ! Для начала выделим память под хранение матриц с коэффициентами 
  ! Лежандра
  call init_evansmol(mxmdm)

  
  pmom = 0.0_dp
  ext = 0.0_dp
  sca = 0.0_dp
  asy = 0.0_dp

  ! Инициализация параметров
  !params%r0 = 0.1_dp
  !params%r1 = 1.0_dp
  !params%npts = 101
  !params%wl = 0.750_dp
  !params%midx = DCMPLX(1.4_dp, 0.01_dp)
  !params%gamma = -3.5_dp
  !params%dens = 300
  !hpbl = 3000.0_dp

  ! Расчет рассеяния
  call miesdist1(params, pmom, ext=ext, sca=sca, asy=asy, ierr=ierr)
  omegaa = sca/ext

  if (ierr/=0) then
    stop 'Ошибка в подпрограмме miesdist1'
  end if
  
  print '(A)', 'C  Параметры распределения:'
  print '(A)', 'C  ========================'
  print '(A, 2F5.2, A)', 'C  Диапазон размеров частиц: [', params%r0, params%r1, ' ]'
  print '(A, F5.2, A)', 'C  Показатель степени: [', params%gamma, ' ]'
  print '(A, F5.2, A)', 'C  Оптическая толща: [', taua, ']'
  print '(A, 2F5.2, A)', 'C  Показатель преломления: [', params%midx, ' ]'
  print '(A, F5.2, A)', 'C  Альбедо поверхноси: [', galbedo, ' ]'

  print '(A)', 'C  Оптические параметры заданного распределения при '
  print '(A)', 'C  условии единичной концентрации'
  print '(A,4A12)','C  ', 'EXT', 'SCA', 'ASY', 'OMEGA'
  print '(A,4E12.4)', 'C  ', ext, sca, asy, omegaa

  ! Матрица Эванса имеет следующие элеметы
  ! P1  P2  0   0  
  ! P2  P5  0   0
  ! 0   0   P3  P4
  ! 0   0  -P4  P6
  !print '(/,5X,A,/)', 'Матрица Эванса:'
  !print '(4A4)', 'P1', 'P2', '0', '0'
  !print '(4A4)', 'P2', 'P5', '0', '0'
  !print '(4A4)', '0', '0', 'P3', 'P4'
  !print '(4A4)', '0', '0', '-P4', 'P6'

  call wiscombe2evans(pmom, evans)

  
  ! Вычисление АОТ через концентрацию частиц и их оптические свойства
  !taua = ext*hpbl*1.0d6*params%dens
  

  !print '(A,F5.3,A,F5.3)', 'tau_m(',params%wl,') = ', tau_m(params%wl)
  !print '(A,F5.3,A,F5.3)', 'tau_a(',params%wl,') = ', taua
    

  !open(200, FILE='atmoslay.lay', status='replace')
  ! Вычитсяем оптические характеристики для слоев
  ! Сохраняем матрицы рассеяния
  print '(A, A7,3A10)', 'C  ','H, km', 'am, 1/km', 'aa, 1/km', 'at, 1/km'
  DO I=1, NLAYS
    HEIGHT(I) = DBLE(NLAYS-I)
    TEMPERATURES(I) = 0.0
    GAS_EXTINCT(I) = 0.0
    extm = ext_m(params%wl)*dexp(-HEIGHT(I)*1000_dp/(h_mol))
    exta = taua/hpbl*dexp(-HEIGHT(I)*1000_dp/hpbl)
    extt = exta+extm
    scat = exta*omegaa+extm
    omegat = scat/extt

    
    evanstot = (evans*exta+Evans_mol*extm)/extt
    
    write(SCAT_FILES(I),'(A,I3.3)') '.scat_file', I
    call write_sca_file(SCAT_FILES(I), extt, scat, evanstot)
    !write(200, '(F6.2,F7.2,E12.5,2X,A1,A13,A1)') HEIGHT(I), 0.0, 0.0, "'",fname,"'"
    print '(A, F7.1, 3E10.2)', 'C  ',HEIGHT(I), extm, exta, extt
  END DO

  !close(200)

  ! deallocate(pmom,  evans, evanstot)
  
  ! Начало радиационного расчета
  MAX_DELTA_TAU = 1.0E-6  ! Шаг по tau
  
  NUMMU = nmu             ! Число отсчетов на 90 градусах
  AZIORDER = 2            ! порядок разложения по азимуту
  SRC_CODE = 1            ! тип источника
  QUAD_TYPE = 'G'         ! тип квадратуры
  DELTAM = 'Y'            ! дельта-м масштабирование
  DIRECT_FLUX = 1.0_dp    ! поток радиации на верхней границе атмосфере
  THETA = 10.0_dp         ! зенитный угол солнца
  DIRECT_MU = DABS(DCOS(0.017453292D0*(THETA)))
  GROUND_TEMP = 0.0_dp    ! температура поверхности
  GROUND_TYPE = 'L'       ! тип поверхности
  GROUND_ALBEDO = galbedo  ! альбедо поверхности
  GROUND_INDEX = DCMPLX(1.0_dp, 0.0_dp)
  SKY_TEMP = 0.0_dp       ! температура неба
  WAVELENGTH = params%wl  ! длина волны
  NUM_LAYERS = NLAYS      ! количество слоев
  NOUTLEVELS = 1          ! количество слоев для вывода
  OUTLEVELS(1) = NLAYS    ! номер слоя

  ! Расчет радиации на нижней границе атмосферы
  CALL RADTRAN (NSTOKES, NUMMU, AZIORDER, MAX_DELTA_TAU,&
                         SRC_CODE, QUAD_TYPE, DELTAM,&
                         DIRECT_FLUX, DIRECT_MU,&
                         GROUND_TEMP, GROUND_TYPE,&
                         GROUND_ALBEDO, GROUND_INDEX,&
                         SKY_TEMP, WAVELENGTH,&
                         NUM_LAYERS, HEIGHT, TEMPERATURES,&
                         GAS_EXTINCT, SCAT_FILES,&
                         NOUTLEVELS, OUTLEVELS,&
                         MU_VALUES, UP_FLUX, DOWN_FLUX,&
                         UP_RAD, DOWN_RAD)

  !OUT_FILE = 'rt3.out'
  NUMAZIMUTHS = numazim
  LAYER_FILE = 'NONE'
  OUTPOL='IQ'
  ! Отображение результатов на экране
  CALL OUTPUT_FILE (NSTOKES, NUMMU, AZIORDER,&
                         SRC_CODE, LAYER_FILE, OUT_FILE,&
                         QUAD_TYPE, DELTAM, DIRECT_FLUX, DIRECT_MU,&
                         GROUND_TEMP, GROUND_TYPE,&
                         GROUND_ALBEDO, GROUND_INDEX,&
                         SKY_TEMP, WAVELENGTH, UNITS, OUTPOL,&
                         NUM_LAYERS, HEIGHT,&
                         NOUTLEVELS, OUTLEVELS, NUMAZIMUTHS,&
                         MU_VALUES, UP_FLUX, DOWN_FLUX,&
                         UP_RAD, DOWN_RAD)

!  CALL OUTPUT_RADS (NSTOKES, NUMMU, AZIORDER,&
!                         SRC_CODE, LAYER_FILE, OUT_FILE,&
!                         QUAD_TYPE, DELTAM, DIRECT_FLUX, DIRECT_MU,&
!                         GROUND_TEMP, GROUND_TYPE,&
!                         GROUND_ALBEDO, GROUND_INDEX,&
!                         SKY_TEMP, WAVELENGTH, UNITS, OUTPOL,&
!                         NUM_LAYERS, HEIGHT,&
!                         NOUTLEVELS, OUTLEVELS, NUMAZIMUTHS,&
!                         MU_VALUES, PHIVAL, UP_FLUX, DOWN_FLUX,&
!                         UP_RAD, DOWN_RAD, UP, DOWN)
  

end subroutine DO_CALC

subroutine DO_CALC1(r0, r1, npts, wl, mre, mim, gamma, dens, hpbl, taua, numazim, galbedo, nmu, OUT_FILE) BIND(C)
  use miev0mod
  use iso_c_binding
  implicit none
  real(c_double), intent(in), value     ::  hpbl, taua, galbedo, r0, r1, wl, mre, mim, gamma, dens
  integer(c_int), intent(in), value     ::  numazim, nmu, npts
  
  character(1), intent(in)              ::  OUT_FILE
  type(SDistParams)                     ::  params
  external DO_CALC

  print *, OUT_FILE
  params%r0 = r0
  params%r1 = r1
  params%npts = npts
  params%wl = wl
  params%midx = DCMPLX(mre, mim)
  params%gamma = gamma
  params%dens = dens
  
  call DO_CALC(params, hpbl, taua, numazim, galbedo, nmu, OUT_FILE)

  
end subroutine DO_CALC1
