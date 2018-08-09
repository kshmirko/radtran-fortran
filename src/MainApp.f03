program MainApp
  use miev0mod
  use rayleigh  
  use aerosol
  use atmos
  implicit none
  
  real(kind=dp) ::  ext, sca, asy, hpbl, taua, extm, exta, extt,&
                    omegat, omegaa, scat
  integer mxmdm, ierr
  real(kind=dp), allocatable  :: pmom(:,:), evans(:,:), evanstot(:,:)
  type(SDistParams)           :: params
  character(len=13)           :: fname
  integer I, J

  INTEGER   MAXV, MAXA, MAXLAY, NLAYS
  PARAMETER (MAXV=128, MAXA=32, NLAYS=40)
  PARAMETER (MAXLAY=200)

  INTEGER   NSTOKES, NUMMU, AZIORDER
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



  
  ! Для начала выделим память под хранение матриц с коэффициентами 
  ! Лежандра
  mxmdm = 60
  call init_evansmol(mxmdm)

  ! allocate memory for matrices
  allocate(pmom(0:mxmdm,4), evans(0:mxmdm, 6), evanstot(0:mxmdm, 6))
  
  pmom = 0.0_dp
  ext = 0.0_dp
  sca = 0.0_dp
  asy = 0.0_dp

  ! Инициализация параметров
  params%r0 = 0.1_dp
  params%r1 = 1.0_dp
  params%npts = 101
  params%wl = 0.750_dp
  params%midx = DCMPLX(1.4_dp, 0.01_dp)
  params%gamma = -3.5_dp
  params%dens = 300

  ! Расчет рассеяния
  call miesdist1(params, pmom, ext=ext, sca=sca, asy=asy, ierr=ierr)
  omegaa = sca/ext

  if (ierr/=0) then
    stop 'Ошибка в подпрограмме miesdist1'
  end if
  
  print '(A)', 'Оптические параметры заданного распределения при '
  print '(A)', 'условии единичной концентрации'
  print '(4A12)', 'EXT', 'SCA', 'ASY', 'OMEGA'
  print '(4E12.4)', ext, sca, asy, sca/ext

  ! Матрица Эванса имеет следующие элеметы
  ! P1  P2  0   0  
  ! P2  P5  0   0
  ! 0   0   P3  P4
  ! 0   0  -P4  P6
  print '(/,5X,A,/)', 'Матрица Эванса:'
  print '(4A4)', 'P1', 'P2', '0', '0'
  print '(4A4)', 'P2', 'P5', '0', '0'
  print '(4A4)', '0', '0', 'P3', 'P4'
  print '(4A4)', '0', '0', '-P4', 'P6'

  call wiscombe2evans(pmom, evans)
  !print '(/,6A10)', 'P1', 'P2', 'P3', 'P4', 'P5', 'P6'
  !DO I=0, mxmdm
  !  print '(6F10.6)', evans(I,:)
  !END DO

  hpbl = 3000.0_dp
  ! Вычисление АОТ через концентрацию частиц и их оптические свойства
  taua = ext*hpbl*1.0d6*params%dens
  

  print '(A,F5.3,A,F5.3)', 'tau_m(',params%wl,') = ', tau_m(params%wl)
  print '(A,F5.3,A,F5.3)', 'tau_a(',params%wl,') = ', taua
    

  open(200, FILE='atmoslay.lay', status='replace')
  ! Вычитсяем оптические характеристики для слоев
  ! Сохраняем матрицы рассеяния
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
    write(200, '(F6.2,F7.2,E12.5,2X,A1,A13,A1)') HEIGHT(I), 0.0, 0.0, "'",fname,"'"
    print '(F7.1, 3E10.2)', HEIGHT(I), extm, exta, extt
  END DO

  close(200)

  deallocate(pmom,  evans, evanstot)
  
  ! Начало радиационного расчета
  MAX_DELTA_TAU = 1.0E-6
  NSTOKES = 2
  NUMMU = 32
  AZIORDER = 2
  SRC_CODE = 1
  QUAD_TYPE = 'G'
  DELTAM = 'Y'
  DIRECT_FLUX = 1.0_dp
  THETA = 10.0_dp
  DIRECT_MU = DABS(DCOS(0.017453292D0*(THETA)))
  GROUND_TEMP = 0.0_dp
  GROUND_TYPE = 'L'
  GROUND_ALBEDO = 0.0_dp
  GROUND_INDEX = DCMPLX(1.0_dp, 0.0_dp)
  SKY_TEMP = 0.0_dp
  WAVELENGTH = params%wl
  NUM_LAYERS = NLAYS
  NOUTLEVELS = 1
  OUTLEVELS(1) = NLAYS

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

  OUT_FILE = 'rt3.out'
  NUMAZIMUTHS = 2
  LAYER_FILE = 'NONE'
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

  

end program