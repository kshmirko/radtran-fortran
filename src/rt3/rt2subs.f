      SUBROUTINE OUTPUT_FILE (NSTOKES, NUMMU, AZIORDER,
     .                    SRC_CODE, LAYER_FILE, OUT_FILE,
     .                    QUAD_TYPE, DELTAM, DIRECT_FLUX, DIRECT_MU,
     .                    GROUND_TEMP, GROUND_TYPE,
     .                    GROUND_ALBEDO, GROUND_INDEX,
     .                    SKY_TEMP, WAVELENGTH, UNITS, OUTPOL,
     .                    NUM_LAYERS, HEIGHT,
     .                    NOUTLEVELS, OUTLEVELS, NUMAZIMUTHS,
     .                    MU_VALUES, UP_FLUX, DOWN_FLUX,
     .                    UP_RAD, DOWN_RAD)
      INTEGER  NSTOKES, NUMMU, NUMAZI, AZIORDER, SRC_CODE, NUM_LAYERS
      INTEGER  NOUTLEVELS, OUTLEVELS(*), NUMAZIMUTHS
      REAL*8   GROUND_TEMP, GROUND_ALBEDO
      REAL*8   SKY_TEMP, WAVELENGTH
      REAL*8   DIRECT_FLUX, DIRECT_MU
      REAL*8   HEIGHT(NUM_LAYERS+1)
      REAL*8   MU_VALUES(NUMMU)
      REAL*8   UP_FLUX(NSTOKES,NOUTLEVELS)
      REAL*8   DOWN_FLUX(NSTOKES,NOUTLEVELS)
      REAL*8   UP_RAD(NSTOKES,NUMMU,AZIORDER+1,NOUTLEVELS)
      REAL*8   DOWN_RAD(NSTOKES,NUMMU,AZIORDER+1,NOUTLEVELS)
      COMPLEX*16  GROUND_INDEX
      CHARACTER*(*) LAYER_FILE, OUT_FILE
      CHARACTER QUAD_TYPE*1, DELTAM*1, UNITS*1, OUTPOL*2, GROUND_TYPE*1
      CHARACTER*32 QUAD_NAME, UNITS_NAME, GROUND_NAME
      CHARACTER*64 FORM1
      INTEGER  I, J, K, L, LI, M, N
      REAL*4   OUT(4), PHI, PHID, PI
      PARAMETER (PI=3.1415926535897932384D0)
      external CONVERT_OUTPUT

      N = NUMMU*(AZIORDER+1)*NOUTLEVELS
      CALL CONVERT_OUTPUT (UNITS, OUTPOL, NSTOKES, N, 
     .                     WAVELENGTH, 0, UP_RAD)
      CALL CONVERT_OUTPUT (UNITS, OUTPOL, NSTOKES, N, 
     .                     WAVELENGTH, 0, DOWN_RAD)
      CALL CONVERT_OUTPUT (UNITS, OUTPOL, NSTOKES, NOUTLEVELS, 
     .                     WAVELENGTH, 1, UP_FLUX)
      CALL CONVERT_OUTPUT (UNITS, OUTPOL, NSTOKES, NOUTLEVELS, 
     .                     WAVELENGTH, 1, DOWN_FLUX)

      NUMAZI = 2*AZIORDER+1
      IF (NSTOKES .LE. 2) NUMAZI = AZIORDER+1
      QUAD_NAME = 'GAUSSIAN'
      IF (QUAD_TYPE .EQ. 'D')  QUAD_NAME = 'DOUBLEGAUSS'
      IF (QUAD_TYPE .EQ. 'L')  QUAD_NAME = 'LOBATTO'
      IF (QUAD_TYPE .EQ. 'E')  QUAD_NAME = 'EXTRA-ANGLES'
      UNITS_NAME = 'WATTS/(M^2 MICRON STER)'
      IF (UNITS .EQ. 'T') UNITS_NAME = 'KELVINS - EBB'
      IF (UNITS .EQ. 'R') UNITS_NAME = 'KELVINS - RJ'
      GROUND_NAME = 'LAMBERTIAN'
      IF (GROUND_TYPE .EQ. 'F')  GROUND_NAME = 'FRESNEL'

      OPEN (UNIT=3, FILE=OUT_FILE, STATUS='UNKNOWN')

C           Output the parameters
      WRITE (3,'(A,I3,A,I3,A,I3,A,I1)')
     .                'C  NUMMU=', NUMMU,  '  NUMAZI=',NUMAZI,
     .                '  AZIORDER=',AZIORDER, '  NSTOKES=',NSTOKES
      WRITE (3,'(A,A32,A,A1)')
     .                'C  LAYER_FILE=',    LAYER_FILE,
     .                '   DELTA-M=',DELTAM
      WRITE (3,'(A,I1,A,A16)')
     .                'C  SRC_CODE=',      SRC_CODE,
     .                '   QUAD_TYPE=',     QUAD_NAME
      IF (SRC_CODE .EQ. 1 .OR. SRC_CODE .EQ. 3) THEN
          WRITE (3,'(A,E11.5,A,F8.6)')
     .                'C  DIRECT_FLUX=',   DIRECT_FLUX,
     .                '   DIRECT_MU=',     DIRECT_MU
      ENDIF
      WRITE (3,'(A,F8.2,A,A16)')
     .                'C  GROUND_TEMP=',   GROUND_TEMP,
     .                '   GROUND_TYPE=',   GROUND_NAME
      IF (GROUND_TYPE(1:1) .EQ. 'F') THEN
          WRITE (3,'(A,2F9.4,A,F8.2)')
     .                'C  GROUND_INDEX=',  GROUND_INDEX,
     .                '   SKY_TEMP=',      SKY_TEMP
      ELSE
          WRITE (3,'(A,F8.5,A,F8.2)')
     .                'C  GROUND_ALBEDO=', GROUND_ALBEDO,
     .                '   SKY_TEMP=',      SKY_TEMP
      ENDIF
      WRITE (3,'(A,E12.6)') 'C  WAVELENGTH=',    WAVELENGTH
      WRITE (3,'(A,A25,A,A2)') 'C  UNITS='     ,    UNITS_NAME,
     .                '   OUTPUT_POLARIZATION=', OUTPOL  


      IF (UNITS(1:1) .EQ. 'T') THEN
          FORM1 = '(F8.3,1X,F5.1,1X,F8.5,4(1X,F7.2),:)'
      ELSE
          FORM1 = '(F8.3,1X,F5.1,1X,F8.5,4(1X,E13.6),:)'
      ENDIF
 
      IF (OUTPOL .EQ. 'VH') THEN
        WRITE (3,'(A)') 
     .    'C    Z      PHI     MU    FLUX/RADIANCE (V,H,U,V)'
      ELSE
        WRITE (3,'(A)') 
     .    'C    Z      PHI     MU    FLUX/RADIANCE (I,Q,U,V)'
      ENDIF
 
      DO L = 1, NOUTLEVELS
        LI = OUTLEVELS(L)
C               Output fluxes at this level
        WRITE (3,FORM1) HEIGHT(LI), 0., -2.0,
     .        (SNGL(UP_FLUX(I,L)),I=1,NSTOKES)
        WRITE (3,FORM1) HEIGHT(LI), 0., +2.0,
     .        (SNGL(DOWN_FLUX(I,L)),I=1,NSTOKES)
 
C               For each azimuth and zenith at this level sum the Fourier
C               azimuth series appropriate for the particular Stokes parameter
C               and output the radiance.
        DO K = 1, NUMAZIMUTHS
          IF (NUMAZIMUTHS .EQ. 1) THEN
            PHID = 0.0
          ELSE
            PHID = 180.0*FLOAT(K-1)/(NUMAZIMUTHS-1)
          ENDIF
          PHI = PI*PHID/180.0
C               Output upwelling radiance: -1 < mu < 0
          DO J = NUMMU, 1, -1
            DO I = 1, NSTOKES
              OUT(I) = 0.0
              DO M = 0, AZIORDER
                IF (I .LE. 2) THEN
                  OUT(I) = OUT(I) + COS(M*PHI)*UP_RAD(I,J,M+1,L)
                ELSE
                  OUT(I) = OUT(I) + SIN(M*PHI)*UP_RAD(I,J,M+1,L)
                ENDIF
              ENDDO
            ENDDO
            WRITE (3,FORM1) HEIGHT(LI), PHID, -MU_VALUES(J),
     .                      (OUT(I),I=1,NSTOKES)
          ENDDO
C               Output downwelling radiance: 0 < mu < 1
          DO J = 1, NUMMU
            DO I = 1, NSTOKES
              OUT(I) = 0.0
              DO M = 0, AZIORDER
                IF (I .LE. 2) THEN
                  OUT(I) = OUT(I) + COS(M*PHI)*DOWN_RAD(I,J,M+1,L)
                ELSE
                  OUT(I) = OUT(I) + SIN(M*PHI)*DOWN_RAD(I,J,M+1,L)
                ENDIF
              ENDDO
            ENDDO
            WRITE (3,FORM1) HEIGHT(LI), PHID,  MU_VALUES(J),
     .                      (OUT(I),I=1,NSTOKES)
          ENDDO
        ENDDO
      ENDDO

      CLOSE (3)

      RETURN
      END





      SUBROUTINE CONVERT_OUTPUT (UNITS, OUTPOL, NSTOKES, NOUT,
     .                           WAVELEN, FLUXCODE, OUTPUT)
C       Converts the output radiance or flux arrays to VH polarization
C     and effective blackbody temperature if desired.  OUTPOL='VH'
C     converts the polarization basis of the first two Stokes parameters
C     to vertical/horizontal polarization.  If UNITS='T' the radiance is
C     converted to effective blackbody brightness temperature, and if
C     UNITS='R' the radiance is converted to Rayleigh-Jeans brightness
C     temperature.  If the output is flux then FLUXCODE=1, and the flux 
C     is divided by pi before converting to brightness temperature.
      INTEGER NSTOKES, NOUT, FLUXCODE
      REAL*8  WAVELEN, OUTPUT(NSTOKES,NOUT)
      CHARACTER UNITS*1, OUTPOL*2
      INTEGER I, J
      REAL*8  IV, IH, RAD, TEMP

      DO J = 1, NOUT      
C           Convert to Vertical and Horizontal polarization if desired
        IF (OUTPOL .EQ. 'VH') THEN
          IV = 0.5*(OUTPUT(1,J) + OUTPUT(2,J))
          IH = 0.5*(OUTPUT(1,J) - OUTPUT(2,J))
          OUTPUT(1,J) = IV
          OUTPUT(2,J) = IH
        ENDIF
C           Convert to brightness temperature
        IF (UNITS .EQ. 'T' .OR. UNITS .EQ. 'R') THEN
          DO I = 1, NSTOKES
            RAD = OUTPUT(I,J)
            IF (OUTPOL .EQ. 'VH' .AND. I .LE. 2)  RAD = 2.0*RAD
            IF (FLUXCODE .EQ. 1)  RAD = RAD/ACOS(-1.0)
            IF (UNITS .EQ. 'R') THEN
              TEMP = RAD * WAVELEN**4 * 1.4388D4/1.1911D8
            ELSE
              IF (RAD .GT. 0.0) THEN
                TEMP = 1.4388D4 /
     .            (WAVELEN*DLOG(1.0+ 1.1911D8/(RAD*WAVELEN**5)))
              ELSE IF (RAD .EQ. 0.0) THEN
                TEMP = 0.0D0
              ELSE
                TEMP = -1.4388D4 /
     .            (WAVELEN*DLOG(1.0+ 1.1911D8/(-RAD*WAVELEN**5)))
              ENDIF
            ENDIF
            OUTPUT(I,J) = TEMP
          ENDDO
        ENDIF
      ENDDO
      RETURN
      END