      PROGRAM MAIN 
!     THIS PROGRAM CALCULATES THE SCATTERING OF A PLANE ELECTROMAGNETIC
!     WAVE BY A LOSSY DIELECTRIC BODY.
!     THE BODY IS ASSUMED TO BE A SPHEROID OF REVOLUTION.
!     THE BODY IS IMMERSED IN A LOSSY DIELECTRIC MEDIUM.                                                                                                         
!     EXP(-JWT) TIME CONVENTION IS USED.                                
!     PERFORM THE NUMERICAL INTEGRATION AND FILL THE A,B,C,D MATRICES   
!     FOR THE OUTER SURFACE AND X,Y MATRICES FOR THE INNER SURFACE.  
!     THE MATRICES ARE STORED IN THE ARRAYS A,B,C,D,X,Y.
!-----------------------------------------------------------------------
      ! declare variables
      integer ib,opt,flag_finished
      real tempt,alfa,beta,inoblate 
      real wavelength 
      DIMENSION TEMP(16,8,4) 
      DIMENSION CLRMTX(25600),CLRTOT(724),RESULT(16,8) 
      DIMENSION BUFF(5) 
      CHARACTER*99 FILENAME1,FILENAME2,FILENAME3,FILENAME4,tablefile 
      CHARACTER*99 FOLDERNAME,INPUTFILE,INFOFILE
      CHARACTER*2  CTEMPT,CWAVEL 
      CHARACTER*1  COPT,CWAVEL1 
      ! The COMMON block
      COMMON DTR,RTD,CPI 
      COMMON /BDYCOM/ DCNR,DCNI,CKPRR,CKPRI,CKR,DCKR,CONK,AOVRB,SIGMA,IB 
      COMMON /CMVCOM/ NM,KMV,CMI(20),CMV,CM2,TWM,PRODM 
      COMMON/VARIAM/QEM,B89 
      COMMON/MTXCOM/NRANK,NRANKI,A(80,80,2),B(80,80,2),CMXNRM(80) 
      COMMON /UVCCOM/ANGINC,ACANS(181,2,2),UANG(181),RTSFCT,DLTANG,NUANG 
      COMMON/ROWCOL/IROW,IROW1,CROW,CROWM,CROW1,ICOL,ICOL1,CCOL,CCOLM,  &
     &CCOL1,CRIJ,CRSSIJ                                                 
      COMMON/BRIN/ML,IJK,HEIGHT 
      COMMON/ENDPNT/EPPS(4),NSECT 
      COMMON/SCATT/SUMM1,SUMM2 
      COMMON/RAVI/IOPT 
      COMMON/VIVEK2/DPART
      EQUIVALENCE (A(1,1,1),CLRMTX(1)),(ACANS(1,1,1),CLRTOT(1)) 
  !============================================================================
   ! Finished the declaration of the variables
   !============================================================================
      ! read the FOLDER CITE 
      read(*,*),FOLDERNAME
      FILENAME1= TRIM(FOLDERNAME)//'/out1_tmat'
      FILENAME2= TRIM(FOLDERNAME)//'/out2_tmat'
      FILENAME3= TRIM(FOLDERNAME)//'/scat_parall'
      FILENAME4= TRIM(FOLDERNAME)//'/scat_perpn'
      INPUTFILE= TRIM(FOLDERNAME)//'/input'
      INFOFILE= TRIM(FOLDERNAME)//'/raininfo.dat'
      ! FILENAME2=FOLDERNAME+'/out2_tmat'
      ! FILENAME3=FOLDERNAME+'/scat_parall'
      ! FILENAME4=FOLDERNAME+'/scat_perpn'
   !***************** SETTING THE CONSTANT PARAMETER *********************
      CPI = ACOS(-1.)
      DTR = CPI/180
      RTD = 180/CPI
      NSECT=2
      !***************** OPEN THE FILES (INPUT/OUTPUT) *********************
      !Read the input file
      OPEN(107,FILE=INPUTFILE)
      OPEN(108,FILE=INFOFILE)
      !open the output file
      OPEN(109,FILE=FILENAME1)
      OPEN(110,FILE=FILENAME2)
      OPEN(111,FILE=FILENAME3)
      OPEN(112,FILE=FILENAME4)
      !============================================================================
      ! Finished opening the files
      !============================================================================
      write(*,*)'!!!!   START MAKE THE FILE FOR THE TMAT   !!!!!'
      !*** READ THE INPUT (SETTING THE RADAR INFORMATION AND SOME VARAIBLES)
      ! write(*,*)'Please input the IB = 8 or 9 (SYMMETRY CODE)...'
      read(107,*) IB
      ! write(*,*) 'Please input the Droplet Temperature in C'
      read(107,*) TEMPT
      ! write(*,*) 'Please input the Wavelength in cm'
      read(107,*) WAVELENGTH
      !display the variable from input
      write(*,*) 'IB = ',IB
      write(*,*) 'TEMPT = ',TEMPT,'C'
      write(*,*) 'WAVELENGTH = ',WAVELENGTH,'cm'
   !***************** READ THE raininfo.dat *********************
   !     SET PROGRAM CONSTANTS.
   !     ML is the number of T-matrices to be calculated.
   !    IOPT is the option for the hydrometer.
      READ(108,*)ML,IOPT
      print*, 'ML = ',ML
      IJK = 0 ! IJK is the number of the which line in raininfo.dat
      iflag_finished=0
      ! print *, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      WRITE(109,"(I3,5X,'NUMBER OF T-matrices')")ML 
        
   20 CALL RDDATA(tempt,wavelength,ALAMD,iflag_finished)
      if (iflag_finished.EQ.1 ) GO TO 87
      ! display the OUTPUT from RDDATA
      ! write(*,*) 'NRANK = ',NRANK
      ! write(*,*) 'NM = ',NM
      ! write(*,*) 'CONK = ',CONK
      ! write(*,*) 'ALAMD = ',ALAMD
      ! write(*,*) 'DPART = ',DPART
      ! write(*,*) 'AOVRB = ',AOVRB
      ! write(*,*) 'SIGMA = ',SIGMA
      ! write(*,*) 'DCNR = ',DCNR
      ! write(*,*) 'DCNI = ',DCNI
      ! write(*,*) 'ANGINC = ',ANGINC
      ! STOP
      IFLAG=0 
!     SET THE NUMBER OF SCATTERING ANGLES (NUANG = 181 MAXIMUM AND NUANG-1
!     SHOULD BE DIVISIBLE INTO 180 BY A WHOLE NUMBER.
      NUANG= 2
      DLTANG=FLOAT(180/(NUANG-1)) 
      UANG(1) = ANGINC
      DO 30 I=2,NUANG 
      UANG(I) = UANG(I-1)+DLTANG
   30 CONTINUE
!     CLEAR THE ACCUMULATING ANSWER REGISTER (USED IN ADDPRC).          
      DO 40 J=1,724 
      CLRTOT(J)=0.0 
   40 END DO 
      SUMM1=0.0 
      SUMM2=0.0 
      RTSFCT=8.0/CONK 
!     SET MULTIPLIER B89 DEPENDENT ON IB VALUE (SYMMETRY INDICATOR).    
      B89 = 1.0 
      IF(IB.EQ.8) B89=2.0 
      BDYFCT=1.0 
!     SET UP A LOOP FOR EACH M VALUE.                                   
      DO 900 IM = 1,NM 
!      print*,'MN=',IM                                                  
!     SET M DEPENDENT VARIABLES.                                        
      CMV = CMI(IM) 
      KMV = CMV 
      CM2 = CMV**2 
      PRODM=1.0 
      IF(KMV.GT.0) GO TO 44 
      EM=1.0 
      GO TO 60 
   44 EM=2.0 
      QUANM = CMV 
      DO 52 IFCT = 1,KMV 
      QUANM=QUANM+1.0 
      PRODM=QUANM*PRODM/2.0 
   52 END DO 
   60 QEM=-2.0/EM 
      TWM=2.0*CMV 
!     INITIALIZE ALL MATRIX AREAS TO ZERO                               
      DO 80 I=1,25600 
      CLRMTX(I)=0.0 
   80 END DO 
!     SET UP A LOOP FOR EACH ROW OF THE MATRICES.                       
      CROW=0.0 
      CROWM = CMV 
      DO 600 IROW = 1,NRANK 
      IROW1 = IROW+NRANK 
      CROW=CROW+1.0 
      CROWM=CROWM+1.0 
      CROW1=CROW+1.0 
!     SET UP A LOOP FOR EACH COLUMN OF THE MATRICES.                    
      CCOL=0.0 
      CCOLM = CMV 
      DO 400 ICOL = 1,NRANK 
      ICOL1 = ICOL+NRANK 
      CCOL=CCOL+1.0 
      CCOLM=CCOLM+1.0 
      CCOL1=CCOL+1.0 
!                                                                       
!     CALCULATE MATRICES A,B ASSOCIATED WITH THE OUTER SURFACE,         
!     FOLLOWING NOTATION BY PETERSON AND STROM, Q1(OUT,RE) IS STORED    
!     IN A, Q1(RE,RE) IS STORED IN   B.                                 
!     ALL MATRICES ARE TRANSPOSED IN THE FOLLOWING CODE.                
!                                                                       
!     PERFORM INTEGRATION USING A SEQUENCE OF 31                        
!     POINT EXTENDED GAUSS-TYPE QUADRATURE FORMULAE.                    
!     RESULT(16,K) CONTAINS THE VALUES OF THE INTEGRALS .               
!     THERE ARE 16 INTEGRATIONS TO                                      
!     BE PERFORMED FOR EACH LOOPING THRO IROW AND ICOL. THESE CORRESPOND
!     TO 4 SUB-MATRIX ELEMENTS FOR EACH OF THE 2 MATRICES (A,B)         
!     AND ASSICIATED REAL AND IMAGINARY PARTS.                          
!                                                                       
      NSECT1=NSECT-1 
      DO 301 J=1,NSECT1 
       JS=J 
        ITH=0 
!      print*,'BEFORE QUAD, K =',K,',NPTS=',NPTS                        
      CALL QUAD(EPPS(J),EPPS(J+1),K,RESULT,NPTS,ITH,JS) 
!      print*,' AFTER QUAD, K =',K,',NPTS=',NPTS                        
      DO 401 I=1,16 
      TEMP(I,K,J)=RESULT(I,K) 
  401 END DO 
  301 END DO 
      DO 501 J=1,NSECT1 
      A(ICOL,IROW1,1)=TEMP(1,K,J)+A(ICOL,IROW1,1) 
      A(ICOL,IROW1,2)=TEMP(2,K,J)+A(ICOL,IROW1,2) 
      B(ICOL,IROW1,1)=TEMP(3,K,J)+B(ICOL,IROW1,1) 
      B(ICOL,IROW1,2)=TEMP(4,K,J)+B(ICOL,IROW1,2) 
      A(ICOL1,IROW,1)=TEMP(5,K,J)+A(ICOL1,IROW,1) 
      A(ICOL1,IROW,2)=TEMP(6,K,J)+A(ICOL1,IROW,2) 
      B(ICOL1,IROW,1)=TEMP(7,K,J)+B(ICOL1,IROW,1) 
      B(ICOL1,IROW,2)=TEMP(8,K,J)+B(ICOL1,IROW,2) 
      A(ICOL1,IROW1,1)=TEMP(9,K,J)+A(ICOL1,IROW1,1) 
      A(ICOL1,IROW1,2)=TEMP(10,K,J)+A(ICOL1,IROW1,2) 
      B(ICOL1,IROW1,1)=TEMP(11,K,J)+B(ICOL1,IROW1,1) 
      B(ICOL1,IROW1,2)=TEMP(12,K,J)+B(ICOL1,IROW1,2) 
      A(ICOL,IROW,1)=TEMP(13,K,J)+A(ICOL,IROW,1) 
      A(ICOL,IROW,2)=TEMP(14,K,J)+A(ICOL,IROW,2) 
      B(ICOL,IROW,1)=TEMP(15,K,J)+B(ICOL,IROW,1) 
      B(ICOL,IROW,2)=TEMP(16,K,J)+B(ICOL,IROW,2) 
  501 END DO 
      IF(IFLAG.EQ.0) WRITE(110,101) NPTS 
  101 FORMAT(/,10X,'NO. OF GAUSS POINTS USED',I5) 
      IFLAG=IFLAG+1 
  400 END DO 
!     CALCULATE THE NORMALIZATION FACTOR (USED IN ADDPRC).              
      CKROW = IROW 
      IF(KMV.GT.0) GO TO 426 
      FCTKI=1.0 
      GO TO 440 
  426 IF(IROW.GE.KMV) GO TO 430 
      CMXNRM(IROW)=1.0 
      GO TO 600 
  430 IBFCT = IROW-KMV+1 
      IEFCT = IROW+KMV 
      FPROD = IBFCT 
      FCTKI=1.0 
      DO 432 LFCT = IBFCT,IEFCT 
      FCTKI = FCTKI*FPROD 
      FPROD=FPROD+1.0 
  432 END DO 
  440 CMXNRM(IROW) = 4.0  *CKROW*(CKROW+1.0  )*FCTKI/(EM*(2.0  *CKROW+1.&
     &0  ))                                                             
  600 END DO 
      NN=2*NRANK 
!     PROCESS COMPUTED MATRICES                                         
!                                                                       
!     CALCULATE T(1,1)=B(INVERSE)*A AND STORE IN B.                     
      CALL PRCSSM(A,B,NRANK,NRANKI) 
      WRITE(109,11) CMI(IM) 
   11  FORMAT(2X,F10.5) 
       DO  7  ICOM1=1,2 
       DO  7  ICOM2=1,NN 
       WRITE(109,9) (B(ICOM2,ICOM3,ICOM1),ICOM3=1,NN) 
    9  FORMAT(2X,8E15.7) 
    7  CONTINUE 
!      CALL PRINTM(X,NN,40)                                             
      CALL ADDPRC(IM) 
  900 END DO 
      GO TO 20 

      !------------------------------------------------------------------
      ! finish the main program
   87 write(*,*) 'Finished the main program!!!!!!!!!!!!'  
      END PROGRAM MAIN
!============================================================================
! FINISHED THE MAIN PROGRAM
!============================================================================

!============================================================================
!  SUBROUTINE LIST
!----------------------------------------------------------------------------

      
!*****************************************************************
! THE SUBROUTINE RDDATA is A PROGRAM TO READ INPUT DATA FOR THE SCATTERING PROGRAM.
! INCLUDE THE SUB: EPSLON
!     OUTPUT: NRANK,NM,CONK,ALAMD,DPART,AOVRB,SIGMA,DCNR,DCNI
      SUBROUTINE RDDATA(tempt,wavelength,ALAMD,iflag_finished)
      COMMON DTR,RTD,CPI
      COMMON/GAUSS/KGAUSS
      COMMON/MTXCOM/NRANK,NRANKI,A(80,80,2),B(80,80,2),CMXNRM(80)
      COMMON/CMVCOM/ NM,KMV,CMI(20),CMV,CM2,TWM,PRODM
      COMMON/BDYCOM/ DCNR,DCNI,CKPRR,CKPRI,CKR,DCKR,CONK,AOVRB,SIGMA,IB
      COMMON/UVCCOM/ ANGINC,ACANS(181,2,2),UANG(181),RTSFCT,DTLANG,NUANG
      COMMON/BRIN/ML,IJK,HEIGHT
      COMMON/ENDPNT/EPPS(4),NSECT
      COMMON/VIVEK2/DPART
      COMMON/RAVI/IOPT
      real tempt,alfa,beta,inoblate
      integer opt,iflag_finished
      DIMENSION EPDEG(4)
      real RMINOR,RMAJOR
      !     CHECK TO SEE IF ALL DATA SETS ARE READ.                           
      IF(IJK.EQ.ML) GO TO 200                                                                    
      !     READ NECESSARY INPUT DATA.                                                                                                          
      !     NM = NUMBER OF M VALUES,NRANK = N VALUE(MATRIX ORDER),  
      !     NSECT = NUMBER OF SECTIONS IN THE BODY,
      !     IB = SYMMETRY CODE  IB = 8 FOR MIRROR SYMMETRY ABOUT THETA = 90 DEGREES
      !     IB = 9 FOR GENERAL SHAPED BODY,ANGINC = ANGLE OF INCIDENCE OF THE INCIDENT WAVE.     
      READ(108,*) NM,NRANK,DPART,AOVRB,ANGINC
      print*,'NM=',NM,',NRANK=',NRANK,',DPART=',DPART,',AOVRB=',AOVRB,',ANGINC=',ANGINC
      ALAMD=WAVELENGTH
      !  ***  Routine returns DCNR and DCNI for a Wavelength, for WATER.      
      CALL EPSLON(ALAMD,TEMPT,DCNR,DCNI)                                                                  
      ! Non-rain branches in the original NASA/Ravi code were controlled by
      ! hand-edited dielectric constants. The Python wrapper now applies the
      ! documented constants for dry ice and low-density dry aggregate tests.
      ! These are executable assumptions, not a full validated hydrometeor
      ! microphysics model.
      IF (IOPT .EQ. 2) THEN
            DCNR=3.168351
            DCNI=0.02492
      ELSEIF (IOPT .EQ. 3) THEN
            DCNR=3.168351
            DCNI=0.02492
      ELSEIF (IOPT .EQ. 8) THEN
            DCNR=1.336817
            DCNI=0.000080937403
      ELSEIF (IOPT .EQ. 9) THEN
            DCNR=1.336817
            DCNI=0.000080937403
      ELSEIF (IOPT .EQ. 10) THEN
            DCNR=1.97194
            DCNI=0.000271361
      ELSEIF (IOPT .EQ. 11) THEN
            DCNR=2.78180
            DCNI=0.000612490
      ENDIF
      IF ( AOVRB .GE. 1. ) AOVRB = 0.999999999 
      IF ( AOVRB .LE. 0.0 ) AOVRB = 0.0000000000001 

      ! CALCULATE THE diameters (RMINOR,RMAJOR)
      RMINOR=(AOVRB**(2./3.))*0.5*DPART 
      RMAJOR=(AOVRB**(-1./3.))*0.5*DPART 
      ! write(*,'(4f15.5)') DPART,AOVRB,RMAJOR,RMINOR 

      ! Calculate the varaibles (NRANKI,TEMPT,ALAMD,DPART,AOVRB,SIGMA,CONK)
      NRANKI = NRANK+1 
      SIGMA=2.*ALAMD/(CPI*100.)                                                                                                                                                
      CONK=CPI*DPART/(ALAMD*(AOVRB**(1./3.))) 

      ! display the information to out1_tmat
      WRITE(109,"(I10,5X,'NRANK')")NRANK 
      WRITE(109,"(I10,5X,'NM')")NM 
      WRITE(109,"(E15.7,5X,'CONKB')")CONK 
      WRITE(109,"(E15.7,5X,'WAVELENGTH IN CM')")ALAMD 
      WRITE(109,"(E15.7,5X,'EQUIVALENT PARTICLE DIAMETER IN CM')")DPART 
      WRITE(109,"(E15.7,5X,'AOVRB')")AOVRB 
      WRITE(109,"(E15.7,5X,'SIGMA')")SIGMA 
      WRITE(109,"(E15.7,5X,'REAL PART OF DIELECTRIC CONSTANT')")DCNR 
      WRITE(109,"(E15.7,5X,'IMAGINARY PART OF DIELECTRIC CONS')")DCNI 
      
      ! display the information to out2_tmat                                                                                                                                                              
      WRITE(110,556) 
      IF(IOPT .EQ. 1) WRITE(110,FMT ='(20X,"RAINDROPS")') 
      IF(IOPT .EQ. 2) WRITE(110,FMT ='(20X,"FROZEN RAINDROPS (FOR SPONGY GRAUPEL)")') 
      IF(IOPT .EQ. 3) WRITE(110,FMT ='(20X,"DRY HAIL")') 
      IF(IOPT .EQ. 4) WRITE(110,FMT ='(20X,"SPONGY HAIL")') 
      IF(IOPT .EQ. 5) WRITE(110,FMT ='(20X,"WET HAIL")') 
      IF(IOPT .EQ. 6) WRITE(110,FMT ='(20X,"NEEDLES CRYSTALS")') 
      IF(IOPT .EQ. 7) WRITE(110,FMT ='(20X,"PLATES CRYSTALS")') 
      IF(IOPT .EQ. 8) WRITE(110,FMT ='(20X,"AGGREGATES SPONGY SNOW")') 
      IF(IOPT .EQ. 9) WRITE(110,FMT ='(20X,"AGGREGATES RHO OF 0.2")') 
      IF(IOPT .EQ. 10) WRITE(110,FMT ='(20X,"AGGREGATES RHO OF 0.5")') 
      IF(IOPT .EQ. 11) WRITE(110,FMT ='(20X,"AGGREGATES RHO OF 0.8")') 
      WRITE(110,FMT ='(8X,F5.3,7X,"DEQ IN CM")')DPART 
      WRITE(110,FMT ='(F15.5,5X,"WAVELENGTH IN CM")')ALAMD 
      WRITE(110,FMT ='(F15.5,5X,"TEMPERATURE IN CENTIGRADE")')TEMPT 
      WRITE(110,FMT ='(8X,I5,7X,"CASES")')NM 
      WRITE(110,FMT ='(8X,I5,7X,"MATRIX RANK")')NRANK 
      WRITE(110,FMT ='(8X,I5,7X,"SECTIONS")')NSECT 
      WRITE(110,FMT ='(8X,I5,7X,"BODY SHAPE")')IB 
      WRITE(110,FMT ='(F12.2,8X,"ANGLE OF INCIDENCE IN DEGREES")')ANGINC                                                                 
      WRITE(110,FMT ='(F15.5,5X,"AOVRB")')AOVRB 
      WRITE(110,FMT ='(F15.5,5X,"CONK")')CONK 
      WRITE(110,FMT ='(F15.5,5X,"SIGMA")')SIGMA 
      WRITE(110,FMT ='(E15.5,5X,"REAL PART OF DIELECTRIC CONST")')DCNR 
      WRITE(110,FMT ='(E15.5,5X,"IMAG PART OF DIELECTRIC CONST")')DCNI 
      !                
                                                             
      !   CMI = M VALUES (ONLY M=1 IS REQUIRED FOR ANGINC=0,BUT 
      !     M=0,1,2, ... IS REQUIRED FOR GENERAL INCIDENCE ANGLE).
      DO   I=1,NM
            CMI(I)=FLOAT(I-1)
      ENDDO
      IF(NM.EQ.1) THEN 
            CMI(1)=1.0
      ENDIF
      PRODM=1.0     
      ! SET PARAMETER IN GAUSS QUADRATURE ROUTINE.                                                                       
      KGAUSS=5 
      CALL CALENP 
      DO I=1,NSECT 
            EPDEG(I)=RTD*EPPS(I) 
      ENDDO
      WRITE(110,148) (EPDEG(I),I=1,NSECT) 
      !     JUMP UP DATA SET COUNTER                                          
      IJK=IJK+1 
      RETURN 

      ! if all data sets are read, then go to 200
  200 WRITE (*,201)
      iflag_finished=1
      RETURN 

      556 FORMAT(//) 
      148 FORMAT(1X,'      END POINTS',8E12.4,/(1X,23X,8E12.4)) 
      201 FORMAT(1X,'ALL DATA SETS HAVE BEEN READ') 
      END  SUBROUTINE RDDATA

      !------------------------------------------------------------------
      !*****************************************************************
 

      SUBROUTINE GENER(XXX,T,ITH,JS) 
            COMMON DTR,RTD,CPI 
            COMMON/THTCOM/THETA,SINTH,COSTH 
            COMMON/MTXCOM/NRANK,NRANKI,A(80,80,2),B(80,80,2),CMXNRM(80) 
            COMMON /BDYCOM/ DCNR,DCNI,CKPRR,CKPRI,CKR,DCKR,CONK,AOVRB,SIGMA,IB 
            COMMON/FNCCOM/PNMLLG(81),BSSLSP(81,31,3),CNEUMN(81,31,3),         &
           &BSLKPR(81,31,3),BSLKPI(81,31,3),CNEUMR(81,31,3),CNEUMI(81,31,3)   
            COMMON/ROWCOL/IROW,IROW1,CROW,CROWM,CROW1,ICOL,ICOL1,CCOL,CCOLM,  &
           &CCOL1,CRIJ,CRSSIJ                                                 
            COMMON /CMVCOM/ NM,KMV,CMI(20),CMV,CM2,TWM,PRODM 
            COMMON/VARIAM/QEM,B89 
            DIMENSION RHANKL(80,2),RBSKPR(80,2),RBESSL(80),HANKLR(81),        &
           &HANKLI(81),RHSKPR(80,2)                                           
            DIMENSION T(4,2,2) 
              XXXR=XXX*RTD 
            DO 1101 I=1,4 
            DO 1101 J=1,2 
            DO 1101 K=1,2 
            T(I,J,K)=0.0 
       1101 CONTINUE 
            SQREAL=CROOTR(DCNR,DCNI) 
            SQIMAG=CROOTI(DCNR,DCNI) 
            QSREAL=CDDVDR(1.0  ,0.0  ,SQREAL,SQIMAG) 
            QSIMAG=CDDVDI(1.0  ,0.0  ,SQREAL,SQIMAG) 
            THETA=XXX 
            COSTH=COS(THETA) 
            SINTH=SIN(THETA) 
            SRMSIN=SINTH 
            TEMPP=CPI-THETA 
            IF ( ABS(TEMPP).LT.1.0E-8) SINTH=0.0 
      !     GENERATE THE LEGENDRE POLYNOMIALS.                                
            CALL GENLGP 
      !     EVALUATE KR AND ITS DERIVATIVE AS A FUNCTION OF THETA.            
        348 CALL GENKR 
      !     GENERATE ALL NECESSARY BESSEL AND NEUMANN FUNCTIONS AND THEIR RATI
            IF((IROW.EQ.1).AND.(ICOL.EQ.1)) CALL GENBSL(ITH,JS) 
      !     CALCULATE K*R FOR OUTER SURFACE.                                  
            CKPRR = SQREAL*CKR 
            CKPRI = SQIMAG*CKR 
            ISWT=1 
            IF((IROW.EQ.1).AND.(ICOL.EQ.1))CALL GENBKR(CKPRR,CKPRI,ISWT,      &
           &ITH,JS)                                                           
                                                                              
            DO 347 K=1,NRANKI 
      !V      write(6,'(a3,i4,4e15.7,/)') 'new',k,bslkpr(k,ith,js),           
      !V     1bslkpi(k,ith,js),cneumr(k,ith,js),cneumi(k,ith,js)              
            HANKLR(K)=BSLKPR(K,ITH,JS)-CNEUMI(K,ITH,JS) 
            HANKLI(K)=BSLKPI(K,ITH,JS)+CNEUMR(K,ITH,JS) 
        347  CONTINUE 
            DO 350 K = 1,NRANK 
            RBESSL(K) = BSSLSP(K,ITH,JS)/BSSLSP(K+1,ITH,JS) 
            RHANKL(K,1)=CDDVDR(BSSLSP(K,ITH,JS),CNEUMN(K,ITH,JS),             &
           &BSSLSP(K+1,ITH,JS),CNEUMN(K+1,ITH,JS))                            
            RHANKL(K,2) = CDDVDI(BSSLSP(K,ITH,JS),CNEUMN(K,ITH,JS),           &
           &BSSLSP(K+1,ITH,JS),CNEUMN(K+1,ITH,JS))                            
            RBSKPR(K,1) = CDDVDR(BSLKPR(K,ITH,JS),BSLKPI(K,ITH,JS),           &
           &BSLKPR(K+1,ITH,JS),BSLKPI(K+1,ITH,JS))                            
            RBSKPR(K,2) = CDDVDI(BSLKPR(K,ITH,JS),BSLKPI(K,ITH,JS),           &
           &BSLKPR(K+1,ITH,JS),BSLKPI(K+1,ITH,JS))                            
            BKR = CDMPYR(SQREAL,SQIMAG,RBSKPR(K,1),RBSKPR(K,2)) 
            BKI = CDMPYI(SQREAL,SQIMAG,RBSKPR(K,1),RBSKPR(K,2)) 
            RBSKPR(K,1) = BKR 
            RBSKPR(K,2) = BKI 
            TEMPNR=BSLKPR(K,ITH,JS)-CNEUMI(K,ITH,JS) 
            TEMPNI=BSLKPI(K,ITH,JS)+CNEUMR(K,ITH,JS) 
            TEMPDR=BSLKPR(K+1,ITH,JS)-CNEUMI(K+1,ITH,JS) 
            TEMPDI=BSLKPI(K+1,ITH,JS)+CNEUMR(K+1,ITH,JS) 
            RHSKPR(K,1)=CDDVDR(TEMPNR,TEMPNI,TEMPDR,TEMPDI) 
            RHSKPR(K,2)=CDDVDI(TEMPNR,TEMPNI,TEMPDR,TEMPDI) 
            HKR=CDMPYR(SQREAL,SQIMAG,RHSKPR(K,1),RHSKPR(K,2)) 
            HKI=CDMPYI(SQREAL,SQIMAG,RHSKPR(K,1),RHSKPR(K,2)) 
            RHSKPR(K,1)=HKR 
            RHSKPR(K,2)=HKI 
        350 END DO 
            BR = RBESSL(IROW) 
            HR = RHANKL(IROW,1) 
            HI = RHANKL(IROW,2) 
      !     CALCULATE FREQUENTLY USED VARIABLE COMBINATIONS FOR USE IN A,B    
      !      MATRICES.                                                        
            CRIJ = CROW+CCOL 
            CRSSIJ = CROW*CCOL 
            CMCRCO = CM2-QEM*CRSSIJ*COSTH**2 
            PNR0C0 = PNMLLG(IROW)*PNMLLG(ICOL) 
            PNR0C1 = PNMLLG(IROW)*PNMLLG(ICOL+1) 
            PNR1C0 = PNMLLG(IROW+1)*PNMLLG(ICOL) 
            PNR1C1 = PNMLLG(IROW+1)*PNMLLG(ICOL+1) 
            B1A = CROW*COSTH*PNR1C1-CROWM*PNR0C1 
            B1B = CCOL*COSTH*PNR1C1-CCOLM*PNR1C0 
            BKR = RBSKPR(ICOL,1) 
            BKI = RBSKPR(ICOL,2) 
            HKR=RHSKPR(ICOL,1) 
            HKI=RHSKPR(ICOL,2) 
            HBKMLR=CDMPYR(BSSLSP(IROW+1,ITH,JS),CNEUMN(IROW+1,ITH,JS),        &
           &BSLKPR(ICOL+1,ITH,JS),BSLKPI(ICOL+1,ITH,JS))                      
            HBKMLI=CDMPYI(BSSLSP(IROW+1,ITH,JS),CNEUMN(IROW+1,ITH,JS),        &
           &BSLKPR(ICOL+1,ITH,JS),BSLKPI(ICOL+1,ITH,JS))                      
            BBKMLR = BSSLSP(IROW+1,ITH,JS)*BSLKPR(ICOL+1,ITH,JS) 
            BBKMLI = BSSLSP(IROW+1,ITH,JS)*BSLKPI(ICOL+1,ITH,JS) 
            HHKMLR=CDMPYR(BSSLSP(IROW+1,ITH,JS),CNEUMN(IROW+1,ITH,JS),        &
           &HANKLR(ICOL+1),HANKLI(ICOL+1))                                    
            HHKMLI=CDMPYI(BSSLSP(IROW+1,ITH,JS),CNEUMN(IROW+1,ITH,JS),        &
           &HANKLR(ICOL+1),HANKLI(ICOL+1))                                    
            BHKMLR=BSSLSP(IROW+1,ITH,JS)*HANKLR(ICOL+1) 
            BHKMLI=BSSLSP(IROW+1,ITH,JS)*HANKLI(ICOL+1) 
            HEPSR = CDMPYR(QSREAL,QSIMAG,HBKMLR,HBKMLI) 
            HEPSI = CDMPYI(QSREAL,QSIMAG,HBKMLR,HBKMLI) 
            BEPSR = CDMPYR(QSREAL,QSIMAG,BBKMLR,BBKMLI) 
            BEPSI = CDMPYI(QSREAL,QSIMAG,BBKMLR,BBKMLI) 
            HHEPSR=CDMPYR(QSREAL,QSIMAG,HHKMLR,HHKMLI) 
            HHEPSI=CDMPYI(QSREAL,QSIMAG,HHKMLR,HHKMLI) 
            BBEPSR=CDMPYR(QSREAL,QSIMAG,BHKMLR,BHKMLI) 
            BBEPSI=CDMPYI(QSREAL,QSIMAG,BHKMLR,BHKMLI) 
             FTHETA=THETA*RTD 
      !      WRITE(6,902) JS,ITH,IROW,ICOL,FTHETA,BKPR2(IROW+1,ITH,JS),       
      !    1 BKPI2(IROW+1,ITH,JS),BSSLSP(IROW+1,ITH,JS)                       
            IF(IB.EQ.9) GO TO 380 
      !     IF IB = 8 (MIRROR SYMMETRY BODY), I=L=0 IF IROW AND ICOL ARE BOTH 
      !     ODD OR BOTH EVEN, J=K=0 IF IROW AND ICOL ARE ODD,EVEN OR EVEN,ODD.
            IF((IROW+ICOL).EQ.((IROW+ICOL)/2)*2) GO TO 392 
      !     TEST FOR M=0 (IF M=0 THE I AND L SUBMATRICES ARE ZERO).           
        380 IF(KMV.EQ.0) GO TO 390 
      !     FILL OUT ELEMENTS FOR EQUIVALENT I-SUBMATRIX POSITION.            
            B1 = B1A+B1B 
            HTBKR = CDMPYR(HR,HI,BKR,BKI) 
            HTBKI = CDMPYI(HR,HI,BKR,BKI) 
            BTBKR = BR*BKR 
            BTBKI = BR*BKI 
            TEMPP=(CROW*CROW1*BKR+CCOL*CCOL1*HR-CRSSIJ*(CRIJ+2.0  )/CKR       &
           &)*DCKR*SINTH                                                      
            SUMAR=TEMPP*PNR1C1 
            SUMR=(CKR*(1.0  +HTBKR)-CCOL*HR-CROW*BKR+CRSSIJ/CKR)*B1*CKR+SUMAR 
            SUMAI = PNR1C1*(CROW*CROW1*BKI+CCOL*CCOL1*HI)*DCKR*SINTH 
            SUMI = (CKR*HTBKI-CCOL*HI-CROW*BKI)*B1*CKR+SUMAI 
            T(1,1,1)=B89*CMV*SRMSIN*CDMPYR(SUMR,SUMI,HBKMLR,HBKMLI) 
            T(1,1,2)=B89*CMV*SRMSIN*CDMPYI(SUMR,SUMI,HBKMLR,HBKMLI) 
            SUMBR = PNR1C1*(CROW*CROW1*BKR+CCOL*CCOL1*BR-CRSSIJ*(CRIJ+2.0D0)/C&
           &KR)*DCKR*SINTH                                                    
            SUMR=(CKR*(1.0  +BTBKR)-CCOL*BR-CROW*BKR+CRSSIJ/CKR)*B1*CKR+SUMBR 
            SUMBI = PNR1C1*CROW*CROW1*BKI*DCKR*SINTH 
            SUMI = (CKR*BTBKI-CROW*BKI)*B1*CKR+SUMBI 
            T(1,2,1)=B89*CMV*SRMSIN*CDMPYR(SUMR,SUMI,BBKMLR,BBKMLI) 
            T(1,2,2)=B89*CMV*SRMSIN*CDMPYI(SUMR,SUMI,BBKMLR,BBKMLI) 
            BTHKR=BR*HKR 
            BTHKI=BR*HKI 
      !     FILL OUT ELEMENTS FOR EQUIVALENT L-SUBMATRIX POSITION.            
            SUMR=(CKR*(DCNR+HTBKR)-CCOL*HR-CROW*BKR+CRSSIJ/CKR)*B1*CKR+SUMAR 
            SUMI=(CKR*(DCNI+HTBKI)-CCOL*HI-CROW*BKI)*B1*CKR+SUMAI 
            T(2,1,1)=-B89*CMV*SRMSIN*CDMPYR(SUMR,SUMI,HEPSR,HEPSI) 
            T(2,1,2)=-B89*CMV*SRMSIN*CDMPYI(SUMR,SUMI,HEPSR,HEPSI) 
            SUMR = (CKR*(DCNR+BTBKR)-CCOL*BR-CROW*BKR+CRSSIJ/CKR)*B1*CKR+SUMBR 
            SUMI = (CKR*(DCNI+BTBKI)-CROW*BKI)*B1*CKR+SUMBI 
            T(2,2,1)=-B89*CMV*SRMSIN*CDMPYR(SUMR,SUMI,BEPSR,BEPSI) 
            T(2,2,2)=-B89*CMV*SRMSIN*CDMPYI(SUMR,SUMI,BEPSR,BEPSI) 
        390 IF (IB.EQ.8) GO TO 400 
      !     FILL OUT ELEMENTS FOR EQIIVALENT J-SUBMATRIX POSITION.            
        392 A12=CMCRCO*PNR1C1+QEM*(CROW*CCOLM*COSTH*PNR1C0+CCOL*CROWM*COSTH*PN&
           &R0C1-CROWM*CCOLM*PNR0C0)                                          
            B1A = CCOL*CCOL1*B1A 
            B1B = CROW*CROW1*B1B 
            B1 = (B1A-B1B)*SINTH 
            DD=-QEM*DCKR 
            CR = CDMPYR(DCNR,DCNI,HR,HI) 
            CI = CDMPYI(DCNR,DCNI,HR,HI) 
            SUMR=(CKR*(BKR-CR)+DCNR*CROW-CCOL)*A12*CKR+(B1A-DCNR*B1B)*SINTH*DD 
            SUMI=(CKR*(BKI-CI)+DCNI*CROW)*A12*CKR-(DCNI*B1B)*SINTH*DD 
            T(3,1,1)=B89*SRMSIN*CDMPYR(SUMR,SUMI,HEPSR,HEPSI) 
            T(3,1,2)=B89*SRMSIN*CDMPYI(SUMR,SUMI,HEPSR,HEPSI) 
            CR = BR*DCNR 
            CI = BR*DCNI 
            SUMR=(CKR*(BKR-CR)+DCNR*CROW-CCOL)*A12*CKR+(B1A-DCNR*B1B)*SINTH*DD 
            SUMI=(CKR*(BKI-CI)+DCNI*CROW)*A12*CKR-(DCNI*B1B)*SINTH*DD 
            T(3,2,1)=B89*SRMSIN*CDMPYR(SUMR,SUMI,BEPSR,BEPSI) 
            T(3,2,2)=B89*SRMSIN*CDMPYI(SUMR,SUMI,BEPSR,BEPSI) 
      !     FILL OUT ELEMENTS FOR EQUIVALENT K-SUBMATRIX.                     
            SUMR = (CKR*(BKR-HR)+CROW-CCOL)*A12*CKR+B1*DD 
            SUMI = (CKR*(BKI-HI))*A12*CKR 
            T(4,1,1)=B89*SRMSIN*CDMPYR(SUMR,SUMI,HBKMLR,HBKMLI) 
            T(4,1,2)=B89*SRMSIN*CDMPYI(SUMR,SUMI,HBKMLR,HBKMLI) 
            SUMR = (CKR*(BKR-BR)+CROW-CCOL)*A12*CKR+B1*DD 
            SUMI = BKI*A12*CKR**2 
            T(4,2,1)=B89*SRMSIN*CDMPYR(SUMR,SUMI,BBKMLR,BBKMLI) 
            T(4,2,2)=B89*SRMSIN*CDMPYI(SUMR,SUMI,BBKMLR,BBKMLI) 
        400 RETURN 
            END                                           
            SUBROUTINE QUAD(A,B,K,RESULT,NPTS,ITH,JS) 
      !     NUMERICAL INTEGRATION USING GAUSS-LEGENDRE METHOD.                
               DOUBLE PRECISION P,D1,D2,D3,D4,D5,D6,D7,D8 
            COMMON/GAUSS/KGAUSS 
            DIMENSION P(381),FUNCT(16,127),FZERO(4,2,2),ACUM(16),TEST(16),    &
           &T(4,2,2),RESULT(16,8)                                             
            DIMENSION D1(54),D2(54),D3(54),D4(54),D5(54),D6(54),D7(54),D8(3) 
            EQUIVALENCE (P(1),D1(1)),(P(55),D2(1)),(P(109),D3(1)),(P(163),D4  &
           &(1)),(P(217),D5(1)),(P(271),D6(1)),(P(325),D7(1)),(P(379),D8(1))  
            DATA D1/                                                          &
           & 7.74596669241483D-01, 5.55555555555557D-01, 8.88888888888889D-01,&
           & 2.68488089868333D-01, 9.60491268708019D-01, 1.04656226026467D-01,&
           & 4.34243749346802D-01, 4.01397414775962D-01, 4.50916538658474D-01,&
           & 1.34415255243784D-01, 5.16032829970798D-02, 2.00628529376989D-01,&
           & 9.93831963212756D-01, 1.70017196299402D-02, 8.88459232872258D-01,&
           & 9.29271953151245D-02, 6.21102946737228D-01, 1.71511909136392D-01,&
           & 2.23386686428967D-01, 2.19156858401588D-01, 2.25510499798206D-01,&
           & 6.72077542959908D-02, 2.58075980961766D-02, 1.00314278611795D-01,&
           & 8.43456573932111D-03, 4.64628932617579D-02, 8.57559200499902D-02,&
           & 1.09578421055925D-01, 9.99098124967666D-01, 2.54478079156187D-03,&
           & 9.81531149553739D-01, 1.64460498543878D-02, 9.29654857429739D-01,&
           & 3.59571033071293D-02, 8.36725938168868D-01, 5.69795094941234D-02,&
           & 7.02496206491528D-01, 7.68796204990037D-02, 5.31319743644374D-01,&
           & 9.36271099812647D-02, 3.31135393257977D-01, 1.05669893580235D-01,&
           & 1.12488943133187D-01, 1.11956873020953D-01, 1.12755256720769D-01,&
           & 3.36038771482077D-02, 1.29038001003512D-02, 5.01571393058995D-02,&
           & 4.21763044155885D-03, 2.32314466399103D-02, 4.28779600250078D-02,&
           & 5.47892105279628D-02, 1.26515655623007D-03, 8.22300795723591D-03/
            DATA D2/                                                          &
           & 1.79785515681282D-02, 2.84897547458336D-02, 3.84398102494556D-02,&
           & 4.68135549906281D-02, 5.28349467901166D-02, 5.59784365104763D-02,&
           & 9.99872888120358D-01, 3.63221481845531D-04, 9.97206259372224D-01,&
           & 2.57904979468569D-03, 9.88684757547428D-01, 6.11550682211726D-03,&
           & 9.72182874748583D-01, 1.04982469096213D-02, 9.46342858373402D-01,&
           & 1.54067504665595D-02, 9.10371156957005D-01, 2.05942339159128D-02,&
           & 8.63907938193691D-01, 2.58696793272147D-02, 8.06940531950218D-01,&
           & 3.10735511116880D-02, 7.39756044352696D-01, 3.60644327807826D-02,&
           & 6.62909660024781D-01, 4.07155101169443D-02, 5.77195710052045D-01,&
           & 4.49145316536321D-02, 4.83618026945841D-01, 4.85643304066732D-02,&
           & 3.83359324198731D-01, 5.15832539520484D-02, 2.77749822021825D-01,&
           & 5.39054993352661D-02, 1.68235251552208D-01, 5.54814043565595D-02,&
           & 5.63443130465928D-02, 5.62776998312542D-02, 5.63776283603847D-02,&
           & 1.68019385741038D-02, 6.45190005017574D-03, 2.50785696529497D-02,&
           & 2.10881524572663D-03, 1.16157233199551D-02, 2.14389800125039D-02,&
           & 2.73946052639814D-02, 6.32607319362634D-04, 4.11150397865470D-03,&
           & 8.98927578406411D-03, 1.42448773729168D-02, 1.92199051247278D-02,&
           & 2.34067774953141D-02, 2.64174733950583D-02, 2.79892182552381D-02/
            DATA D3/                                                          &
           & 1.80739564445388D-04, 1.28952408261042D-03, 3.05775341017553D-03,&
           & 5.24912345480885D-03, 7.70337523327974D-03, 1.02971169579564D-02,&
           & 1.29348396636074D-02, 1.55367755558440D-02, 1.80322163903913D-02,&
           & 2.03577550584721D-02, 2.24572658268161D-02, 2.42821652033366D-02,&
           & 2.57916269760242D-02, 2.69527496676331D-02, 2.77407021782797D-02,&
           & 2.81388499156271D-02, 9.99982430354891D-01, 5.05360952078625D-05,&
           & 9.99598799671912D-01, 3.77746646326985D-04, 9.98316635318407D-01,&
           & 9.38369848542380D-04, 9.95724104698407D-01, 1.68114286542147D-03,&
           & 9.91495721178104D-01, 2.56876494379402D-03, 9.85371499598521D-01,&
           & 3.57289278351730D-03, 9.77141514639705D-01, 4.67105037211432D-03,&
           & 9.66637851558417D-01, 5.84344987583563D-03, 9.53730006425761D-01,&
           & 7.07248999543356D-03, 9.38320397779592D-01, 8.34283875396818D-03,&
           & 9.20340025470011D-01, 9.64117772970252D-03, 8.99744899776941D-01,&
           & 1.09557333878379D-02, 8.76513414484705D-01, 1.22758305600827D-02,&
           & 8.50644494768350D-01, 1.35915710097655D-02, 8.22156254364980D-01,&
           & 1.48936416648152D-02, 7.91084933799848D-01, 1.61732187295777D-02,&
           & 7.57483966380512D-01, 1.74219301594641D-02, 7.21423085370098D-01,&
           & 1.86318482561388D-02, 6.82987431091078D-01, 1.97954950480975D-02/
            DATA D4/                                                          &
           & 6.42276642509760D-01, 2.09058514458120D-02, 5.99403930242243D-01,&
           & 2.19563663053178D-02, 5.54495132631931D-01, 2.29409642293877D-02,&
           & 5.07687757533716D-01, 2.38540521060385D-02, 4.59130011989833D-01,&
           & 2.46905247444876D-02, 4.08979821229888D-01, 2.54457699654648D-02,&
           & 3.57403837831532D-01, 2.61156733767061D-02, 3.04576441556714D-01,&
           & 2.66966229274503D-02, 2.50678730303482D-01, 2.71855132296248D-02,&
           & 1.95897502711100D-01, 2.75797495664819D-02, 1.40424233152560D-01,&
           & 2.78772514766137D-02, 8.44540400837110D-02, 2.80764557938172D-02,&
           & 2.81846489497457D-02, 2.81763190330167D-02, 2.81888141801924D-02,&
           & 8.40096928705192D-03, 3.22595002508787D-03, 1.25392848264749D-02,&
           & 1.05440762286332D-03, 5.80786165997757D-03, 1.07194900062519D-02,&
           & 1.36973026319907D-02, 3.16303660822264D-04, 2.05575198932735D-03,&
           & 4.49463789203206D-03, 7.12243868645840D-03, 9.60995256236391D-03,&
           & 1.17033887476570D-02, 1.32087366975291D-02, 1.39946091276191D-02,&
           & 9.03727346587510D-05, 6.44762041305726D-04, 1.52887670508776D-03,&
           & 2.62456172740443D-03, 3.85168761663987D-03, 5.14855847897819D-03,&
           & 6.46741983180368D-03, 7.76838777792199D-03, 9.01610819519566D-03,&
           & 1.01788775292361D-02, 1.12286329134080D-02, 1.21410826016683D-02/
            DATA D5/                                                          &
           & 1.28958134880121D-02, 1.34763748338165D-02, 1.38703510891399D-02,&
           & 1.40694249578135D-02, 2.51578703842806D-05, 1.88873264506505D-04,&
           & 4.69184924247851D-04, 8.40571432710723D-04, 1.28438247189701D-03,&
           & 1.78644639175865D-03, 2.33552518605716D-03, 2.92172493791781D-03,&
           & 3.53624499771678D-03, 4.17141937698409D-03, 4.82058886485126D-03,&
           & 5.47786669391895D-03, 6.13791528004137D-03, 6.79578550488277D-03,&
           & 7.44682083240758D-03, 8.08660936478883D-03, 8.71096507973207D-03,&
           & 9.31592412806942D-03, 9.89774752404876D-03, 1.04529257229060D-02,&
           & 1.09781831526589D-02, 1.14704821146939D-02, 1.19270260530193D-02,&
           & 1.23452623722438D-02, 1.27228849827324D-02, 1.30578366883530D-02,&
           & 1.33483114637252D-02, 1.35927566148124D-02, 1.37898747832410D-02,&
           & 1.39386257383068D-02, 1.40382278969086D-02, 1.40881595165083D-02,&
           & 9.99997596379750D-01, 6.93793643241083D-06, 9.99943996207055D-01,&
           & 5.32752936697805D-05, 9.99760490924434D-01, 1.35754910949228D-04,&
           & 9.99380338025023D-01, 2.49212400482998D-04, 9.98745614468096D-01,&
           & 3.89745284473282D-04, 9.97805354495956D-01, 5.54295314930373D-04,&
           & 9.96514145914890D-01, 7.40282804244503D-04, 9.94831502800622D-01,&
           & 9.45361516858527D-04, 9.92721344282788D-01, 1.16748411742996D-03/
            DATA D6/                                                          &
           & 9.90151370400771D-01, 1.40490799565515D-03, 9.87092527954033D-01,&
           & 1.65611272815445D-03, 9.83518657578632D-01, 1.91971297101387D-03,&
           & 9.79406281670862D-01, 2.19440692536384D-03, 9.74734459752401D-01,&
           & 2.47895822665757D-03, 9.69484659502459D-01, 2.77219576459345D-03,&
           & 9.63640621569812D-01, 3.07301843470258D-03, 9.57188216109859D-01,&
           & 3.38039799108691D-03, 9.50115297521293D-01, 3.69337791702565D-03,&
           & 9.42411565191083D-01, 4.01106872407503D-03, 9.34068436157727D-01,&
           & 4.33264096809299D-03, 9.25078932907077D-01, 4.65731729975685D-03,&
           & 9.15437587155765D-01, 4.98436456476553D-03, 9.05140358813263D-01,&
           & 5.31308660518706D-03, 8.94184568335557D-01, 5.64281810138445D-03,&
           & 8.82568840247341D-01, 5.97291956550816D-03, 8.70293055548114D-01,&
           & 6.30277344908575D-03, 8.57358310886234D-01, 6.63178124290190D-03,&
           & 8.43766882672707D-01, 6.95936140939044D-03, 8.29522194637402D-01,&
           & 7.28494798055382D-03, 8.14628787655138D-01, 7.60798966571904D-03,&
           & 7.99092290960843D-01, 7.92794933429486D-03, 7.82919394118284D-01,&
           & 8.24430376303287D-03, 7.66117819303759D-01, 8.55654356130769D-03,&
           & 7.48696293616938D-01, 8.86417320948252D-03, 7.30664521242183D-01,&
           & 9.16671116356077D-03, 7.12033155362253D-01, 9.46368999383007D-03/
            DATA D7/                                                          &
           & 6.92813769779114D-01, 9.75465653631741D-03, 6.73018830230419D-01,&
           & 1.00391720440569D-02, 6.52661665410019D-01, 1.03168123309476D-02,&
           & 6.31756437711193D-01, 1.05871679048852D-02, 6.10318113715188D-01,&
           & 1.08498440893373D-02, 5.88362434447664D-01, 1.11044611340069D-02,&
           & 5.65905885423653D-01, 1.13506543159806D-02, 5.42965666498311D-01,&
           & 1.15880740330440D-02, 5.19559661537457D-01, 1.18163858908302D-02,&
           & 4.95706407918762D-01, 1.20352707852796D-02, 4.71425065871658D-01,&
           & 1.22444249816120D-02, 4.46735387662029D-01, 1.24435601907140D-02,&
           & 4.21657686626164D-01, 1.26324036435421D-02, 3.96212806057616D-01,&
           & 1.28106981638774D-02, 3.70422087950079D-01, 1.29782022395374D-02,&
           & 3.44307341599437D-01, 1.31346900919602D-02, 3.17890812068477D-01,&
           & 1.32799517439305D-02, 2.91195148518247D-01, 1.34137930851101D-02,&
           & 2.64243372410927D-01, 1.35360359349562D-02, 2.37058845589829D-01,&
           & 1.36465181025713D-02, 2.09665238243181D-01, 1.37450934430019D-02,&
           & 1.82086496759252D-01, 1.38316319095064D-02, 1.54346811481378D-01,&
           & 1.39060196013255D-02, 1.26470584372302D-01, 1.39681588065169D-02,&
           & 9.84823965981194D-02, 1.40179680394566D-02, 7.04069760428552D-02,&
           & 1.40553820726499D-02, 4.22691647653637D-02, 1.40803519625536D-02/
            DATA D8/                                                          &
           & 1.40938864107825D-02, 1.40928450691604D-02, 1.40944070900962D-02/
            IF(A.EQ.B)GO TO 107 
            SUM=(B+A)/2.0 
            DIFF=(B-A)/2.0 
      !     ONE POINT FORMULA.                                                
      !     SET UP VARIABLE COMBINATIONS FOR USE IN EVALUATION OF INTEGRANDS  
             ITH=ITH+1 
            CALL GENER(SUM,T,ITH,JS) 
            DO 1010 II=1,4 
            DO 1020 JJ=1,2 
            DO 1030 KK=1,2 
      !     JJ=1,2 CORRESPONDS TO MATRICES A AND B . II=1,4                   
      !     CORRESPONDS TO FILLING OUT EQ. I,L,J,K POSITOONS OR SUB-MATRICES  
      !     KK=1,2 CORRESPONDS TO REAL AND IMAGINARY PARTS.                   
            FZERO(II,JJ,KK)=T(II,JJ,KK) 
            L=(II-1)*4+(JJ-1)*2+KK 
            RESULT(L,1)=2.0  *T(II,JJ,KK)*DIFF 
       1030 END DO 
       1020 END DO 
       1010 END DO 
            I=0 
            IOLD=0 
            INEW=1 
            K=2 
            DO 1040 N=1,16 
            ACUM(N)=0.0 
       1040  CONTINUE 
            GO TO 103 
        101 CONTINUE 
            IF(K.EQ.KGAUSS) GO TO 105 
            K=K+1 
            DO 1050 N=1,16 
            ACUM(N)=0.0 
       1050 END DO 
      !     CONTRIBUTION FROM FUNCTION VALUES ALREADY COMPUTED.               
            DO 102 J=1,IOLD 
            I=I+1 
            DO 1060 N=1,16 
            ACUM(N)=ACUM(N)+P(I)*FUNCT(N,J) 
       1060 END DO 
        102 END DO 
      !     CONTRIBUTION FROM NEW VALUES.                                     
        103 CONTINUE 
            IOLD=IOLD+INEW 
            DO 104 J=INEW,IOLD 
            I=I+1 
            X=P(I)*DIFF 
            TEMP1=SUM+X 
              ITH=ITH+1 
            CALL GENER(TEMP1,T,ITH,JS) 
            DO 1070 II=1,4 
            DO 1080 JJ=1,2 
            DO 1090 KK=1,2 
            L=(II-1)*4+(JJ-1)*2+KK 
      !     L GOES FROM 1 TO 16.                                              
            FUNCT(L,J)=T(II,JJ,KK) 
       1090 END DO 
       1080 END DO 
       1070 END DO 
            TEMP2=SUM-X 
              ITH=ITH+1 
            CALL GENER(TEMP2,T,ITH,JS) 
            DO 1071 II=1,4 
            DO 1081 JJ=1,2 
            DO 1091 KK=1,2 
            L=(II-1)*4+(JJ-1)*2+KK 
            FUNCT(L,J)=FUNCT(L,J)+T(II,JJ,KK) 
       1091 END DO 
       1081 END DO 
       1071 END DO 
            I=I+1 
            DO 1100 N=1,16 
            ACUM(N)=ACUM(N)+P(I)*FUNCT(N,J) 
       1100 END DO 
        104 END DO 
            INEW=IOLD+1 
            I=I+1 
            DO 1200 II=1,4 
            DO 1300 JJ=1,2 
            DO 1400 KK=1,2 
            N=(II-1)*4+(JJ-1)*2+KK 
            RESULT(N,K)=(ACUM(N)+P(I)*FZERO(II,JJ,KK))*DIFF 
       1400 END DO 
       1300 END DO 
       1200 END DO 
            GO TO 101 
        105 CONTINUE 
      !     NORMAL TERMINATION.                                               
        106 CONTINUE 
            NPTS=INEW+IOLD 
            RETURN 
      !     TRIVIAL CASE                                                      
        107 CONTINUE 
            K=2 
            DO 1600 M=1,16 
            DO 1700 N=1,2 
            RESULT(M,N)=0.0 
       1700 END DO 
       1600 END DO 
            NPTS=0 
            RETURN 
            END                                           
            SUBROUTINE GENLGP 
      !     A ROUTINE TO GENERATE LEGENDRE POLYNOMIALS.                       
      !     THE INDEX ON THE FUNCTION IS INCREMENTED BY ONE.                  
            COMMON DTR,RTD,CPI 
            COMMON/MTXCOM/NRANK,NRANKI,A(80,80,2),B(80,80,2),CMXNRM(80) 
            COMMON/FNCCOM/PNMLLG(81),BSSLSP(81,31,3),CNEUMN(81,31,3),         &
           &BSLKPR(81,31,3),BSLKPI(81,31,3),CNEUMR(81,31,3),CNEUMI(81,31,3)   
            COMMON /CMVCOM/ NM,KMV,CMI(20),CMV,CM2,TWM,PRODM 
            COMMON /THTCOM/ THETA,SINTH,COSTH 
            DTWM=TWM+1.0 
      !     THIS IS SPECIAL CASE WHEN THETA EQUALS CPIAND M=0.                
      !     WHEN THETA EQUALS CPI ALL INTEGRANDS ARE 0 AND ANY VALUES CAN BE  
      !     PUT IN PNMLLG(41).HERE WE HAVE PUT THEM EQUAL TO 0.               
            IF((SINTH.EQ.0.0  ).AND.(KMV.EQ.0)) GO TO 6 
      !     AT THIS POINT THETA LIES STRICTLY BETWEEN 0 AND CPI.              
            IF(THETA)16,4,16 
          4 IF(KMV-1)6,12,6 
          6 DO 8 ILG = 1,NRANKI 
            PNMLLG(ILG)=0.0 
          8 END DO 
            GO TO 88 
         12 PNMLLG(1)=0.0 
            PNMLLG(2)=1.0 
            PLA=1.0 
            GO TO 48 
         16 IF(KMV)20,20,40 
      !     THE SPECIAL CASE WHEN M = 0.                                      
         20 PLA=1.0/SINTH 
            PLB = COSTH*PLA 
            PNMLLG(1) = PLA 
            PNMLLG(2) = PLB 
            IBEG = 3 
            GO TO 60 
      !     GENERAL CASE FOR M NOT EQUAL TO 0.                                
         40 DO 44 ILG = 1,KMV 
            PNMLLG(ILG)=0.0 
         44 END DO 
            IF((SINTH.EQ.0.0  ).AND.(KMV.EQ.1)) GO TO 1001 
            PLA = PRODM*SINTH**(KMV-1) 
            GO TO 1002 
       1001 PLA=0.0 
       1002 CONTINUE 
            PNMLLG(KMV+1) = PLA 
         48 PLB = DTWM*COSTH*PLA 
            PNMLLG(KMV+2) = PLB 
            IBEG = KMV+3 
      !     DO RECURSION FORMULA FOR ALL REMAINING LEGENDRE POLYNOMIALS.      
         60 CNMUL = IBEG+IBEG-3 
            CNM=2.0 
            CNMM = DTWM 
            DO 80 ILGR = IBEG,NRANKI 
            PLC = (CNMUL*COSTH*PLB-CNMM*PLA)/CNM 
            PNMLLG(ILGR) = PLC 
            PLA = PLB 
            PLB = PLC 
            CNMUL=CNMUL+2.0 
            CNM=CNM+1.0 
            CNMM=CNMM+1.0 
         80 END DO 
         88 RETURN 
            END                                           
            FUNCTION CDABX(A,B) 
            IF(A) 4,22,4 
          4 IF(B) 8,30,8 
          8 E = AMAX1(A,B) 
            F = AMIN1(A,B) 
            G = F/E 
            CDABX= ABS(E)*SQRT(1.0  +G*G) 
            RETURN 
         22 IF(B) 28,26,28 
         26 CDABX=0.0 
            RETURN 
         28 CDABX=ABS(B) 
            RETURN 
         30 CDABX=ABS(A) 
            RETURN 
            END                                           
            FUNCTION CDMPYR(A,B,C,D) 
            CDMPYR = A*C-B*D 
            RETURN 
            END                                           
            FUNCTION CDMPYI(A,B,C,D) 
            CDMPYI = B*C+A*D 
            RETURN 
            END                                           
            FUNCTION CDDVDR(A,B,C,D) 
            E = C*C+D*D 
            F=A*C+B*D 
            CDDVDR=F/E 
            RETURN 
            END                                           
            FUNCTION CDDVDI(A,B,C,D) 
            E = C*C+D*D 
            F=B*C-A*D 
            CDDVDI=F/E 
            RETURN 
            END                                           
            SUBROUTINE GENBSL(ITH,JS) 
      !     GENERATE BESSEL AND NEUMANN FUNCTIONS FOR REAL ARGUMENTS.         
      !     THE INDEX ON THE FUNCTION IS INCREMENTED BY ONE.                  
            COMMON DTR,RTD,CPI 
            COMMON/MTXCOM/NRANK,NRANKI,A(80,80,2),B(80,80,2),CMXNRM(80) 
            COMMON /THTCOM/ THETA,SINTH,COSTH 
            COMMON /BDYCOM/ DCNR,DCNI,CKPRR,CKPRI,CKR,DCKR,CONK,AOVRB,SIGMA,IB 
            COMMON/FNCCOM/PNMLLG(81),BSSLSP(81,31,3),CNEUMN(81,31,3),         &
           &BSLKPR(81,31,3),BSLKPI(81,31,3),CNEUMR(81,31,3),CNEUMI(81,31,3)   
      !     SET UP A LOOP TO GET 2 SUCCESSIVE BESSEL FUNCTIONS                
            NVAL=NRANK-1 
            PCKR=CKR 
            DO 40 I=1,4 
            CALL BESSEL(NVAL,PCKR,ANSWR,IERROR) 
            IF(IERROR) 20,20,32 
         20 ANSA=ANSWR 
            NVAL=NVAL+1 
            CALL BESSEL(NVAL,PCKR,ANSWR,IERROR) 
            IF(IERROR) 24,24,28 
         24 ANSB=ANSWR 
            GO TO 60 
         28 NVAL=NVAL-1 
         32 NVAL=NVAL+NRANK 
         40 END DO 
      !     PROGRAM UNABLE TO GENERATE BESSEL FUNCTIONS                       
            WRITE(6,1001) 
       1001 FORMAT(///,5X,'UNABLE TO GENERATE BESSEL FUNCTIONS',//) 
      !     SET UP FOR PROPER RECURSION OF THE BESSEL FUNCTONS                
         60 IF(NVAL-NRANK)100,100,64 
         64 IEND=NVAL-NRANK 
            CONN=2*(NVAL-1)+1.0 
            DO 72 IP=1,IEND 
            ANSC=CONN*ANSA/PCKR-ANSB 
            CONN=CONN-2.0 
            ANSB=ANSA 
            ANSA=ANSC 
         72  CONTINUE 
      !     PROGRAM IS READY TO RECURSE DOWNWARD INTO BESSEL FUNCTION         
        100 BSSLSP(NRANKI,ITH,JS)=ANSB 
            BSSLSP(NRANKI-1,ITH,JS)=ANSA 
            CONN=    (FLOAT(NRANK+NRANK-1)) 
            IE=NRANKI-2 
            JE=IE 
            DO 120 JB=1,JE 
            ANSC=CONN*ANSA/PCKR-ANSB 
            BSSLSP(IE,ITH,JS)=ANSC 
            ANSB=ANSA 
            ANSA=ANSC 
            IE=IE-1 
            CONN=CONN-2.0 
        120 END DO 
      !     GENERATE NEUMANN FUNCTIONS                                        
            CSKRX= COS(PCKR)/PCKR 
            SNKRX= SIN(PCKR)/PCKR 
            CKR2=PCKR**2 
            CMULN=3.0 
            SNSA=-CSKRX 
            SNSB=-CSKRX/PCKR-SNKRX 
            CNEUMN(1,ITH,JS)=SNSA 
            CNEUMN(2,ITH,JS)=SNSB 
            DO 280 I=3,NRANKI 
            SNSC=CMULN*SNSB/PCKR-SNSA 
            CNEUMN(I,ITH,JS)=SNSC 
            SNSA=SNSB 
            SNSB=SNSC 
            CMULN=CMULN+2.0 
        280 END DO 
      !     PERFORM WRONSKIAN TEST ON ORDERS 0 AND 1 AND ORDERS NRANK-1 AND NR
            QUANBT=ABS(CKR2*(BSSLSP(2,ITH,JS)*CNEUMN(1,ITH,JS)-               &
           &BSSLSP(1,ITH,JS)*CNEUMN(2,ITH,JS))-1.0)                           
            QUANNT=ABS(CKR2*(BSSLSP(NRANKI,ITH,JS)*CNEUMN(NRANK,ITH,JS)-      &
           &BSSLSP(NRANK,ITH,JS)*CNEUMN(NRANKI,ITH,JS))-1.0  )                
            IF(QUANBT-1.0E-10)360,352,352 
        352 THTPRT = RTD*THETA 
            WRITE(6,356) THTPRT,PCKRR,QUANBT,QUANNT 
        356 FORMAT(/,10X,'THETA=',F9.4,'KR=',F10.4,'BESSEL TEST=',E12.5,      &
           &'NEUMANN TEST=',E12.5)                                            
            GO TO 362 
        360 IF(QUANNT-1.0E-10)362,352,352 
        362 RETURN 
            END                                           
            SUBROUTINE BESSEL(NORDER,ARGMNT,ANSWR,IERROR) 
            IERROR=0 
            N=NORDER 
            X=ARGMNT 
            CN=N 
            SUM=1.0 
            APR=1.0 
            TOPR=-0.5*X*X 
            CI=1.0 
            CNI=FLOAT(2*N) +3.0 
            CNI=FLOAT(2*N) +3.0 
            DO 60 I=1,100 
            ACR=TOPR*APR/(CI*CNI) 
            SUM=SUM+ACR 
            IF( ABS(ACR/SUM)-1.0E-20) 100,100,40 
         40 APR=ACR 
            CI=CI+1.0 
            CNI=CNI+2.0 
         60 END DO 
            IERROR=1 
            GO TO 200 
      !     THE SERIES HAS CONVERGED                                          
        100 PROD=    (FLOAT(2*N))+1.0 
            FACT=1.0D0 
            IF(N) 160,160,120 
        120  DO 140 IFCT=1,N 
            FACT=FACT*X/PROD 
            PROD=PROD-2.0 
        140 END DO 
        160 ANSWR=FACT*SUM 
        200 RETURN 
            END                                           
            FUNCTION CROOTR(A,B) 
            DMAG=(A*A+B*B)**0.25 
            ANGLE=0.5*ATAN2(B,A) 
            CROOTR=DMAG*COS(ANGLE) 
            RETURN 
            END                                           
            FUNCTION CROOTI(A,B) 
            DMAG=(A*A+B*B)**0.25 
            ANGLE=0.5*ATAN2(B,A) 
            CROOTI=DMAG*SIN(ANGLE) 
            RETURN 
            END                                           
            SUBROUTINE PRCSSM(A,B,NR,NRI) 
      !     A ROUTINE TO SOLVE THE EQUATION T = (A-INVERSE)*B  ( ALL MATRICES 
      !     ARE TRANSPOSED) USING GAUSS-JORDAN ELIMINATION.                   
            DIMENSION A(80,80,2),B(80,80,2) 
            DIMENSION AIJMAX(2),ARAT(2) 
            EQUIVALENCE (L,FL),(K,FK) 
            N = 2*NR 
      !     START REDUCTION OF THE A MATRIX.                                  
            DO 80 I = 1,N 
      !     SEARCH FOR THE MAXIMUM ELEMENT IN THE ITH ROW OF THE A-MATRIX.    
            AIJMAX(1) = A(I,1,1) 
            AIJMAX(2) = A(I,1,2) 
            JMAX = 1 
            DO 10 J = 2,N 
            IF(CDABX(A(I,J,1),A(I,J,2)).LE.CDABX(AIJMAX(1),AIJMAX(2))) GOTO 10 
            AIJMAX(1) = A(I,J,1) 
            AIJMAX(2) = A(I,J,2) 
            JMAX = J 
         10 END DO 
      !     IF AIJMAX IS ZERO ( AS IT WILL BE FOR ANY ROW (OR COLUMN) WHERE TH
      !     INDEX M IS .GT. THE INDEX N, I.E., THE LEGENDRE FUNCTIONS FORCE TH
      !     MATRIX ELEMENTS TO ZERO),THEN THE MATRIX IS SINGULAR SO SOLVE THE 
      !     REDUCED MATRIX (ORDER = 2*(NRANK-M)).                             
            IF(CDABX(AIJMAX(1),AIJMAX(2)).GT.0.0  ) GO TO 20 
            JMAX = I 
            GO TO 75 
      !     NORMALIZE THE ITH ROW BY AIJMAX (JMAX ELEMENT OF THE ITH ROW).    
         20 DO 30 J = 1,N 
            T1 = A(I,J,1) 
            T2 = A(I,J,2) 
            A(I,J,1) = CDDVDR(T1,T2,AIJMAX(1),AIJMAX(2)) 
            A(I,J,2) = CDDVDI(T1,T2,AIJMAX(1),AIJMAX(2)) 
      !     NORMALIZE THE ITH ROW OF B.                                       
            T1 = B(I,J,1) 
            T2 = B(I,J,2) 
            B(I,J,1) = CDDVDR(T1,T2,AIJMAX(1),AIJMAX(2)) 
            B(I,J,2) = CDDVDI(T1,T2,AIJMAX(1),AIJMAX(2)) 
         30 END DO 
      !     USE ROW TRANSFORMATIONS TO GET ZEROS ABOVE AND BELOW THE JMAX     
      !     ELEMENT OF THE ITH ROW OF A.  APPLY SAME ROW TRANSFORMATIONS      
      !     TO THE B MATRIX.                                                  
            DO 70 K = 1,N 
            IF(K.EQ.I) GO TO 70 
            ARAT(1) = -A(K,JMAX,1) 
            ARAT(2) = -A(K,JMAX,2) 
            DO 50 J = 1,N 
            IF(CDABX(A(I,J,1),A(I,J,2)).LE.0.0  ) GO TO 50 
            A(K,J,1) = CDMPYR(ARAT(1),ARAT(2),A(I,J,1),A(I,J,2))+A(K,J,1) 
            A(K,J,2) = CDMPYI(ARAT(1),ARAT(2),A(I,J,1),A(I,J,2))+A(K,J,2) 
         50 END DO 
            A(K,JMAX,1)=0.0 
            A(K,JMAX,2)=0.0 
            DO 60 J=1,N 
            IF(CDABX(B(I,J,1),B(I,J,2)).LE.0.0  ) GO TO 60 
            B(K,J,1) = CDMPYR(ARAT(1),ARAT(2),B(I,J,1),B(I,J,2))+B(K,J,1) 
            B(K,J,2) = CDMPYI(ARAT(1),ARAT(2),B(I,J,1),B(I,J,2))+B(K,J,2) 
         60 END DO 
         70 END DO 
      !     STORE ROW COUNTER (I) IN TOP ELEMENT OF JMAX COLUMN.  THUS,       
      !     THE TOP ROW OF A WILL CONTAIN THE LOCATION OF THE PIVOT           
      !     (UNITY) ELEMENT OF EACH COLUMN (AFTER REDUCTION).                 
         75 L = I 
      !     STORE THE INTEGER I IN THE TOP ROW OF A.                          
            A(1,JMAX,1) = FL 
         80 END DO 
      !     THE REDUCTION OF A IS COMPLETE.  PERFORM ROW INTERCHANGES         
      !     AS INDICATED IN THE FIRST ROW OF A.                               
            DO 120 I = 1,N 
            K=I 
      !     PUT THE INTEGER VALUE IN A INTO K.                                
         90 FK = A(1,K,1) 
            IF(K-I) 90,120,100 
      !     IF K(1,I) IS LESS THAN I, THEN THAT ROW HAS ALREADY BEEN          
      !     INVOLVED IN AN INTERCHANGE, AND WE USE K(1,K) UNTIL WE GET        
      !     A VALUE OF K GREATER THAN I (CORRESPONDING TO A ROW STORED        
      !     BELOW THE ITH ROW).                                               
        100 DO 110 J=1,N 
            ARAT(1) = B(I,J,1) 
            ARAT(2) = B(I,J,2) 
            B(I,J,1) = B(K,J,1) 
            B(I,J,2) = B(K,J,2) 
            B(K,J,1) = ARAT(1) 
            B(K,J,2) = ARAT(2) 
        110 END DO 
        120 END DO 
            RETURN 
            END                                           
            SUBROUTINE ADDPRC(IM)
      !     A ROUTINE TO OBTAIN THE SCATTERED FIELD COEFFICIENTS AND CALCULATE
      !     THE DIFFERENTIAL SCATTERING CROSS SECTION IN THE AZIMUTHAL PLANE. 
            INTEGER IM
            COMPLEX CI,CIM 
            COMMON DTR,RTD,CPI 
            COMMON/MTXCOM/NRANK,NRANKI,A(80,80,2),TMAT(80,80,2),CMXNRM(80) 
            COMMON/FNCCOM/PNMLLG(81),BSSLSP(81,31,3),CNEUMN(81,31,3),         &
           &BSLKPR(81,31,3),BSLKPI(81,31,3),CNEUMR(81,31,3),CNEUMI(81,31,3)   
            COMMON /BDYCOM/ DCNR,DCNI,CKPRR,CKPRI,CKR,DCKR,CONK,AOVRB,SIGMA,IB 
            COMMON /CMVCOM/ NM,KMV,CMI(20),CMV,CM2,TWM,PRODM 
            COMMON /THTCOM/ THETA,SINTH,COSTH 
            COMMON /UVCCOM/ANGINC,ACANS(181,2,2),UANG(181),RTSFCT,DLTANG,NUANG 
            COMMON/SCATT/SUMM1,SUMM2 
             COMMON/VIVEK2/TWOA 
            DIMENSION ZXOLD(181),ZYOLD(181),AB1(80,2),AB2(80,2),FG1(80,2),FG2(&
           &80,2),FGANS(181,2,2)                                              
            LOGICAL TEST 
            DATA TEST/.TRUE./ 
            CI = (0.0,1.0) 
      !     GENERATE THE LEGENDRE FUNCTIONS FOR THE INCIDENT ANGLE.           
            IF(ANGINC) 15,5,15 
          5 COSTH=1.0 
         10  SINTH=0.0 
            THETA=0.0 
            GO TO 30 
         15 IF(ANGINC-180.0  ) 25,20,25 
         20 COSTH=-1.0 
            GO TO 10 
         25 THETA = DTR*ANGINC 
            SINTH=SIN(THETA) 
            COSTH=COS(THETA) 
         30 CALL GENLGP 
      !                                                                       
      !     GENERATE THE INCIDENT FIELD COEFFICIENTS -- AB1 = THETA POLARIZATI
      !     AND AB2 = PHI POLARIZATION.                                       
      !                                                                       
            CN=0.0 
            DO 35 N=1,NRANK 
            NP = N+NRANK 
            CN=CN+1.0 
            N1 = N+1 
            CI1R=REAL(CI**N) 
            CI1I=AIMAG(CI**N) 
            CI2R=REAL(CI**N1) 
            CI2I=AIMAG(CI**N1) 
            P1 = CN*COSTH*PNMLLG(N1)-(CN+CMV)*PNMLLG(N) 
            P2 = CMV*PNMLLG(N1) 
            AB1(N,1) = -CI1R*P2 
            AB1(N,2) = -CI1I*P2 
            AB1(NP,1) = CI2R*P1 
            AB1(NP,2) = CI2I*P1 
            AB2(N,1) = CI1R*P1 
            AB2(N,2) = CI1I*P1 
            AB2(NP,1) = -CI2R*P2 
            AB2(NP,2) = -CI2I*P2 
         35 END DO 
      !                                                                       
      !     THE SCATTERED FIELD COEFFICIENTS = THE ROW VECTOR OF INCIDENT FIEL
      !     COEFFICIENTS TIMES THE T-TRANSPOSED MATRIX.                       
      !                                                                       
            NR2 = 2*NRANK 
            DO 45 J = 1,NR2 
            S1R=0.0 
            S1I=0.0 
            S2R=0.0 
            S2I=0.0 
            DO 40 I = 1,NR2 
            S1R = S1R+CDMPYR(AB1(I,1),AB1(I,2),TMAT(I,J,1),TMAT(I,J,2)) 
            S1I = S1I+CDMPYI(AB1(I,1),AB1(I,2),TMAT(I,J,1),TMAT(I,J,2)) 
            S2R = S2R+CDMPYR(AB2(I,1),AB2(I,2),TMAT(I,J,1),TMAT(I,J,2)) 
            S2I = S2I+CDMPYI(AB2(I,1),AB2(I,2),TMAT(I,J,1),TMAT(I,J,2)) 
         40 END DO 
            FG1(J,1) = S1R 
            FG1(J,2) = S1I 
            FG2(J,1) = S2R 
            FG2(J,2) = S2I 
         45 END DO 
      !     CALCULATE SCATTERING COSSECTIONS NORMALIZED FOR PARALLEL AND      
      !     PERPN. POL.                                                       
            SUM1=0.0 
            SUM2=0.0 
            DO 1001 I=1,NRANK 
            II=I+NRANK 
            TEMP1=FG1(I,1)**2+FG1(I,2)**2+FG1(II,1)**2+FG1(II,2)**2 
            TEMP2=FG2(I,1)**2+FG2(I,2)**2+FG2(II,1)**2+FG2(II,2)**2 
            SUM1=SUM1+TEMP1/CMXNRM(I) 
            SUM2=SUM2+TEMP2/CMXNRM(I) 
       1001 END DO 
      !     NORMALIZE SCATTERING CROSSECTIONS.                                
            SUM1=(RTSFCT*2.0/CONK)*SUM1 
            SUM2=(RTSFCT*2.0/CONK)*SUM2 
      !         NORMALIZE W.R.T. EQ. SPHERICAL DIA.                           
                CNORM=AOVRB**(-2./3.) 
                SUM1=SUM1*CNORM 
                SUM2=SUM2*CNORM 
      !     ACCUMULATE RESULTS FOR EACH M VALUE                               
            SUMM1=SUM1+SUMM1 
            SUMM2=SUM2+SUMM2 
      !                                                                       
      !     EVALUATE THE SCATTERED FIELD AT EACH SCATTERING ANGLE.            
      !                                                                       
            DO 170 IU = 1,NUANG 
      !     GENERATE THE LEGENDRE MULTIPLIERS.                                
            IF(UANG(IU)) 95,85,95 
         85 COSTH=1.0 
         90 SINTH=0.0 
            THETA=0.0 
            GO TO 110 
         95 IF(UANG(IU)-180.0  ) 105,100,105 
        100 COSTH=-1.0 
            GO TO 90 
        105 THETA = DTR*UANG(IU) 
            SINTH=SIN(THETA) 
            COSTH=COS(THETA) 
        110 CALL GENLGP 
            FGANS(IU,1,1)=0.0 
            FGANS(IU,1,2)=0.0 
            FGANS(IU,2,1)=0.0 
            FGANS(IU,2,2)=0.0 
            CN=0.0 
            DO 160 N = 1,NRANK 
            NP = N+NRANK 
            N1 = N+1 
            CN=CN+1.0 
            P1 = CN*COSTH*PNMLLG(N1)-(CN+CMV)*PNMLLG(N) 
            P2 = CMV*PNMLLG(N1) 
            CIM = (-CI)**N1 
            CIR=REAL(CIM) 
            CII=AIMAG(CIM) 
            F1R = FG1(N,1)*P2 
            F1I = FG1(N,2)*P2 
            G1R = -FG1(NP,2)*P1 
            G1I = FG1(NP,1)*P1 
            FGANS(IU,1,1) = FGANS(IU,1,1)+CDMPYR(CIR,CII,F1R+G1R,F1I+G1I)/CMXN&
           &RM(N)                                                             
            FGANS(IU,1,2) = FGANS(IU,1,2)+CDMPYI(CIR,CII,F1R+G1R,F1I+G1I)/CMXN&
           &RM(N)                                                             
            F2R = FG2(N,1)*P1 
            F2I = FG2(N,2)*P1 
            G2R = -FG2(NP,2)*P2 
            G2I = FG2(NP,1)*P2 
            FGANS(IU,2,1) = FGANS(IU,2,1)-CDMPYR(CIR,CII,F2R+G2R,F2I+G2I)/CMXN&
           &RM(N)                                                             
            FGANS(IU,2,2) = FGANS(IU,2,2)-CDMPYI(CIR,CII,F2R+G2R,F2I+G2I)/CMXN&
           &RM(N)                                                             
        160 END DO 
      !                                                                       
      !     THE NORMALIZED DIFF.SCAT.CROSS SECT. IS GIVEN BY ((8/KA)*FGANS)**2
      !     SCALE FGANS TO CALCULATE DIFF. SCAT. CROSS SECT. (RTSFCT = 8/KA)  
      !                                                                       
            FGANS(IU,1,1) = RTSFCT*FGANS(IU,1,1) 
            FGANS(IU,1,2) = RTSFCT*FGANS(IU,1,2) 
            FGANS(IU,2,1) = RTSFCT*FGANS(IU,2,1) 
            FGANS(IU,2,2) = RTSFCT*FGANS(IU,2,2) 
        170 END DO 
      !     ACCUMULATE THE RESULTS FOR EACH M VALUE.                          
            DO 172 IUP = 1,NUANG 
            ACANS(IUP,1,1) = ACANS(IUP,1,1)+FGANS(IUP,1,1) 
            ACANS(IUP,1,2) = ACANS(IUP,1,2)+FGANS(IUP,1,2) 
            ACANS(IUP,2,1) = ACANS(IUP,2,1)+FGANS(IUP,2,1) 
            ACANS(IUP,2,2) = ACANS(IUP,2,2)+FGANS(IUP,2,2) 
             WRITE(110,202) UANG(IUP),ACANS(IUP,1,1),ACANS(IUP,1,2) 
        202  FORMAT(2X,'PARRALEL',F7.3,2X,2E15.7) 
             WRITE(110,203) UANG(IUP),ACANS(IUP,2,1),ACANS(IUP,2,2) 
        203  FORMAT(2X,'PERPENDICULAR',F7.3,2X,2E15.7) 
        172 END DO 
      !     CALCULATE THE EXTINCTION CROSSECTIONS.                            
            EXTPP=ACANS(1,1,2)*RTSFCT/4.0 
            EXTPER=ACANS(1,2,2)*RTSFCT/4.0 
      !        NORMALIZE W.R.T. EQ. SPHERICAL DIA.                            
               EXTPP=EXTPP*CNORM 
               EXTPER=EXTPER*CNORM 
      !     CALCULATE FORWARD AND BACKWARD AMPLITUDE IN FAR ZONE.             
      !     SIGMA EQUALS 4.0/K                                                
            FORRP=SIGMA*ACANS(1,1,1)/RTSFCT 
            FORIP=SIGMA*ACANS(1,1,2)/RTSFCT 
            FORPE=SIGMA*ACANS(1,2,1)/RTSFCT 
            FORIPE=SIGMA*ACANS(1,2,2)/RTSFCT 
            BORRP=SIGMA*ACANS(NUANG,1,1)/RTSFCT 
            BORIP=SIGMA*ACANS(NUANG,1,2)/RTSFCT 
            BORRPE=SIGMA*ACANS(NUANG,2,1)/RTSFCT 
            BORIPE=SIGMA*ACANS(NUANG,2,2)/RTSFCT 
      !     CALCULATE NORMALIZED RADAR CROSSECTIONS FOR BOTH POLARIZATIONS    
            XHOR=ACANS(NUANG,1,1)**2+ACANS(NUANG,1,2)**2 
            YVER=ACANS(NUANG,2,1)**2+ACANS(NUANG,2,2)**2 
      !         NORMALIZE W.R.T. EQ. SPHERICAL DIA.                           
                XHOR=XHOR*CNORM 
                YVER=YVER*CNORM 
      !     PRINT THE SCATTERING RESULTS.                                     
      !     WRITE(110,175) KMV,ANGINC                                           
      ! 175 FORMAT(///,35X,35H********** ACCUMULATED SUMS FOR M =,I3,11H *****
      !    1*****/1H0,42X,20HANGLE OF INCIDENCE =,F6.2,8H DEGREES/1H0,40X,37HD
      !    2IFFERENTIAL SCATTERING CROSS SECTION/1H0,50X,16HAZIMUTHAL  PLANE/1
      !    3H0,9X,9HANGLE,12X,21HVERTICAL POLARIZATION,23X,23HHORIZONTAL POLAR
      !    4IZATION//)                                                        
      !                                                                       
            WRITE(110,555) 
        555 FORMAT(//) 
            WRITE(110,175) KMV,ANGINC 
        175 FORMAT(1H ,35X,35H********** ACCUMULATED SUMS FOR M =,I3,11H *****&
           &*****/1H ,42X,20HANGLE OF INCIDENCE =,F6.2,8H DEGREES)            
      !                                                                       
            WRITE(110,2001) 
      !2001  FORMAT(////,7X,'SCAT',9X,'EXTN',10X,'RADAR',10X,'RE(FOR)',10X,   
      !    1'IM(FOR)',8X,'RE(BACK)',8X,'IM(BACK)',//)                         
       2001 FORMAT(//,10X,'SCAT',10X,'EXTN',10X,'RADAR',10X,'RE(FOR)',8X,     &
           &'IM(FOR)',8X,'RE(BACK)',7X,'IM(BACK)',/)                          
            WRITE(110,2002) SUMM1,EXTPP,XHOR,FORRP,FORIP,BORRP,BORIP 

      !2002  FORMAT(1X,7(3X,E12.5),/,'PARRALLEL POL')                         
            WRITE(110,2003) SUMM2,EXTPER,YVER,FORPE,FORIPE,BORRPE,BORIPE 
            IF (IM.EQ.NM) THEN
                  WRITE(112,2005) SUMM1,EXTPP,XHOR,FORRP,FORIP,BORRP,BORIP 
                  WRITE(111,2005) SUMM2,EXTPER,YVER,FORPE,FORIPE,BORRPE,BORIPE 
            ENDIF

      !2003  FORMAT(1X,7(3X,E12.5),/,'PERPN POL') 
      2002 FORMAT(1X,7(3X,E12.5),10X,'PARALLEL POL')                    
      2003 FORMAT(1X,7(3X,E12.5),10X,'PERPN POL') 
      2005 FORMAT(1X,7(3X,E12.5),10X) 

      !                                                                       
      !                                                                       
      !                                                                       
      !     IF(KMV.EQ.(NM-1)) THEN                                            
      !     WRITE(8,2018)SUMM1,EXTPP,XHOR                                     
      !     WRITE(8,2018)SUMM2,EXTPER,YVER                                    
      !     WRITE(8,2020)FORRP,FORIP,BORRP,BORIP                              
      !     WRITE(8,2020)FORPE,FORIPE,BORRPE,BORIPE                           
      !     ENDIF                                                             
      !                                                                       
       2018 FORMAT(3E15.7) 
       2020 FORMAT(4E15.7) 
            RETURN 
            END                                           
            SUBROUTINE PRINTM(P,N,ND) 
      !     A ROUTINE TO PRINT A MATRIX.                                      
            DIMENSION P(ND,ND,2) 
            DO 20 K=1,2 
            DO 10 I=1,N 
            WRITE(6,100) I,(P(I,J,K),J=1,N) 
         10 END DO 
         20 END DO 
        100 FORMAT (/,' ROWI3',2X,1P8E15.6,/,(10X,8E15.6)) 
            RETURN 
            END                                           
            SUBROUTINE GENKR 
              COMMON DTR,RTD,CPI 
      !     CALCULATE CKR AND DCKR AS A FUNCTION OF THETA FOR A OBLATE SPHEROI
            COMMON /BDYCOM/ DCNR,DCNI,CKPRR,CKPRI,CKR,DCKR,CONK,AOVRB,SIGMA,IB 
            COMMON /THTCOM/ THETA,SINTH,COSTH 
            BOVRA=1.0/AOVRB 
            QB = 1.000/SQRT((BOVRA*COSTH)**2+SINTH**2) 
            CKR = CONK*QB 
            DCKR = CONK*COSTH*SINTH*(BOVRA**2-1.000)*QB**3 
      !       THET1=THETA*RTD                                                 
         10   FORMAT(2X,3(E15.7,4X)) 
            RETURN 
            END                                           
             SUBROUTINE CALENP 
             COMMON DTR,RTD,CPI 
             COMMON/ENDPNT/EPPS(4),NSECT 
             COMMON/BDYCOM/DCNR,DCNI,CKPRR,CKPRI,CKR,DCKR,CONK,AOVRB,SIGMA,IB 
             EPPS(1)=0. 
              EPPS(2)=CPI/2. 
             RETURN 
            END                                           
              SUBROUTINE EPSLON(ALAM,TEMP,REPS,AIMEPS) 
      !     SET UP CONSTANTS.                                                 
            SIGMA=12.5664E08 
            CPI=3.141592 
            EINF=5.27137+.0216474*TEMP-.00131198*(TEMP**2) 
            ALPHA=-16.8129/(TEMP+273.)+.0609265 
            BRIN=2513.98D0/(TEMP+273.000) 
            ALS=.00033836*EXP(BRIN) 
            T=TEMP-25.0 
            EPS=78.54*(1.-4.579E-03*T+1.19E-05*(T**2)-2.8E-08*(T**3)) 
      !     SET UP COMMON PARAMETERS.                                         
            SLAM=(ALS/ALAM)**(1.0-ALPHA) 
            SINAL=SIN(ALPHA*CPI/2.000) 
            COSAL=COS(ALPHA*CPI/2.000) 
      !     CALCULATE RE(EPSLON)                                              
            REPS=(EPS-EINF)*(1.0+SLAM*SINAL) 
            X=1.0+2.0*SLAM*SINAL+SLAM**2 
            REPS=REPS/X 
            REPS=REPS+EINF 
      !     CALCULATE IM(EPSLON)                                              
            AIMEPS=(EPS-EINF)*SLAM*COSAL 
            AIMEPS=(AIMEPS/X)+SIGMA*ALAM/18.8496E10 
            RETURN 
            END                                           
             SUBROUTINE SPICE(ALAM,TEMP,DCNR,DCNI,DCNR2,DCNI2) 
             CALL EPSLON(ALAM,TEMP,DCNR,DCNI) 
             WRITE(110,2) DCNR,DCNI 
          2  FORMAT(1X,'DIELECTRIC CONS.OF WATER=',E15.8,'+J',E15.8) 
             EMX=REPS 
             EMY=-AIMEPS 
      !      ASSIGN THE INCLUSION EPSILON.                                    
      !      VALUE TAKEN FROM TIURI,TRANS. IEEE,GE-20,PAGE:52                 
      !       COMPUTE THE REAL PART OF EPSILON.                               
              DCNR2=1.508**3 
      !       COMPUTE THE IMAGINARY PART OF EPSILON.                          
              DCNI2=(0.34*0.001)/((1.0-0.417)**2.0) 
             WRITE(110,30) DCNR2,DCNI2 
         30  FORMAT(2X,'DIEL. CONST. ICE= ',2(E15.7,4X)) 
              RETURN 
            END                                           
            SUBROUTINE GENBKR(AR,AI,ISWT,ITH,JS) 
      !     4 TH JANUARY 1989                                                 
      !                                                                       
      !                                                                       
             COMMON/MTXCOM/NRANK,N,A(80,80,2),B(80,80,2),CMXNRM(80) 
             COMMON/FNCCOM/PNMLLG(81),BSSLSP(81,31,3),CNEUMN(81,31,3),        &
           & BSLKPR(81,31,3),BSLKPI(81,31,3),CNEUMR(81,31,3),CNEUMI(81,31,3)  
      !        DIMENSION SJR(2000),SJI(2000),SYR(100),SYI(100),SHR(1),SHI(1)  
              DIMENSION SJR(2000),SJI(2000),SYR(100),SYI(100),SHR(2),SHI(2) 
      !       LOGICAL LT,LF                                                   
      !       write(6,'(a6,2e15.7,/)') 'newarg',ar,ai                         
              NAK=2 
      !     WRITE (5,104)                                                     
        104 FORMAT ('SPHERICAL BESSEL FUNCTIONS FOR COMPLEX ARGUMENT') 
            ZERO=0.0D0 
            ONE=1.0D0 
            TWO=2.0D0 
            THREE=3.0D0 
            IZ=0 
            DR=AR*AR-AI*AI 
              DI=TWO*AR*AI 
             CC=TWO 
             EPS=1.0D-16 
             WUNR=ONE 
             WUNI=ZERO 
              CALL DVDD(WUNR,WUNI,DR,DI,T1,T2) 
      !         WRITE(6,*) 'AFT CSP DV1'                                      
              SRARG=SQRT(AR*AR+AI*AI) 
              IF(SRARG.GT.0.5D0)GO TO 29 
              NP=N+1 
              CALL MLTD(AR,AI,AR,AI,ZR,ZI) 
              ZR=ZR/TWO 
              ZI=ZI/TWO 
              FDNM=THREE 
              HDN=ONE 
              HDNM=ONE 
              HDNI=ZERO 
              DO 14 I=1,NP 
               NN=I-1 
              EN=NN 
      !       CALCULATE..                                                     
              IF(NN-1)2,6,3 
          6   FNR=AR/THREE 
              FNI=AI/THREE 
              GO TO 5 
          2   FNR=ONE 
              FNI=ZERO 
              GO TO 5 
          3   CALL MLTD(FNR,FNI,AR,AI,FNR,FNI) 
              FDNM=FDNM+TWO 
              FNR=FNR/FDNM 
              FNI=FNI/FDNM 
          5   ANSR=ONE 
              ANSI=ZERO 
              PANSR=ONE 
              PANSI=ZERO 
              TRM=-ONE 
              TIM=ZERO 
              AB=ONE 
              BA=THREE 
          7   GNU=AB*(TWO*EN+BA) 
              ZRS=-ZR/GNU 
              ZIS=-ZI/GNU 
              CALL MLTD(TRM,TIM,ZRS,ZIS,TRM,TIM) 
              ANSR=ANSR-TRM 
              ANSI=ANSI-TIM 
              IF(ANSR.EQ.ZERO)GO TO 15 
              IF (ANSI.EQ.ZERO)GO TO 16 
            IF (ABS((PANSR-ANSR)/ANSR).LE.EPS.AND.ABS((PANSI-ANSI)/ANSI)      &
           &    .LE.EPS)GO TO 8                                               
              GO TO 17 
         15   IF (ABS((PANSI-ANSI)/ANSI).LE.EPS)GO TO 8 
              GO TO 17 
         16   IF (ABS((PANSR-ANSR)/ANSR).LE.EPS)GO TO 8 
         17   PANSR=ANSR 
              PANSI=ANSI 
              AB=AB+ONE 
              BA=BA+TWO 
              GO  TO 7 
          8   CALL MLTD(FNR,FNI,ANSR,ANSI,SJR(I),SJI(I)) 
      !        CALCULATE...                                                   
              IF(NN-1)4,10,9 
          4   GDR=-ONE 
              GDI=ZERO 
              CALL DVDD(GDR,GDI,AR,AI,HR,HI) 
      !        WRITE(6,*) 'AFT CSP DV2'                                       
              GO TO 11 
         10   HDR=AR 
              HDI=AI 
          9   CALL MLTD(HDR,HDI,AR,AI,HDR,HDI) 
              HDNM=HDNM*HDN 
              HDN=HDN+TWO 
              CALL DVDD(HDNM,HDNI,HDR,HDI,HR,HI) 
      !          WRITE(6,*) 'AFT CSP DV3'                                     
              HR=-HR 
              HI=-HI 
         11   ALSR=ONE 
              ALSI=ZERO 
              PALSR=ONE 
              PALSI=ZERO 
              TRN=-ONE 
              TIN=ZERO 
              AC=ONE 
              CA=ONE 
         12   HNU=AC*(CA-TWO*EN) 
              XRS=-ZR/HNU 
              XIS=-ZI/HNU 
              CALL MLTD(TRN,TIN,XRS,XIS,TRN,TIN) 
              ALSR=ALSR-TRN 
              ALSI=ALSI-TIN 
              IF(ALSR.EQ.ZERO)GO TO 18 
              IF(ALSI.EQ.ZERO)GO TO 19 
            IF(ABS((PALSR-ALSR)/ALSR).LE.EPS.AND.ABS((PALSI-ALSI)/ALSI)       &
           &      .LE.EPS)GO TO 13                                            
              GO TO 20 
         18   IF(ABS((PALSI-ALSI)/ALSI).LE.EPS)GO TO 13 
             GO TO 20 
         19    IF (ABS((PALSR-ALSR)/ALSR).LE.EPS)GO TO 13 
         20    PALSR=ALSR 
               PALSI=ALSI 
               AC=AC+ONE 
               CA=CA+TWO 
               GO TO 12 
         13    CALL MLTD(HR,HI,ALSR,ALSI,SYR(I),SYI(I)) 
      !                                                                       
      !                                                                       
      !                                                                       
      !                                                                       
      !                                                                       
               IF(NAK.EQ.2)GO TO 14 
               IF(NAK.EQ.5)GO TO 50 
               YRR=SYR(I) 
               YII=SYI(I) 
               IF(AI.LT.ZERO)GO TO 51 
               SYR(I)=SJR(I)-YII 
               SYI(I)=SJI(I)+YRR 
               GOTO 14 
         51   SYR(I)=SJR(I)+YII 
              SYI(I)=SJI(I)-YRR 
              GO TO 14 
         50     IF(AI.LT.ZERO)GO TO 48 
                SHR(I)=SJR(I)-SYI(I) 
                SHI(I)=SJI(I)+SYR(I) 
                GO TO 14 
         48     SHR(I)=SJR(I)+SYI(I) 
              SHI(I)=SJI(I)-SYR(I) 
         14   CONTINUE 
              DO 545 KA=1,N 
              BSLKPR(KA,ITH,JS)=SJR(KA) 
              BSLKPI(KA,ITH,JS)=SJI(KA) 
              CNEUMR(KA,ITH,JS)=SYR(KA) 
              CNEUMI(KA,ITH,JS)=SYI(KA) 
        545   CONTINUE 
              RETURN 
         29   DSN=SIN(AR) 
              DCS=COS(AR) 
              EXYL=EXP(AI) 
              EXYS=EXP(-AI) 
              XSN=AR*DSN 
              XCO=AR*DCS 
              YSN=AI*DSN 
              YCO=AI*DCS 
              ZXY=AR*AR+AI*AI 
              TZXY=TWO*ZXY 
              SJZRL=(XSN+YCO)/TZXY 
              SJZRS=(XSN-YCO)/TZXY 
              SJZIL=(XCO-YSN)/TZXY 
              SJZIS=(-XCO-YSN)/TZXY 
              SYZRL=EXYL*(-SJZIL) 
              SYZRS=EXYS*SJZIS 
              SYZIL=EXYL*SJZRL 
              SYZIS=EXYS*(-SJZRS) 
              SJZRL=EXYL*SJZRL 
              SJZRS=EXYS*SJZRS 
              SJZIL=EXYL*SJZIL 
              SJZIS=EXYS*SJZIS 
              SJR(1)=SJZRL+SJZRS 
              SJI(1)=SJZIL+SJZIS 
             SJR(2)=((AR*SJZRL+AI*SJZIL)/ZXY+SYZRL)+                          &
           &     ((AR*SJZRS+AI*SJZIS)/ZXY+SYZRS)                              
             SJI(2)=((-AI*SJZRL+AR*SJZIL)/ZXY+SYZIL)+                         &
           &      ((-AI*SJZRS+AR*SJZIS)/ZXY+SYZIS)                            
             NHO=0 
             IF(ABS(AI).LT.5.0D0)GO TO 43 
      !      .....                                                            
             NHO=1 
             YEX=EXP(-ABS(AI)) 
             ANUR=YEX*DSN 
             ANUI=YEX*DCS 
             IF(AI.GE.ZERO)ANUI=-ANUI 
             CALL  DVDD(ANUR,ANUI,AR,AI,HRZ,HIZ) 
      !       WRITE(6,*) 'AFT CSP DV4'                                        
             CALL  MLTD(AR,AI,AR,AI,ZSR,ZSI) 
             CALL  DVDD(ANUR,ANUI,ZSR,ZSI,HRW,HIW) 
      !          WRITE(6,*) 'AFT CSP DV5'                                     
             IF(AI)38,39,39 
         38   ANUR=-ANUR 
              GO TO 40 
         39    ANUI=-ANUI 
         40   CALL  DVDD(ANUI,ANUR,AR,AI,HOA,HOB) 
      !          WRITE(6,*) 'AFT CSP DV6'                                     
      ! ....                                                                  
      ! ....                                                                  
      ! ....                                                                  
               IF(NAK.LT.5)GO TO 54 
               SHR(1)=HRZ 
               SHI(1)=HIZ 
               SHR(2)=HRW-HOA 
               SHI(2)=HIW-HOB 
             GO TO 55 
         54  IF (NAK.EQ.2)GO TO 56 
             SYR(1)=HRZ 
             SYI(1)=HIZ 
              SYR(2)=HRW-HOA 
              SYI(2)=HIW-HOB 
              GO TO 36 
         56   HRW=HRW-HOA 
              HIW=HIW-HOB 
              SYR(1)=-SJI(1)+HIZ 
              SYI(1)=SJR(1)-HRZ 
              SYR(2)=-SJI(2)+HIW 
              SYI(2)=SJR(2)-HRW 
              GO TO 57 
         55   SYR(1)=-SJI(1)+SHI(1) 
              SYI(1)=SJR(1)-SHR(1) 
              SYR(2)=-SJI(2)+SHI(2) 
              SYI(2)=SJR(2)-SHR(2) 
         57   IF(AI.GE.ZERO)GO TO 36 
              SYR(1)=-SYR(1) 
              SYI(1)=-SYI(1) 
              SYR(2)=-SYR(2) 
              SYI(2)=-SYI(2) 
             GO TO 36 
      ! ....                                                                  
         43   SYR(1)=SYZRL+SYZRS 
              SYI(1)=SYZIL+SYZIS 
              SYR(2)=((AR*SYZRL+AI*SYZIL)/ZXY-SJZRL)+                         &
           &      ((AR*SYZRS+AI*SYZIS)/ZXY-SJZRS)                             
            SYI(2)=((-AI*SYZRL+AR*SYZIL)/ZXY-SJZIL)+                          &
           &    ((-AI*SYZRS+AR*SYZIS)/ZXY-SJZIS)                              
      !     ...                                                               
      !     ...                                                               
      !     ...                                                               
         42   IF(NAK.EQ.2)GO TO 36 
              IF(NAK.EQ.5)GO TO 52 
              YRZ=SYR(1) 
              YIZ=SYI(1) 
              YRW=SYR(2) 
              YIW=SYI(2) 
              IF(AI.LT.ZERO)GO TO 53 
              SYR(1)=SJR(1)-YIZ 
              SYI(1)=SJI(1)+YRZ 
              SYR(2)=SJR(2)-YIW 
              SYI(2)=SJI(2)+YRW 
              GO TO 36 
         53   SYR(1)=SJR(1)+YIZ 
              SYI(1)=SJI(1)-YRZ 
              SYR(2)=SJR(2)+YIW 
              SYI(2)=SJI(2)-YRW 
              GO TO 36 
         52   IF(AI.LT.ZERO)GO TO 41 
              SHR(1)=SJR(1)-SYI(1) 
              SHI(1)=SJI(1)+SYR(1) 
              SHR(2)=SJR(2)-SYI(2) 
              SHI(2)=SJI(2)+SYR(2) 
              GO TO 36 
         41   SHR(1)=SJR(1)+SYI(1) 
              SHI(1)=SJI(1)-SYR(1) 
              SHR(2)=SJR(2)+SYI(2) 
              SHI(2)=SJI(2)-SYR(2) 
         36   IF(N.GT.1)GOTO 111 
              DO 551 KB=1,N 
              BSLKPR(KB,ITH,JS)=SJR(KB) 
              BSLKPI(KB,ITH,JS)=SJI(KB) 
              CNEUMR(KB,ITH,JS)=SYR(KB) 
              CNEUMI(KB,ITH,JS)=SYI(KB) 
        551   CONTINUE 
              RETURN 
      !       .....                                                           
        111   M=N+1 
      !       .....                                                           
              NN=SRARG+30 
              IF ((N+30).GT.NN)NN=N+30 
              GDR=SJR(2) 
              GDI=SJI(2) 
         30   SJR(NN)=ZERO 
              SJI(NN)=ZERO 
              SJR(NN-1)=1.0D-20 
              SJI(NN-1)=1.0D-20 
              NM=NN-2 
              DO 31 K=2,NM 
              KK=NN-K 
              CALL DVDD(SJR(KK+1),SJI(KK+1),AR,AI,SJR(KK),SJI(KK)) 
      !             WRITE (6,*) 'AFT CSP DV7'                                 
      !       CALL ERRSET(72,LT,LF,LF,LF,)                                    
      !       CALL ERRSNS                                                     
              SJR(KK)=(CC*KK+ONE)*SJR(KK)-SJR(KK+2) 
      !       CALL ERRSNS(NUM,,,)                                             
      !       IF(NUM.EQ.72)GO TO 24                                           
      !       CALL ERRSNS                                                     
              SJI(KK)=(CC*KK+ONE)*SJI(KK)-SJI(KK+2) 
      !     CALL ERRSET(72,LF,LF,LF,LT,)                                      
      !     CALL ERRSNS(NUM,,,)                                               
      !     IF(NUM.EQ.72)GO TO 24                                             
         31 END DO 
            CALL DVDD(GDR,GDI,SJR(2),SJI(2),RAR,RAI) 
      !          WRITE (6,*) 'AFT CSP DV8'                                    
      ! ...                                                                   
      ! ...                                                                   
            IF(RAR.NE.ZERO.AND.RAI.NE.ZERO)GO TO 67 
            IF(ABS(SJR(2)).LT.ABS(SJI(2)))GO TO 68 
            IF(RAR.NE.ZERO)GO TO 69 
            IF(GDR.EQ.ZERO.AND.SJI(2).EQ.ZERO)GO TO 69 
            GO TO 24 
         69   IF(RAI.NE.ZERO)GO TO 67 
              IF(GDI.EQ.ZERO.AND.SJI(2).EQ.ZERO)GO TO 67 
              GO TO 24 
         68   IF(RAR.NE.ZERO)GO TO 70 
              IF(GDI.EQ.ZERO.AND.SJR(2).EQ.ZERO)GO TO 70 
              GO TO 24 
         70  IF(RAI.NE.ZERO)GO TO  67 
               IF(GDR.EQ.ZERO.AND.SJR(2).EQ.ZERO)GO TO 67 
               GO TO 24 
         67    DO 32 K=3,M 
               TR=SJR(K) 
               TI=SJI(K) 
         32    CALL MLTD(TR,TI,RAR,RAI,SJR(K),SJI(K)) 
               SJR(2)=GDR 
               SJI(2)=GDI 
      ! ..                                                                    
               IF(NHO.EQ.1)GO TO 44 
      !                                                                       
      !                                                                       
      !                                                                       
               IF(NAK.EQ.3)GO TO 66 
         22    DO 23 K=3,M 
             CALL DVDD(SYR(K-1),SYI(K-1),AR,AI,SYR(K),SYI(K)) 
      !          WRITE (6,*) 'AFT CSP DV9'                                    
             SYR(K)=(CC*K-THREE)*SYR(K)-SYR(K-2) 
             SYI(K)=(CC*K-THREE)*SYI(K)-SYI(K-2) 
              IF(NAK.EQ.2)GO TO 23 
               IF(AI.LT.ZERO)GO TO 45 
         47    SHR(K)=SJR(K)-SYI(K) 
               SHI(K)=SJI(K)+SYR(K) 
              GO TO 23 
         45   SHR(K)=SJR(K)+SYI(K) 
              SHI(K)=SJI(K)-SYR(K) 
         23   CONTINUE 
              DO 555 KC=1,N 
              BSLKPR(KC,ITH,JS)=SJR(KC) 
              BSLKPI(KC,ITH,JS)=SJI(KC) 
              CNEUMR(KC,ITH,JS)=SYR(KC) 
              CNEUMI(KC,ITH,JS)=SYI(KC) 
        555   CONTINUE 
              RETURN 
         66   DO 60 K=3,M 
              CALL DVDD(YRW,YIW,AR,AI,YRT,YIT) 
      !            WRITE(6,*) 'AFT CSP DV10'                                  
              YRT=(CC*K-THREE)*YRT-YRZ 
               YIT=(CC*K-THREE)*YIT-YIZ 
               IF(AI.LT.ZERO)GO TO 58 
              SYR(K)=SJR(K)-YIT 
              SYI(K)=SJI(K)-YRT 
              GO TO 59 
         58   SYR(K)=SJR(K)+YIT 
              SYI(K)=SJI(K)-YRT 
         59  YRZ=YRW 
             YIZ=YIW 
             YRW=YRT 
             YIW=YIT 
         60  CONTINUE 
              DO 565 KD=1,N 
              BSLKPR(KD,ITH,JS)=SJR(KD) 
              BSLKPI(KD,ITH,JS)=SJI(KD) 
              CNEUMR(KD,ITH,JS)=SYR(KD) 
              CNEUMI(KD,ITH,JS)=SYI(KD) 
        565   CONTINUE 
             RETURN 
      !                                                                       
      !                                                                       
      !                                                                       
         44  IF(NAK.NE.5)GO TO 61 
             DO 46 K=3,M 
             CALL DVDD(SHR(K-1),SHI(K-1),AR,AI,SHR(K),SHI(K)) 
      !           WRITE (6,*) 'AFT CSP DV11'                                  
             SHR(K)=(CC*K-THREE)*SHR(K)-SHR(K-2) 
              SHI(K)=(CC*K-THREE)*SHI(K)-SHI(K-2) 
              SYR(K)=-SJI(K)+SHI(K) 
              SYI(K)=SJR(K)-SHR(K) 
              IF(AI.GE.ZERO)GO TO 46 
              SYR(K)=-SYR(K) 
              SYI(K)=-SYI(K) 
         46     CONTINUE 
              DO 575 KE=1,N 
              BSLKPR(KE,ITH,JS)=SJR(KE) 
              BSLKPI(KE,ITH,JS)=SJI(KE) 
              CNEUMR(KE,ITH,JS)=SYR(KE) 
              CNEUMI(KE,ITH,JS)=SYI(KE) 
        575   CONTINUE 
                RETURN 
         61   IF(NAK.EQ.3)GO TO 62 
              DO 63 K=3 ,M 
              CALL DVDD(HRW,HIW,AR,AI,HRT,HIT) 
      !            WRITE(6,*) 'AFT CSP DV12'                                  
              HRT=(CC*K-THREE)*HRT-HRZ 
              HIT=(CC*K-THREE)*HIT-HIZ 
            SYR(K)=-SJI(K)+HIT 
            SYI(K)=SJR(K)-HRT 
             IF(AI.GE.ZERO)GO TO 64 
              SYR(K)=-SYR(K) 
              SYI(K)=-SYI(K) 
         64   HRZ=HRW 
              HIZ=HIW 
              HRW=HRT 
              HIW=HIT 
         63   CONTINUE 
              DO 585 KF=1,N 
              BSLKPR(KF,ITH,JS)=SJR(KF) 
              BSLKPI(KF,ITH,JS)=SJI(KF) 
              CNEUMR(KF,ITH,JS)=SYR(KF) 
              CNEUMI(KF,ITH,JS)=SYI(KF) 
        585   CONTINUE 
              RETURN 
         62   DO 65 K=3,M 
              CALL DVDD(SYR(K-1),SYI(K-1),AR,AI,SYR(K),SYI(K)) 
      !           WRITE (6,*) 'AFT CSP DV13'                                  
              SYR(K)=(CC*K-THREE)*SYR(K)-SYR(K-2) 
              SYI(K)=(CC*K-THREE)*SYI(K)-SYI(K-2) 
         65   CONTINUE 
              DO 595 KH=1,N 
              BSLKPR(KH,ITH,JS)=SJR(KH) 
              BSLKPI(KH,ITH,JS)=SJI(KH) 
              CNEUMR(KH,ITH,JS)=SYR(KH) 
              CNEUMI(KH,ITH,JS)=SYI(KH) 
        595   CONTINUE 
              RETURN 
         24   NN=NN-1 
              WRITE(6,26)NN 
         26   FORMAT (1X,'NN REDUCED TO',I6) 
              IZ=IZ+1 
              IF(IZ.LE.25)GOTO 211 
              DO 591 KG=1,N 
              BSLKPR(KG,ITH,JS)=SJR(KG) 
              BSLKPI(KG,ITH,JS)=SJI(KG) 
              CNEUMR(KG,ITH,JS)=SYR(KG) 
              CNEUMI(KG,ITH,JS)=SYI(KG) 
        591   CONTINUE 
              RETURN 
        211   GO TO 30 
            END                                           
      !                                                                       
              SUBROUTINE DVDD(XA,YA,XB,YB,XC,YC) 
      !                                                                       
      !                                                                       
      !       IMPLICIT REAL*8(A-H,O-Z)                                        
      !      LOGICAL LT,LF                                                    
      !      DATA LT/.TRUE./,LF/.FALSE./                                      
             ZERO=0.0D0 
             IF(XB.NE.ZERO.OR.YB.NE.ZERO)GO TO 3 
             WRITE(6,100) 
        100  FORMAT('BOTH REAL AND IMAGINARY PARTS OF DENOMINATOR ARE ZERO') 
             RETURN 
      !        WRITE(6,500)                                                   
      !        FORMAT('ERROR 1 IN DVDD')                                      
      !      CALL ERRSET(73,LT,LF,LF,LF,)                                     
      !      CALL ERRSET(74,LT,LF,LF,LF,)                                     
          3  DENOM=XB*XB+YB*YB 
      !         IF(DENOM.EQ.ZERO) RETURN                                      
             IF(DENOM.EQ.ZERO)GO TO 1 
             XX=(XA/DENOM)*XB+(YA/DENOM)*YB 
             IF(XX.EQ.ZERO)GO TO 1 
             YC=YA*(XB/DENOM)-(XA/DENOM)*YB 
             IF(YC.EQ.ZERO)GO TO 1 
             XC=XX 
              RETURN 
      !            WRITE(6,501)                                               
      !        FORMAT('ERROR 3 IN DVDD')                                      
      !       CALL ERRSET(73,LF,LF,LF,LT,)                                    
      !       CALL ERRSET(74,LF,LF,LF,LT,)                                    
          1      WRITE(6,*) 'OVERFLOW,DENOM=0' 
             IF(ABS(XB).LT.ABS(YB))GO TO 2 
          8  DC=YB/XB 
              AC=XA/XB 
              BC=YA/XB 
      !       CALL ERRSET(74,LT,LF,LF,LF,)                                    
              DENOM=1.0D0+DC*DC 
      !                                                                       
             XC=(AC+BC*DC)/DENOM 
                                                                              
             YC=(BC-AC*DC)/DENOM 
      !      CALL ERRSET(74,LF,LF,LF,LT,)                                     
             RETURN 
          2 AD=XA/YB 
             CD=XB/YB 
             BD=YA/YB 
      !      CALL ERRSET(74,LT,LF,LF,LF,)                                     
             DENOM=1.0D0+CD*CD 
      !                                                                       
             XC=(BD+AD*CD)/DENOM 
             YC=(-AD+BD*CD)/DENOM 
      !      CALL ERRSET(74,LF,LF,LF,LT,)                                     
             RETURN 
            END                                           
      !     $                                                                 
      !****************************************** RAVI, NASA/MSFC. 12/17/92   
      !                                                                       
             SUBROUTINE MLTD(XA,YA,XB,YB,XC,YC) 
      !      IMPLICIT REAL*8 (A-H,O-Z)                                        
             XX=XA*XB-YA*YB 
            YC=XA*YB+YA*XB 
             XC=XX 
             RETURN 
            END                                           
      
