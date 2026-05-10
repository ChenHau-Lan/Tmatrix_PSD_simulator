                                                                        
       PROGRAM MUELLR 
!                                                                       
!      !!! LAST MODIFICATION: RAVI: NASA/MSFC: 1/5/93 !!!               
!                                                                       
!                                                                       
!                            *** NOTE ***                               
!      !!! FOR SINGLE PARTICLE SCATTERING CHARACTERISTICS ONLY !!!      
!                                                 NASA/MSFC: RAVI: 12/21
!      !!! CODE MODIFIED TO HANDLE NRANK <= 40 AND MODES <= 20 !!!      
!                                                 NASA/MSFC: RAVI: 12/28
!                                                                       
! ... ................................................................. 
!     CALCULATE MUELLER MATRIX (AVERAGED STOKES MATRIX) FOR A DISTRIBUTI
!     OF SCATTERERS WITH A PRESCRIBED ORIENTATION DISTRIBUTION.         
!                                                                       
!     ALSO CALCULATES THE RADAR OBSERVABLES FROM THE MUELLER MATRICES   
!     FOR THE ENSEMBLE OF PARTICLES.                                    
!                                                 NASA/MSFC: RAVI: 12/21
! ... ................................................................. 
! *** CAVEAT: MAX PARMATERS THIS PROGRAM CAN HANDLE DUE                 
! ***         TO MEMORY LIMITATIONS ARE : NDSD = 180                    
! ***                                     NELANG = 10                   
!                                                            RAVI: 12/30
! ... ................................................................. 
!                                                                       
      PARAMETER (NMU=4) 
!                                                                       
      INTEGER FLAG, DISTYP, COUNT, NTHETA, NPHI, CHECK, ROW, NANG 
      INTEGER IZ, IR, IC 
      INTEGER DTYP1, DTYP2, DTYP3, DTYP4, DTYP5, DTYP6, DTYP7, DTYP8,   &
     &        DTYP9, DTYP10, DTYP11, DTYP12, DTYP13, DTYP14             
!                                                            
      character*99 MULLERFILE, TMATFILE, FOLDERNAME, OUTPUTFILE,        &
      & OUTPUTFILE2, OUTPUTFILE3       
      REAL LOWER, UPPER, THETAM, SUM, THETA, SHOP, PNORM, ELVANG 
      REAL ANGINT, OLDTHT, MEAN, SIGMA,  TCL, TCU, PROB 
      REAL LAMBDA, EASTRT, EASTOP, EAINC, GP, LP, KAPPA, THETA0, PHI0 
      REAL DSTEP, DMIN, DMAX 
      REAL FP, PHNORM, ABI, ABQ, ABV, ABH 
!                                                                       
      CHARACTER*20 FNAME 
      COMPLEX TMAT, ACANS 
      LOGICAL RDROP, FRDROP, DHAIL, SPHAIL, WHAIL, DGRAUP, SPGRAUP,     &
     &        WGRAUP, NCRYST, PCRYST, SPAGG, AGG2, AGG5, AGG8           
!                                                                       
      COMMON DEGRAD, RADDEG, PI 
      COMMON /MTXCOM/ NRANK,NRANKI,TMAT(20,80,80),CMXNRM(80) 
      COMMON /FNCCOM/ PNMLLG(81) 
      COMMON /CMVCOM/ NM,KMV,CMI(20),CMV,CM2,TWM,PRODM 
      COMMON /UVCCOM/ ANGINT,ANGINP,UTHETA,UPHI,RTSFCT,IP 
      COMMON /ANGCOM/ TANG(91),PANG(91),COSB(91),SINB(91),              &
     &                ACANS(91,2,2),SCV,SCH,EXV,EXH                     
      COMMON /MAIN01/ ABV, ABH, CONK, S(91,4,4), SS(91,4,4), LAMBDA,    &
     &EXM(4,4),EXMP(4,4),EXMPT(4,4),SCHPT,SCVPT,EXHPT,EXVPT,PHNORM      
      COMMON /MAIN02/ NTHETA, TASC(20), TWT(20), NPHI, PASC(20), PWT(20) 
      COMMON /MAIN03/ DISTYP, THETAM, PNORM, MEAN, SIGMA, TCL, TCU 
      COMMON /MAIN04/ NANG, NUANG, OLDTHT, CHECK 
      COMMON /MAIN05/ FNAME 
      COMMON /MAIN06/ EASTRT, EASTOP, EAINC, ELVANG, KAPPA, THETA0, PHI0 
      COMMON/INT/QMU(NMU),QWT(NMU),QMUE(2*NMU),QWTE(2*NMU) 
      COMMON /DSD1/ DMIN, DMAX, DSTEP 
      COMMON /DSD2/ DMIN1, DMIN2, DMIN3, DMIN4, DMIN5, DMIN6, DMIN7,    &
     &       DMIN8, DMIN9, DMIN10, DMIN11, DMIN12, DMIN13, DMIN14,      &
     &              DMAX1, DMAX2, DMAX3, DMAX4, DMAX5, DMAX6, DMAX7,    &
     &       DMAX8, DMAX9, DMAX10, DMAX11, DMAX12, DMAX13, DMAX14,      &
     &              DSTEP1, DSTEP2, DSTEP3, DSTEP4, DSTEP5, DSTEP6,     &
     &              DSTEP7, DSTEP8, DSTEP9, DSTEP10, DSTEP11, DSTEP12,  &
     &              DSTEP13, DSTEP14                                    
      COMMON /DSD3/ DPAR1, DPAR2, DPAR3 
      COMMON /DSD4/ S1PAR1, S1PAR2, S1PAR3, S2PAR1, S2PAR2, S2PAR3,     &
     &                 S3PAR1, S3PAR2, S3PAR3, S4PAR1, S4PAR2, S4PAR3,  &
     &                 S5PAR1, S5PAR2, S5PAR3, S6PAR1, S6PAR2, S6PAR3,  &
     &                 S7PAR1, S7PAR2, S7PAR3, S8PAR1, S8PAR2, S8PAR3,  &
     &                 S9PAR1, S9PAR2, S9PAR3, S10PAR1, S10PAR2,        &
     &                 S10PAR3, S11PAR1, S11PAR2, S11PAR3, S12PAR1,     &
     &                 S12PAR2, S12PAR3, S13PAR1, S13PAR2, S13PAR3,     &
     &                 S14PAR1, S14PAR2, S14PAR3                        
      COMMON /DSD5/ RDROP, FRDROP, DHAIL, SPHAIL, WHAIL, DGRAUP,        &
     & SPGRAUP, WGRAUP, NCRYST, PCRYST, SPAGG, AGG2, AGG5, AGG8         
      COMMON /DSD6/ DEQ 
      COMMON WX(0:1000) 
!                                                                       
      DIMENSION IA(2),SP(91,4,4) 
! ... ................................................................. 
! read the filename
      read(*,*),FOLDERNAME
      MULLERFILE= TRIM(FOLDERNAME)//'/mueller.inp'
      TMATFILE= TRIM(FOLDERNAME)//'/out1_tmat'
      OUTPUTFILE= TRIM(FOLDERNAME)//'/out1_musingle'
      OUTPUTFILE2= TRIM(FOLDERNAME)//'/out1_Smatrix'
      OUTPUTFILE3= TRIM(FOLDERNAME)//'/out1_Kmatrix'
!                                                                       
! *** OPEN INPUT DATA FILE                                              
      OPEN(12,FILE=MULLERFILE) 
! *** OPEN OUTPUT RESULT FILE                                           
      OPEN(7,FILE=OUTPUTFILE) 
      OPEN(108,FILE=OUTPUTFILE2)
      OPEN(109,FILE=OUTPUTFILE3)


!                                                                       
! *** THE DIAGNOSTIC FILE IS UNIT # 9. THE WRITE(9,*) AND WRITE (9,*)   
! *** HAVE BEEN COMMENTED AT THE PRESENT TIME.                          
! *** THE WRITE(3,*) AND WRITE (3,*) TO UNIT # 3 HAVE ALSO BEEN         
! *** COMMENTED AT THE PRESENT TIME.                                    
! ***                                                    RAVI: 12/30/91 
! ... ................................................................. 
! ----------------------------------------------------------------------                                                    
! *** READ INPUT FROM mueller.inp                                       
! *** ALL DIAMETERS IN UNITS OF CMS.                                    
      READ(12,*)CHECK 
      READ(12,*)EASTRT 
      READ(12,*)EASTOP 
      READ(12,*)EAINC 
      READ(12,*)NTHETA 
      READ(12,*)NPHI 
      READ(12,*)NDSD 
      READ(12,*)NPART 
      READ(12,*)RDROP,DMIN1,DMAX1,DSTEP1,DTYP1,S1PAR1,S1PAR2,S1PAR3 
      READ(12,*)FRDROP,DMIN2,DMAX2,DSTEP2,DTYP2,S2PAR1,S2PAR2,S2PAR3 
      READ(12,*)DHAIL,DMIN3,DMAX3,DSTEP3,DTYP3,S3PAR1,S3PAR2,S3PAR3 
      READ(12,*)SPHAIL,DMIN4,DMAX4,DSTEP4,DTYP4,S4PAR1,S4PAR2,S4PAR3 
      READ(12,*)WHAIL,DMIN5,DMAX5,DSTEP5,DTYP5,S5PAR1,S5PAR2,S5PAR3 
      READ(12,*)DGRAUP,DMIN6,DMAX6,DSTEP6,DTYP6,S6PAR1,S6PAR2,S6PAR3 
      READ(12,*)SPGRAUP,DMIN7,DMAX7,DSTEP7,DTYP7,S7PAR1,S7PAR2,S7PAR3 
      READ(12,*)WGRAUP,DMIN8,DMAX8,DSTEP8,DTYP8,S8PAR1,S8PAR2,S8PAR3 
      READ(12,*)NCRYST,DMIN9,DMAX9,DSTEP9,DTYP9,S9PAR1,S9PAR2,S9PAR3 
      READ(12,*)PCRYST,DMIN10,DMAX10,DSTEP10,DTYP10,S10PAR1,S10PAR2,    &
     &          S10PAR3                                                 
      READ(12,*)SPAGG,DMIN11,DMAX11,DSTEP11,DTYP11,S11PAR1,S11PAR2,     &
     &          S11PAR3                                                 
      READ(12,*)AGG2,DMIN12,DMAX12,DSTEP12,DTYP12,S12PAR1,S12PAR2,      &
     &          S12PAR3                                                 
      READ(12,*)AGG5,DMIN13,DMAX13,DSTEP13,DTYP13,S13PAR1,S13PAR2,      &
     &          S13PAR3                                                 
      READ(12,*)AGG8,DMIN14,DMAX14,DSTEP14,DTYP14,S14PAR1,S14PAR2,      &
     &          S14PAR3                                                 
! ... ................................................................. 
!                                                                       
! *** NELANG IS THE # OF ELEVATION ANGLES                               
      NELANG = NINT(((EASTOP-EASTRT)/EAINC)) + 1 
      PRINT*,'*** # OF ELEVATION ANGLES = ',NELANG 
!                                                                       
! ... ................................................................. 
!                                                                       
!     WRITE THE HEADER FOR THE OUTPUT FILE                              
      WRITE(7,FMT ='(I22,8X,"# OF PARTICLE TYPES")')NPART 
      WRITE(7,FMT ='(F22.2,8X,"START ELEVATION ANGLE DEG")')EASTRT 
      WRITE(7,FMT ='(F22.2,8X,"STOP ELEVATION ANGLE DEG")')EASTOP 
      WRITE(7,FMT ='(F22.2,8X,"INCREMENTAL ANGLE DEG")')EAINC 
      WRITE(7,FMT ='(3X,"PARTICLE TYPE CONSIDERED")') 
!                                                                       
  555 CONTINUE 
!                                                                       
! ... ................................................................. 
      IF (RDROP) THEN 
        OPEN(2,FILE=TMATFILE) 
        DISTYP = DTYP1 
        PRINT*,'*** PROCESSING BEING DONE FOR RAINDROPS' 
        WRITE(7,FMT ='(4X,"RAINDROPS")') 
      ELSEIF (FRDROP) THEN 
        OPEN(2,FILE=TMATFILE) 
        DISTYP = DTYP2 
        PRINT*,'*** PROCESSING BEING DONE FOR FROZEN RAINDROPS' 
        WRITE(7,FMT ='(4X,"FROZEN RAINDROPS")') 
      ELSEIF (DHAIL) THEN 
        OPEN(2,FILE=TMATFILE) 
        DISTYP = DTYP3 
        PRINT*,'*** PROCESSING BEING DONE FOR HAIL ( ICE )' 
        WRITE(7,FMT ='(4X,"DRY HAIL")') 
      ELSEIF (SPHAIL) THEN 
        OPEN(2,FILE=TMATFILE) 
        DISTYP = DTYP4 
        PRINT*,'*** PROCESSING BEING DONE FOR HAIL ( SPONGY )' 
        WRITE(7,FMT ='(4X,"SPONGY HAIL")') 
      ELSEIF (WHAIL) THEN 
        OPEN(2,FILE=TMATFILE) 
        DISTYP = DTYP5 
        PRINT*,'*** PROCESSING BEING DONE FOR HAIL ( WET )' 
        WRITE(7,FMT ='(4X,"WET HAIL")') 
      ELSEIF (DGRAUP) THEN 
        OPEN(2,FILE=TMATFILE) 
        DISTYP = DTYP6 
        PRINT*,'*** PROCESSING BEING DONE FOR GRAUPEL ( ICE )' 
        WRITE(7,FMT ='(4X,"DRY GRAUPEL")') 
      ELSEIF (SPGRAUP) THEN 
        OPEN(2,FILE=TMATFILE) 
        DISTYP = DTYP7 
        PRINT*,'*** PROCESSING BEING DONE FOR GRAUPEL ( SPONGY )' 
        WRITE(7,FMT ='(4X,"SPONGY GRAUPEL")') 
      ELSEIF (WGRAUP) THEN 
        OPEN(2,FILE=TMATFILE) 
        DISTYP = DTYP8 
        PRINT*,'*** PROCESSING BEING DONE FOR GRAUPEL ( WET )' 
        WRITE(7,FMT ='(4X,"WET GRAUPEL")') 
      ELSEIF (NCRYST) THEN 
        OPEN(2,FILE=TMATFILE) 
        DISTYP = DTYP9 
        PRINT*,'*** PROCESSING BEING DONE FOR CRYSTALS ( NEEDLES )' 
        WRITE(7,FMT ='(4X,"NEEDLES CRYSTALS")') 
      ELSEIF (PCRYST) THEN 
        OPEN(2,FILE=TMATFILE) 
        DISTYP = DTYP10 
        PRINT*,'*** PROCESSING BEING DONE FOR CRYSTALS ( PLATES )' 
        WRITE(7,FMT ='(4X,"PLATES CRYSTALS")') 
      ELSEIF (SPAGG) THEN 
        OPEN(2,FILE=TMATFILE) 
        DISTYP = DTYP11 
        PRINT*,'*** PROCESSING BEING DONE FOR AGGREGATES ( SPONGY )' 
        WRITE(7,FMT ='(4X,"SNOW: SPONGY AGGREGATES")') 
!       WRITE(7,FMT ='(4X,"SPONGY AGGREGATES")')                        
      ELSEIF (AGG2) THEN 
        OPEN(2,FILE=TMATFILE) 
        DISTYP = DTYP12 
        PRINT*,'*** PROCESSING BEING DONE FOR AGGREGATES ( RHO=0.2 )' 
        WRITE(7,FMT ='(4X,"SNOW AGGREGATES RHO=0.2")') 
      ELSEIF (AGG5) THEN 
        OPEN(2,FILE=TMATFILE) 
        DISTYP = DTYP13 
        PRINT*,'*** PROCESSING BEING DONE FOR AGGREGATES ( RHO=0.5 )' 
        WRITE(7,FMT ='(4X,"AGGREGATES RHO=0.5")') 
      ELSEIF (AGG8) THEN 
        OPEN(2,FILE=TMATFILE) 
        DISTYP = DTYP14 
        PRINT*,'*** PROCESSING BEING DONE FOR AGGREGATES ( RHO=0.8 )' 
        WRITE(7,FMT ='(4X,"AGGREGATES RHO=0.8")') 
      ENDIF 
! ... ................................................................. 
!                                                                       
! *** THIS ROUTINE SETS PARAMETERS THAT ARE REQUIRED FOR PROCESSING     
! *** DIFFERENT TYPES ( SPECIES ) OF PARTICLES.          RAVI: 12/26/91 
      CALL SETPAR(DISTYP) 
      ! OUTPUT: DMIN, DMAX, DSTEP, DPAR1, DPAR2, DPAR3                  
! ... ................................................................. 
!                                                                       
! *** NDS IS NUMBER OF INPUT T-MATRICES, THE FIRST LINE IN THE          
! *** INPUT T-MATRIX FILE (tcb_*** RESIDING ON THE CRAY MSS)            
      READ(2,*) NDS 
! *** RECALCULATE NDS BASED ON THE DMIN, DMAX AND DSTEP REQUIRED        
! *** FOR EACH SPECIES.                                RAVI: 12/20/91   
      ! NDS = ((DMAX-DMIN)/DSTEP) + 1 
      PRINT*,'*** THE # OF T-MATRICES TO BE PROCESSED ARE ',NDS 
!                                                                       
! ... ................................................................. 
! *** THE 6001 LOOP CALCULATES MUELLER MATRICES CORRESPONDING           
! *** TO EACH T-MATRIX.                                                 
!                                                                       
      DO 6001 NDSET=1,NDS 
!                                                                       
!                                                                       
! *** IANG COUNTS THE DIFFERENT ELEVATION ANGLES.                       
! *** THE REQUIRED INCREMENTS ARE SPECIFIED IN THE INPUT FILE           
      IANG = 1 
!                                                                       
! ... ................................................................. 
! *** INITIALIZE CONSTANTS AND READ IN T-MATRIX.                        
! *** NOTE: INIT WILL NOT RETURN TO THIS PROGRAM UNLESS THE PARTICLE DIA
! ***       IN THE INPUT T-MATRIX FILE IS: DMIN <= DEQ <= DMAX          
! ***                                                    RAVI: 12/30/91 
      CALL INIT 
      !GET OUTPUT: TMAT & NRANK, NM, CONK, LAMBDA, DEQ, AOVERB, EPR, EPI
      IF(NDSET .EQ. 1) THEN 
        WRITE(7,FMT ='(E22.13,8X,"WAVELENGTH CM")')LAMBDA 
      ENDIF 
! ... ................................................................. 
!                                                                       
! *** GET USER INPUT. THIS ROUTINE SELECTS THE SPECIFIED PARTICLE       
! *** ORIENTATION DISTRIBUTION. REQUIRED INPUT TO THE ROUTINE IS        
! *** AVAILABLE IN THE MAIN PROGRAM AND IS PASSED THRU' COMMON BLKS.    
! ***                                                    RAVI: 12/30/91 
      CALL USERIN 
      ! OUTPUT: DISTYP, THETAM, PNORM, MEAN, SIGMA, TCL, TCU            
! ... ................................................................. 
!                                                                       
! *** SET UP LOOP OVER ELEVATION ANGLES                                 
! *** RESET ELVANG TO EASTRT FOR EACH T-MATRIX.          RAVI: 12/23/91 
      IF (IANG .EQ. 1 ) THEN 
        ELVANG = EASTRT 
      ENDIF 
!                                                                       
! ... ................................................................. 
! *** THIS LOOPS OVER THE ELEVATION ANGLES                              
   27 CONTINUE 
!                                                                       
! 25     WRITE (*,907) ELVANG                                           
!                                                                       
!        INITIALIZE SOME VARIABLES                                      
         SCVPT = 0.0 
         SCHPT = 0.0 
         EXVPT = 0.0 
         EXHPT = 0.0 
                                                                        
         DO 100 IZ = 1,NANG 
            DO 75 IR = 1,4 
               DO 50 IC = 1,4 
                  SS(IZ,IR,IC) = 0.0 
   50          CONTINUE 
   75       CONTINUE 
  100    CONTINUE 
                                                                        
          DO  110  IR =1,4 
              DO  110 IC = 1,4 
  110             EXMPT(IR,IC) = 0. 
                                                                        
  ! FINISHED INITIALIZING VARIABLES                                     
  !---------------------------------------------------------------------
                                                                        
!        OUTER LOOP:  INTEGRATE OVER THETA                              
         DO 600 ITH = 1,NTHETA 
!           CALCULATE THETA OF CURRENT PARTICLE ORIENTATION             
            ANGINT = 2.0*ATAN(SQRT(1.0-TASC(ITH)**2)/                   &
     &               (1.0 + TASC(ITH)))*RADDEG                          
            OLDTHT = ANGINT 
!            WRITE (*,900) ANGINT                                       
            IF (CHECK .EQ. 1) THEN 
!              WRITE (9,901)                                            
!              WRITE (9,902)                                            
            ENDIF 
                                                                        
!           ZERO AVERAGE MUELLER-MATRIX ACCUMULATOR                     
            DO 200 IZ = 1,NANG 
               DO 150 IR = 1,4 
                  DO 125 IC = 1,4 
                     S(IZ,IR,IC) = 0.0 
  125             CONTINUE 
  150          CONTINUE 
  200       CONTINUE 
                                                                        
            DO  115  IR = 1,4 
                DO  115  IC = 1,4 
  115       EXMP(IR,IC) = 0. 
                                                                        
                                                                        
            SCVP = 0.0 
            SCHP = 0.0 
            EXVP = 0.0 
            EXHP = 0.0 
                                                                        
!           INNER LOOP:  INTEGRATE OVER PHI                             
            DO 300 IPHI = 1,NPHI 
!              CALCULATE PHI OF CURRENT PARTICLE ORIENTATION; RECOVER   
!              UNROTATED THETA OF PARTICLE ORIENTATION                  
               ANGINP = PASC(IPHI) 
               OLDPHI = ANGINP 
               ANGINT = OLDTHT 
                                                                        
!              CALCULATE SCATTERED FIELD FOR CURRENT THETA, PHI ORIENTAT
!              OF PARTICLE                                              
               CALL SCFLD(ITH,IPHI) 
               ! GET OUTPUT ACANS, SCV, SCH, EXV, EXH                   
                                                                        
!              ACCUMULATE TOTAL SCATTERING IN MUELLER-MATRIX ACCUMULATOR
               CALL MMCALC(IPHI,SCHP, SCVP, EXHP, EXVP) 
               ! GET F(4,4), S(91,4,4), SCHP, SCVP, EXHP, EXVP          
  300       CONTINUE 
                                                                        
!           ACCUMULATE PHI-AVERAGED CONTRIBUTION TO CURRENT THETA TERM: 
!           GET ORIENTATION PROBABILITY DENSITY FUNCTION, TAKING CARE   
!           TO USE THE ORIGINAL (NOT ROTATED) THETA AND PHI.            
            IF (DISTYP .EQ. 0) THEN 
!              P(THETA) = 0.5; P(PHI) = 1.0/360.0 (BOTH RANDOM)         
               PROB = 0.5/360.0 
            ELSE IF (DISTYP .EQ. 1) THEN 
!              SIMPLE HARMONIC OSCILLATOR MODEL                         
               PROB = SHOP(OLDTHT, MEAN, THETAM, OLDPHI, PNORM) 
!              ADJUSTMENTS IF THETA IS IN AN AMBIGUOUS ZONE             
                                                                        
               IF (OLDTHT .LT. TCU) THEN 
                  PROB = PROB +                                         &
     &            SHOP(-OLDTHT,MEAN,THETAM,(OLDPHI+180.0),PNORM)        
               ELSE IF (OLDTHT .GT. TCL) THEN 
                  PROB = PROB +                                         &
     &            SHOP((360.0-OLDTHT),MEAN,THETAM,(OLDPHI+180.0),PNORM) 
               ENDIF 
            ELSE IF (DISTYP .EQ. 2) THEN 
!              2-SIGMA-TRUNCATED GAUSSIAN MODEL                         
               PROB = GP(OLDTHT, MEAN, SIGMA, OLDPHI, PNORM) 
!              ADJUSTMENTS IF THETA IS IN AN AMBIGUOUS ZONE             
               IF (OLDTHT .LT. TCU) THEN 
                  PROB = PROB +                                         &
     &            GP(-OLDTHT,MEAN,SIGMA,(OLDPHI+180.0),PNORM)           
               ELSE IF (OLDTHT .GT. TCL) THEN 
                  PROB = PROB +                                         &
     &            GP((360.0-OLDTHT),MEAN,SIGMA,(OLDPHI+180.0),PNORM)    
               ENDIF 
            ELSE IF (DISTYP .EQ. 3) THEN 
!              LANGEVIN UNIMODAL DISTRIBUTION                           
               PROB = LP(KAPPA, OLDTHT, OLDPHI, THETA0, PHI0, PNORM) 
            ELSE IF (DISTYP .EQ. 4) THEN 
!              FISHER (MEAN THETA = 0) DISTRIBUTION                     
               PROB = FP(KAPPA, OLDTHT, OLDPHI, PNORM) 
            ENDIF 
!           IF (CHECK .EQ. 1) WRITE (9,903) PROB, PNORM                 
                                                                        
            DO 500 IS = 1,NANG 
               DO 400 IR = 1,4 
                  DO 350 IC = 1,4 
                     SS(IS,IR,IC) = SS(IS,IR,IC) +                      &
     &                              PROB*TWT(ITH)*S(IS,IR,IC)           
  350             CONTINUE 
  400          CONTINUE 
  500       CONTINUE 
               DO 550 IR = 1,4 
                  DO 550 IC = 1,4 
                     EXMPT(IR,IC) = EXMPT(IR,IC) +                      &
     &                              PROB*TWT(ITH)*EXMP(IR,IC)           
  550             CONTINUE 
            SCHPT = SCHPT + PROB*TWT(ITH)*SCHP 
            EXHPT = EXHPT + PROB*TWT(ITH)*EXHP 
            SCVPT = SCVPT + PROB*TWT(ITH)*SCVP 
            EXVPT = EXVPT + PROB*TWT(ITH)*EXVP 
  600    CONTINUE 
                                                                        
!        SCATTERING, EXTINCTION AND ABSORPTION CROSS SECTIONS IN MM^2   
         WVNMBR = 2.*PI/LAMBDA 
         SCHPT = SCHPT*16.0*PI/(WVNMBR*WVNMBR) 
         EXHPT = EXHPT*16.0*PI/(WVNMBR*WVNMBR) 
         SCVPT = SCVPT*16.0*PI/(WVNMBR*WVNMBR) 
         EXVPT = EXVPT*16.0*PI/(WVNMBR*WVNMBR) 
         ABH = EXHPT - SCHPT 
         ABV = EXVPT - SCVPT 
         IF (ABH .LT. 1.0E-08) ABH = 0.0 
         IF (ABV .LT. 1.0E-08) ABV = 0.0 
                                                                        
!        NORMALIZATION OF PHASE MATRIX ELEMENTS                         
                                                                        
                                                                        
         DO 800 K = 1,NANG 
            DO 700 IR = 1,4 
               DO 650 IC = 1,4 
                  SS(K,IR,IC) = SS(K,IR,IC)/((SCHPT+SCVPT)*.5) 
  650          CONTINUE 
  700       CONTINUE 
  800    CONTINUE 
                                                                        
         ABI = ABV + ABH 
         ABQ = ABV - ABH 
!        WRITE (3,8001)ELVANG,ABI,ABQ                                   
 8001    FORMAT(F5.2,2(2X,E15.7,1X)) 
              DO 825 J = 1,NANG 
            IF(J.LE.(NANG/2)) THEN 
!U             QMUE(J)=QMU(J)                                           
!U             QWTE(J)=QWT(J)                                           
             QMUE(J)=QMU(NMU+1-J) 
             QWTE(J)=QWT(NMU+1-J) 
             ELSE 
             QMUE(J)=-QMU(NANG+1-J) 
             QWTE(J)=QWT(NANG+1-J) 
!U              QMUE(J)=-QMU(J-NMU)                                     
!U              QWTE(J)=QWT(J-NMU)                                      
             ENDIF 
!        WRITE(6,*) ACOS(QMUE(J))*180/3.1415927                         
!        WRITE(3,*) ACOS(QMUE(J))*180/3.1415927                         
  825    CONTINUE 
                                                                        
!        PHASE FUNCTION NORMALIZATION                                   
         PHNORM = 0. 
         DO 830  J = 1, NANG 
            PHNORM = PHNORM + SS(J,1,1)*QWTE(J) 
  830    CONTINUE 
                                                                        
!        MULTIPLY BY 2*PI FOR AZIMUTHAL INTEGRATION                     
                                                                        
         PHNORM = PHNORM*2.*PI 
                                                                        
!        DUMP MUELLER MATRIX TO OUTPUT FILE                             
!V         DO 875 J = 1,NANG                                            
!V            WRITE (3,905) RADDEG*ACOS(QMUE(J))                        
!V            DO 850 ROW = 1,4                                          
!V               WRITE (3,906) (SS(J,ROW,L), L = 1,4)                   
!V 850        CONTINUE                                                  
!V 875     CONTINUE                                                     
!                                                                       
!        PRINT *,'ELVANG',ELVANG,'QMUE(NMU+1)',QMUE(NMU+1)              
! ... ................................................................. 
!      CALL ROUTINE TO COMPUTE RADAR OBSERVABLES FOR BACKSCATTER        
       CALL RMP(5) 
! ... ................................................................. 
!                                                                       
!        INCREMENT ELEVATION ANGLE; CHECK FOR LIMIT                     
         PRINT*,'ELEVATION ANGLE = ',ELVANG 
         ELVANG = ELVANG + EAINC 
         IANG = IANG + 1 
!                                                                       
         IF (ELVANG .LE. EASTOP) GO TO 27 
! ... ................................................................. 
!                                                                       
 6001   CONTINUE 
!                                                                       
!                                                                       
      write(*,*) 'Finished the main program!!!!!!!!!!!!'  
      ! STOP 0 
!900  FORMAT (1X,'   INTEGRATING OVER THETA = ',F7.3)                   
  901 FORMAT (/,1X,'THETA',5X,'PHI',4X,'ORIGINAL',2X,'ORIGINAL',2X,     &
     &       'ROTATED',2X,'ROTATED')                                    
  902 FORMAT ('COUNTER',2X,'COUNTER',3X,'THETA',6X,'PHI',6X,            &
     &        'THETA',5X,'PHI',/)                                       
  903 FORMAT (1X,'(ORIENTATION PROBABILITY DENSITY FUNCTION = ',F8.5,   &
     &        ';',/,2X,'NORMALIZATION FACTOR = ',E12.5,')',/)           
  904 FORMAT(' SCAT XSECT = ',1PE15.5,/,'  ABS XSECT = ',1PE15.5,/,     &
     &       '  EXT XSECT = ',1PE15.5,/)                                
  905 FORMAT (1X,F6.2) 
  906 FORMAT (1X,4E12.4) 
!907  FORMAT (1X,'ELEVATION ANGLE: ',F4.1)                              
                                                                        
      END                                           
                                                                        
!__________________________________________________________________     
      SUBROUTINE INIT 
!     DATA INITIALIZATION                                               
! *** THERE ARE SOME CHANGES INTRODUCED HERE       *** RAVI: 12/30/91   
!                                                                       
      INTEGER NM, NRANK, NR2, I, NRANKI 
      REAL DEGRAD, RADDEG, PI, CONK 
      REAL SCAT, EXT, LAMBDA, MODE 
      REAL DEQ, AOVERB, EPR, EPI, TEMP 
      REAL TREAL(80,80), TIMAG(80,80) 
      REAL DMIN, DMAX, DSTEP 
      COMPLEX TMAT, ACANS, A1, A2, A3, A4 
      CHARACTER*20 FNAME 
                                                                        
      COMMON DEGRAD, RADDEG, PI 
      COMMON /MTXCOM/ NRANK,NRANKI,TMAT(20,80,80),CMXNRM(80) 
      COMMON /CMVCOM/ NM,KMV,CMI(20),CMV,CM2,TWM,PRODM 
      COMMON /MAIN01/ ABV, ABH, CONK, S(91,4,4), SS(91,4,4), LAMBDA,    &
     &EXM(4,4),EXMP(4,4),EXMPT(4,4),SCHPT,SCVPT,EXHPT,EXVPT,PHNORM      
      COMMON /DSD1/ DMIN, DMAX, DSTEP 
      COMMON /DSD6/ DEQ 
      COMMON /DIEL/ EPR 
                                                                        
      DEGRAD = 0.017453292519943 
      RADDEG = 57.295779513 
      PI = 3.1415926535898 
!                                                                       
! *** THIS IS TO ONLY PROCESS THE REQUIRED T-MATRICES ( DMIN TO DMAX )  
! ***                                                  RAVI: 12/30/91   
  100 CONTINUE 
!                                                                       
                                                                        
!     READ IN BODY PARAMETERS                                           
       READ(2,*,END=101)NRANK 
       READ(2,*)NM 
       READ(2,*)CONK 
       READ(2,*)LAMBDA 
       READ(2,*)DEQ 
       READ(2,*)AOVERB 
      READ(2,*)SIGMA 
      READ(2,*)EPR 
      READ(2,*)EPI 
!     WRITE(3,9001)NRANK                                                
 9001 FORMAT('C  ',I5,5X,'NRANK') 
!     WRITE(3,9002)NM                                                   
 9002 FORMAT('C  ',I5,5X,'MODES') 
!     WRITE(3,9003)CONK                                                 
 9003 FORMAT('C  ',E15.8,5X,'KB') 
!     WRITE(3,9004)LAMBDA                                               
 9004 FORMAT('C  ',E15.8,5X,'WAVELENGTH IN CM') 
!     WRITE(3,9005)DEQ                                                  
 9005 FORMAT('C  ',E15.8,5X,'EQUIVALENT DIAMETER IN CM') 
!     WRITE(3,9006)AOVERB                                               
 9006 FORMAT('C  ',E15.8,5X,'AXIS RATIO') 
!     WRITE(3,9007)EPR                                                  
 9007 FORMAT('C  ',E15.8,5X,'REAL PART OF DIELECTRIC CONSTANT') 
!     WRITE(3,9008)EPI                                                  
 9008 FORMAT('C  ',E15.8,5X,'IMAG PART OF DIELECTRIC CONSTANT') 
!     WRITE(3,'(A)') 'C  ANGLE          ABSORP (I, Q)'                  
!     WRITE(3,'(A)') 'C                 EXTINCTION (4 X 4)'             
!     WRITE(3,'(A)') 'C                 FWDSCATTER (4 X 4)'             
!     WRITE(3,'(A)') 'C                 BACKSCATTER (4 X 4)'            
                                                                        
!     WRITE(9,9010)NRANK                                                
 9010 FORMAT('C  ',I5,5X,'NRANK') 
!     WRITE(9,9011)NM                                                   
 9011 FORMAT('C  ',I5,5X,'MODES') 
!     WRITE(9,9012)CONK                                                 
 9012 FORMAT('C  ',E15.8,5X,'KB') 
!     WRITE(9,9013)LAMBDA                                               
 9013 FORMAT('C  ',E15.8,5X,'WAVELENGTH IN CM') 
!     WRITE(9,9014)DEQ                                                  
 9014 FORMAT('C  ',E15.8,5X,'EQUIVALENT DIAMETER IN CM') 
!     WRITE(9,9015)AOVERB                                               
 9015 FORMAT('C  ',E15.8,5X,'AXIS RATIO') 
!     WRITE(9,9016)EPR                                                  
 9016 FORMAT('C  ',E15.8,5X,'REAL PART OF DIELECTRIC CONSTANT') 
!     WRITE(9,9017)EPI                                                  
 9017 FORMAT('C  ',E15.8,5X,'IMAG PART OF DIELECTRIC CONSTANT') 
      IF ((NRANK .GT. 40) .OR. (NM .GT. 20)) THEN 
         WRITE (6,*) 'T-MATRIX IS TOO BIG; RECOMPILE SOURCE CODE.' 
         STOP 1 
      ENDIF 
                                                                        
!     READ IN THE MATRIX                                                
      DO 800 MCTR = 1,NM 
         READ (2,*) MODE 
!        READ IN REAL VALUES FIRST:                                     
!V         READ (2,902) ((TREAL(I,J), J = 1,2*NRANK), I = 1,2*NRANK)    
!                                                                       
! ***  CAVEAT: THE T-MATRIX READ IN HERE IS THE TRANSPOSED MATRIX.      
! ***          THE T-MATRIX PROGRAM WRITES INTO OUPUT THE TRANSPOSED    
! ***          T-MATRIX. (i.e. f^ = a^T^, where ^=transpose)            
! ***          IN ROUTINE ADDPRC IN THIS PROGRAM, THE SCATTERED FIELD   
! ***          COEFFS. ARE CALCULATED AS f = Ta.         RAVI: 1/9/92   
         DO  20 I =1, 2*NRANK 
   20      READ(2,902) (TREAL(I,J), J=1, 2*NRANK) 
!        READ IN IMAGINARY VALUES NEXT:                                 
         DO  30  I =1, 2*NRANK 
   30    READ(2,902) (TIMAG(I,J), J=1, 2*NRANK) 
!        PLACE INTO COMPLEX T-MATRIX ARRAY                              
         DO 700 I = 1,2*NRANK 
            DO 650 J = 1,2*NRANK 
               TMAT(MCTR,I,J) = CMPLX(TREAL(I,J),TIMAG(I,J)) 
  650       CONTINUE 
  700    CONTINUE 
  800 END DO 
                                                                        
!     REWIND(2)                                                         
      DO 10 I = 1,NM 
         CMI(I) = FLOAT(I-1) 
   10 END DO 
      NR2 = 2*NRANK 
      NRANKI = NRANK+1 
!                                                                       
!                                                                       
! *** THIS IS TO PROCESS ONLY THE REQUIRED T-MATRICES ( DMIN TO DMAX )  
! ***                                                  RAVI: 12/30/91   
      IF ( DEQ .GE. DMIN .AND. DEQ .LE. DMAX ) THEN 
        PRINT*,'*** THE PARTICLE DIAMETER BEING PROCESSED IS ',DEQ 
        RETURN 
      ELSE 
        GO TO 100 
      ENDIF 
  101 RETURN 
!                                                                       
  900 FORMAT (A20) 
  901 FORMAT (1X,A20) 
  902 FORMAT(2X,8E15.7) 
      END                                           
!     END OF SUBROUTINE INIT                                            
                                                                        
!_________________________________________________________________      
                                                                        
      SUBROUTINE USERIN 
      PARAMETER (NMU=4) 
!     GET USER INPUT                                                    
! *** THIS ROUTINE HAS BEEN MODIFIED TO DO THE PROCESSING REQUIRED      
! *** AT THE PRESENT TIME                                RAVI: 12/20/91 
                                                                        
      INTEGER NTHETA, NPHI, DISTYP, COUNT, NANG, NUANG, IZ, IR 
      INTEGER CHECK 
      REAL DEGRAD, RADDEG, PI 
      REAL ELVANG, LOWER, UPPER, MEAN, THETAM, TWT, TASC, PWT, PASC 
      REAL SIGMA, SUM, THETA, P, SHOP, GP, LP, PNORM, SS, OLDTHT 
      REAL TCL, TCU, LAMBDA, EASTRT, EASTOP, EAINC, KAPPA, THETA0, PHI0 
      REAL DPAR1, DPAR2, DPAR3 
      CHARACTER*20 FNAME 
      CHARACTER*1 QDTYPE 
      CHARACTER*12 QUADNAM 
                                                                        
      COMMON DEGRAD, RADDEG, PI 
      COMMON /MAIN01/ ABV, ABH, CONK, S(91,4,4), SS(91,4,4), LAMBDA,    &
     &EXM(4,4),EXMP(4,4),EXMPT(4,4),SCHPT,SCVPT,EXHPT,EXVPT,PHNORM      
      COMMON /MAIN02/ NTHETA, TASC(20), TWT(20), NPHI, PASC(20), PWT(20) 
      COMMON /MAIN03/ DISTYP, THETAM, PNORM, MEAN, SIGMA, TCL, TCU 
      COMMON /MAIN04/ NANG, NUANG, OLDTHT, CHECK 
      COMMON /MAIN05/ FNAME 
      COMMON /MAIN06/ EASTRT, EASTOP, EAINC, ELVANG, KAPPA, THETA0, PHI0 
      COMMON/INT/QMU(NMU),QWT(NMU),QMUE(2*NMU),QWTE(2*NMU) 
      COMMON /DSD3/ DPAR1, DPAR2, DPAR3 
!                                                                       
      IF (CHECK .EQ. 1) THEN 
!        WRITE (9,905) EASTRT, EASTOP, EAINC                            
!        WRITE (9,908) NTHETA, NPHI                                     
      ENDIF 
!                                                                       
      IF (DISTYP .EQ. 0) THEN 
!        RANDOM DISTRIBUTION                                            
         LOWER = -1.0 
         UPPER = 1.0 
         PNORM = 1.0 
                                                                        
!        OUTPUT TO DIAGNOSTIC FILE                                      
         IF (CHECK .EQ. 1) THEN 
!           WRITE (9,909)                                               
!           WRITE (9,907) LOWER, UPPER                                  
         ENDIF 
      ELSE IF (DISTYP .EQ. 1) THEN 
!        SIMPLE HARMONIC OSCILLATOR:                                    
         CALL SHOTYP (MEAN,THETAM,PNORM,LOWER,UPPER,CHECK,TCL,TCU) 
      ELSE IF (DISTYP .EQ. 2) THEN 
!        GAUSSIAN:                                                      
         CALL GTYPE (MEAN,SIGMA,PNORM,LOWER,UPPER,CHECK,TCL,TCU) 
      ELSE IF (DISTYP .EQ. 3) THEN 
!        LANGEVIN:                                                      
         LOWER = -1.0 
         UPPER = 1.0 
         CALL LTYPE (KAPPA,THETA0,PHI0,PNORM,CHECK) 
      ELSE IF (DISTYP .EQ. 4) THEN 
!        FISHER:                                                        
         CALL FTYPE (KAPPA,PNORM,LOWER,UPPER,CHECK) 
      ENDIF 
                                                                        
!                                                                       
!                                                                       
         QDTYPE='L' 
!                                                                       
!                                                                       
!      SET THE NAME OF THE QUADRATURE ACCORDING TO SELECTION            
!                                                                       
        IF(QDTYPE.EQ.'G') THEN 
        QUADNAM='GAUSS-LEGEND' 
        ELSEIF(QDTYPE.EQ.'D') THEN 
        QUADNAM='DOUBLE-GAUSS' 
        ELSEIF(QDTYPE.EQ.'L') THEN 
        QUADNAM='LOBATTO-QUAD' 
        ELSE 
        QUADNAM='USER SPECIFY' 
        ENDIF 
!                                                                       
!                                                                       
!                                                                       
!                                                                       
                                                                        
!           MAKE THE DESIRED QUADRATURE ABSCISSAS AND WEIGHTS           
      IF (QDTYPE(1:1) .EQ. 'G') THEN 
          CALL GLQUAD                                                   &
     &                                                                  &
     &                                                                  &
     &                                                                  &
     &                                                                  &
     &                                                                  &
     &                       (NMU, QMU, QWT)                            
      ELSE IF (QDTYPE(1:1) .EQ. 'D') THEN 
          CALL DGQUAD                                                   &
     &                       (NMU, QMU, QWT)                            
        ELSE IF (QDTYPE(1:1) .EQ. 'L') THEN 
          CALL LBQUAD                                                   &
     &                       (NMU, QMU, QWT)                            
      ELSE 
!         WRITE (6,*) ' ENTER QUADRATURE MU VALUES:'                    
          DO 50 I = 1, NMU 
!             WRITE (*,'(1X,A,I2,A)') 'MU VALUE (', I, ') : '           
              READ (4,*) QMU(I) 
   50     CONTINUE 
          CALL QDWTS (NMU, QMU, QWT) 
      ENDIF 
!                                                                       
!                                                                       
!                                                                       
!                                                                       
      NANG = 2*NMU 
!     WRITE (6,901) NANG                                                
      NUANG = NANG 
      CALL GAUSS(TWT,TASC,NTHETA,LOWER,UPPER) 
      CALL GAUSS(PWT,PASC,NPHI,0.0,360.0) 
!     WRITE (6,*) '     '                                               
                                                                        
      RETURN 
  900 FORMAT (1X,F7.2) 
  901 FORMAT (1X,I2) 
  903 FORMAT (A20) 
  904 FORMAT (1X,A20) 
  905 FORMAT (1X,'ELEVATION ANGLE START, STOP, INCREMENT = ',           &
     &        2(F7.2,', '),F7.2,'.')                                    
  906 FORMAT (1X,'THETA DISTRIBUTION TYPE = ',I1,/,' (0 = RANDOM;',     &
     &      /,'  1 = SIMPLE HARMONIC OSCILLATOR;',/,'  2 = GAUSSIAN)')  
  907 FORMAT (1X,'COS(THETA) RANGES FROM ',F7.2,' TO ',F7.2) 
  908 FORMAT (1X,'INTEGRATION ORDERS:  THETA -- ',I2,'; PHI -- ',I2) 
  909 FORMAT (/,1X,'RANDOM DISTRIBUTION IN THETA:') 
      END                                           
!     END OF SUBROUTINE USERIN                                          
!_________________________________________________________________      
                                                                        
      SUBROUTINE SHOTYP (MEAN,THETAM,PNORM,LOWER,UPPER,CHECK,TCL,TCU) 
!     GET INPUT ON SIMPLE HARMONIC OSCILLATOR DISTRIBUTION              
                                                                        
      COMMON DEGRAD, RADDEG, PI 
      INTEGER CHECK, PCOUNT, TCOUNT 
      REAL DEGRAD, RADDEG, PI 
      REAL LOWER, UPPER, MEAN, THETAM, PNORM, PROB, SHOP 
      REAL TWT(20), TASC(20), PWT(20), PASC(20), THETA, PHI 
      REAL TCL, TCU, PSUM, TSUM 
!                                                                       
      REAL DPAR1, DPAR2, DPAR3 
      COMMON /DSD3/ DPAR1, DPAR2, DPAR3 
!                                                                       
!     GET MEAN, MAXIMUM AMPLITUDE OF THE OSCILLATION                    
!     WRITE (6,'(\A)') 'ENTER MEAN THETA OF OSCILLATION:'               
!     TYPE *, 'ENTER MEAN THETA OF OSCILLATION:'                        
! *** READ (4,*) MEAN                                                   
! *** THIS IS BEING READ IN THE MAIN PROGRAM         *** RAVI: 12/26/91 
      MEAN = DPAR1 
!     WRITE (6,900) MEAN                                                
!                                                                       
!     TYPE *, 'ENTER AMPLITUDE OF OSCILLATION (DEGREES): '              
! *** READ (4,*) THETAM                                                 
! *** THIS IS BEING READ IN THE MAIN PROGRAM         *** RAVI: 12/26/91 
      THETAM = DPAR2 
!     WRITE (6,900) THETAM                                              
                                                                        
!     WRITE(3,9020)MEAN,THETAM                                          
 9020 FORMAT('C   ',2(E15.7,1X),' SIMPLE HARMONIC ',                    &
     & 'MEAN = AMPLITUDE= ')                                            
!     NORMALIZE THE DISTRIBUTION OVER GIVEN THETA RANGE                 
      TCL = MEAN + THETAM 
      IF (TCL .GT. 180.0) THEN 
         TCL = 360.0 - TCL 
         LOWER = -1.00 
      ELSE 
         TCL = 180.0 
         LOWER = COS(DEGRAD*(MEAN + THETAM)) 
      ENDIF 
                                                                        
      TCU = MEAN - THETAM 
      IF (TCU .LT. 0.0) THEN 
         TCU = -TCU 
         UPPER = 1.00 
      ELSE 
         TCU = 0.00 
         UPPER = COS(DEGRAD*(MEAN - THETAM)) 
      ENDIF 
                                                                        
      CALL GAUSS (TWT, TASC, 20, LOWER, UPPER) 
      CALL GAUSS(PWT, PASC, 20, 0.0, 360.0) 
      PSUM = 0.0 
      DO 200 PCOUNT = 1,20 
         PHI = PASC(PCOUNT) 
         TSUM = 0.0 
         DO 100 TCOUNT = 1,20 
            THETA = 2.0*ATAN(SQRT(1.0-TASC(TCOUNT)**2)/                 &
     &              (1.0 + TASC(TCOUNT)))*RADDEG                        
                                                                        
            PROB = SHOP(THETA,MEAN,THETAM,PHI,1.00) 
!           ADJUSTMENTS IF THETA IS IN AN AMBIGUOUS ZONE                
            IF (THETA .LT. TCU) THEN 
               PROB = PROB +                                            &
     &                SHOP(-THETA,MEAN,THETAM,(PHI+180.0),1.00)         
            ELSE IF (THETA .GT. TCL) THEN 
               PROB = PROB +                                            &
     &                SHOP((360.0-THETA),MEAN,THETAM,(PHI+180.0),1.00)  
            ENDIF 
            TSUM = TSUM + PROB*TWT(TCOUNT) 
  100    CONTINUE 
         PSUM = PSUM + TSUM*PWT(PCOUNT) 
  200 END DO 
      PNORM = PSUM 
                                                                        
!     MAKE SURE NORMALIZATION HAS BEEN DONE PROPERLY                    
      PSUM = 0.0 
      DO 400 PCOUNT = 1,20 
         PHI = PASC(PCOUNT) 
         TSUM = 0.0 
         DO 300 TCOUNT = 1,20 
            THETA = 2.0*ATAN(SQRT(1.0-TASC(TCOUNT)**2)/                 &
     &              (1.0 + TASC(TCOUNT)))*RADDEG                        
                                                                        
            PROB = SHOP(THETA,MEAN,THETAM,PHI,PNORM) 
!           ADJUSTMENTS IF THETA IS IN AN AMBIGUOUS ZONE                
            IF (THETA .LT. TCU) THEN 
               PROB = PROB +                                            &
     &                SHOP(-THETA,MEAN,THETAM,(PHI+180.0),PNORM)        
            ELSE IF (THETA .GT. TCL) THEN 
               PROB = PROB +                                            &
     &                SHOP((360.0-THETA),MEAN,THETAM,(PHI+180.0),PNORM) 
            ENDIF 
            TSUM = TSUM + PROB*TWT(TCOUNT) 
  300    CONTINUE 
         PSUM = PSUM + TSUM*PWT(PCOUNT) 
  400 END DO 
                                                                        
!     OUTPUT TO DIAGNOSTIC FILE                                         
      IF (CHECK .EQ. 1) THEN 
!        WRITE (9,901)                                                  
!        WRITE (9,902) MEAN, THETAM                                     
!        WRITE (9,903) PSUM                                             
      ENDIF 
                                                                        
      RETURN 
  900 FORMAT (1X,F7.2) 
  901 FORMAT (/,1X,'SIMPLE HARMONIC OSCILLATOR DISTRIBUTION IN THETA:') 
  902 FORMAT (1X,'MEAN THETA = ',F7.2,'; AMPLITUDE OF OSCILLATION = ',  &
     &        F7.2)                                                     
  903 FORMAT (1X,'NORMALIZATION OF P(THETA,PHI) YIELDS ',F6.4) 
      END                                           
!     END OF SUBROUTINE SHOTYP()                                        
!_________________________________________________________________      
                                                                        
      SUBROUTINE GTYPE (MEAN,SIGMA,PNORM,LOWER,UPPER,CHECK,TCL,TCU) 
!     GET USER INPUT ON GAUSSIAN DISTRIBUTION FUNCTION                  
                                                                        
      COMMON DEGRAD, RADDEG, PI 
      INTEGER CHECK, PCOUNT, TCOUNT 
      REAL DEGRAD, RADDEG, PI 
      REAL LOWER, UPPER, MEAN, THETAM, PNORM, PROB, GP 
      REAL TWT(20), TASC(20), PWT(20), PASC(20), THETA, PHI 
      REAL TCL, TCU, PSUM, TSUM 
!                                                                       
      REAL DPAR1, DPAR2, DPAR3 
      COMMON /DSD3/ DPAR1, DPAR2, DPAR3 
!                                                                       
!     GET MEAN, SIGMA OF THE OSCILLATION                                
!     WRITE (6,*) 'ENTER MEAN THETA OF OSCILLATION:'                    
! *** READ (4,*) MEAN                                                   
! *** THIS IS BEING READ IN THE MAIN PROGRAM         *** RAVI:12/26/91  
      MEAN = DPAR1 
!     WRITE (6,900) MEAN                                                
!     PRINT*, 'MEAN = ',MEAN                                            
!     TYPE *, 'ENTER SIGMA OF GAUSSIAN (DEGREES): '                     
! *** READ (4,*) SIGMA                                                  
! *** THIS IS BEING READ IN THE MAIN PROGRAM         *** RAVI:12/26/91  
      SIGMA = DPAR2 
!     PRINT*, 'SIGMA = ',SIGMA                                          
!     WRITE (6,900) SIGMA                                               
                                                                        
!     WRITE(3,9030)MEAN,SIGMA                                           
 9030 FORMAT('C   ',2(E15.7,1X),' GAUSSIAN  ', 'MEAN = SIGMA= ') 
                                                                        
!     NORMALIZE THE DISTRIBUTION OVER GIVEN THETA RANGE                 
      TCL = MEAN + 2.0*SIGMA 
      IF (TCL .GT. 180.0) THEN 
         TCL = 360.0 - TCL 
         LOWER = -1.00 
      ELSE 
         TCL = 180.0 
         LOWER = COS(DEGRAD*(MEAN + 2.0*SIGMA)) 
      ENDIF 
                                                                        
      TCU = MEAN - 2.0*SIGMA 
      IF (TCU .LT. 0.0) THEN 
         TCU = -TCU 
         UPPER = 1.00 
      ELSE 
         TCU = 0.00 
         UPPER = COS(DEGRAD*(MEAN - 2.0*SIGMA)) 
      ENDIF 
                                                                        
      CALL GAUSS (TWT, TASC, 20, LOWER, UPPER) 
      CALL GAUSS(PWT, PASC, 20, 0.0, 360.0) 
!     NORMALIZE THE DISTRIBUTION OVER GIVEN THETA RANGE                 
      PSUM = 0.0 
      DO 200 PCOUNT = 1,20 
         PHI = PASC(PCOUNT) 
         TSUM = 0.0 
         DO 100 TCOUNT = 1,20 
            THETA = 2.0*ATAN(SQRT(1.0-TASC(TCOUNT)**2)/                 &
     &              (1.0 + TASC(TCOUNT)))*RADDEG                        
                                                                        
            PROB = GP(THETA,MEAN,SIGMA,PHI,1.00) 
!           ADJUSTMENTS IF THETA IS IN AN AMBIGUOUS ZONE                
            IF (THETA .LT. TCU) THEN 
               PROB = PROB +                                            &
     &                GP(-THETA,MEAN,SIGMA,(PHI+180.0),1.00)            
            ELSE IF (THETA .GT. TCL) THEN 
               PROB = PROB +                                            &
     &                GP((360.0-THETA),MEAN,SIGMA,(PHI+180.0),1.00)     
            ENDIF 
            TSUM = TSUM + PROB*TWT(TCOUNT) 
  100    CONTINUE 
         PSUM = PSUM + TSUM*PWT(PCOUNT) 
  200 END DO 
      PNORM = PSUM 
                                                                        
!     MAKE SURE NORMALIZATION HAS BEEN DONE PROPERLY                    
      PSUM = 0.0 
      DO 400 PCOUNT = 1,20 
         PHI = PASC(PCOUNT) 
         TSUM = 0.0 
         DO 300 TCOUNT = 1,20 
            THETA = 2.0*ATAN(SQRT(1.0-TASC(TCOUNT)**2)/                 &
     &              (1.0 + TASC(TCOUNT)))*RADDEG                        
                                                                        
            PROB = GP(THETA,MEAN,SIGMA,PHI,PNORM) 
!           ADJUSTMENTS IF THETA IS IN AN AMBIGUOUS ZONE                
            IF (THETA .LT. TCU) THEN 
               PROB = PROB +                                            &
     &                GP(-THETA,MEAN,SIGMA,(PHI+180.0),PNORM)           
            ELSE IF (THETA .GT. TCL) THEN 
               PROB = PROB +                                            &
     &                GP((360.0-THETA),MEAN,SIGMA,(PHI+180.0),PNORM)    
            ENDIF 
            TSUM = TSUM + PROB*TWT(TCOUNT) 
  300    CONTINUE 
         PSUM = PSUM + TSUM*PWT(PCOUNT) 
  400 END DO 
                                                                        
!     OUTPUT TO DIAGNOSTIC FILE                                         
      IF (CHECK .EQ. 1) THEN 
!        WRITE (9,901)                                                  
!        WRITE (9,902) MEAN, SIGMA                                      
!        WRITE (9,903) PNORM                                            
!        WRITE (9,904) PSUM                                             
      ENDIF 
                                                                        
      RETURN 
  900 FORMAT (1X,F7.2) 
  901 FORMAT (/,1X,'GAUSSIAN DISTRIBUTION IN THETA:') 
  902 FORMAT (1X,'MEAN THETA = ',F7.2,'; SIGMA = ',F7.2) 
  903 FORMAT (1X,'AREA UNDER SURFACE = ',F6.3) 
  904 FORMAT (1X,'NORMALIZATION OF P(THETA,PHI) YIELDS ',F6.4) 
      END                                           
!     END OF SUBROUTINE GTYPE()                                         
!_________________________________________________________________      
                                                                        
      SUBROUTINE LTYPE (KAPPA,THETA0,PHI0,PNORM,CHECK) 
!     GET USER INPUT ON LANGEVIN DISTRIBUTION FUNCTION                  
                                                                        
      COMMON DEGRAD, RADDEG, PI 
      INTEGER CHECK, PCOUNT, TCOUNT 
      REAL DEGRAD, RADDEG, PI 
      REAL LOWER, UPPER, KAPPA, THETA0, PHI0, PNORM, PROB, LP 
      REAL TWT(20), TASC(20), PWT(20), PASC(20), THETA, PHI 
      REAL PSUM, TSUM 
!                                                                       
      REAL DPAR1, DPAR2, DPAR3 
      COMMON /DSD3/ DPAR1, DPAR2, DPAR3 
!                                                                       
!     GET KAPPA, THETA0, AND PHI0 OF THE OSCILLATION                    
!     WRITE (6,*) 'ENTER KAPPA:'                                        
! *** READ (4,*) KAPPA                                                  
! *** THIS IS BEING READ IN THE MAIN PROGRAM         *** RAVI: 12/26/91 
      KAPPA = DPAR1 
!     WRITE (6,*) 'ENTER THETA0 OF DISTRIBUTION:'                       
! *** READ (4,*) THETA0                                                 
! *** THIS IS BEING READ IN THE MAIN PROGRAM         *** RAVI: 12/26/91 
      THETA0 = DPAR2 
!     WRITE (6,*) 'ENTER PHI0 OF DISTRIBUTION:'                         
! *** READ (4,*) PHI0                                                   
! *** THIS IS BEING READ IN THE MAIN PROGRAM         *** RAVI: 12/26/91 
      PHI0 = DPAR3 
                                                                        
!     WRITE(3,9040)KAPPA,THETA0,PHI0                                    
 9040 FORMAT('C   ',3(E15.7,1X),' LANGEVIN ','KAPPA = THETA0 =  PHI0=') 
                                                                        
!     NORMALIZE THE DISTRIBUTION OVER GIVEN THETA RANGE                 
      LOWER = -1.0 
      UPPER = 1.0 
                                                                        
      CALL GAUSS (TWT, TASC, 20, LOWER, UPPER) 
      CALL GAUSS(PWT, PASC, 20, 0.0, 360.0) 
!     NORMALIZE THE DISTRIBUTION OVER GIVEN THETA RANGE                 
      PSUM = 0.0 
      DO 200 PCOUNT = 1,20 
         PHI = PASC(PCOUNT) 
         TSUM = 0.0 
         DO 100 TCOUNT = 1,20 
            THETA = 2.0*ATAN(SQRT(1.0-TASC(TCOUNT)**2)/                 &
     &              (1.0 + TASC(TCOUNT)))*RADDEG                        
                                                                        
            PROB = LP(KAPPA, THETA, PHI, THETA0, PHI0, 1.00) 
            TSUM = TSUM + PROB*TWT(TCOUNT) 
  100    CONTINUE 
         PSUM = PSUM + TSUM*PWT(PCOUNT) 
  200 END DO 
      PNORM = PSUM 
                                                                        
!     MAKE SURE NORMALIZATION HAS BEEN DONE PROPERLY                    
      PSUM = 0.0 
      DO 400 PCOUNT = 1,20 
         PHI = PASC(PCOUNT) 
         TSUM = 0.0 
         DO 300 TCOUNT = 1,20 
            THETA = 2.0*ATAN(SQRT(1.0-TASC(TCOUNT)**2)/                 &
     &              (1.0 + TASC(TCOUNT)))*RADDEG                        
                                                                        
            PROB = LP(KAPPA, THETA, PHI, THETA0, PHI0, PNORM) 
            TSUM = TSUM + PROB*TWT(TCOUNT) 
  300    CONTINUE 
         PSUM = PSUM + TSUM*PWT(PCOUNT) 
  400 END DO 
                                                                        
!     OUTPUT TO DIAGNOSTIC FILE                                         
      IF (CHECK .EQ. 1) THEN 
!        WRITE (9,901)                                                  
!        WRITE (9,902) KAPPA, THETA0, PHI0                              
!        WRITE (9,903) PNORM                                            
!        WRITE (9,904) PSUM                                             
      ENDIF 
                                                                        
      RETURN 
  900 FORMAT (1X,F7.2) 
  901 FORMAT (/,1X,'LANGEVIN UNIMODAL DISTRIBUTION:') 
  902 FORMAT (1X,'KAPPA = ',F4.1,'; THETA0 = ',F7.2,                    &
     &        '; PHI0 = ',F7.2)                                         
  903 FORMAT (1X,'AREA UNDER SURFACE = ',E12.5) 
  904 FORMAT (1X,'NORMALIZATION OF P(THETA,PHI) YIELDS ',F6.4) 
      END                                           
!     END OF SUBROUTINE LTYPE()                                         
!_________________________________________________________________      
                                                                        
      SUBROUTINE FTYPE (KAPPA,PNORM,LOWER,UPPER,CHECK) 
!     GET USER INPUT ON FISHER DISTRIBUTION FUNCTION                    
                                                                        
      COMMON DEGRAD, RADDEG, PI 
      INTEGER CHECK 
      REAL DEGRAD, RADDEG, PI 
      REAL KAPPA, PNORM, FP, UPPER, LOWER, THETA0, THETA1, THETA 
      REAL FMAX, TARGET, F 
!                                                                       
      REAL DPAR1, DPAR2, DPAR3 
      COMMON /DSD3/ DPAR1, DPAR2, DPAR3 
!                                                                       
!     GET KAPPA OF THE DISTRIBUTION                                     
!     WRITE (6,*) 'ENTER KAPPA:'                                        
! *** READ (4,*) KAPPA                                                  
! *** THIS IS BEING READ IN THE MAIN PROGRAM     *** RAVI: 12/26/91     
      KAPPA = DPAR1 
!     WRITE(3,9051)KAPPA                                                
 9051 FORMAT('C   ',E15.7,1X,' FISHER ','KAPPA = ') 
                                                                        
!     GET LIMITS OF DISTRIBUTION FOR INTEGRATION PURPOSES               
      UPPER = 1.00 
      IF (KAPPA .GT. 0.0) THEN 
         THETA0 = (180.0/PI)*ACOS((-1.0 + SQRT(1.0 +                    &
     &            4.0*KAPPA*KAPPA))/(2.0*KAPPA))                        
         THETA1 = 180.0 
         FMAX = FP(KAPPA,THETA0,0.0,1.0) 
         TARGET = 0.01*FMAX 
!        REPEAT                                                         
  100       THETA = (THETA0 + THETA1)/2.0 
            F = FP(KAPPA,THETA,0.0,1.0) 
            IF (F .GT. TARGET) THEN 
               THETA0 = THETA 
            ELSE 
               THETA1 = THETA 
            ENDIF 
            IF (ABS((F - TARGET)/TARGET) .GT. 0.01) GO TO 100 
!        UNTIL (F - TARGET) < 0.01                                      
         LOWER = COS(DEGRAD*THETA) 
      ELSE 
         THETA = 180.0 
         LOWER = -1.00 
      ENDIF 
                                                                        
!     NORMALIZE THE DISTRIBUTION OVER GIVEN THETA RANGE                 
      IF (KAPPA .EQ. 0.00) THEN 
!        LIMITING CASE FOR UNIFORM DISTRIBUTION:                        
         PNORM = 0.5 
      ELSE 
!        GENERAL CASE:                                                  
         PNORM = KAPPA/(EXP(KAPPA) - EXP(-KAPPA)) 
      ENDIF 
                                                                        
!     OUTPUT TO DIAGNOSTIC FILE                                         
      IF (CHECK .EQ. 1) THEN 
!        WRITE (9,901)                                                  
!        WRITE (9,902) KAPPA                                            
!        WRITE (9,903) THETA                                            
      ENDIF 
                                                                        
      RETURN 
  900 FORMAT (1X,F7.2) 
  901 FORMAT (/,1X,'FISHER DISTRIBUTION:') 
  902 FORMAT (1X,'KAPPA = ',F4.1) 
  903 FORMAT (1X,'THETA RANGES FROM 0.0 TO ',F6.2,' DEGREES') 
      END                                           
!     END OF SUBROUTINE FTYPE()                                         
!_________________________________________________________________      
                                                                        
      SUBROUTINE SCFLD (TCOUNT, PCOUNT) 
!     CALCULATE SCATTERED FIELD FOR CURRENT PARTICLE ORIENTATION        
                                                                        
      INTEGER I, IZ, IP, NUANG, IM, NMODES, KMV, IFCT, CHECK 
      INTEGER TCOUNT, PCOUNT, NRANK 
      REAL ELVANG, ANGINT, ANGINP, SCH, SCV, UTHETA, UPHI, TANG, PANG 
      REAL COSB, SINB, CMI, CMV, CM2, PRODM, QUANM, TWM, CMXNRM 
      REAL OLDTHT, MEAN, SIGMA, TEMPTH, TEMPPH 
      REAL CONK, LAMBDA, DEQ, AOVERB, EPR, EPI, TEMP 
      REAL TCL, TCU, MODE, WVNMBR, EASTRT, EASTOP, EAINC 
      REAL DEGRAD, RADDEG, PI,KAPPA,THETAO,PHIO 
      COMPLEX ACANS, TMAT 
                                                                        
      COMMON DEGRAD, RADDEG, PI 
      COMMON /MTXCOM/ NRANK,NRANKI,TMAT(20,80,80),CMXNRM(80) 
      COMMON /FNCCOM/ PNMLLG(81) 
      COMMON /CMVCOM/ NMODES,KMV,CMI(20),CMV,CM2,TWM,PRODM 
      COMMON /UVCCOM/ ANGINT,ANGINP,UTHETA,UPHI,RTSFCT,IP 
      COMMON /ANGCOM/ TANG(91),PANG(91),COSB(91),SINB(91),              &
     &                ACANS(91,2,2),SCV,SCH,EXV,EXH                     
      COMMON /MAIN01/ ABV, ABH, CONK, S(91,4,4), SS(91,4,4), LAMBDA,    &
     &EXM(4,4),EXMP(4,4),EXMPT(4,4),SCHPT,SCVPT,EXHPT,EXVPT,PHNORM      
      COMMON /MAIN02/ NTHETA, TASC(20), TWT(20), NPHI, PASC(20), PWT(20) 
      COMMON /MAIN03/ DISTYP, THETAM, PNORM, MEAN, SIGMA, TCL, TCU 
      COMMON /MAIN04/ NANG, NUANG, OLDTHT, CHECK 
!J    COMMON /MAIN06/ EASTRT, EASTOP, EAINC, ELVANG                     
      COMMON /MAIN06/ EASTRT, EASTOP, EAINC, ELVANG, KAPPA, THETA0, PHI0 
                                                                        
!     ADJUST ORIENTATION VECTOR FOR +Z-INCIDENT WAVE                    
      TEMPTH = ANGINT 
      TEMPPH = ANGINP 
      CALL ROTATE(ELVANG,ANGINT,ANGINP) 
      IF (CHECK .EQ. 1) THEN 
!        WRITE (9,900) TCOUNT, PCOUNT, TEMPTH, TEMPPH, ANGINT, ANGINP   
      ENDIF 
                                                                        
!     ZERO SCATTERED-FIELD ARRAY                                        
      DO 200 I=1,NUANG 
         DO 100 IZ = 1,2 
            ACANS(I,1,IZ) = 0.0 
            ACANS(I,2,IZ) = 0.0 
  100    CONTINUE 
  200 END DO 
      SCH = 0.0 
                                                                        
!     SET POLARIZATION FOR FIRST LOOP THROUGH                           
      UTHETA = 0.0 
      UPHI = 1.0 
      IP = 1 
                                                                        
!     GENERATE TRANSFORMATION ANGLES BETWEEN LAB AND BODY FRAMES        
      CALL GENANG(NUANG,TANG,PANG,COSB,SINB) 
                                                                        
!     REPEAT                                                            
                                                                        
!        ONE CALL TO ADDPRC FOR EACH MODE:                              
  300    DO 500 IM = 1,NMODES 
            CMV = CMI(IM) 
            KMV = CMV 
            CM2 = CMV**2 
            PRODM = 1.0 
            IF (KMV .NE. 0) THEN 
               QUANM = CMV 
               DO 400 IFCT = 1,KMV 
                  QUANM = QUANM + 1.0 
                  PRODM = QUANM*PRODM/2.0 
  400          CONTINUE 
            ENDIF 
            TWM = 2.0*CMV 
                                                                        
            CALL ADDPRC 
  500    CONTINUE 
                                                                        
         IF (IP .EQ. 1) THEN 
!           SET UP FOR NEXT RUN, ORTHOGONAL POLARIZATION                
            SCV = 0.0 
                                                                        
            UTHETA = 1.0 
            UPHI = 0.0 
            IP = 2 
            GO TO 300 
         ENDIF 
!     UNTIL IP = 2                                                      
                                                                        
!     CALCULATE WAVENUMBER AND NORMALIZING FACTOR                       
      WVNMBR = 2.0*PI/LAMBDA 
      DO 700 I=1,NUANG 
         DO 600 IZ = 1,2 
            ACANS(I,1,IZ) = (4.0/WVNMBR)*ACANS(I,1,IZ) 
            ACANS(I,2,IZ) = (4.0/WVNMBR)*ACANS(I,2,IZ) 
  600    CONTINUE 
  700 END DO 
                                                                        
!JH   PRINT SCATTERING 2X2 MATRICIES HERE                               
                                                                        
      RETURN 
  900 FORMAT (2X,I2,7X,I2,6X,F6.2,4X,F6.2,3X,F6.2,3X,F6.2) 
      END                                           
!     END OF SUBROUTINE SCFLD                                           
!_________________________________________________________________      
                                                                        
      SUBROUTINE MMCALC(IPHI,SCHP,SCVP,EXHP, EXVP) 
!     ACCUMULATE CURRENT CONTRIBUTION TO MUELLER-MATRIX ACCUMULATOR     
                                                                        
      INTEGER J, NANG, IPHI, CHECK, ROW, COL 
      REAL T1, T2, T3, T4, S, PWT,  MMATRX(4,4) 
      REAL SCH, SCV, EXH, EXV, SCHP, SCVP, EXHP, EXVP 
      REAL OLDTHT, MEAN, SIGMA, TCL, TCU, LAMBDA 
      REAL EASTRT, EASTOP, EAINC, KAPPA, THETA0, PHI0 
      COMPLEX A1, A2, A3, A4, ACANS, F(2,2) 
                                                                        
      COMMON /ANGCOM/ TANG(91),PANG(91),COSB(91),SINB(91),              &
     &                ACANS(91,2,2),SCV,SCH,EXV,EXH                     
      COMMON /MAIN01/ ABV, ABH, CONK, S(91,4,4), SS(91,4,4), LAMBDA,    &
     &EXM(4,4),EXMP(4,4),EXMPT(4,4),SCHPT,SCVPT,EXHPT,EXVPT,PHNORM      
      COMMON /MAIN02/ NTHETA, TASC(20), TWT(20), NPHI, PASC(20), PWT(20) 
      COMMON /MAIN03/ DISTYP, THETAM, PNORM, MEAN, SIGMA, TCL, TCU 
      COMMON /MAIN04/ NANG, NUANG, OLDTHT, CHECK 
!     COMMON /MAIN06/ EASTRT, EASTOP, EAINC, ELVANG                     
      COMMON /MAIN06/ EASTRT, EASTOP, EAINC, ELVANG, KAPPA, THETA0, PHI0 
                                                                        
      DO 300 J = 1,NANG 
!        PREPARE TO PASS SCATTERING-AMPLITUDE MATRIX TO "MMAT"          
         F(1,1) = ACANS(J,1,1) 
         F(1,2) = ACANS(J,1,2) 
         F(2,1) = ACANS(J,2,1) 
         F(2,2) = ACANS(J,2,2) 
         IF(J.EQ.5) THEN 
!        WRITE(6,*)'SCATTERING M.'                                      
!        WRITE(6,*) -F(2,1),-F(2,2)                                     
!        WRITE(6,*) -F(1,1),-F(1,2)                                     
         ENDIF 
                                                                        
!        CALL "MMAT"; GET MUELLER MATRIX BACK                           
         CALL MMAT(F,MMATRX,LAMBDA) 
                                                                        
!        WEIGHT SUMMATION TERMS APPROPRIATELY -- TOTAL OF 16            
!        (NOTE:  NORMALIZATION OF PHI INTEGRATION IS TAKEN CARE         
!        OF BY THE VARIABLE "PNORM" IN THE THETA INTEGRATION.)          
         DO 200 ROW = 1,4 
            DO 100 COL = 1,4 
               S(J,ROW,COL) = S(J,ROW,COL) +                            &
     &                        PWT(IPHI)*MMATRX(ROW,COL)                 
  100       CONTINUE 
  200    CONTINUE 
  300 END DO 
                                                                        
!     PREPARE TO PASS SCATTERING AMPLITUDES TO EXTINCTION MATRIX        
                                                                        
      F(1,1) = ACANS(1,1,2) 
      F(1,2) = ACANS(1,1,1) 
      F(2,1) = ACANS(1,2,2) 
      F(2,2) = ACANS(1,2,1) 
                                                                        
      CALL EMAT(F, EXM, LAMBDA) 
                                                                        
         DO 500 ROW = 1,4 
            DO 400 COL = 1,4 
               EXMP(ROW,COL) = EXMP(ROW,COL) +                          &
     &                        PWT(IPHI)*EXM(ROW,COL)                    
  400       CONTINUE 
  500    CONTINUE 
                                                                        
      SCHP = SCHP + PWT(IPHI)*SCH 
      EXHP = EXHP + PWT(IPHI)*EXH 
      SCVP = SCVP + PWT(IPHI)*SCV 
      EXVP = EXVP + PWT(IPHI)*EXV 
                                                                        
      RETURN 
      END                                           
!     END OF SUBROUTINE MMCALC                                          
!_________________________________________________________________      
                                                                        
      SUBROUTINE MMAT (F, M, LAMBDA) 
!     CALCULATE (I,Q,U,V) MUELLER MATRIX FROM SCATTERING-AMPLITUDE MATRI
!     THE SCATTERING-AMPLITUDE MATRIX IS DEFINED IN TERMS OF THE E-FIELD
!     COMPONENTS PARALLEL AND PERPENDICULAR TO THE PLANE OF SCATTERING. 
                                                                        
!     REFERENCE:                                                        
!        LIGHT SCATTERING BY SMALL PARTICLES                            
!        H. C. VAN DE HULST                                             
!        DOVER PUBLICATIONS, INC. (NEW YORK), 1981;                     
!        PP. 34, 43-44, 46.                                             
!        (NOTE VAN DE HULST'S NOTATION FOR THE M AND A MATRICES; THE    
!        COMPONENTS OF THE F MATRIX WILL BE COPIED INTO THE APPROPRIATE 
!        A-MATRIX COMPONENTS FIRST TO AVOID CONFUSION.)                 
                                                                        
!     W. M. ADAMS     CSU -- FT. COLLINS, CO     31 JULY 1988           
                                                                        
!     INPUT TO SUBROUTINE:                                              
!        F  --  COMPLEX 2-BY-2 SCATTERING-AMPLITUDE MATRIX              
!        LAMBDA  --  WAVELENGTH IN MM                                   
!     RETURNED BY SUBROUTINE:                                           
!        M  --  REAL 4-BY-4 MUELLER MATRIX                              
                                                                        
!     LOCAL VARIABLES (SEE REFERENCE, PAGE 44 FOR NOTATION              
!     USED IN LOCAL VARIABLES):                                         
!        A1   --  ELEMENT OF TRANSFORMATION MATRIX A = F(2,1)           
!        A2   --  ELEMENT OF TRANSFORMATION MATRIX A = F(1,2)           
!        A3   --  ELEMENT OF TRANSFORMATION MATRIX A = -F(1,1)          
!        A4   --  ELEMENT OF TRANSFORMATION MATRIX A = -F(2,2)          
!        M1   --  SQUARED MAGNITUDE OF A1                               
!        M2   --  SQUARED MAGNITUDE OF A2                               
!        M3   --  SQUARED MAGNITUDE OF A3                               
!        M4   --  SQUARED MAGNITUDE OF A4                               
!        S21  --  (A2*CONJ(A1) + A1*CONJ(A2))/2                         
!        S23  --  (A2*CONJ(A3) + A3*CONJ(A2))/2                         
!        S24  --  (A2*CONJ(A4) + A4*CONJ(A2))/2                         
!        S31  --  (A3*CONJ(A1) + A1*CONJ(A3))/2                         
!        S34  --  (A3*CONJ(A4) + A4*CONJ(A3))/2                         
!        S41  --  (A4*CONJ(A1) + A1*CONJ(A4))/2                         
!        D21  --  (I/2)*(A2*CONJ(A1) - A1*CONJ(A2))                     
!        D23  --  (I/2)*(A2*CONJ(A3) - A3*CONJ(A2))                     
!        D24  --  (I/2)*(A2*CONJ(A4) - A4*CONJ(A2))                     
!        D31  --  (I/2)*(A3*CONJ(A1) - A1*CONJ(A3))                     
!        D34  --  (I/2)*(A3*CONJ(A4) - A4*CONJ(A3))                     
!        D41  --  (I/2)*(A4*CONJ(A1) - A1*CONJ(A4))                     
                                                                        
!     VERY IMPORTANT NOTES:                                             
!     ->  VAN DE HULST USES L (PARALLEL) AND R (PERPENDICULAR) NOTATION:
!         VECTOR R CROSS VECTOR L = DIRECTION OF PROPAGATION.           
!     ->  THIS PROGRAM USES UNIT VECTORS A(THETA) AND A(PHI); A(THETA)  
!         CROSS A(PHI) = DIRECTION OF PROPAGATION.                      
!     ->  VECTOR R AND A(PHI) POINT IN OPPOSITE DIRECTIONS, NECESSITATIN
!         A SIGN CHANGE (WHEN CONVERTING TO VAN DE HULST'S NOTATION) EAC
!         TIME THE PERPENDICULAR POLARIZATION IS INVOLVED.  SEE COPY INT
!         VAN DE HULST'S "A" MATRIX, BELOW.                             
                                                                        
      COMPLEX F(2,2) 
      REAL M(4,4) 
      COMPLEX A1, A2, A3, A4 
      REAL M1, M2, M3, M4 
      REAL S21, S23, S24, S31, S34, S41 
      REAL D21, D23, D24, D31, D34, D41 
      INTEGER ROW, COL 
      REAL   LAMBDA, PI 
      DATA  PI/ 3.14159265/ 
                                                                        
!     SET TRANSFORMATION MATRIX VALUES                                  
                                                                        
!     F(1,1):  PARALLEL COMPONENT, PERPENDICULAR INCIDENT POLARIZATION  
!                (VAN DE HULST'S A3 -- ONE SIGN CHANGE)                 
      A3 = -F(1,1) 
                                                                        
!     F(1,2):  PARALLEL COMPONENT, PARALLEL INCIDENT POLARIZATION       
!                (VAN DE HULST'S A2 -- NO SIGN CHANGES)                 
      A2 = F(1,2) 
                                                                        
!     F(2,1):  PERPENDICULAR COMPONENT, PERPENDICULAR INCIDENT POLARIZAT
!                (VAN DE HULST'S A1 -- TWO SIGN CHANGES WHICH CANCEL)   
      A1 = F(2,1) 
                                                                        
!     F(2,2):  PERPENDICULAR COMPONENT, PARALLEL INCIDENT POLARIZATION  
!                (VAN DE HULST'S A4 -- ONE SIGN CHANGE)                 
      A4 = -F(2,2) 
                                                                        
!     CALCULATE CONSTANTS TO BE USED IN MUELLER MATRIX CALCULATION --   
!     -> M TERMS:                                                       
         M1 = CABS(A1)*CABS(A1) 
         M2 = CABS(A2)*CABS(A2) 
         M3 = CABS(A3)*CABS(A3) 
         M4 = CABS(A4)*CABS(A4) 
!     -> S TERMS:                                                       
         S21 = REAL(A2)*REAL(A1) + AIMAG(A2)*AIMAG(A1) 
         S23 = REAL(A2)*REAL(A3) + AIMAG(A2)*AIMAG(A3) 
         S24 = REAL(A2)*REAL(A4) + AIMAG(A2)*AIMAG(A4) 
         S31 = REAL(A3)*REAL(A1) + AIMAG(A3)*AIMAG(A1) 
         S34 = REAL(A3)*REAL(A4) + AIMAG(A3)*AIMAG(A4) 
         S41 = REAL(A4)*REAL(A1) + AIMAG(A4)*AIMAG(A1) 
!     -> D TERMS:                                                       
         D21 = REAL(A2)*AIMAG(A1) - AIMAG(A2)*REAL(A1) 
         D23 = REAL(A2)*AIMAG(A3) - AIMAG(A2)*REAL(A3) 
         D24 = REAL(A2)*AIMAG(A4) - AIMAG(A2)*REAL(A4) 
         D31 = REAL(A3)*AIMAG(A1) - AIMAG(A3)*REAL(A1) 
         D34 = REAL(A3)*AIMAG(A4) - AIMAG(A3)*REAL(A4) 
         D41 = REAL(A4)*AIMAG(A1) - AIMAG(A4)*REAL(A1) 
                                                                        
!     CALCULATE MUELLER MATRIX PARAMETERS --                            
!     -> ROW 1:                                                         
         M(1,1) = (M2 + M3 + M4 + M1)/2.0 
         M(1,2) = (M2 - M3 + M4 - M1)/2.0 
         M(1,3) = S23 + S41 
         M(1,4) = -D23 - D41 
!     -> ROW 2:                                                         
         M(2,1) = (M2 + M3 - M4 - M1)/2.0 
         M(2,2) = (M2 - M3 - M4 + M1)/2.0 
         M(2,3) = S23 - S41 
         M(2,4) = -D23 + D41 
!     -> ROW 3:                                                         
         M(3,1) = S24 + S31 
         M(3,2) = S24 - S31 
         M(3,3) = S21 + S34 
         M(3,4) = -D21 + D34 
!     -> ROW 4:                                                         
         M(4,1) = D24 + D31 
         M(4,2) = D24 - D31 
         M(4,3) = D21 + D34 
         M(4,4) = S21 - S34 
                                                                        
!     MULTIPLY MUELLER-MATRIX TERMS BY APPROPRIATE SCALING FACTOR       
      RETURN 
      END                                           
!     END OF SUBROUTINE MMAT                                            
!_________________________________________________________________      
                                                                        
      SUBROUTINE EMAT(F, EX, LAMDA) 
!                                                                       
!     CALCULATE EXTINCTION MATRIX FROM FORWARD AMPLITUDE MATRIX 'F'     
!                                                                       
!     RETURNS 4-BY-4 REAL EXTINCTION MATRIX IN THE ARRAY 'EX'           
!                                                                       
      COMPLEX F(2,2) 
      REAL EX(4,4), LAMDA, PI, FACT 
                                                                        
      PI = 3.14159265 
      EX(1,1) = LAMDA*AIMAG(F(1,1)+F(2,2)) 
      EX(2,2) = EX(1,1) 
      EX(3,3) = EX(1,1) 
      EX(4,4) = EX(1,1) 
                                                                        
      EX(1,2) = LAMDA*AIMAG(F(1,1)-F(2,2)) 
      EX(2,1) = EX(1,2) 
                                                                        
      EX(3,4) = LAMDA*REAL(F(2,2)-F(1,1)) 
      EX(4,3) = -EX(3,4) 
                                                                        
      EX(1,3) = LAMDA*AIMAG(F(1,2)+F(2,1)) 
      EX(3,1) = EX(1,3) 
                                                                        
      EX(1,4) = LAMDA*REAL(F(1,2)-F(2,1)) 
      EX(4,1) = EX(1,4) 
                                                                        
      EX(2,3) = LAMDA*AIMAG(F(1,2)-F(2,1)) 
      EX(3,2) = -EX(2,3) 
                                                                        
      EX(2,4) = LAMDA*REAL(F(1,2)+F(2,1)) 
      EX(4,2) = -EX(2,4) 
      RETURN 
      END                                           
                                                                        
                                                                        
!_________________________________________________________________      
      SUBROUTINE ROTATE (ELVANG, THETAP, PHIP) 
!     IN THIS PROGRAM, THE INCIDENT WAVE TRAVELS IN THE PLUS-Z DIRECTION
!     THIS ROUTINE, ASSUMING A MINUS-X-DIRECTED INCIDENT WAVE AND AN    
!     ELEVATION ANGLE, ROTATES THE PARTICLE ORIENTATION VECTOR SO THAT  
!     THE PLUS-Z-DIRECTED WAVE "SEES" THE SAME ORIENTATION AS THE USER- 
!     SPECIFIED WAVE.                                                   
                                                                        
      REAL DEGRAD, RADDEG, PI, ELVANG, THETAP, PHIP, MAGXY 
      REAL ELVRAD, THERAD, PHIRAD, X, Y, Z, XPRIME, YPRIME, ZPRIME 
                                                                        
      DATA DEGRAD, RADDEG, PI/0.0174532925,57.295779513,3.14159265/ 
                                                                        
!     CONVERT TO RADIANS                                                
      ELVRAD = DEGRAD*ELVANG 
      THERAD = DEGRAD*THETAP 
      PHIRAD = DEGRAD*PHIP 
                                                                        
!     CALCULATE CARTESIAN COMPONENTS OF UNROTATED ORIENTATION VECTOR    
      X = SIN(THERAD)*COS(PHIRAD) 
      Y = SIN(THERAD)*SIN(PHIRAD) 
      Z = COS(THERAD) 
                                                                        
!     MULTIPLY ORIENTATION VECTOR BY ROTATION MATRIX                    
      XPRIME = SIN(ELVRAD)*X + COS(ELVRAD)*Z 
      YPRIME = Y 
      ZPRIME = -COS(ELVRAD)*X + SIN(ELVRAD)*Z 
      IF (ABS(XPRIME) .LT. 1.0E-03) XPRIME = 0.0 
      IF (ABS(YPRIME) .LT. 1.0E-03) YPRIME = 0.0 
      IF (ABS(ZPRIME) .LT. 1.0E-03) ZPRIME = 0.0 
                                                                        
!     CONVERT ROTATED CARTESIAN COMPONENTS INTO ANGLES THETA AND PHI    
!     THETA:                                                            
      MAGXY = SQRT(XPRIME*XPRIME + YPRIME*YPRIME) 
      IF (ABS(ZPRIME) .GT. 0.0) THEN 
         THETAP = ATAN(MAGXY/ZPRIME) 
         IF (THETAP .LT. 0.0) THETAP = THETAP + PI 
         IF ((THETAP .EQ. 0.0) .AND. (ZPRIME .LT. 0.0)) THETAP = PI 
      ELSE IF (ABS(ZPRIME) .EQ. 0.0) THEN 
         IF (ABS(MAGXY) .GT. 0.0) THEN 
            THETAP = PI/2.0 
         ELSE 
            THETAP = 0.0 
         ENDIF 
      ENDIF 
                                                                        
!     PHI:                                                              
      IF (YPRIME .EQ. 0.0) THEN 
         IF (XPRIME .LT. 0.0) THEN 
            PHIP = PI 
         ELSE 
            PHIP = 0.0 
         ENDIF 
      ELSE 
         PHIP = 2.0*ATAN(YPRIME/(XPRIME + MAGXY)) 
         IF (PHIP .LT. 0.0) PHIP = PHIP + 2.0*PI 
      ENDIF 
                                                                        
      THETAP = RADDEG*THETAP 
      PHIP = RADDEG*PHIP 
                                                                        
      RETURN 
      END                                           
!     END OF SUBROUTINE ROTATE                                          
!_________________________________________________________________      
                                                                        
      SUBROUTINE GAUSS (WT,ASC,N,AA,BB) 
                                                                        
      DIMENSION WT(N),ASC(N) 
                                                                        
      DATA PI,CONST,TOL/3.14159265358979,.148678816357,1.0E-07/ 
      DATA C1,C2,C3,C4/.125,-.0807291666,.2460286458,-1.824438767/ 
                                                                        
      IF (N .EQ. 1) THEN 
         ASC(1) = 0.5773502692 
         WT(1) = 1.0 
      ELSE 
         DN = N 
         NDIV2 = N/2 
         NP1 = N + 1 
         NNP1 = N*NP1 
         APPFCT = 1./SQRT((N + 0.5)**2 + CONST) 
         CON1 = 0.5*(BB - AA) 
         CON2 = 0.5*(BB + AA) 
         DO 120 K = 1,NDIV2 
            B = (K - 0.25)*PI 
            BISQ = 1.0/(B*B) 
            BFROOT = B*(1.0 + BISQ*(C1 + BISQ*(C2 + BISQ*(C3            &
     &               + C4*BISQ))))                                      
            X = COS(APPFCT*BFROOT) 
!           REPEAT                                                      
  100          PM2 = 1.0 
               PM1 = X 
               DO 110 IN = 2,N 
                  P = ((2*IN - 1)*X*PM1 - (IN - 1)*PM2)/IN 
                  PM2 = PM1 
                  PM1 = P 
  110          CONTINUE 
               PM1 = PM2 
               AUX = 1.0/(1.0 - X*X) 
               DER1P = DN*(PM1 - X*P)*AUX 
               DER2P = (2.0*X*DER1P - NNP1*P)*AUX 
               RATIO = P/DER1P 
               XI = X - RATIO*(1.0 + RATIO*DER2P/(2.0*DER1P)) 
               IF ((ABS(XI - X) - TOL) .GT. 0.0) THEN 
                  X = XI 
                  GO TO 100 
               ENDIF 
!           UNTIL (ABS(XI - X) - TOL) <= 0.0                            
            ASC(K) = -X 
            WT(K) = 2.0*(1.0 - X*X)/(DN*PM1)**2 
            ASC(NP1-K) = -ASC(K) 
            WT(NP1-K) = WT(K) 
  120    CONTINUE 
         IF (MOD(N,2) .NE. 0) THEN 
            ASC(NDIV2 + 1) = 0.0 
            NM1 = N-1 
            NM2 = N-2 
            PROD = DN 
            DO 130 K = 1,NM2,2 
               PROD = FLOAT(NM1 - K)/FLOAT(N - K)*PROD 
  130       CONTINUE 
            WT(NDIV2 + 1) = 2.0/PROD**2 
         ENDIF 
         DO 140 K = 1,N 
            ASC(K) = CON1*ASC(K) + CON2 
            WT(K) = CON1*WT(K) 
  140    CONTINUE 
      ENDIF 
      RETURN 
      END                                           
                                                                        
!_________________________________________________________________      
      SUBROUTINE GENLGP 
      COMPLEX TMAT 
      COMMON /MTXCOM/ NRANK,NRANKI,TMAT(20,80,80),CMXNRM(80) 
      COMMON /FNCCOM/ PNMLLG(81) 
      COMMON /CMVCOM/ NM,KMV,CMI(20),CMV,CM2,TWM,PRODM 
      COMMON /THTCOM/ THETA,SINTH,COSTH 
      DTWM = TWM+1.0 
      IF(THETA.NE.0.0) GO TO 16 
      IF(KMV.EQ.1) GO TO 12 
      DO 8 ILG = 1,NRANKI 
      PNMLLG(ILG) = 0.0 
    8 END DO 
      GO TO 88 
   12 PNMLLG(1) = 0.0 
      PNMLLG(2) = 1.0 
      PLA = 1.0 
      GO TO 48 
   16 IF(KMV.GT.0) GO TO 40 
!     THE SPECIAL CASE WHEN M = 0.                                      
      PLA = 1.0/SINTH 
      PLB = COSTH*PLA 
      PNMLLG(1) = PLA 
      PNMLLG(2) = PLB 
      IBEG = 3 
      GO TO 60 
!     GENERAL CASE FOR M NOT EQUAL TO 0.                                
   40 DO 44 ILG = 1,KMV 
      PNMLLG(ILG) = 0.0 
   44 END DO 
      PLA = PRODM*SINTH**(KMV-1) 
      PNMLLG(KMV+1) = PLA 
   48 PLB = DTWM*COSTH*PLA 
      PNMLLG(KMV+2) = PLB 
      IBEG = KMV+3 
!     DO RECURSION FORMULA FOR ALL REMAINING LEGENDRE POLYNOMIALS.      
   60 CNMUL = IBEG+IBEG-3 
      CNM = 2.0 
      CNMM = DTWM 
      DO 80 ILGR = IBEG,NRANKI 
      PLC = (CNMUL*COSTH*PLB-CNMM*PLA)/CNM 
      PNMLLG(ILGR) = PLC 
      PLA = PLB 
      PLB = PLC 
      CNMUL = CNMUL+2.0 
      CNM = CNM+1.0 
      CNMM = CNMM+1.0 
   80 END DO 
   88 RETURN 
      END                                           
                                                                        
!_________________________________________________________________      
      SUBROUTINE GENANG (N,TSP,PSP,COSBTA,SINBTA) 
!     A ROUTINE TO TRANSFORM SCATTERING ANGLES FROM THE LAB FRAME TO THE
!     BODY FRAME  AND TO GENERATE THE TRANSFORMATION MATRIX FOR THE     
!     POLARIZATION VECTORS                                              
                                                                        
      PARAMETER(NMU=4) 
      INTEGER I, CHECK 
      REAL DEGRAD, RADDEG, PI 
      REAL CHI, ALPHA, SINCHI, COSCHI, SINALF, COSALF 
      REAL THETAS, SINTS, COSTS, COSPS, XPRIME, YPRIME, ZPRIME, MAGXY 
      REAL TSP, PSP, COSBTA, SINBTA, PHISP, SINPSP, COSPSP, OLDTHT 
      REAL EASTRT, EASTOP, EAINC, KAPPA, THETA0, PHI0 
                                                                        
      COMMON DEGRAD,RADDEG,PI 
      COMMON/CMVCOM/NM,KMV,CMI(20),CMV,CM2,TWM,PRODM 
      COMMON/UVCCOM/ANGINT,ANGINP,UTHETA,UPHI,RTSFCT,IP 
      COMMON /FNCCOM/ PNMLLG(81) 
      COMMON /MAIN04/ NANG, NUANG, OLDTHT, CHECK 
      COMMON /INT/QMU(NMU),QWT(NMU),QMUE(2*NMU),QWTE(2*NMU) 
      DIMENSION TSP(N),PSP(N),COSBTA(N),SINBTA(N) 
!J    COMMON /MAIN06/ EASTRT, EASTOP, EAINC, ELVANG                     
      COMMON /MAIN06/ EASTRT, EASTOP, EAINC, ELVANG, KAPPA, THETA0, PHI0 
                                                                        
!     CALCULATE CHI AND ALPHA (WANG'S THESIS, EQNS. 72 AND 73)          
      CHI = ANGINP - 180.0 
      ALPHA = -ANGINT 
                                                                        
!     CALCULATE SINES AND COSINES OF CHI AND ALPHA                      
      SINCHI = SIN(DEGRAD*CHI) 
      COSCHI = COS(DEGRAD*CHI) 
      SINALF = SIN(DEGRAD*ALPHA) 
      COSALF = COS(DEGRAD*ALPHA) 
                                                                        
!     CALCULATE THETA (SCATTERED) PRIME -- "TSP" -- AND PHI             
!     (SCATTERED) PRIME -- "PSP" -- FOR ALL LAB-FRAME THETA             
!     (SCATTERED) -- "THETAS" -- ANGLES                                 
      DO 100 I = 1,N 
!        THETA (SCATTERED) AND ITS SINE AND COSINE:                     
            IF(I.LE.(N/2)) THEN 
!U            QMUE(I)=QMU(I)                                            
             QMUE(I)=QMU(NMU+1-I) 
             ELSE 
!U             QMUE(I)=-QMU(I-NMU)                                      
            QMUE(I)=-QMU(N+1-I) 
             ENDIF 
         THETAS = ACOS(QMUE(I)) 
         SINTS = SIN(THETAS) 
         COSTS = COS(THETAS) 
                                                                        
!        PHI (SCATTERED) AND ITS COSINE (ACTUALLY ONLY NEED COSINE OF PH
!        (SCATTERED)).  THIS WAS ADDED TO IMPOSE THE REQUIREMENT THAT   
!        PHI (SCATTERED) BE EQUAL TO ZERO EXCEPT AT THETA (SCATTERED)   
!        EQUAL TO PI, WHERE PHI (SCATTERED) IS REQUIRED TO EQUAL PI.    
!        NOTE THAT THIS CHANGE MAKES NO DIFFERENCE IN THIS LOOP, SINCE  
!        THE DIFFERENT COS(PS) IS MULTIPLIED BY SIN(PI), WHICH IS ZERO  
!        ANYWAY.  REAL DIFFERENCE COMES IN THE SINE BETA/ COSINE BETA   
!        LOOP AT END OF THIS SUBROUTINE.                                
         IF (I .LT. N) THEN 
!           ALL OTHER ANGLES OTHER THAN BACKSCATTER:                    
            COSPS = 1.00 
         ELSE IF (I .EQ. N) THEN 
!           BACKSCATTER:  PHI (SCATTERED) = PI                          
            COSPS = -1.00 
         ENDIF 
                                                                        
!        X', Y' AND Z' (WANG, EQNS. 84 MODIFIED FOR PHI (SCATTERED))    
         XPRIME = SINTS*COSALF*COSCHI*COSPS - COSTS*SINALF 
         YPRIME = -SINTS*SINCHI*COSPS 
         ZPRIME = SINTS*SINALF*COSCHI*COSPS + COSTS*COSALF 
                                                                        
!        THETA (SCATTERED) PRIME (WANG, EQN. 85A)                       
         MAGXY = SQRT(XPRIME*XPRIME + YPRIME*YPRIME) 
         IF (ABS(ZPRIME) .GT. 0.0) THEN 
            TSP(I) = ATAN(MAGXY/ZPRIME) 
            IF (TSP(I) .LT. 0.0) TSP(I) = TSP(I) + PI 
            IF ((TSP(I) .EQ. 0.0) .AND. (ZPRIME .LT. 0.0)) TSP(I) = PI 
         ELSE IF (ABS(ZPRIME) .EQ. 0.0) THEN 
            IF (ABS(MAGXY) .GT. 0.0) THEN 
               TSP(I) = PI/2.0 
            ELSE 
               TSP(I) = 0.0 
            ENDIF 
         ENDIF 
                                                                        
!        PHI (SCATTERED) PRIME (WANG, EQN. 85B -- NOTE ERROR IN         
!        EQUATION:  SHOULD BE ATAN(Y'/X'))                              
         IF (ABS(XPRIME) .GT. 0.0) THEN 
            IF (ABS(XPRIME + MAGXY) .LT. 1.0E-03) THEN 
               PSP(I) = PI 
            ELSE 
               PSP(I) = 2.0*ATAN(YPRIME/(XPRIME + MAGXY)) 
            ENDIF 
            IF (PSP(I) .LT. 0.0) PSP(I) = PSP(I) + 2.0*PI 
         ELSE IF (ABS(XPRIME) .EQ. 0.0) THEN 
            IF (YPRIME .LT. 0.0) THEN 
               PSP(I) = 1.5*PI 
            ELSE IF (YPRIME .EQ. 0.0) THEN 
               PSP(I) = 0.0 
            ELSE IF (YPRIME .GT. 0.0) THEN 
               PSP(I) = 0.5*PI 
            ENDIF 
         ENDIF 
                                                                        
         TSP(I) = RADDEG*TSP(I) 
         PSP(I) = RADDEG*PSP(I) 
  100 END DO 
                                                                        
!     CALCULATE THE ELEMENTS OF THE TRANSFORMATION                      
!     MATRIX BETWEEN LAB AND BODY FRAMES FOR THE                        
!     DIRECTION OF POLARIZATION                                         
      DO 110 I = 1,N 
!        PHI (SCATTERED) PRIME AND SINE AND COSINE                      
         PHISP = DEGRAD*PSP(I) 
         SINPSP = SIN(PHISP) 
         COSPSP = COS(PHISP) 
                                                                        
!        THETA (SCATTERED) AND SINE AND COSINE; COSINE OF PHI (SCATTERED
         THETAS = ACOS(QMUE(I)) 
         SINTS = SIN(THETAS) 
         COSTS = COS(THETAS) 
!        SEE COMMENTS IN "DO 100" LOOP ABOUT PHI (SCATTERED)            
         IF (I .LT. N) THEN 
!           ALL OTHER ANGLES OTHER THAN BACKSCATTER:                    
            COSPS = 1.00 
         ELSE IF (I .EQ. N) THEN 
!           BACKSCATTER:  PHI (SCATTERED) = PI                          
            COSPS = -1.00 
         ENDIF 
                                                                        
!        COSINE BETA AND SINE BETA (WANG, EQNS. 93 -- NOTE COUPLE OF    
!        ERRORS IN 93B)                                                 
         COSBTA(I) = -SINPSP*COSALF*SINCHI*COSPS + COSPSP*COSCHI*COSPS 
         SINBTA(I) = SINPSP*COSALF*COSCHI*COSTS*COSPS                   &
     &               + SINALF*SINTS*SINPSP + SINCHI*COSTS*COSPSP*COSPS  
  110 END DO 
                                                                        
      RETURN 
      END                                           
!     END OF SUBROUTINE GENANG                                          
!_________________________________________________________________      
                                                                        
                                                                        
      SUBROUTINE ADDPRC 
!     A ROUTINE TO OBTAIN THE SCATTERED FIELD COEFFICIENTS AND CALCULATE
!     THE DIFFERENTIAL SCATTERING CROSS SECTION.                        
                                                                        
      INTEGER CHECK 
      REAL ELVANG, OLDTHT, EASTRT, EASTOP, EAINC, KAPPA, THETA0, PHI0 
      COMPLEX TMAT,AB1(80),AB2(80),CI,C1,C2,S1,S2,CIM,ACANS 
      COMPLEX FG1(80),FG2(80) 
                                                                        
      COMMON DEGRAD,RADDEG,PI 
      COMMON /MTXCOM/ NRANK,NRANKI,TMAT(20,80,80),CMXNRM(80) 
      COMMON /FNCCOM/ PNMLLG(81) 
      COMMON /CMVCOM/ NM,KMV,CMI(20),CMV,CM2,TWM,PRODM 
      COMMON /THTCOM/ THETA,SINTH,COSTH 
      COMMON /UVCCOM/ ANGINT,ANGINP,UTHETA,UPHI,RTSFCT,IP 
      COMMON /ANGCOM/ TANG(91),PANG(91),COSB(91),SINB(91),              &
     &                ACANS(91,2,2),SCV,SCH,EXV,EXH                     
      COMMON /MAIN04/ NANG, NUANG, OLDTHT, CHECK 
!J    COMMON /MAIN06/ EASTRT, EASTOP, EAINC, ELVANG                     
      COMMON /MAIN06/ EASTRT, EASTOP, EAINC, ELVANG, KAPPA, THETA0, PHI0 
                                                                        
      NR2 = 2*NRANK 
      CI = (0.0,1.0) 
!     TRANSFORM THE PROBLEM FROM THE LAB FRAME TO THE BODY FRAME        
!     TANG(1) = INCIDENT ANGLE IN THE BODY SYSTEM                       
!     GENERATE THE LEGENDRE FUNCTIONS FOR THE INCIDENT ANGLE.           
      THETA=DEGRAD*TANG(1) 
      COSTH = COS(THETA) 
      SINTH = SIN(THETA) 
      CALL GENLGP 
                                                                        
!     TRANSFORM THE INCIDENT FIELD POLARIZATION VECTOR FROM THE LAB FRAM
!     THE BODY FRAME:                                                   
!     UU1 = COMPONENT IN THETA DIRECTION; UU2 = COMPONENT IN PHI DIRECTI
      UU1 = COSB(1)*UTHETA + SINB(1)*UPHI 
      UU2 = -SINB(1)*UTHETA + COSB(1)*UPHI 
                                                                        
!     GENERATE THE INCIDENT FIELD COEFFICIENTS -- AB1 = THETA POLARIZATI
!     AND AB2 = PHI POLARIZATION.                                       
      CN = 0.0 
      DO 35 N = 1,NRANK 
         NP = N + NRANK 
         CN = CN + 1.0 
         N1 = N + 1 
         C1 = CI**N 
         C2 = CI**N1 
         P1 = CN*COSTH*PNMLLG(N1) - (CN+CMV)*PNMLLG(N) 
         P2 = CMV*PNMLLG(N1) 
         AB1(N) = -C1*P2*UU1 
         AB1(NP) = C2*P1*UU1 
         AB2(N) = C1*P1*UU2 
         AB2(NP) = -C2*P2*UU2 
   35 END DO 
                                                                        
!     THE SCATTERED FIELD COEFFICIENTS = THE TRANSITION MATRIX TIMES THE
!     INCIDENT FIELD COEFFICIENTS. (f = Ta, *** refer routine INIT ***) 
!                                                        RAVI: 1/9/92   
      DO 45 I = 1,NR2 
         S1 = 0.0 
         S2 = 0.0 
         DO 40 J = 1,NR2 
            S1 = S1 + TMAT(KMV+1,J,I)*AB1(J) 
            S2 = S2 + TMAT(KMV+1,J,I)*AB2(J) 
   40    CONTINUE 
         FG1(I) = S1 
         FG2(I) = S2 
   45 END DO 
                                                                        
!     CALCULATE THE NORMALIZATION FACTOR                                
      IF (KMV .EQ. 0) THEN 
         EM = 1.0 
      ELSE 
         EM = 2.0 
      ENDIF 
      DO 450 IROW = 1,NRANK 
      CKROW = IROW 
      IF(KMV.GT.0) GO TO 426 
      FCTKI = 1.0 
      GO TO 440 
  426 IF(IROW.GE.KMV) GO TO 430 
      CMXNRM(IROW) = 1.0 
      GO TO 450 
  430 IBFCT = IROW-KMV+1 
      IEFCT = IROW+KMV 
      FPROD = IBFCT 
      FCTKI = 1.0 
      DO 432 LFCT = IBFCT,IEFCT 
      FCTKI = FCTKI*FPROD 
      FPROD = FPROD+1.0 
  432 END DO 
  440 CMXNRM(IROW) = 4.0*CKROW*(CKROW+1.0)*FCTKI/(EM*(2.0*CKROW+1.0)) 
!      WRITE(6,"(' CMXNRM', I5,E15.7)") IROW,CMXNRM(IROW)               
  450 END DO 
                                                                        
!     EVALUATE THE SCATTERED FIELD AT EACH SCATTERING ANGLE.            
      DO 170 IU = 1,NUANG 
!        GENERATE THE LEGENDRE MULTIPLIERS.                             
         THETA = DEGRAD*TANG(IU) 
         SINTH = SIN(THETA) 
         COSTH = COS(THETA) 
         IF (TANG(IU) .EQ. 180.0) THEN 
            THETA = 0.0 
            COSTH = -1.0 
         ENDIF 
         CALL GENLGP 
         PHI = CMV*PANG(IU)*DEGRAD 
         SINPHI = SIN(PHI) 
         COSPHI = COS(PHI) 
         S1 = 0.0 
         S2 = 0.0 
         CN = 0.0 
         DO 160 N = 1,NRANK 
            NP = N + NRANK 
            IF ((IU .EQ. 1) .AND. (IP .EQ. 1)) THEN 
               SCH = SCH + (CABS(FG1(N))**2 + CABS(FG1(NP))**2          &
     &              + CABS(FG2(N))**2 + CABS(FG2(NP))**2)/CMXNRM(N)     
            ENDIF 
            IF ((IU .EQ. 1) .AND. (IP .EQ. 2)) THEN 
               SCV = SCV + (CABS(FG1(N))**2 + CABS(FG1(NP))**2          &
     &              + CABS(FG2(N))**2 + CABS(FG2(NP))**2)/CMXNRM(N)     
            ENDIF 
            N1 = N + 1 
            CN = CN + 1.0 
            P1 = CN*COSTH*PNMLLG(N1) - (CN + CMV)*PNMLLG(N) 
            P2 = CMV*PNMLLG(N1) 
            AA = SINPHI*P1 
            BB = COSPHI*P1 
            CC = SINPHI*P2 
            DD = COSPHI*P2 
            CIM = (-CI)**N1 
                                                                        
!           SOLVE FOR THE THETA-POLARIZED SCATTERED FIELD               
!           IN THE BODY FRAME                                           
            S1 = S1 + CIM*(FG1(N)*DD + CI*FG1(NP)*BB - FG2(N)*CC        &
     &           - CI*FG2(NP)*AA)/CMXNRM(N)                             
                                                                        
!           SOLVE FOR THE PHI-POLARIZED SCATTERING FIELD                
!           IN THE BODY FRAME                                           
            S2 = S2 - CIM*(FG1(N)*AA + CI*FG1(NP)*CC + FG2(N)*BB        &
     &           + CI*FG2(NP)*DD)/CMXNRM(N)                             
  160    CONTINUE 
                                                                        
         ACANS(IU,1,IP) = ACANS(IU,1,IP) + S1 
         ACANS(IU,2,IP) = ACANS(IU,2,IP) + S2 
  170 END DO 
                                                                        
!     TRANSFORM THE SCATTERED FIELD (ACANS) BACK TO THE LAB FRAME       
      IF (KMV .EQ. (NM-1)) THEN 
         DO 180 KUP = 1,NUANG 
            S1 = ACANS(KUP,1,IP) 
            S2 = ACANS(KUP,2,IP) 
            ACANS(KUP,1,IP) = COSB(KUP)*S1 - SINB(KUP)*S2 
            ACANS(KUP,2,IP) = SINB(KUP)*S1 + COSB(KUP)*S2 
  180    CONTINUE 
         IF (IP .EQ. 1) EXH = AIMAG(ACANS(1,2,1)) 
         IF (IP .EQ. 2) EXV = AIMAG(ACANS(1,1,2)) 
      ENDIF 
      RETURN 
      END                                           
                                                                        
!_____________________________________________________________________  
      REAL FUNCTION SHOP (THETA, THBAR, THMAX, PHI, NORM) 
!     SIMPLE HARMONIC OSCILLATOR ORIENTATION PROBABILITY DENSITY FUNCTIO
!     THETA, THBAR, THETA-MAX (THMAX) AND PHI ARE ASSUMED GIVEN IN DEGRE
!     "NORM" IS CALCULATED IN SUBROUTINE "USERIN" AND                   
!     NORMALIZES THE FUNCTION TO UNITY OVER THE RANGE                   
!     OF THETA AND PHI SPECIFIED.                                       
!     NOTE THAT P(THETA,PHI) IS ASSUMED SEPARABLE IN THETA AND PHI.     
                                                                        
      REAL PI, DEGRAD 
      REAL THETA, THBAR, THMAX, CONST, RATIO, PHI, NORM 
      REAL PTHETA, PPHI 
      DATA PI, DEGRAD /3.141592654, 0.0174532925/ 
                                                                        
      CONST = 1.0/(2.0*PI*PI*DEGRAD*THMAX) 
      RATIO = (THETA - THBAR)/THMAX 
      PTHETA = CONST/SQRT(1.0 - RATIO*RATIO) 
      PPHI = 1.0/360.0 
      SHOP = PTHETA*PPHI/NORM 
                                                                        
      RETURN 
      END                                           
!     END OF FUNCTION SHOP()                                            
                                                                        
      REAL FUNCTION GP (THETA, THBAR, SIGMA, PHI, NORM) 
!     GAUSSIAN ORIENTATION PROBABILITY DENSITY FUNCTION                 
!     THETA, THBAR, SIGMA AND PHI ARE ASSUMED GIVEN IN DEGREES          
!     "NORM" IS CALCULATED IN SUBROUTINE "USERIN" AND                   
!     NORMALIZES THE FUNCTION TO UNITY OVER THE RANGE                   
!     OF THETA AND PHI SPECIFIED.                                       
!     NOTE THAT P(THETA,PHI) IS ASSUMED SEPARABLE IN THETA AND PHI.     
                                                                        
      REAL PI, DEGRAD 
      REAL THETA, THBAR, SIGMA, CONST, RATIO, PHI, NORM 
      REAL PTHETA, PPHI 
      DATA PI, DEGRAD /3.141592654, 0.0174532925/ 
                                                                        
      CONST = 1.0/(SQRT(2.0*PI)*DEGRAD*SIGMA) 
      RATIO = (THETA - THBAR)/(SQRT(2.0)*SIGMA) 
      PTHETA = CONST*EXP(-RATIO*RATIO) 
      PPHI = 1.0/360.0 
      GP = PTHETA*PPHI/NORM 
                                                                        
      RETURN 
      END                                           
!     END OF FUNCTION GP()                                              
                                                                        
!_________________________________________________________________      
      REAL FUNCTION LP (KAPPA, THETA, PHI, THETA0, PHI0, NORM) 
!     LANGEVIN ORIENTATION PROBABILITY DISTRIBUTION                     
!     THETA, PHI, THETA0, AND PHI0 ARE ASSUMED GIVEN IN DEGREES         
!     "NORM" IS CALCULATED IN SUBROUTINE "USERIN" AND                   
!     NORMALIZES THE FUNCTION TO UNITY.                                 
                                                                        
      REAL KAPPA, THETA, PHI, THETA0, PHI0, NORM 
      REAL DEGRAD, TERM1, TERM2 
      DATA DEGRAD /0.0174532925/ 
                                                                        
      TERM1 = COS(DEGRAD*THETA0)*COS(DEGRAD*THETA) 
      TERM2 = SIN(DEGRAD*THETA0)*SIN(DEGRAD*THETA)*                     &
     &        COS(DEGRAD*(PHI - PHI0))                                  
      LP = SIN(DEGRAD*THETA)*EXP(KAPPA*(TERM1 + TERM2))/(NORM*360.) 
                                                                        
      RETURN 
      END                                           
!     END OF FUNCTION LP()                                              
!_________________________________________________________________      
                                                                        
      REAL FUNCTION FP (KAPPA, THETA, PHI, NORM) 
!     LANGEVIN ORIENTATION PROBABILITY DISTRIBUTION                     
!     (ALSO FISHER DISTRIBUTION):  MEAN THETA = ZERO.                   
!     THETA AND PHI ARE ASSUMED GIVEN IN DEGREES.                       
!     "NORM" IS CALCULATED IN SUBROUTINE "USERIN" AND                   
!     NORMALIZES THE FUNCTION TO UNITY.  FP IS ASSUMED                  
!     SEPARABLE IN THETA AND PHI.                                       
                                                                        
      REAL KAPPA, THETA, PHI, PTHETA, PPHI, NORM, DEGRAD 
      DATA DEGRAD /0.0174532925/ 
                                                                        
      PTHETA = EXP(KAPPA*COS(DEGRAD*THETA))*SIN(DEGRAD*THETA) 
      PPHI = 1.0/360.0 
      FP = NORM*PTHETA*PPHI 
                                                                        
      RETURN 
      END                                           
!     END OF FUNCTION FP()                                              
!                                                                       
!                                                                       
!                                                                       
!_________________________________________________________________      
      SUBROUTINE QDWTS (NUM, ASBCIS, WTS) 
      INTEGER NUM, I 
      REAL  ASBCIS(1), WTS(1) 
      REAL  X(100), ALPHAS(100) 
!                                                                       
      DO 100 I = 1, NUM 
          X(I) = ASBCIS(I)**2 
  100 END DO 
      DO 120 I = 1, NUM 
        ALPHAS(I) = 1.00/(2*I-1) 
  120 END DO 
      CALL VANDER (X, WTS, ALPHAS, NUM) 
!                                                                       
      RETURN 
      END                                           
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!_________________________________________________________________      
      SUBROUTINE GLQUAD                                                 &
     &                          (NUM, ASBCIS, WTS)                      
!        GENERATES THE ABSCISSAS AND WEIGHTS FOR AN EVEN 2*NUM POINT    
!      GAUSS-LEGENDRE QUADRATURE.  ONLY THE NUM POSITIVE POINTS ARE RETU
      INTEGER  NUM 
      REAL   ASBCIS(1), WTS(1) 
      INTEGER  N, I, J, K 
      REAL   Z, Z1, P1, P2, P3, PP, EPS 
      PARAMETER (EPS=3.0E-14) 
!                                                                       
      N=2*NUM 
      DO 12 I=1,NUM 
        Z=COS(3.14159265400*(I-.2500)/(N+.500)) 
        K=0 
    1   CONTINUE 
          P1=1.00 
          P2=0.00 
          DO 11 J=1,N 
            P3=P2 
            P2=P1 
            P1=((2.00*J-1.00)*Z*P2-(J-1.00)*P3)/J 
   11     CONTINUE 
          PP=N*(Z*P1-P2)/(Z*Z-1.00) 
          Z1=Z 
          Z=Z1-P1/PP 
          K=K+1 
        IF (ABS(Z-Z1).GT.EPS .AND. K.LT.10) GO TO 1 
        ASBCIS(NUM+1-I)=Z 
        WTS(NUM+1-I)=2.00/((1.00-Z*Z)*PP*PP) 
   12 END DO 
!                                                                       
      RETURN 
      END                                           
!                                                                       
!                                                                       
!                                                                       
!_________________________________________________________________      
      SUBROUTINE DGQUAD                                                 &
     &                          (NUM, ASBCIS, WTS)                      
!        GENERATES THE ABSCISSAS AND WEIGHTS FOR AN NUM POINT           
!      GAUSS-LEGENDRE QUADRATURE BETWEEN 0 AND 1.                       
      INTEGER  NUM 
      REAL   ASBCIS(1), WTS(1) 
      INTEGER  N, I, J, M, K 
      REAL   Z, Z1, P1, P2, P3, PP, EPS 
      PARAMETER (EPS=3.00-14) 
!                                                                       
      N=NUM 
      M=(N+1)/2 
      DO 12 I=1,M 
        Z=COS(3.14159265400*(I-.2500)/(N+.500)) 
        K=0 
    1   CONTINUE 
          P1=1.00 
          P2=0.00 
          DO 11 J=1,N 
            P3=P2 
            P2=P1 
            P1=((2.00*J-1.00)*Z*P2-(J-1.00)*P3)/J 
   11     CONTINUE 
          PP=N*(Z*P1-P2)/(Z*Z-1.00) 
          Z1=Z 
          Z=Z1-P1/PP 
          K=K+1 
        IF(ABS(Z-Z1).GT.EPS .AND. K.LT.10) GO TO 1 
        ASBCIS(I)=.500-.500*Z 
        ASBCIS(N+1-I)=.500+.500*Z 
        WTS(I)=1.00/((1.00-Z*Z)*PP*PP) 
        WTS(N+1-I)=WTS(I) 
   12 END DO 
!                                                                       
      RETURN 
      END                                           
!                                                                       
!                                                                       
!                                                                       
!_________________________________________________________________      
      SUBROUTINE LBQUAD (NUM, ASBCIS, WTS) 
!        GENERATES THE ABSCISSAS AND WEIGHTS FOR AN EVEN 2*NUM POINT    
!      LOBATTO QUADRATURE.  ONLY THE NUM POSITIVE POINTS ARE RETURNED.  
      INTEGER  NUM 
      REAL   ASBCIS(1), WTS(1) 
      INTEGER  N, N1, I, J, K 
      REAL   Z, Z1, P1, P2, P3, PP, PPP, CI, EPS 
      PARAMETER (EPS=3.00-14) 
!                                                                       
      N=2*NUM 
      N1=N-1 
      CI = 0.50 
      IF (MOD(N,2) .EQ. 1)  CI = 1.00 
      DO 12 I=1,NUM-1 
        Z=SIN(3.14159265400*(I-CI)/(N-.5)) 
        K=0 
    1   CONTINUE 
          P1=1.00 
          P2=0.00 
          DO 11 J=1,N1 
            P3=P2 
            P2=P1 
            P1=((2.00*J-1.00)*Z*P2-(J-1.00)*P3)/J 
   11     CONTINUE 
          PP=N1*(Z*P1-P2)/(Z*Z-1.00) 
          PPP=(2.00*Z*PP-N1*(N1+1)*P1)/(1.00-Z*Z) 
          Z1=Z 
          Z=Z1-PP/PPP 
          K=K+1 
        IF(ABS(Z-Z1).GT.EPS .AND. K.LT.10) GO TO 1 
        ASBCIS(I)=Z 
        WTS(I)=2.00/(N*(N-1)*P1*P1) 
   12 END DO 
      ASBCIS(NUM)=1.00 
      WTS(NUM)=2.00/(N*(N-1)) 
!                                                                       
      RETURN 
      END                                           
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!_________________________________________________________________      
      SUBROUTINE VANDER(X,W,Q,N) 
!     IMPLICIT REAL  (A-H,O-Z)                                          
      PARAMETER (NMAX=100,ZERO=0.000,ONE=1.000) 
      DIMENSION X(N),W(N),Q(N),C(NMAX) 
      IF(N.EQ.1)THEN 
        W(1)=Q(1) 
      ELSE 
        DO 11 I=1,N 
          C(I)=ZERO 
   11   CONTINUE 
        C(N)=-X(1) 
        DO 13 I=2,N 
          XX=-X(I) 
          DO 12 J=N+1-I,N-1 
            C(J)=C(J)+XX*C(J+1) 
   12     CONTINUE 
          C(N)=C(N)+XX 
   13   CONTINUE 
        DO 15 I=1,N 
          XX=X(I) 
          T=ONE 
          B=ONE 
          S=Q(N) 
          K=N 
          DO 14 J=2,N 
            K1=K-1 
            B=C(K)+XX*B 
            S=S+Q(K1)*B 
            T=XX*T+B 
            K=K1 
   14     CONTINUE 
          W(I)=S/T 
   15   CONTINUE 
      ENDIF 
      RETURN 
      END                                           
!                                                                       
!                                                                       
!                                                                       
                                                                        
!_________________________________________________________________      
! ... ................................................................. 
                                                                        
      SUBROUTINE RMP(NANG) 
!     CALCULATE RADAR MEASURABLE PARAMETERS                             
!     FROM MUELLER MATRIX (AT BACKSCATTER)                              
!                                                                       
                                                                        
      INTEGER NANG, CTR1, CTR2 
      REAL SCAT, EXT, CONK, S, SS, ZDPT, ZDP 
      REAL ZHH, ZDR, LDR, NUM, DENOM, MAGS11, MAGS12, MAGS22 
      REAL PHI11, PHI12, PHI22, B, PHI, THETA, I(4) 
      REAL IS, QS, US, VS, P, SIGMA, DEPTH, UMAGSQ 
      REAL RHOMAG, RHOFAZ, RADDEG, ARG, PI, DELHV 
      REAL MAG1, MAG2, RHORE, RHOIM, LAMBDA, NFAC, MAGKSQ, MRHO 
      REAL DPHA, DATT, ATTH, ATTV 
      COMPLEX S11, S12, S22, ROOT, RHO(4), A, C, U, RHOHV 
                                                                        
                                                                        
      COMMON /MAIN01/ ABV, ABH, CONK, S(91,4,4), SS(91,4,4), LAMBDA,    &
     &EXM(4,4),EXMP(4,4),EXMPT(4,4),SCHPT,SCVPT,EXHPT,EXVPT,PHNORM      
      COMMON /DIEL/ EPR 
      COMMON /DSD6/ DEQ 
      DATA PI, RADDEG /3.141592654, 57.295779513/ 
!                                                                       
!                                                                       
!     AS IN RADAR CASE MAGINITUDE OF |K|^2 IS TAKEN TO BE .92           
      MAGKSQ = .92 
      NFAC =1.E6*(LAMBDA**4)/((PI**5)*(MAGKSQ)) 
!                                                                       
!     TO GET BACK BACKSCATTER CROSS-SECTION MUTIPLY (SCAT*4*PI)         
! *** THIS IS BEING DONE IN THE MAIN PROGRAM    
      ! print*,'((SCHPT+SCVPT)*.5)=',((SCHPT+SCVPT)*.5)
      NFAC = NFAC*((SCHPT+SCVPT)*.5)*4.*PI 
      ! print*,'NFAC=',NFAC*((SCHPT+SCVPT)*.5)*4.*PI 
      ! NFAC= NFAC*4.*PI   
      ! print*,'NFAC=',NFAC
      ! STOP                        

!RAVI NFAC = NFAC*4.*PI                                                 
!     PARTKSQ = (((EPR-1.)/(EPR+2.))**2.)                               
!                                                                       
!     ABSOLUTE REFLECTIVITY IN LOG(MM^6/M^3)                            
      ZHH = 0.5*NFAC*(SS(NANG,1,1) - SS(NANG,1,2) - SS(NANG,2,1) +      &
     &           SS(NANG,2,2))                                          
      IF (ZHH .GT. 0.0) THEN 
         ZHH = 10.0*ALOG10(ZHH) 
      ELSE 
         ZHH = -500.0 
      ENDIF 
                                                                        
!     DIFFERENTIAL REFLECTIVITY                                         
      NUM = SS(NANG,1,1) - SS(NANG,1,2) - SS(NANG,2,1) + SS(NANG,2,2) 
      DENOM = SS(NANG,1,1) + SS(NANG,1,2) + SS(NANG,2,1) + SS(NANG,2,2) 
                                                                        
! *** THIS IS TO ACCOUNT FOR ZDPT=0. AT AN ELEVATION ANGLE = 90DEG      
! ***                                                   RAVI  12/23/91. 
      ZDPT = NUM - DENOM 
      IF(ZDPT .GT. 0.) THEN 
         ZDP = 10.*ALOG10(.5*NFAC*ZDPT) 
      ELSEIF(ZDPT .LT. 0.) THEN 
         ZDP = -10.*ALOG10(.5*NFAC*(-ZDPT)) 
      ELSE 
         ZDP = 0. 
      ENDIF 
                                                                        
      ZDR = 10.0*ALOG10(NUM/DENOM) 
                                                                        
!     LINEAR DEPOLARIZATION RATIO                                       
      NUM = SS(NANG,1,1) - SS(NANG,1,2) + SS(NANG,2,1) - SS(NANG,2,2) 
      DENOM = SS(NANG,1,1) - SS(NANG,1,2) - SS(NANG,2,1) + SS(NANG,2,2) 
      LDR = NUM/DENOM 
      IF (LDR .GT. 0.0) THEN 
         LDR = 10.0*ALOG10(LDR) 
      ELSE 
         LDR = -500.0 
      ENDIF 
                                                                        
!     RHO                                                               
      MAG1 = SS(NANG,1,1) - SS(NANG,1,2) - SS(NANG,2,1) + SS(NANG,2,2) 
      MAG2 = SS(NANG,1,1) + SS(NANG,1,2) + SS(NANG,2,1) + SS(NANG,2,2) 
      DENOM = SQRT(MAG1)*SQRT(MAG2) 
                                                                        
!     -1 IS INTRODUCED IN RHORE BECAUSE OF S(3,3), S(4,4) ARE -VE IN    
!     SCATTERING MATRIX. S(3,3) AND S(4,4) SHOULD BE +VE .              
                                                                        
      RHORE =-1.*((SS(NANG,3,3) + SS(NANG,4,4)))/DENOM 
      RHOIM = -(SS(NANG,4,3) - SS(NANG,3,4))/DENOM 
      RHOHV = CMPLX(RHORE,RHOIM) 
      DELHV = ATAN(AIMAG(RHOHV)/REAL(RHOHV))*RADDEG 
      MRHO = CABS(RHOHV) 
                                                                        
!     CO-POL AND CROSS-POL NULLS                                        
!     (POLARIZATION RATIO, SPREAD, DEPTH):                              
                                                                        
!        CALCULATE SCATTERING-MATRIX ELEMENTS                           
         MAGS11 = SQRT(0.5*(SS(NANG,1,1) + SS(NANG,1,2) +               &
     &                      SS(NANG,2,1) + SS(NANG,2,2)))               
         MAGS12 = SQRT(0.5*(SS(NANG,1,1) - SS(NANG,1,2) +               &
     &                      SS(NANG,2,1) - SS(NANG,2,2)))               
         MAGS22 = SQRT(0.5*(SS(NANG,1,1) - SS(NANG,1,2) -               &
     &                      SS(NANG,2,1) + SS(NANG,2,2)))               
         PHI12 = 0.0 
         PHI11 = PHI12 - ATAN2((SS(NANG,1,4) + SS(NANG,2,4)),           &
     &                         (SS(NANG,1,3) + SS(NANG,2,3)))           
         PHI22 = PHI12 - ATAN2((SS(NANG,4,1) - SS(NANG,4,2)),           &
     &                         (SS(NANG,3,1) - SS(NANG,3,2)))           
         S11 = CMPLX(MAGS11*COS(PHI11),MAGS11*SIN(PHI11)) 
         S12 = CMPLX(MAGS12*COS(PHI12),MAGS12*SIN(PHI12)) 
         S22 = CMPLX(MAGS22*COS(PHI22),MAGS22*SIN(PHI22)) 
                                                                        
!        CALCULATE MEAN VALUES OF CO-POL, CROSS-POL NULLS               
!        CO-POL:                                                        
            ROOT = CSQRT(S12*S12 - S11*S22) 
            RHO(1) = (-S12 + ROOT)/S22 
            RHO(2) = (-S12 - ROOT)/S22 
!        CROSS-POL:                                                     
            A = S22*CONJG(S12) + CONJG(S11)*S12 
            B = MAGS22*MAGS22 - MAGS11*MAGS11 
            C = -CONJG(A) 
            ROOT = CSQRT(B*B - 4.0*A*C) 
            IF (MAGS12 .EQ. 0.00) THEN 
!              DEGENERATE CASE (LINEAR, NOT QUADRATIC)                  
               RHO(3) = (0.0,0.0) 
               RHO(4) = RHO(3) 
            ELSE 
               RHO(3) = (-B + ROOT)/(2.0*A) 
               RHO(4) = (-B - ROOT)/(2.0*A) 
            ENDIF 
                                                                        
!         WRITE (3,904)                                                 
         DO 200 CTR1 = 1,4 
                                                                        
                                                                        
!           NORMALIZED INPUT STOKES VECTOR PARAMETERS                   
            IF (RHO(CTR1) .EQ. (0.0,1.0)) THEN 
               PHI = -PI/2.0 
               THETA = 0.0 
            ELSE 
               U = (1.0 - (0.0,1.0)*RHO(CTR1))/                         &
     &             (1.0 + (0.0,1.0)*RHO(CTR1))                          
               IF (CABS(U) .EQ. 0.0) THEN 
                  PHI = 0.0 
               ELSE 
                  PHI = -ATAN2(AIMAG(U),REAL(U)) 
               ENDIF 
               UMAGSQ = REAL(U)*REAL(U) + AIMAG(U)*AIMAG(U) 
               THETA = ACOS((UMAGSQ - 1.0)/                             &
     &                      (UMAGSQ + 1.0))                             
            ENDIF 
            I(1) = 1.0 
            I(2) = SIN(THETA)*COS(PHI) 
            I(3) = SIN(THETA)*SIN(PHI) 
            I(4) = -COS(THETA) 
                                                                        
!           OUTPUT (SCATTERED-WAVE) STOKES VECTOR                       
            IS = 0.0 
            QS = 0.0 
            US = 0.0 
            VS = 0.0 
            DO 100 CTR2 = 1,4 
               IS = IS + SS(NANG,1,CTR2)*I(CTR2) 
               QS = QS + SS(NANG,2,CTR2)*I(CTR2) 
               US = US + SS(NANG,3,CTR2)*I(CTR2) 
               VS = VS + SS(NANG,4,CTR2)*I(CTR2) 
  100       CONTINUE 
                                                                        
!           DEGREE OF POLARIZATION FOR THIS NULL                        
            P = SQRT((QS*QS + US*US + VS*VS)/(IS*IS)) 
                                                                        
!           SPREAD, DEPTH                                               
            ARG = 1.0 - P*P 
            IF (ABS(ARG) .LT. 1.0E-03) THEN 
               SIGMA = 0.0 
            ELSE 
               SIGMA = 0.5*SQRT(1.0 - P*P) 
            ENDIF 
            DEPTH = 0.25*(1.0 - P*P) 
                                                                        
!           OUTPUT TO FILE                                              
            RHOMAG = SQRT(REAL(RHO(CTR1))*REAL(RHO(CTR1)) +             &
     &                    AIMAG(RHO(CTR1))*AIMAG(RHO(CTR1)))            
            IF (CABS(RHO(CTR1)) .EQ. 0.0) THEN 
               RHOFAZ = 0.0 
            ELSE 
               RHOFAZ = RADDEG*ATAN2(AIMAG(RHO(CTR1)),REAL(RHO(CTR1))) 
            ENDIF 
            IF (CTR1 .LE. 2) THEN 
!               WRITE (3,905) RHOMAG, RHOFAZ, SIGMA, DEPTH              
!V               CALL POLE(RHO(CTR1))                                   
            ELSE 
!               WRITE (3,906) RHOMAG, RHOFAZ, SIGMA, DEPTH              
!V               CALL POLE(RHO(CTR1))                                   
            ENDIF 
  200    CONTINUE 
                                                                        
         DO  300 IR = 1, 4 
         DO  300 IC = 1, 4 
         SS(1,IR,IC) = SS(1,IR,IC)*((SCHPT+SCVPT)*.5) 
         SS(5,IR,IC) = SS(5,IR,IC)*((SCHPT+SCVPT)*.5) 
  300    CONTINUE 
!                                                                       
! ***    DPHA NEEDS TO MULTIPLIED BY 0.1 *** RAVI 12/20/91 ***          
! RAVI   DPHA = RADDEG*EXMPT(3,4)                                       
!v         DPHA = RADDEG*EXMPT(3,4)*0.1                                 
!v      see the explantion below which was added on 3/15/94 for the fact
!v      in DPHA, ATTH, ATTV, DATT parameters.                           
!         DPHA = 2.*RADDEG*EXMPT(3,4)*0.1  ::: compare to mueller progra
!        cwy 20050817                                                   
         DPHA = RADDEG*EXMPT(3,4)*0.1 
! ***    THESE ARE FOR ATTENUATION (dB/Km) AT H, V -POLS.               
! ***    (1,1)=MHH+MVV  (1,2)=MVV-MHH                                   
!v                                                                      
!v       Attenuation (i.e. extinction) is 4*pi/K times the forward ampli
!v       according to optical theorm.                                   
!v       EXMPT matrices are made of M_vv, M_hh and M_hv terms           
!v       which are for 2_by_2 field extinction matrices for             
!v       field quantities Eq.(2) and (3) Page 139 Tsangs book.          
!v         ATTH = 0.4343*0.5*(EXMPT(1,1)-EXMPT(1,2))                    
!v         ATTV = 0.4343*0.5*(EXMPT(1,1)+EXMPT(1,2))                    
!v         DATT = -.4343*EXMPT(1,2)                                     
         ATTH = 0.4343*(EXMPT(1,1)-EXMPT(1,2)) 
         ATTV = 0.4343*(EXMPT(1,1)+EXMPT(1,2)) 
         DATT = -0.4343*2.*EXMPT(1,2) 
                                                                        
!        WRITE(3,'(X,4E14.6,X,E14.6,X,A)') (EXMPT(1,J),J=1,4),          
         WRITE(109,'(X,4E14.6)') (EXMPT(1,J),J=1,4) 
!    .   PHNORM, 'PHNORM'                                               
!        WRITE(3,'(X,4E14.6,E14.6,X,A)') (EXMPT(2,J),J=1,4),            
         WRITE(109,'(X,4E14.6)') (EXMPT(2,J),J=1,4) 
!    .   DPHA,'DPH DEG'                                                 
!        WRITE(3,'(X,4E14.6,X,F7.2,X,F5.2,X,A,X,A)') (EXMPT(3,J),J=1,4),
         WRITE(109,'(X,4E14.6)') (EXMPT(3,J),J=1,4) 
!    .   ZDP,MRHO,'RHO','ZDP DB'                                        
!        WRITE(3,'(X,4E14.6,X,E14.6,X,A)') (EXMPT(4,J),J=1,4),          
         WRITE(109,'(X,4E14.6)') (EXMPT(4,J),J=1,4) 
!    .   DELHV,'DHV DEG'                                                
!        WRITE(3,'(X,4E14.6,X,E14.6,X,A)') (SS(1,1,J),J=1,4),           
!    .   EXHPT,'EXH CM^2'                                               
!        WRITE(3,'(X,4E14.6,X,E14.6,X,A)') (SS(1,2,J),J=1,4),           
!    .   EXVPT,'EXV CM^2'                                               
!        WRITE(3,'(X,4E14.6,X,E14.6,X,A)') (SS(1,3,J),J=1,4),           
!    .   SCHPT,'SCH CM^2'                                               
!        WRITE(3,'(X,4E14.6,X,E14.6,X,A)') (SS(1,4,J),J=1,4),           
!    .   SCVPT,'SCV CM^2'                                               
                                                                        
!        WRITE(3,'(X,4E14.6,2X,E14.6,X,A)') (SS(5,1,J),J=1,4),          
         WRITE(108,'(X,4E14.6)') (SS(5,1,J),J=1,4) 
!    .   ZHH,'ZHH DBZ'                                                  
!        WRITE(3,'(X,4E14.6,2X,E14.6,X,A)') (SS(5,2,J),J=1,4),          
         WRITE(108,'(X,4E14.6)') (SS(5,2,J),J=1,4) 
!    .   ZDR,'ZDR DB'                                                   
!        WRITE(3,'(X,4E14.6,X,E14.6,X,A)') (SS(5,3,J),J=1,4),           
         WRITE(108,'(X,4E14.6)') (SS(5,3,J),J=1,4) 
!    .   LDR,'LDR DB'                                                   
!        WRITE(3,'(X,4E14.6,X,E14.6,X,A)') (SS(5,4,J),J=1,4),           
!    .   DATT,'DATT DB'                                                 
         WRITE(108,'(X,4E14.6)') (SS(5,4,J),J=1,4) 
                                                                        
!        WRITE (3,*) '          '                                       
!                                                                       
!                                                                       
!  ***   THESE STATEMENTS HAVE BEEN MADE BY RAVI ON 7/18/91. ***        
!                                                                       
      WRITE(7,*) 
      WRITE(7,FMT ='(E22.13,8X,"PARTICLE DIAMETER IN CM")')DEQ 
      WRITE(7,FMT ='(E22.13,8X,"REFLECTIVITY DBZ")')ZHH 
      WRITE(7,FMT ='(E22.13,8X,"DIFFERENTIAL REFLECTIVITY DB")')ZDR 
      WRITE(7,FMT ='(E22.13,8X,"LINEAR DEPOL RATIO DB")')LDR 
      WRITE(7,FMT ='(E22.13,8X,"DIFFERENTIAL PHASE DEG/KM")')DPHA 
      WRITE(7,FMT ='(E22.13,8X,"DIFFERENCE REFLECTIVITY DB")')ZDP 
      WRITE(7,FMT ='(E22.13,8X,"HV CORRELATION")')MRHO 
      WRITE(7,FMT ='(E22.13,8X,"PHASESHIFT UPON BACKSCATTER DEG")')DELHV 
      WRITE(7,FMT ='(E22.13,8X,"ATTENUATION HPOL DB/KM")')ATTH 
      WRITE(7,FMT ='(E22.13,8X,"ATTENUATION VPOL DB/KM")')ATTV 
      WRITE(7,FMT ='(E22.13,8X,"DIFFERENTIAL ATTENUATION DB/KM")')DATT 
!                                                                       
! *** END OF MODIFICATIONS BY RAVI ON 7/19/81. ***                      
!                                                                       
!                                                                       
                                                                        
      RETURN 
  900 FORMAT (1X,'ZHH = ',F6.1,' DB') 
  901 FORMAT (1X,'ZDR = ',F6.1,' DB') 
  902 FORMAT (1X,'LDR = ',F6.1,' DB') 
  903 FORMAT (1X,'RHO (H,V) = ',2(F6.2,2X)) 
  904 FORMAT (/,20X,'NULL',9X,'SPREAD    DEPTH',/,15X,'(MAG)   (PHASE)') 
  905 FORMAT (1X,'CO-POL:     ',F8.3,2X,F7.2,2(3X,F6.3)) 
  906 FORMAT (1X,'CROSS-POL:  ',F8.3,2X,F7.2,2(3X,F6.3)) 
      END                                           
!     END OF SUBROUTINE RMP                                             
!                                                                       
!_________________________________________________________________      
!                                                                       
      SUBROUTINE POLE(RHO) 
!                                                                       
!     CALCULATE AXIS RATIO AND TILT FROM COMPLEX RHO                    
!                                                                       
!                                                                       
      REAL TAU, PHI, MAGUSQ, R 
      COMPLEX RHO, U, J 
      DATA PI, RADDEG /3.141592654, 57.295779513/ 
                                                                        
      J = CMPLX(0.,1.) 
      U = (1.-J*RHO)/(1.+J*RHO) 
      MAGUSQ = CABS(U)*CABS(U) 
      TAU = 45.-.5*RADDEG*ACOS((MAGUSQ-1.)/(MAGUSQ+1.)) 
      R = 1./(TAN(TAU/RADDEG)) 
      PHI=-.5*RADDEG*ATAN(AIMAG(U)/REAL(U)) 
!      WRITE (3,"('TAU, R, PHI', 3(E15.7,2X))") TAU, R, PHI             
      RETURN 
      END                                           
!                                                                       
!     END OF SUBROUTINE POLE                                            
!                                                                       
! ... ................................................................. 
!                                                                       
      SUBROUTINE SETPAR(DISTYP) 
!                                                                       
! *** THIS ROUTINE SETS PARAMETERS THAT ARE REQUIRED FOR PROCESSING     
! *** DIFFERENT TYPES ( SPECIES ) OF PARTICLES.          RAVI: 12/26/91 
!                                                                       
      INTEGER DISTYP 
      LOGICAL RDROP, FRDROP, DHAIL, SPHAIL, WHAIL, DGRAUP, SPGRAUP,     &
     &        WGRAUP, NCRYST, PCRYST, SPAGG, AGG2, AGG5, AGG8           
      COMMON /DSD1/ DMIN, DMAX, DSTEP 
      COMMON /DSD2/ DMIN1, DMIN2, DMIN3, DMIN4, DMIN5, DMIN6, DMIN7,    &
     &       DMIN8, DMIN9, DMIN10, DMIN11, DMIN12, DMIN13, DMIN14,      &
     &              DMAX1, DMAX2, DMAX3, DMAX4, DMAX5, DMAX6, DMAX7,    &
     &       DMAX8, DMAX9, DMAX10, DMAX11, DMAX12, DMAX13, DMAX14,      &
     &              DSTEP1, DSTEP2, DSTEP3, DSTEP4, DSTEP5, DSTEP6,     &
     &              DSTEP7, DSTEP8, DSTEP9, DSTEP10, DSTEP11, DSTEP12,  &
     &              DSTEP13, DSTEP14                                    
      COMMON /DSD3/ DPAR1, DPAR2, DPAR3 
      COMMON /DSD4/ S1PAR1, S1PAR2, S1PAR3, S2PAR1, S2PAR2, S2PAR3,     &
     &                 S3PAR1, S3PAR2, S3PAR3, S4PAR1, S4PAR2, S4PAR3,  &
     &                 S5PAR1, S5PAR2, S5PAR3, S6PAR1, S6PAR2, S6PAR3,  &
     &                 S7PAR1, S7PAR2, S7PAR3, S8PAR1, S8PAR2, S8PAR3,  &
     &                 S9PAR1, S9PAR2, S9PAR3, S10PAR1, S10PAR2,        &
     &                 S10PAR3, S11PAR1, S11PAR2, S11PAR3, S12PAR1,     &
     &                 S12PAR2, S12PAR3, S13PAR1, S13PAR2, S13PAR3,     &
     &                 S14PAR1, S14PAR2, S14PAR3                        
      COMMON /DSD5/ RDROP, FRDROP, DHAIL, SPHAIL, WHAIL, DGRAUP,        &
     & SPGRAUP, WGRAUP, NCRYST, PCRYST, SPAGG, AGG2, AGG5, AGG8         
!                                                                       
! ... ................................................................. 
      IF (RDROP) THEN 
        DMIN = DMIN1 
        DMAX = DMAX1 
        DSTEP = DSTEP1 
        DPAR1 = S1PAR1 
        DPAR2 = S1PAR2 
        DPAR3 = S1PAR3 
        RDROP = .FALSE. 
!      print*,'cwy',DMIN,DMAX,DSTEP,DPAR1,DPAR2,DPAR3                   
      ELSEIF (FRDROP) THEN 
        DMIN = DMIN2 
        DMAX = DMAX2 
        DSTEP = DSTEP2 
        DPAR1 = S2PAR1 
        DPAR2 = S2PAR2 
        DPAR3 = S2PAR3 
        FRDROP = .FALSE. 
      ELSEIF (DHAIL) THEN 
        DMIN = DMIN3 
        DMAX = DMAX3 
        DSTEP = DSTEP3 
        DPAR1 = S3PAR1 
        DPAR2 = S3PAR2 
        DPAR3 = S3PAR3 
        DHAIL = .FALSE. 
      ELSEIF (SPHAIL) THEN 
        DMIN = DMIN4 
        DMAX = DMAX4 
        DSTEP = DSTEP4 
        DPAR1 = S4PAR1 
        DPAR2 = S4PAR2 
        DPAR3 = S4PAR3 
        SPHAIL = .FALSE. 
      ELSEIF (WHAIL) THEN 
        DMIN = DMIN5 
        DMAX = DMAX5 
        DSTEP = DSTEP5 
        DPAR1 = S5PAR1 
        DPAR2 = S5PAR2 
        DPAR3 = S5PAR3 
        WHAIL = .FALSE. 
      ELSEIF (DGRAUP) THEN 
        DMIN = DMIN6 
        DMAX = DMAX6 
        DSTEP = DSTEP6 
        DPAR1 = S6PAR1 
        DPAR2 = S6PAR2 
        DPAR3 = S6PAR3 
        DGRAUP = .FALSE. 
      ELSEIF (SPGRAUP) THEN 
        DMIN = DMIN7 
        DMAX = DMAX7 
        DSTEP = DSTEP7 
        DPAR1 = S7PAR1 
        DPAR2 = S7PAR2 
        DPAR3 = S7PAR3 
        SPGRAUP = .FALSE. 
      ELSEIF (WGRAUP) THEN 
        DMIN = DMIN8 
        DMAX = DMAX8 
        DSTEP = DSTEP8 
        DPAR1 = S8PAR1 
        DPAR2 = S8PAR2 
        DPAR3 = S8PAR3 
        WGRAUP = .FALSE. 
      ELSEIF (NCRYST) THEN 
        DMIN = DMIN9 
        DMAX = DMAX9 
        DSTEP = DSTEP9 
        DPAR1 = S9PAR1 
        DPAR2 = S9PAR2 
        DPAR3 = S9PAR3 
        NCRYST = .FALSE. 
      ELSEIF (PCRYST) THEN 
        DMIN = DMIN10 
        DMAX = DMAX10 
        DSTEP = DSTEP10 
        DPAR1 = S10PAR1 
        DPAR2 = S10PAR2 
        DPAR3 = S10PAR3 
        PCRYST = .FALSE. 
      ELSEIF (SPAGG) THEN 
        DMIN = DMIN11 
        DMAX = DMAX11 
        DSTEP = DSTEP11 
        DPAR1 = S11PAR1 
        DPAR2 = S11PAR2 
        DPAR3 = S11PAR3 
        SPAGG = .FALSE. 
      ELSEIF (AGG2) THEN 
        DMIN = DMIN12 
        DMAX = DMAX12 
        DSTEP = DSTEP12 
        DPAR1 = S12PAR1 
        DPAR2 = S12PAR2 
        DPAR3 = S12PAR3 
        AGG2 = .FALSE. 
      ELSEIF (AGG5) THEN 
        DMIN = DMIN13 
        DMAX = DMAX13 
        DSTEP = DSTEP13 
        DPAR1 = S13PAR1 
        DPAR2 = S13PAR2 
        DPAR3 = S13PAR3 
        AGG5 = .FALSE. 
      ELSEIF (AGG8) THEN 
        DMIN = DMIN14 
        DMAX = DMAX14 
        DSTEP = DSTEP14 
        DPAR1 = S14PAR1 
        DPAR2 = S14PAR2 
        DPAR3 = S14PAR3 
        AGG8 = .FALSE. 
      ENDIF 
!                                                                       
      WRITE(7,FMT ='(F22.2,8X,"MIN PARTICLE DIAMETER CMS")')DMIN 
      WRITE(7,FMT ='(F22.2,8X,"MAX PARTICLE DIAMETER CMS")')DMAX 
      WRITE(7,FMT ='(F22.3,8X," INCREMENTAL DIAMETER CMS")')DSTEP 
      IF ( DISTYP .EQ. 0 ) THEN 
        WRITE(7,FMT ='(10X,"RANDOM PDF FOR ORIENTATION")') 
      ELSEIF ( DISTYP .EQ. 1 ) THEN 
        WRITE(7,FMT ='(10X,"S HARMONIC OSCILLATION FOR ORIENTATION")') 
        WRITE(7,FMT ='(F22.2,8X,"MEAN OF OSCILLATION DEG")')DPAR1 
        WRITE(7,FMT ='(F22.2,8X,"AMPLITUDE OF OSCILLATION DEG")')DPAR2 
      ELSEIF ( DISTYP .EQ. 2 ) THEN 
        WRITE(7,FMT ='(10X,"GAUSSIAN PDF FOR ORIENTATION")') 
        WRITE(7,FMT ='(F22.2,8X,"MEAN OF DISTRIBUTION DEG")')DPAR1 
        WRITE(7,FMT ='(F22.2,8X,"SIGMA OF DISTRIBUTION DEG")')DPAR2 
      ELSEIF ( DISTYP .EQ. 3 ) THEN 
        WRITE(7,FMT ='(10X,"LANGEVIN PDF FOR ORIENTATION")') 
        WRITE(7,FMT ='(F22.2,8X,"KAPPA OF DISTRIBUTION DEG")')DPAR1 
        WRITE(7,FMT ='(F22.2,8X,"THETA0 OF DISTRIBUTION DEG")')DPAR2 
        WRITE(7,FMT ='(F22.2,8X,"PHI0 OF DISTRIBUTION DEG")')DPAR3 
      ELSEIF ( DISTYP .EQ. 4 ) THEN 
        WRITE(7,FMT ='(10X,"FISHER PDF FOR ORIENTATION")') 
        WRITE(7,FMT ='(F22.2,8X,"KAPPA OF DISTRIBUTION DEG")')DPAR1 
      ENDIF 
! ... ................................................................. 
!                                                                       
      RETURN 
      END                                           
