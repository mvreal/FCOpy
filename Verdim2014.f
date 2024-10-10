      PROGRAM VERDIM 
      COMMON E,FYD,FCD,NP,XP,YP,XG,YG,LX,LY,AREA,JOX,JOY,JOXY, 
     * NA,MAX,MAY,NB,XB,YB,PERC,NRC,IL(5),GS,GC,AS 
     * ,EPSC2,EPSCU,A1,A2,CONST 
      REAL*8 EPSC2(5),EPSCU(5),A1(5),A2(5),CONST,FCK 
      REAL*8 NN,Y,X,X2,X3,X4,XY,X2Y,DDX 
      REAL*8 E,FYD(100),FCD(5),XP(100),YP(100),XG,YG,AREA,XB(100), 
     * YB(100),PERC(100) 
      REAL*8 JOX,JOY,JOXY,LX,LY,NA,MAX,MAY 
      REAL*8 JX,JY,JXY,AS,DX,DY,GC,GS,SOX,SOY,SX,SY,XMAX,XMIN,YMAX,YMIN 
      CHARACTER ARQ*12 
      INTEGER OP 
   10 WRITE(*,11) 
   11 FORMAT(////,'>>>>>> QUAL O NOME DO ARQUIVO DE DADOS ?') 
      READ(*,'(A12)',ERR=10) ARQ 
      IF(ARQ.EQ.' ') STOP 
      IR=1 
      IW=0 
      OPEN(UNIT=1,FILE=ARQ,STATUS='OLD',FORM='FORMATTED', 
     *ACCESS='SEQUENTIAL') 
      READ(IR,21) OP 
      READ(IR,20) GC,GS 
      READ(IR,21) NP 
      DO 30 J1=1,NP 
         READ(IR,20) XP(J1),YP(J1) 
   30 CONTINUE 
      READ(IR,21) NRC 
      DO 40 J1=1,NRC 
         READ(IR,24) FCD(J1),IL(J1) 
         FCD(J1)=FCD(J1)/GC 
   40 CONTINUE 
      READ(IR,22) NB,E,AS 
      DO 50 J1=1,NB 
         READ(IR,23) XB(J1),YB(J1),FYD(J1),PERC(J1) 
         FYD(J1)=FYD(J1)/GS 
   50 CONTINUE 
      READ(IR,20) NA,MAX,MAY 
      CONST=210000.D0/E 
      DO 55 J1=1,NRC 
          FCK = FCD(J1)*GC*CONST 
          IF(FCK.LE.50) THEN 
              EPSC2(J1)=0.002D0 
              EPSCU(J1)=0.0035D0 
              A1(J1)=1000.D0 
              A2(J1)=250000.D0 
          ELSE 
              EPSC2(J1)=0.002D0+0.000085d0*(FCK-50.D0)**0.53 
              EPSCU(J1)=0.0026D0+0.035D0*((90-FCK)/100)**4 
              NN=1.4D0+23.4D0*((90-FCK)/100)**4 
              DDX=EPSC2(J1)/1000.D0 
              X=0.D0 
              X2=0.D0 
              X3=0.D0 
              X4=0.D0 
              XY=0.D0 
              X2Y=0.D0 
              DO 51 J2=1,1000 
                  Y=1.D0-(1.D0-X/EPSC2(J1))**NN 
                  X2=X2+X*X 
                  X3=X3+X*X*X 
                  X4=X4+X*X*X*X 
                  XY=XY+X*Y 
                  X2Y=X2Y+X*X*Y 
                  X=X+DDX 
   51         CONTINUE               
              A1(J1)=(X2Y-X4*XY/X3)/(X3-X2*X4/X3) 
              A2(J1)=-(X2Y-X3*A1(J1))/X4 
          END IF 
   55 CONTINUE              
      XMAX=-1.E11 
      YMAX=-1.E11 
      XMIN=1E11 
      YMIN=1E11 
      DO 60 J1=1,NP 
         IF(XP(J1).GT.XMAX) XMAX=XP(J1) 
         IF(XP(J1).LT.XMIN) XMIN=XP(J1) 
         IF(YP(J1).GT.YMAX) YMAX=YP(J1) 
         IF(YP(J1).LT.YMIN) YMIN=YP(J1) 
   60 CONTINUE 
      LX=XMAX-XMIN 
      LY=YMAX-YMIN 
      NP1=NP-1 
      AREA=0 
      SX=0 
      SY=0 
      JX=0 
      JY=0 
      JXY=0 
      DO 70 J1=1,NP1 
         DX=XP(J1+1)-XP(J1) 
         DY=YP(J1+1)-YP(J1) 
         AREA=AREA+(XP(J1)+DX/2.)*DY 
         SX=SX+(XP(J1)*(YP(J1)+DY/2.)+DX*(YP(J1)/2.+DY/3.))*DY 
         SY=SY+(XP(J1)*(XP(J1)+DX)+DX*DX/3.)*DY/2. 
         JX=JX+(XP(J1)*(YP(J1)*(DY+YP(J1))+DY*DY/3.)+DX*(YP(J1)* 
     .    (YP(J1)/2.+DY/1.5)+DY*DY/4.))*DY 
         JY=JY+(DX**3/4.+XP(J1)*(DX*DX+XP(J1)*(1.5*DX+XP(J1))))*DY/3. 
         JXY=JXY+(XP(J1)*(XP(J1)*(YP(J1)+DY/2.)+DX*(YP(J1)+DY/1.5))+ 
     .    DX*DX*(YP(J1)/3.+DY/4.))*DY/2. 
   70 CONTINUE 
      XG=SY/AREA 
      YG=SX/AREA 
      SOX=SX-YG*AREA 
      SOY=SY-XG*AREA 
      JOX=JX-AREA*YG*YG 
      JOY=JY-AREA*XG*XG 
      JOXY=JXY-XG*YG*AREA 
      DO 80 J1=1,NP 
         XP(J1)=XP(J1)-XG 
         YP(J1)=YP(J1)-YG 
   80 CONTINUE 
      DO 90 J1=1,NB 
         XB(J1)=XB(J1)-XG 
         YB(J1)=YB(J1)-YG 
   90 CONTINUE 
      CALL AJUSTL(OP) 
      CLOSE(1) 
      GO TO 10 
   20 FORMAT(8F10.0) 
   21 FORMAT(8I10) 
   22 FORMAT(I10,7F10.0) 
   23 FORMAT(4F10.0) 
   24 FORMAT(F10.0,9I10) 
      END 
 
      SUBROUTINE AJUSTL(OP) 
      INTEGER OP 
      REAL*8 NR,MRX,MRY,LAM,LAMMIN,NRMIN,MRXMIN,MRYMIN 
      REAL*8 R(3,2),RT(3,3),DP(3) 
      COMMON E,FYD,FCD,NP,XP,YP,XG,YG,LX,LY,AREA,JOX,JOY,JOXY, 
     *   NA,MAX,MAY,NB,XB,YB,PERC,NRC,IL(5),GS,GC,AS 
     *   ,EPSC2,EPSCU,A1,A2 
      REAL*8 EPSC2(5),EPSCU(5),A1(5),A2(5),CONST 
      REAL*8 E,FYD(100),FCD(5),XP(100),YP(100),XG,YG,AREA,XB(100), 
     *   YB(100),PERC(100) 
      REAL*8 JOX,JOY,JOXY,LX,LY,NA,MAX,MAY 
      REAL*8 ALFA0,ALPH,ALPG,AS1,AS2,AS3,AS,B,BAS,C,CA,CA0,CAR,X, 
     *   EPSI,EPSS,FS,GC,GS,PJX,PJY,SA,SA0,SS,TOL,TOLE,PI,PI2,GRAUS 
      DATA PI,PI2,GRAUS/3.1415926535897932385,1.5707963267948966192, 
     *   57.29577951308232/ 
      TOLMIN=1.D0 
      K=0 
      IW=0 
      TOLE=1.D-8 
      BAS=MAX*MAX+MAY*MAY+NA*NA 
      LAM=1.D0 
      ALFA0=0 
      IF(JOX.EQ.JOY.AND.ABS(JOXY).GT.1E-5) THEN 
         ALFA0=PI2 
      ELSE 
         IF(JOX.NE.JOY) ALFA0=ATAN(-2*JOXY/(JOX-JOY))/2. 
      ENDIF 
      CA0=COS(ALFA0) 
      SA0=SIN(ALFA0) 
      CA0=CA0*CA0 
      SA0=SA0*SA0 
      SS=JOXY*SIN(2*ALFA0) 
      PJX=JOX*CA0+JOY*SA0-SS 
      PJY=JOY*CA0+JOX*SA0+SS 
      ALPH=0 
      IF(MAX.EQ.0) THEN 
         ALPH=-DSIGN(PI2,MAY) 
      ELSE 
         ALPH=ATAN(MAY*PJX/(MAX*PJY)) 
         IF(MAX.GT.0) ALPH=ALPH+PI 
      ENDIF 
      ALFA0=ALPH+ALFA0 
      UU=0.5D0 
  300 UU=(PI+UU)**5 
      UU=UU-INT(UU) 
      K0=0 
      X=(LX+LY)*UU 
      IF(OP.EQ.1) THEN 
         AS1=ABS(MAX)/(0.4*LY*FYD(1)) 
         AS2=ABS(MAY)/(0.4*LX*FYD(1)) 
         IF(NA.GT.0) THEN 
            AS3=NA/FYD(1) 
         ELSE 
            AS3=DMAX1(0D0,(NA-FCD(1)*AREA)/FYD(1)) 
         ENDIF 
         AS=AS1+AS2+AS3 
      ENDIF 
      ALPH=ALFA0 
  890 CALL ESFOR(E,FYD,FCD,NP,XP,YP,NB,XB,YB,PERC,X,ALPH,AS,B,C, 
     *EPSS,EPSI,NR,MRX,MRY,R,NRC,IL,EPSC2,EPSCU,A1,A2) 
      DP(1)=MAX-LAM*MRX 
      DP(2)=MAY-LAM*MRY 
      DP(3)=NA-LAM*NR 
      TOL=SQRT((DP(1)**2+DP(2)**2+DP(3)**2)/BAS) 
      IF(TOL.LE.TOLE) GO TO 900 
      K=K+1 
      K0=K0+1 
      IF(K0.LE.50) GO TO 301 
      IF(TOL.LT.TOLMIN) THEN 
        TOLMIN=TOL 
        MRXMIN=MRX 
        MRYMIN=MRY 
        NRMIN=NR 
        ASMIN=AS 
        EPSSMIN=EPSS 
        EPSIMIN=EPSI 
        ALPHMIN=ALPH 
        LAMMIN=LAM 
      END IF 
      GO TO 300 
  301 CA=COS(ALPH) 
      SA=SIN(ALPH) 
      RT(1,1)=LAM*(R(1,1)*CA-R(2,1)*SA) 
      RT(1,2)=LAM*(-MRY) 
      RT(2,1)=LAM*(R(1,1)*SA+R(2,1)*CA) 
      RT(2,2)=LAM*MRX 
      RT(3,1)=LAM*R(3,1) 
      RT(3,2)=0 
      IF(OP.EQ.2) THEN 
         RT(1,3)=MRX 
         RT(2,3)=MRY 
         RT(3,3)=NR 
      ELSE 
         RT(1,3)=R(1,2)*CA-R(2,2)*SA 
         RT(2,3)=R(1,2)*SA+R(2,2)*CA 
         RT(3,3)=R(3,2) 
      ENDIF 
      CALL PIVO(RT,DP,IVER) 
      IF(IVER.EQ.1) GO TO 300 
      X=X+DP(1) 
      IF(OP.EQ.2) THEN 
         LAM=LAM+DP(3) 
      ELSE 
         AS1=AS+DP(3) 
               AS=AS1 
      ENDIF 
      ALPH=ALPH+DP(2) 
      IF(ABS(ALPH).GT.2*PI) ALPH=SIGN(MOD(ALPH,2*PI),ALPH) 
      IF(K.LT.10000) GO TO 890 
      WRITE(*,*) ">>>> NAO CONVERGIU",TOLMIN 
      MRX=MRXMIN 
      MRY=MRYMIN 
      NR=NRMIN 
      AS=ASMIN 
      EPSS=EPSSMIN 
      EPSI=EPSIMIN 
      ALPH=ALPHMIN 
      LAM=LAMMIN 
  900 WRITE(IW,500) 
      WRITE(IW,500) 
      WRITE(IW,501) 
      IF(OP.EQ.2) THEN 
         WRITE(IW,503) 
      ELSE 
         WRITE(IW,502) 
      ENDIF 
      WRITE(IW,501) 
      WRITE(IW,500) 
      WRITE(IW,501) 
      WRITE(IW,504) 
      WRITE(IW,505) MAX,MAY,NA 
      WRITE(IW,506) MRX,MRY,NR 
      WRITE(IW,501) 
      WRITE(IW,500) 
      WRITE(IW,501) 
      IF(OP.EQ.2) THEN 
         FS=1./LAM 
         WRITE(IW,511) FS 
      ELSE 
         WRITE(IW,507) AS 
      ENDIF 
      WRITE(IW,508) EPSS 
      WRITE(IW,509) EPSI 
      IF(ABS(ALPH).GT.2*PI) ALPH=SIGN(MOD(ALPH,2*PI),ALPH) 
      ALPG=GRAUS*ALPH 
      WRITE(IW,510) ALPG 
      WRITE(IW,501) 
	  WRITE(IW,500) 
      WRITE(IW,500) 
      WRITE(*,512) 
      READ(*,513) CAR 
      RETURN 
  500 FORMAT(1X,78('*'))     
  501 FORMAT(1X,'**',74X,'**') 
  502 FORMAT(1X,'**     DIMENSIONAMENTO DE SECAO DE CONCRETO ARMADO A SO 
     *LICITACOES NORMAIS   **') 
  503 FORMAT(1X,'**      VERIFICACAO DE SECAO DE CONCRETO ARMADO A SOLIC 
     *ITACOES NORMAIS      **') 
  504 FORMAT(1X,'**                                           MX         
     *  MY          N      **') 
  505 FORMAT(1X,'**    ESFORCOS ATUANTES DE CALCULO:   ',3E12.4,'  **') 
  506 FORMAT(1X,'**    ESFORCOS RESISTENTES DE CALCULO:',3E12.4,'  **') 
  507 FORMAT(1X,'**',12X,'AREA TOTAL DE ARMADURA:               ', 
     *E15.4,9X,'**') 
  508 FORMAT(1X,'**',12X,'DEFORMACAO NA FIBRA SUPERIOR DA SECAO:', 
     *E15.4,9X,'**') 
  509 FORMAT(1X,'**',12X,'DEFORMACAO NA FIBRA INFERIOR DA SECAO:', 
     *E15.4,9X,'**') 
  510 FORMAT(1X,'**',12X,'INCLINACAO DA LINHA NEUTRA:           ', 
     *E15.4,9X,'**') 
  511 FORMAT(1X,'**',12X'RESERVA:                              ', 
     *E15.4,9X,'**') 
  512 FORMAT(/,25X,' TECLE <ENTER> PARA CONTINUAR',/) 
  513 FORMAT(A1) 
      END 
 
      SUBROUTINE PIVO(A,B,IVER) 
      REAL*8 A(3,3),B(3) 
      INTEGER II(18),IVER 
      REAL*8 AUX1,AUX2,DUM1,DUM2,DUM3,DUM4,DUM5 
      DATA II/1,2,3,1,3,2,2,1,3,2,3,1,3,1,2,3,2,1/ 
      IVER=0 
      DO 10 J1=1,6 
         JJ=3*J1 
         I=II(JJ-2) 
         J=II(JJ-1) 
         K=II(JJ) 
         IF(A(I,I).EQ.0) GO TO 10 
         AUX1=A(J,I)/A(I,I) 
         AUX2=A(K,I)/A(I,I) 
         DUM1=A(K,K)-A(I,K)*AUX2 
         DUM2=A(K,J)-A(I,J)*AUX2 
         DUM3=B(K)-B(I)*AUX2 
         AUX2=A(J,J)-A(I,J)*AUX1 
         IF(AUX2.EQ.0) GO TO 10 
         DUM4=(B(J)-B(I)*AUX1)/AUX2 
         DUM5=(A(J,K)-A(I,K)*AUX1)/AUX2 
         AUX1=DUM1-DUM2*DUM5 
         IF(AUX1.EQ.0) GO TO 10 
         B(K)=(DUM3-DUM2*DUM4)/AUX1 
         B(J)=DUM4-DUM5*B(3) 
         B(I)=(B(I)-B(J)*A(I,J)-B(K)*A(I,K))/A(I,I) 
         GO TO 20 
   10 CONTINUE 
      IVER=1 
   20 RETURN 
      END 
  

      SUBROUTINE ESFOR(E,FYD,FC,NP,XP,YP,NB,XB,YB,PERC,X,ALFA,AS,B,C, 
     *EPSS,EPSI,NRZ,MRX,MRY,R,NRC,IL,EPSC2,EPSCU,A1,A2) 
      REAL*8 NRZTI,NRZT,MRKS,MRET,NRZ,MRX,MRY,KS1I,KS2I,KS1II,KS2II 
      REAL*8 XP(*),YP(*),XB(*),YB(*),PERC(*),KSP(100),ETP(100), 
     *     KSB(100),ETB(100),R(3,2),FYD(*),FC(*),EPSP(100) 
      REAL*8 EPSC2(*),EPSCU(*),A1(*),A2(*) 
      REAL*8 X,DD,ALFA,AS,B,BLX,C,CA,CLX,DUM1,DUM2,E,EPS0,EPSB,EPS1, 
     *EPSI,EPSS,ET,ET01,ET12,ET1I,ET1II,ET2I,ET2II,FCD,SA,SIG,TETIA, 
     *TETIC,TETSA,TETSC,TKSIA,TKSIC,TKSSA,TKSSC 
      INTEGER IL(*) 
      CA=COS(ALFA) 
      SA=SIN(ALFA) 
      TETSC=0 
      TETIC=0 
      DO 10 J1=1,NP 
         KSP(J1)=XP(J1)*CA+YP(J1)*SA 
         ETP(J1)=-XP(J1)*SA+YP(J1)*CA 
         IF(ETP(J1).GT.TETSC) THEN 
            TETSC=ETP(J1) 
            TKSSC=KSP(J1) 
         ENDIF 
         IF(ETP(J1).LT.TETIC) THEN 
            TETIC=ETP(J1) 
            TKSIC=KSP(J1) 
         ENDIF 
   10 CONTINUE 
      TETSA=0 
      TETIA=0 
      DO 20 J1=1,NB 
         KSB(J1)=XB(J1)*CA+YB(J1)*SA 
         ETB(J1)=-XB(J1)*SA+YB(J1)*CA 
         IF(ETB(J1).GT.TETSA) THEN 
            TETSA=ETB(J1) 
            TKSSA=KSB(J1) 
         ENDIF 
         IF(ETB(J1).LT.TETIA) THEN 
            TETIA=ETB(J1) 
            TKSIA=KSB(J1) 
         ENDIF 
   20 CONTINUE 
      H=TETSC-TETIC 
      DD=TETSC-TETIA 
      X23=EPSCU(1)/(0.01D0+EPSCU(1))*DD 
      IF(X.LT.X23) THEN 
         B=-0.01D0/(DD-X) 
         C=B*(X-TETSC) 
         BLX=-0.01D0/(DD-X)**2 
         CLX=(DD-TETSC)*BLX 
         EPSI=0.01D0 
         EPSS=-0.01D0*X/(DD-X) 
      ELSE 
         IF(X.LT.H) THEN 
            B=-EPSCU(1)/X 
            C=-EPSCU(1)-B*TETSC 
            BLX=EPSCU(1)/X**2 
            CLX=-TETSC*BLX 
            EPSS=-EPSCU(1) 
            EPSI=DMAX1(EPSCU(1)*(DD-X)/X,0.D0) 
         ELSE 
            IF(X.GT.1.D150) THEN 
               B=0 
	           C=-EPSC2(1) 
               BLX=1.D-100 
               CLX=1.D-100 
               EPSS=-EPSC2(1) 
               EPSI=-EPSC2(1) 
               GO TO 29 
            END IF 
            B=-EPSC2(1)/(X-(EPSCU(1)-EPSC2(1))/EPSCU(1)*H) 
            C=B*(X-TETSC) 
            BLX=EPSC2(1)/(X-(EPSCU(1)-EPSC2(1))/EPSCU(1)*H)**2 
            CLX=(TETSC-(EPSCU(1)-EPSC2(1))/EPSCU(1)*H)*BLX 
            EPSS=-EPSC2(J1)*X/(X-(EPSCU(1)-EPSC2(1))/EPSCU(1)*H) 
            EPSI=-EPSC2(J1)*(X-H)/(X-(EPSCU(1)-EPSC2(1))/EPSCU(1)*H) 
         END IF 
      END IF 
   29 CONTINUE    
      DO 30 J1=1,NP 
         EPSP(J1)=B*ETP(J1)+C 
   30 CONTINUE 
      NRZT=0 
      MRKS=0 
      MRET=0 
      DO 50 J1=1,3 
         DO 50 J2=1,2 
            R(J1,J2)=0 
   50 CONTINUE 
      DO 60 J1=1,NB 
         EPSB=B*ETB(J1)+C 
         CALL ACO(E,EPSB,FYD(J1),SIG,ET) 
         DUM1=PERC(J1)*AS*ET*(BLX*ETB(J1)+CLX) 
         DUM2=PERC(J1)*SIG 
         NRZTI=AS*DUM2 
         NRZT=NRZT+NRZTI 
         MRKS=MRKS+NRZTI*ETB(J1) 
         MRET=MRET-NRZTI*KSB(J1) 
         R(1,1)=R(1,1)+DUM1*ETB(J1) 
         R(1,2)=R(1,2)+DUM2*ETB(J1) 
         R(2,1)=R(2,1)-DUM1*KSB(J1) 
         R(2,2)=R(2,2)-DUM2*KSB(J1) 
         R(3,1)=R(3,1)+DUM1 
         R(3,2)=R(3,2)+DUM2 
   60 CONTINUE 
      IF(ABS(EPSS-EPSI).LE.1E-10) THEN 
         IF(EPSS.GE.0) GO TO 100 
         DO 70 J1=1,NRC 
            FCD=FC(J1) 
            IF(J1.EQ.1) THEN 
               NP1=1 
            ELSE 
               NP1=IL(J1-1) 
            ENDIF 
            NP2=IL(J1)-1 
            CALL CENTRA(NP1,NP2,FCD,B,C,EPSS,KSP,ETP,NRZT,MRKS,MRET, 
     *         BLX,CLX,R,EPSC2(J1),A1(J1),A2(J1)) 
   70    CONTINUE 
      ELSE 
         IF(EPSS.GE.0.AND.EPSI.GE.0) GO TO 100 
         ET01=-C/B 
         DO 80 J1=1,NRC 
            ET12=(-EPSC2(J1)-C)/B 
            FCD=FC(J1) 
            IF(J1.EQ.1) THEN 
               NP1=1 
            ELSE 
               NP1=IL(J1-1) 
            ENDIF 
            NP2=IL(J1)-1 
            DO 80 J2=NP1,NP2 
               EPS0=EPSP(J2) 
               EPS1=EPSP(J2+1) 
               IF(EPS0.EQ.EPS1) GO TO 80 
               IF(EPS0.GE.0.AND.EPS1.GE.0) GO TO 80 
               CALL DIFER(J2,ET01,ET12,KSP,ETP,EPS0,EPS1,KS1I,ET1I,KS2I, 
     *            ET2I,KS1II,ET1II,KS2II,ET2II,EPSC2(J1)) 
               CALL REGI(FCD,B,C,KS1I,ET1I,KS2I,ET2I,NRZT,MRKS,MRET, 
     *            BLX,CLX,R,A1(J1),A2(J1)) 
               CALL REGII(FCD,KS1II,ET1II,KS2II,ET2II,NRZT,MRKS,MRET) 
   80    CONTINUE 
      ENDIF 
  100 NRZ=NRZT 
      MRX=MRKS*CA-MRET*SA 
      MRY=MRKS*SA+MRET*CA 
      RETURN 
      END 
 
 
      SUBROUTINE ACO(E,EPSB,FYD,SIG,ET) 
      REAL*8 EPSB,FYD 
      REAL*8 A,B,C,DUM1,E,EPS1,EPS2,ET,SIG 
      EPS2=FYD/E 
      IF(ABS(EPSB).LE.EPS2) THEN 
         SIG=E*EPSB 
         ET=E 
      ELSE 
         SIG=SIGN(FYD,EPSB) 
         ET=0 
      ENDIF 
      RETURN 
      END 
 
 
      SUBROUTINE DIFER(I,ET01,ET12,KSP,ETP,EPS0,EPS1,KS1I,ET1I,KS2I, 
     *   ET2I,KS1II,ET1II,KS2II,ET2II,EPSC2) 
      INTEGER T01,T12 
      REAL*8 ETP(*),KSP(*),KS1I,KS2I,KS1II,KS2II,KS01,KS12 
      REAL*8 DET,DET01,DET12,DKSDET,DUM1,DUM2,EPS0,EPS1,ET01,ET12,ET1I, 
     * ET1II,ET2I,ET2II,EPSC2 
      T01=0 
      T12=0 
      KS1I=0 
      ET1I=0 
      KS2I=0 
      ET2I=0 
      KS1II=0 
      ET1II=0 
      KS2II=0 
      ET2II=0 
      I2=I+1 
      DET=ETP(I2)-ETP(I) 
      DKSDET=(KSP(I2)-KSP(I))/DET 
      DUM1=ET01-ETP(I) 
      DUM2=ET12-ETP(I) 
	  KS01=KSP(I)+DUM1*DKSDET 
      KS12=KSP(I)+DUM2*DKSDET 
      DET01=DUM1/DET 
      DET12=DUM2/DET 
      IF(DET01.GT.0.AND.DET01.LT.1) T01=1 
      IF(DET12.GT.0.AND.DET12.LT.1) T12=1 
      IF(EPS0.LT.EPS1) THEN 
         T01=-T01 
         T12=-T12 
      ENDIF 
      IF(T01.EQ.0.AND.T12.EQ.0) THEN 
         IF(EPS0.LT.0) THEN 
            IF(EPS0.GT.-EPSC2) THEN 
               KS1I=KSP(I) 
               ET1I=ETP(I) 
               KS2I=KSP(I2) 
               ET2I=ETP(I2) 
            ELSE 
               KS1II=KSP(I) 
               ET1II=ETP(I) 
               KS2II=KSP(I2) 
               ET2II=ETP(I2) 
            ENDIF 
         ENDIF 
      ELSE 
         IF(T01.EQ.1) THEN 
            KS1I=KS01 
            ET1I=ET01 
            IF(T12.EQ.1) THEN 
               KS2I=KS12 
               ET2I=ET12 
               KS1II=KS12 
               ET1II=ET12 
               KS2II=KSP(I2) 
               ET2II=ETP(I2) 
            ELSE 
               KS2I=KSP(I2) 
               ET2I=ETP(I2) 
            ENDIF 
         ELSE 
            IF(T01.EQ.-1) THEN 
               KS2I=KS01 
               ET2I=ET01 
               IF(T12.EQ.-1) THEN 
                  KS1I=KS12 
                  ET1I=ET12 
                  KS2II=KS12 
                  ET2II=ET12 
                  KS1II=KSP(I) 
                  ET1II=ETP(I) 
               ELSE 
                  KS1I=KSP(I) 
                  ET1I=ETP(I) 
               ENDIF 
            ELSE 
               IF(T12.EQ.1) THEN 
                  KS1I=KSP(I) 
                  ET1I=ETP(I) 
                  KS2I=KS12 
                  ET2I=ET12 
                  KS1II=KS12 
                  ET1II=ET12 
                  KS2II=KSP(I2) 
                  ET2II=ETP(I2) 
               ELSE 
                  KS1I=KS12 
                  ET1I=ET12 
                  KS2I=KSP(I2) 
                  ET2I=ETP(I2) 
                  KS1II=KSP(I) 
                  ET1II=ETP(I) 
                  KS2II=KS12 
                  ET2II=ET12 
               ENDIF 
            ENDIF 
         ENDIF 
      ENDIF 
      RETURN 
      END 
 
 
      SUBROUTINE CENTRA(NP1,NP2,FCD,B,C,EPSS,KSP,ETP,NRZT,MRKS,MRET, 
     *   BLX,CLX,R,EPSC2,A1,A2) 
      REAL*8 A1,A2 
      REAL*8 KSP(*),ETP(*),R(3,2),NRZT,MRKS,MRET 
      REAL*8 B,BLX,C,CLX,EPSS,FCD,EPSC2 
      IF(EPSS.GE.0) RETURN 
      IF(EPSS.LT.-EPSC2) GO TO 20 
      DO 10 J1=NP1,NP2 
         J2=J1+1 
         CALL REGI(FCD,B,C,KSP(J1),ETP(J1),KSP(J2),ETP(J2),NRZT,MRKS, 
     *      MRET,BLX,CLX,R,A1,A2) 
   10 CONTINUE 
      RETURN 
   20 DO 30 J1=NP1,NP2 
         J2=J1+1 
         CALL REGII(FCD,KSP(J1),ETP(J1),KSP(J2),ETP(J2),NRZT,MRKS,MRET) 
   30 CONTINUE 
      RETURN 
      END 
 
 
      SUBROUTINE REGI(FCD,B,C,KS1,ET1,KS2,ET2,NRZT,MRKS,MRET,BLX,CLX,R 
     *,A1,A2) 
      REAL*8 A1,A2 
      REAL*8 R(3,2),KS1,KS2,NRZT,MRKS,MRET 
      REAL*8 B,BLX,BLE,BR,C,CLX,CLE,D0,D1,D2,DET,DET1,DET2,DET3,DKS, 
     * DKS2,E0,E1,E2,ET1,ET2,FCD,G00,G01,G02,G03,G10,G11,G12 
      IF(KS1.EQ.0.AND.ET1.EQ.0.AND.KS2.EQ.0.AND.ET2.EQ.0) RETURN 
      BLE=2.D0*A2*B 
      CLE=2.D0*A2*C+A1 
      D0=C*A1+A2*C*C 
      D1=B*CLE 
      D2=A2*B*B 
      E0=CLE*CLX 
      E1=BLE*CLX+CLE*BLX 
      E2=BLE*BLX 
      BR=.85*FCD 
      DKS=KS2-KS1 
      DET=ET2-ET1 
      DET1=DET/2. 
      DET2=DET*DET 
      DET3=DET2*DET 
      DKS2=DKS*DKS 
      G00=(KS1+DKS/2.)*DET 
      G01=(KS1*(ET1+DET1)+DKS*(ET1/2.+DET/3.))*DET 
      G02=(KS1*(ET1*(DET+ET1)+DET2/3.)+DKS*(ET1*(ET1/2+DET/1.5)+ 
     *   DET2/4.))*DET 
      G03=(KS1*(ET1*(DET2+ET1*(1.5*DET+ET1))+DET3/4.)+DKS*(ET1* 
     *   (0.75*DET2+ET1*(DET+ET1/2.))+DET3/5.))*DET 
      G10=(KS1*(KS1+DKS)+DKS2/3)*DET1 
      G11=(KS1*(KS1*(ET1+DET1)+DKS*(ET1+DET/1.5))+DKS2*(ET1/3.+DET/4.))* 
     *   DET1 
      G12=(KS1*(KS1*(ET1*(ET1+DET)+DET2/3.)+DKS*(ET1*(ET1+DET/0.75)+ 
     *   DET2/2.))+DKS2*(ET1*(ET1/3.+DET1)+DET2/5.))*DET1 
      NRZT=NRZT+BR*(D0*G00+D1*G01+D2*G02) 
      MRKS=MRKS+BR*(D0*G01+D1*G02+D2*G03) 
      MRET=MRET-BR*(D0*G10+D1*G11+D2*G12) 
      R(1,1)=R(1,1)+BR*(E0*G01+E1*G02+E2*G03) 
      R(2,1)=R(2,1)-BR*(E0*G10+E1*G11+E2*G12) 
      R(3,1)=R(3,1)+BR*(E0*G00+E1*G01+E2*G02) 
      RETURN 
      END 
 
 
      SUBROUTINE REGII(FCD,KS1,ET1,KS2,ET2,NRZT,MRKS,MRET) 
      REAL*8 KS1,KS2,NRZT,MRKS,MRET 
      REAL*8 DET,DKS,ET1,ET2,FC,FCD,G00,G01,G10 
      IF(KS1.EQ.0.AND.ET1.EQ.0.AND.KS2.EQ.0.AND.ET2.EQ.0) RETURN 
      DKS=KS2-KS1 
      DET=ET2-ET1 
      G00=(KS1+DKS/2.)*DET 
      G01=(KS1*(ET1+DET/2.)+DKS*(ET1/2.+DET/3.))*DET 
      G10=(KS1*(KS1+DKS)+DKS*DKS/3.)*DET/2. 
      FC=0.85*FCD 
      NRZT=NRZT-FC*G00 
      MRKS=MRKS-FC*G01 
      MRET=MRET+FC*G10 
      RETURN 
      END 
 
 
 
 
 