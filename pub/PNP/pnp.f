C-------------------Beginning of PNP code-------------------------------
C-----------------------------------------------------------------------
C  PROGRAM TO SOLVE PNP EQUATIONS & THEIR APPLICATION IN ION CHANNELS
C 
C  Author : Duan P. Chen
C           dchen@rpslmc.edu
C           Department of Biophysics and Physiology
C                      & Physiology                          
C           Rush Medical College
C           Rush Presbyterian St.-Luke's Medical Ctr.
C           1750 W. Harrison Street., JS1298
C           Chicago, IL 60612, USA
C
C  METHOD : GUMMEL NON-LINEAR BLOCK ITERATION
C  
C  NOTE:    This program reads in the protein charge distribution 
C           from a file, which is specified by the user in PCFILE.
C           The format of PCFILE is two columns of real numbers:
C             the first column is x-along the channel, scaled by
C             the length of the channel;
C             the second column is the channel protein charge distribution.
C             The x-coordinate of the first row must be 0.0, and 
C                 that of the last row must be 1.0.  Furthremore, 
C                 all mesh points are equally spaced.
C               
C            
C           Inputs supply VOLT1, VOLT2, VSTEP, and the program
C            calculates NVSTEP.
C           It calculates the potential, and concentration profiles 
C                             if NVSTEP =1; 
C           It calculates the IV-curves, if NVSTEP >1.
C
C          MODPSE :  1, The modified Poisson Equation is used
C                    Barcilon, V., D.~P. Chen, and R.~S. Eisenberg 1992
C                    SIAM J. Appl. Math. 52:1405--1425;
C                    else Poisson Equation is used.
C          N      :  (N+2) total mesh points, UNIFORM
C                    mesh must be uniform
C          NION   :  # OF KINDS OF IONS IN THE BATHES
C          DION   :  DIFFUSION CONSTANTS in CM**2/SEC.
C          CL, CR :  ION CONCENTRATION IN THE BATHES, LEFT AND RIGHT
C                    in units of MOLEs/LITER
C          EA     :  DIELECTRIC CONSTANT FOR CHANNEL(PORE)
C          EM     :  DIELECTRIC CONSTANT FOR MEMBRANE  
C                         in units of E0
C          RADIUS :  RADIUS OF THE CHANNEL IN ANGSTROMS
C          WIDTH  :  WIDTH  OF THE CHANNEL IN ANGSTROMS
C          VOLT1  :  BEGINNING VOLTAGE (INPUT IN VOLTS)
C          VOLT2  :  FINAL VOLTAGE     (INPUT IN VOLTS)
C          VSTEP  :  VOLTAGE STEP      (INPUT IN VOLTS)
C-----------------------------------------------------------------------
      PARAMETER(MION=3,MXGRID=9002,MAXIV=500)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(0:MXGRID),  PCT(0:MXGRID)
      DIMENSION U(0:MXGRID), V(0:MXGRID), P(0:MXGRID)
      DIMENSION WA(MXGRID),  WB(MXGRID),  WC(MXGRID),  WR(MXGRID)
      DIMENSION CL(MION), CR(MION), DION(MION), JDATE(3)
      DOUBLE PRECISION LAMDA2, TOL
      INTEGER N, MODPSE , ITCPU
      LOGICAL   FIRTIME
      CHARACTER IVFILE*80
      CHARACTER PFILE*80,  CFILE*80, PCFILE*80
      EXTERNAL  SEMIBC, SEMISOL , NPSOLJEX
      DATA TKB /0.025875D0 /
      DATA TKBmV /0.025875D3 /
      DATA NION, Z1, Z2 /2, 1.d0, -1.d0/
C----------------------------------------------------------------------
C     ITCPU = MCLOCK(0)
C----------------------------------------------------------------------
C
C   READ IN THE NECESSARY PARAMETERS: input
C
C----------------------------------------------------------------------
C  geometry and dielectric constants
C----------------------------------------------------------------------
      READ(5,*) MODPSE, RADIUS, WIDTH , EM , EA 
C
C  diffusion constants
C
      READ(5,*) ( DION(I) , I = 1 , NION )
C
C  intrinsic concentration 
C
      READ(5,*) ( CL(I) , I = 1 , NION) 
      READ(5,*) ( CR(I) , I = 1 , NION) 
C
C  applied voltage and steps 
C
      READ(5,*) VOLT1, VOLT2, VSTEP
C
      READ(5,*) TOL, MAXITER
C
C  read in number of mesh-points and
C     protein charge distribution filename
C
      READ(5,*)  N, PCFILE
      N1     =  N + 1
C
C  read in output file names
C  if VOLT1 = VOLT2 then potential and concentraion are written
C   to PFILE and CFILE 
C  if VOLT2 > VOLT1 then IV curve is written to IVFILE
C
      READ(5,*) IVFILE, PFILE, CFILE
C
C reads in  protein charge distribution 
C
      OPEN(UNIT=4, FILE=PCFILE,STATUS='UNKNOWN')
      DO I = 0, N1
       READ(4,*) X(I), PCT(I)
      ENDDO
      CLOSE(UNIT=4)
      XSTEP = X(1)
C-----------------------------------------------------------------------
      ALPHA = RADIUS/WIDTH
      ER    = EM/EA
      EP    = ER/( -ALPHA*ALPHA*DLOG(ALPHA) )
      TKBMV1 = 1.D0 / TKBMV
      TKB1 = 1.D0 / TKB
C-----------------------------------------------------------------------
C     CALL GETHOST(HOSTNAME)
C     CALL IDATE(JDATE)
C     WRITE(6,1) JDATE, HOSTNAME
C  1    FORMAT(1x,' Solve PNP equations on ',3(i4,1h/), ' on',5x,a20 )
      WRITE(6,1) 
1     FORMAT(1x,' Solve PNP equations')
C
      IF( MODPSE.EQ.1 ) then
      WRITE(6,*) '  INCLUDING induced charge term'
      ELSE
      WRITE(6,*) '  NOT INCLUDING induced charge term'
      ENDIF
C-----------------------------------------------------------------------
      WRITE(6,3) NION, RADIUS, WIDTH
3     FORMAT(//1X,'NION, RADIUS, WIDTH=',I5,1p,2E15.5,' angstroms')
      WRITE(6,5) EM, EA
5     FORMAT(1X,'EM, EA =', 2F10.5)
      WRITE(6,13) ( DION(I) , I = 1 , NION )
13    FORMAT(/1X,'THE DIFFUSION CONSTANT IN cm*cm/s are:'
     + ,1X,1P,5E10.2)
15    format(/1x,'IV_CURVE FOR VOLT1, VOLT2, VSTEP, NVSTEP = '/1P,
     & 3E12.2,i7)
C-----------------------------------------------------------------------
C  set up NVSTEP value
C-----------------------------------------------------------------------
        NVSTEP = INT( (VOLT2 - VOLT1)/VSTEP ) 
        TEMP   = VOLT2 - VOLT1- DFLOAT(NVSTEP-1)*VSTEP
        TEMP   = DABS( TEMP )
c       DO 160 WHILE ( TEMP.GT.1.0D-8 )
c       DO WHILE ( TEMP.GT.1.0D-8 )
        NVSTEP = NVSTEP + 1
c       TEMP   = VOLT1 + DFLOAT(NVSTEP-1) * VSTEP
c       TEMP   = DABS( VOLT2 - TEMP )
c       END DO
        WRITE(6,15) VOLT1, VOLT2, VSTEP , NVSTEP
C-----------------------------------------------------------------------
      WRITE(6,19) N, TOL
19    FORMAT(/1x,'Mesh point exclude two end points, TOL',I9,1p,E12.2 )
C-----------------------------------------------------------------------
      WRITE(6,9) (   CL(I) , I = 1 , NION )
9     FORMAT(/1X,' LEFT-INSIDE  ION CONCENTRATION IN mM/L ='
     &,3P,3F12.3)
      WRITE(6,11) (   CR(I) , I = 1 , NION )
11    FORMAT(1X,'RIGHT-OUTSIDE ION CONCENTRATION IN mM/L ='
     +,3P,3F12.3/)
C----------------------------------------------------------------------
C     get the scaling parameters
C----------------------------------------------------------------------
      VLI = CL(1)
      ULI = CL(2)
      VRI = CR(1)
      URI = CR(2)
      CSCALE= DMAX1(ULI, URI, VLI, VRI)
      IF (CSCALE.LT.1.D-10) THEN
          CHARGEL= PCT(0)
          CHARGER= PCT(N1)
          CSCALE= DMAX1(CHARGEL,CHARGER) 
      ENDIF
      IF (CSCALE.LT.1.D-10) CSCALE=1.D0 
      CSCALE1 = 1.D0/CSCALE
      CSCALEM = 1.D3*CSCALE
      YMDA  = 4.206D0 * CSCALE * WIDTH * WIDTH  / EA
      LAMDA2= 1.D0/YMDA
      EP2   = 2.0D0*EP*LAMDA2
      W     = 3.03115182D+6 * CSCALE * RADIUS * RADIUS /WIDTH
C----------------------------------------------------------------------
      DO I = 0, N1
      PCT(I) = CSCALE1 * PCT(I)
      ENDDO
      UCI =  CSCALE1*ULI 
      VCI =  CSCALE1*VLI 
      CHARGEL= PCT(0)
      CALL SEMIBC(UCI, VCI, CHARGEL, PL, UL, VL )
      UCI =  CSCALE1*URI
      VCI =  CSCALE1*VRI
      CHARGER= PCT(N1)
      CALL SEMIBC(UCI, VCI, CHARGER, PR, UR, VR )
C
      WRITE(6,147) UL, VL, PL, CHARGEL, PL*TKBMV 
     &            ,UR, VR, PR, CHARGER, PR*TKBMV 
c     WRITE(6,148) UL*CSCALE, VL*CSCALE, UR*CSCALE, VR*CSCALE
147    FORMAT(/1X,'Modified bc according to Donnan potential are'/1x,1p,
     &'UL, VL, PL, CHARGEL, PL in mV'/5e12.3/
     & 1x, 'UR, VR, PR, CHARGER, PR in mV'/5e12.3/)
148    FORMAT(1X,'Effective concentrations are:'/1x, 'left  ',1p,2e12.5
     & /1x, 'right ',2e12.5, ' M/liter negative positive '/)
C-----------------------------------------------------------------------
C
C   BELOW CALCULATE THE I-V CURVE
C
C----------------------------------------------------------------------
299   DELTA  = VOLT1*TKB1
      DSTEP  = VSTEP*TKB1
      OPEN(UNIT=4, FILE=IVFILE,STATUS='UNKNOWN')
      WRITE(6,*)'Below lists: Vappl (mV), Current (pA), Conductance(pS)
     & (I/V), FLUX1, FLUX2 (ions/sec)'
C----------------------------------------------------------------------
      DO 500 IV = 1, NVSTEP	
C
C initial guess
C
      U(0)  = UL  
      U(N1) = UR  
      V(0)  = VL
      V(N1) = VR
      VA    =  PL - PR + DELTA 
      P(0)  =  VA
      P(N1) =  0.D0
      PSTEP =  VA /DFLOAT(N1)
      DO I = 1, N
         P(I) = DFLOAT(N1-I) * PSTEP
      ENDDO 
C----------------------------------------------------------------------
c
c  modpse = 1 the modified Poisson eq is used
c
C----------------------------------------------------------------------
400   WA(1) = EP2
      IF ( MODPSE.NE.1 ) WA(1)=0.D0
      CALL SEMISOL(N, X, PCT, LAMDA2, TOL, MAXITER, U, V, P,
     &      WA, WB, WC, WR )
      ITER1 = INT(WR(2))
      IF ( IABS(ITER1 - MAXITER ).LT.1 ) THEN
       IF ( FIRTIME ) THEN
        WRITE(6,*) 'FIRST TIME WARNING '
        WRITE(6,201)
        DO I = 1, N
            P(I) = 0.D0
        ENDDO
        FIRTIME = .FALSE.
        GOTO 400
       ENDIF
      WRITE(6,201) 
      WRITE(6,415) 
      WRITE(6,202) 
      ENDIF
201   FORMAT(1x,'WARNING WARNING WARNING --- NOT CONVERGING')
202   FORMAT(1x,'Increase the number of mesh points or reduce tol')
415   FORMAT(1x,'Error in semisol is: ',1p,e13.4,' iter= ',f5.0)
C-----------------------------------------------------------------------
      CALL NPSOLJEX(N, XSTEP, VL, VR, Z1, P, FLUX1)
      CALL NPSOLJEX(N, XSTEP, UL, UR, Z2, P, FLUX2)
      FLUX1 = W*FLUX1*DION(1)
      FLUX2 = W*FLUX2*DION(2)
      TFLUX = FLUX1*Z1 + FLUX2*Z2
      IF ( DABS(DELTA).GT.1.D-10 ) THEN
         GAMA  = TFLUX*TKB1/DELTA
      ENDIF
      FLUX1 = 6.25D6*FLUX1
      FLUX2 = 6.25D6*FLUX2
C-----------------------------------------------------------------------
C
C conductance in pS
C
C-----------------------------------------------------------------------
      WRITE(6,555) DELTA*TKBMV, TFLUX, GAMA, FLUX1, FLUX2
      WRITE(4,555) DELTA*TKBMV, TFLUX, GAMA, FLUX1, FLUX2
C-----------------------------------------------------------------------
C  store total current and voltage for cal reversal
C----------------------------------------------------------------------
      DELTA  = DELTA + DSTEP	
500   CONTINUE		
      CLOSE(UNIT=4)
C----------------------------------------------------------------------
C----------------------------------------------------------------------
      IF (NVSTEP.EQ.1) THEN
        WRITE(6,412)  PFILE , CFILE
        OPEN(UNIT=3, FILE=PFILE,STATUS='UNKNOWN')
        OPEN(UNIT=4, FILE=CFILE,STATUS='UNKNOWN')
        WRITE(3,552) ( X(I)*WIDTH, P(I)*TKBMV, I = 0, N1)
        WRITE(4,553) ( X(I)*WIDTH, V(I)*CSCALE, U(I)*CSCALE,I = 0, N1)
      ELSE 
        WRITE(6,155)  IVFILE
      ENDIF
C----------------------------------------------------------------------
C----------------------------------------------------------------------
155   FORMAT(/'IV_curve written to FILE: ',a15, 'in the above format',
     & ' and units.')
412   FORMAT(/'x vs Potential in file:     ',a10, ' Angstroms vs mV' /
     &        'x vs Concentration in file: ',a10, ' Angstroms vs M/l')
552   format(1x,1p,2e13.4)
553   format(1x,1p,3e13.4)
555   format(1x,1p,5e13.4)
C----------------------------------------------------------------------
C     ITCPU = ( MCLOCK(0) - ITCPU )/100
C     WRITE(6,601) ITCPU
C 601   FORMAT(/1x, ' Total cpu in seconds ', i9 )
      WRITE(6,601) 
601   FORMAT(/1x, ' Finished' )
      STOP
      END
C----------------------------------------------------------------------
C----------------------------------------------------------------------
C***********************************************************************
C
C  to solve for the POSITIVE ION concentration 
C  according to Nernst Plack eq if potential is know
C
C  input :   Mesh equal spacing xstep 
C            P(I), I = 0, 1, ... , N+1
C            VL, VR stored in
C               V(0), V(N+1) respectively
C
C  output:   V(I), I = 0, 1, ... , N+1
C            V IS CONCENTRATION OF POSITIVE ION
C
C***********************************************************************
      SUBROUTINE NPSVEX(N, V, P, WR)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION V(0:*), P(0:*), wr(*)
      INTRINSIC DSQRT, DLOG, DEXP
C
         n1   = n + 1
         d    = dexp(p(1))
         wr(1) =  dexp(p(0)) + d   
         do i = 2, n
         b  =  d
         d    = dexp(p(i))
         wr(i) = wr(i-1) + ( b + d )
         enddo
         vd   = wr(n) +  d + dexp(p(n1))  
         vd1  =  1.d0/vd
c
         do  i = 1, n
         a     = dexp(p(0)-p(i))
         b     = dexp(p(n1)-p(i))
         v(i)  = vd1 * ( v(0) * (vd-wr(i)) * a  + v(n1) * wr(i) * b )
         enddo
      RETURN
      END
C***********************************************************************
C
C  to solve for the NEGATIVE ION concentration 
C  according to Nernst Plack eq if potential is know
C
C  input :   Mesh equal spacing xstep 
C            P(I), I = 0, 1, ... , N+1
C            UL, UR stored in
C               U(0), u(N+1) respectively
c
C  output:   U(I), I = 0, 1, ... , N+1
C            U IS CONCENTRATION OF NEGATIVE ION
C
C***********************************************************************
      SUBROUTINE NPSUEX(N, V, P, WR)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION V(0:*), P(0:*), wr(*)
      INTRINSIC DSQRT, DLOG, DEXP
C
         n1   = n + 1
         d    = dexp(-p(1))
         wr(1) =  dexp(-p(0)) + d   
         do i = 2, n
         b  =  d
         d    = dexp(-p(i))
         wr(i) = wr(i-1) + ( b + d )
         enddo
         vd   = wr(n) +  d + dexp(-p(n1))  
         vd1  =  1.d0/vd
c
         do  i = 1, n
         a     = dexp(-p(0)+p(i))
         b     = dexp(-p(n1)+p(i))
         v(i)  = vd1 * ( v(0) * (vd-wr(i)) * a  + v(n1) * wr(i) * b )
         enddo
      RETURN
      END
C***********************************************************************
C
C  to solve for the FLUX according to Nernst Plack eq
C     if potential AND concentration are known
C
C  input :   Mesh point x(i), and
C            P(I), I = 0, 1, ... , N+1
C            vl, vr 
c            IZ  , the valence 1 or -1
C
C  output:   flux 
C
C***********************************************************************
      SUBROUTINE NPSOLJEX(N, xstep, VL, VR, z, P, FLUX)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION P(0:*)
      DOUBLE PRECISION vl, vr, z, flux, XSTEP
      INTRINSIC   DEXP
C
C     solve for flux
C
         flux = 0.d0
         d    = dexp(z*p(0))
         do i = 0, n
         b = d
         d = dexp(z*p(i+1))
         flux = flux + ( b + d )
         enddo
         FLUX = FLUX * XSTEP
         FLUX = 2.d0*( VL*DEXP(Z*P(0)) - VR*DEXP(Z*P(N+1)) )/FLUX
c
      RETURN
      END
C***********************************************************************
C   eq form is
c   
c  - lamda2 * phi " = p - n + N
C
C  to solve the semiconductor eq by non-linear block iteration method
c   
C***********************************************************************
c  Nernst-Plack eq. is integrated exactly by the integral eq.
C***********************************************************************
c
C
C  input :   Mesh point x(i), i = 1, N
C            tolerance   tol
C            UL, UR, VL, VR, PL, PR stored in
C               U(0), U(N+1), V(0), V(N+1) , P(0), P(N+1) respectively
C            maximum # of iterations maxiter
C            lamda2 is lamda * lamda
C
C  output:   U(I), V(I), P(I), I = 0, 1, ... , N+1
C            U IS CONCENTRATION OF NEGATIVE ION
C            V IS CONCENTRATION OF POSITIVE ION
C            P IS THE ELECTRICAL POTENTIAL IN UNIT OF e/kT
C
C***********************************************************************
      SUBROUTINE SEMISOL(N, x, PCT, LAMDA, TOL, MAXITER, U, V, P,
     &   WA, WB, WC, WR )
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION U(0:*), V(0:*), P(0:*), pct(0:*), x(0:*)
      DIMENSION wa(*), wb(*), wc(*), wr(*)
      DOUBLE PRECISION LAMDA2, LAMDA
      INTRINSIC DSQRT, DLOG, DEXP
      EXTERNAL  F, MYTRIDAG, NPSVEX , NPSUEX
C
C   EXTRACTE EP2 FROM WA
C
      EP2 = WA(1)
      DELTA = p(0) - p(n+1)
      XSTEP1= 1.d0/X(1)
      LAMDA2 = (LAMDA*XSTEP1)*XSTEP1
      DO 900 ITER = 1, MAXITER
c
c  solve for u and v , z(u)=-1 , z(v)=1
C
            CALL NPSVEX(N, V, P, WR)
            CALL NPSUEX(N, U, P, WR)
C
C     use updated u & v to solve new p
C     alpha = rh(i) * DP(i),  denoted  a2
C          alpha is definded in the document of discretization
C     beta  = rh(i-1) * dp , denoted a3
C solve for p
C
      DU2  = U(1) - U(0)
      DV2  = V(1) - V(0)
      DP2  = P(1) - P(0)
      FDP2 = F(DP2)
      FDP2N= F(-DP2)
      DO 300 I = 1, N
      WA(I)  =  LAMDA2 
      WB(I)  = - ( 2.D0*LAMDA2 + U(I) + V(I) ) 
      WC(I)  =  LAMDA2 
C
C   redefine a2 , a3
C
      DU = DU2
      DV = DV2
      DP = DP2
      FDP= FDP2
      FDPN= FDP2N
      DU2  = U(I+1) - U(I)
      DV2  = V(I+1) - V(I)
      DP2  = P(I+1) - P(I)
      FDP2 = F(DP2)
      FDP2N= F(-DP2)
c
      A4 = x(i)
      wr(i) = - ( u(i) + v(i) ) * p(i) + u(i) - v(i) - pct(i)
      wr(i) = wr(i) + du2*fdp2 - du*fdpn - dv2*fdp2n + dv*fdp
      wr(i) = wr(i) + ep2 * ( p(i) - ( 1.d0 - a4 ) * delta )
300   continue
      wr(N) = wr(N) - wc(n) * p(N+1)
      wr(1) = wr(1) - wa(1) * p(0)
c
c
c    wa subdiagonal
c    wb diagonal
c    wc superdiagonal
c    wr right-hand side, output is solution
c
c
c     CALL DGTSL(N, WA, WB, WC, WR, IEER )
      CALL MYTRIDAG(N, WA, WB, WC, WR, IEER)
        IF ( IEER .NE. 0 ) THEN
        WRITE(6,311) ITER, IEER
311     FORMAT(1X,' LINPACK BOBMED PROGRAM STOPPED '/
     &       1X,' IN SOLVING new P AT ITER=', I9, '  BOBMED AT ROW',I4)
        STOP
        ENDIF
c
c  CALCULATE ERROR AND COMPARED TO TOL & CHECK CONVERGENCE
c
      ERR = 0.D0
      sca = 0.d0
      DO 400 I = 1, N
      ERR = ERR + ( P(I) - WR(I) )*( P(I) - WR(I) )
c
C  update P
c
      P(I) = WR(I)
      sca  = sca + dabs(p(i))
400   CONTINUE
c
      if ( sca.gt.dfloat(n) ) then
        err = err/sca 
      else
        err = err/dfloat(n)
      endif
c
      IF ( ERR.LT.TOL ) THEN
      WR(1) = ERR
      WR(2) = DFLOAT(ITER)
      GOTO 999
      ENDIF
900   CONTINUE
      WR(1) = ERR
      WR(2) = DFLOAT(maxITER)
999   RETURN
      END
c
c
c
C     (1.0d0 - dexp(x4)).eq.1.0d0  x4 = -100.0
C     (dexp(x1)-1.0d0).eq.-1.0d0   x1 = -100.0
C     dexp(-x5).ne.0.0d0           x5 = 1280.0
      FUNCTION F(X)
      DOUBLE PRECISION F, X, D, DEXP
      INTRINSIC DEXP
      IF ( DABS(X).LT.1.0D-4 ) THEN
         F = 1.0D0 + X*0.25d0 + X*X*0.05D0
         F = F/(1.0D0+0.5D0*X+X*X*0.16666666666667D0)
         F = F *  0.16666666666667D0
      ELSE IF ( X.GT.400.0D0 ) THEN
         F = 1.0D0/(X*X)
      ELSE
         D = DEXP(X)
         F = D - 1.0D0 - X - 0.5D0 * X * X
         F = F/(X*X)
         F = F/(D-1.0D0)
      ENDIF
      RETURN
      END
C
      SUBROUTINE MYTRIDAG(N,A,B,C,R,IEER)
C
C  N IS ORDER OF MATRIX
C  A IS SUB-DIAGONAL   A(2) through A(N)
C  B IS DIAGONAL       B(1) through B(N)
C  C IS SUPER-DIAGONAL C(1) through C(N-1)
C  R IS RIGHT HAND SIDE, the solusion on output
C  IEER IS ERROR REPORT IF NON-ZERO
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(N),B(N),C(N),R(N)
      INTEGER J1X, J2X
      IEER = 0
      IF(B(1).EQ.0.) THEN
      IEER = 1
      GOTO 13
      ENDIF
      BET=B(1)
      R(1)=R(1)/BET
      DO 11 J=2,N
        C(J-1)=C(J-1)/BET
        BET=B(J)-A(J)*C(J-1)
        IF(BET.EQ.0.)  THEN
        IEER = J
        GOTO 13
        ENDIF
        R(J)=(R(J)-A(J)*R(J-1))/BET
11    CONTINUE
      J2X = N - 1
C     J1X = IAND(MAX0(N - 1,0),1)
      J1X = mod(MAX0(N - 1,0),2)
      IF (J1X .EQ. 1) R(J2X) = R(J2X) - C(J2X)*R(J2X+1)
      DO 12 J = N - J1X - 1, 1, -2
         R(J) = R(J) - C(J)*R(J+1)
         R(J-1) = R(J-1) - C(J-1)*R(J)
   12 CONTINUE
13    RETURN
      END
C-----------------------------------------------------------------------
C  
C  To calculate the modified boundary condition according to
C  1. flux zero;
C  2. zero net charge;
C  input :   CI, CHARGE
C  output:   UP, VP, P
C            UP IS CONCENTRATION OF NEGATIVE ION
C            P  IS THE ELECTRICAL POTENTIAL IN UNIT OF e/kT
C
      SUBROUTINE SEMIBC(CIu, CIv, CHARGE, P, UP, VP )
      DOUBLE PRECISION CIu, CIv, CHARGE, P, UP, VP, DSQRT, DLOG
      INTRINSIC DSQRT, DLOG
      IF ( DABS(CHARGE).LT.1.D-10 ) THEN
        UP = CIU
        VP = CIV
        P  = 0.D0
      ELSE
        P =  DSQRT ( CHARGE*CHARGE + 4.0D0*CIu*CIv )
        UP   = 0.5D0*( P + CHARGE )
        VP   = 0.5D0*( P - CHARGE )
        P = UP / CIu
        P    = DLOG (P)
      ENDIF
      RETURN
      END
C-----------------End of PNP code---------------------------------------
