       SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     & RPL,DDSDDT,DRPLDE,DRPLDT,
     & STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     & NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     & CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C
       INCLUDE 'ABA_PARAM.INC'
C
       CHARACTER*80 CMNAME
C      
       DIMENSION STRESS(NTENS),STATEV(NSTATV),
     & DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     & STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     & PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
C
       DIMENSION n_cosines(42,3),n_weight(42),epsel_v(6),epsE_v(6),
     & epsin_v(6),energy(5),rho_o(42),epsel_o(3,3),epsE_o(3,3),
     & epsin_o(3,3),deps(3,3),eps_o(3,3),sigma_o(3,3),rho_ini(42),
     & zero_T(3,3),matDz0(3,3,3,3),matS0(3,3,3,3),matDz(3,3,3,3),
     & matS(3,3,3,3),dsig(3,3),sigmaTri(3,3),Yd(42),
     & epsin(3,3),epsE(3,3),epsel(3,3),matDz_2(6,6),
     & Depsin(3,3),sigEner(3,3),DepsE(3,3),Depsel(3,3),drho_i(42),
     & sigmaT_i(3,3),R_i(3,3),depsin_i(3,3),ddepsin(3,3),ddrho(42),
     & depsin_f(3,3),sigmaT_f(3,3),drho_f(42),depsin_check(3,3),
     & dg_dsig(3,3),R_f(3,3),sigmaT(3,3),drho(42),rho(42)
       INTEGER flag0,flag1,activated(42),detec,A_K,activated_1(42)	
       DOUBLE PRECISION E0,nu0,N_V0,a_0,Kc,Ec,Alpha,Ko,Eo,c0,c1,
     & n_cosines,n_weight,epsel_v,epsE_v,epsin_v,energy,rho_o,
     & epsel_o,epsE_o,epsin_o,deps,eps_o,sigma_o,rho_ini,
     & zero_T,matDz0,matS0,matDz,matS,dsig,sigmaTri,Yd,
     & epsin,epsE,epsel,Depsin,sigEner,DepsE,Depsel,drho_i,
     & sigmaT_i,R_i,depsin_i,ddepsin,ddrho,
     & depsin_f,sigmaT_f,drho_f,depsin_check,
     & dg_dsig,R_f,sigmaT,drho,rho,matDz_2,T0,T1,T2
C
       COMMON/parameters/E0,nu0,N_V0,a_0,Kc,Ec,Alpha,Ko,Eo,c0,c1
       data zero,half,one,two,theta /0.0D0,0.5D0,1.0D0,2.0D0,0.5D0/
       data beta/33.2699078510/ 
C
       parameter(FTOL0=1.D-8,PI=3.1415926526D0,ITER=100,nwrite=1,
     & FTOL1=1.D-20)
C     
C      DIMENSION  T0,T1,T2
C
       T0=0.707106781186548
       T1= cos(PI/4)*cos(pi/2-beta*PI/180)
       T2= cos(beta*PI/180)
       n_cosines(1,1:3)=(/1.0D0,0.0D0,0.0D0/)
       n_cosines(2,1:3)=(/0.0D0,1,0.0D0/)
       n_cosines(3,1:3)=(/0.0D0,0.0D0,1/)
       n_cosines(4,1:3)=(/T0,T0,0.0D0/)
       n_cosines(5,1:3)=(/T0,-T0,0.0D0/)
       n_cosines(6,1:3)=(/T0,0.0D0,T0/)
       n_cosines(7,1:3)=(/T0,0.0D0,-T0/)
       n_cosines(8,1:3)=(/0.0D0,T0,T0/)
       n_cosines(9,1:3)=(/0.0D0,T0,-T0/)
       n_cosines(10,1:3)=(/T1,T1,T2/)
       n_cosines(11,1:3)=(/T1,T1,-T2/)
       n_cosines(12,1:3)=(/T1,-T1,T2/)
       n_cosines(13,1:3)=(/T1,-T1,-T2/)
       n_cosines(14,1:3)=(/T1,T2,T1/)
       n_cosines(15,1:3)=(/T1,T2,-T1/)
       n_cosines(16,1:3)=(/T1,-T2,T1/)
       n_cosines(17,1:3)=(/T1,-T2,-T1/)
       n_cosines(18,1:3)=(/T2,T1,T1/)
       n_cosines(19,1:3)=(/T2,T1,-T1/)
       n_cosines(20,1:3)=(/T2,-T1,T1/)
       n_cosines(21,1:3)=(/T2,-T1,-T1/) 
       n_cosines(22:42,1:3)=-n_cosines(1:21,1:3)
       n_weight(1:3)=(/0.0265214244093,0.0265214244093,0.0265214244093/)
       DO I=4,9
           n_weight(I)=0.0199301476312
       END DO
       DO I=10,21
           n_weight(I)=0.0250712367487
       END DO
       n_weight(22:42)=n_weight(1:21)
C
C     MATERIAL PARAMETERS-(obtain its value from material input)
C
       E0 = PROPS(1)
       nu0 = PROPS(2)
       N_V0 = PROPS(3)
       a_0 = PROPS(4)
       Kc = PROPS(5)
       Ec = PROPS(6)
       Alpha = PROPS(7)
       Ko = PROPS(8)
       Eo = PROPS(9)
       c0=(16/3)*(1-nu0**2)/E0
       c1=(32/3)*(1-nu0**2)/(2-nu0)/E0
C=============================================================== 
C      DEFAULT NTENS = 6 FOR 3D
C      NDI = 3 NUMBER OF DIRECT STRESS COMPONENS
C      NSHR= 3 NUMBER OF ENGINEERING SHEAR STRESS COMPONENTS
C ---  Get state variables from STATEV(vector)
       rho_o(1:42)=STATEV(1:42)
       epsel_v(1:6)=STATEV(43:48)
       epsE_v(1:6)  = STATEV(49:54)
       epsin_v(1:6) = STATEV(55:60)
       energy(1:5)=STATEV(61:65)
C ---  TRANSFORM VECTOR VALUE INTO 2D TENSOR
       CALL MAT1_MAT2(NTENS,epsel_v,epsel_o,half)
       CALL MAT1_MAT2(NTENS,epsE_v,epsE_o,half)
       CALL MAT1_MAT2(NTENS,epsin_v,epsin_o,half)
       CALL MAT1_MAT2(NTENS,DSTRAN,deps,half)
       CALL MAT1_MAT2(NTENS,STRAN,eps_o,half)
       CALL MAT1_MAT2(NTENS,STRESS,sigma_o,one)
C
C ---  CALCULATION OF THE EFFECTIVE FOURTH-ORDER ELASTIC STIFFNESS TENSOR
C
       DO I=1,42
           rho_ini(I)=N_V0*a_0**3
       END DO
       IF (TIME(2).EQ.0) THEN
           rho_o(1:42)=rho_ini(1:42)
       END IF
       DO I=1,3
         DO J=1,3
          zero_T(I,J)=0.0D0
         END DO
       END DO
C ---  CALCULATION OF THE EFFECTIVE FOURTH-ORDER ELASTIC STIFFNESS TENSOR
       CALL matDO1(rho_ini,n_cosines,n_weight,zero_T,matDz0,matS0)
       CALL matDO1(rho_o,n_cosines,n_weight,sigma_o,matDz,matS)
C      ELASTICAL TRIAL
       CALL Aijkl_Bkl(matDz,deps,dsig)
C      WRITE(7,*),'DEPS=',DEPS
C      WRITE(7,*),'DSTRAN=',STRAN
C      WRITE(7,*),'DSIG=',DSIG
       CALL Aij_plus_Bij(sigma_o,dsig,sigmaTri)  
       CALL fdDP(rho_o,n_cosines,n_weight,sigmaTri,Yd,flag0,activated) 
C      WRITE(7,*),'SIGMATRI=',SIGMATRI
C      DAMAGE OR NOT?
C
       IF (flag0.EQ.0) THEN
C       NO DAMAGE GENERATED, ELASTIC INCREMENT
         DO I=1,3
           DO J=1,3
             sigmaT(I,J)=sigmaTri(I,J)
             epsin(I,J)=epsin_o(I,J)
           END DO
         END DO
         CALL Aijkl_Bkl(matS0,sigmaT,epsel)
         DO I=1,3
           DO J=1,3
             DepsE(I,J)=Deps(I,J)
             epsE(I,J)=epsE_o(I,J)+Deps(I,J)
             Depsin(I,J)=0.0D0
             sigEner(I,J)=sigma_o(I,J)+theta*dsig(I,J)
             Depsel(I,J)=epsel(I,J)-epsel_o(I,J)
           END DO
         END DO
         rho(1:42)=rho_o(1:42)
C - - -  ENERGY CALCULATION
         CALL Aij_Bij(sigEner,Deps,energy(1))
         CALL Aij_Bij(sigEner,DepsE,energy(2))
         CALL Aij_Bij(sigEner,Depsel,energy(3))
         CALL Aij_Bij(sigEner,Depsin,energy(4))
         energy(5)=0.0D0
         
       ELSEIF (flag0.NE.0) THEN

C     DMAGE OCCURS,IRREVERSIBLE INCREMENT
         INC=1
         flag1=0
10       DO WHILE ( flag1 .NE. 1)
           DO I=1,42
             drho_i(I)=0.0D0
           END DO
           DO I=1,3
             DO J=1,3
               sigmaT_i(I,J)=sigmaTri(I,J)
               R_i(I,J)=0.0D0
               depsin_i(I,J)=0.0D0
             END DO 
           END DO
           DO WHILE (INC.LE.ITER)
C Determine the size of coefficient array 
             A_K=0
             Do I=1,42
               A_K=A_K+activated(I)
             END DO
             CALL fd_lam(rho_o,drho_i,sigmaT_i,R_i,n_cosines, 
     &  n_weight,ddepsin,ddrho,dsig,activated,A_K)
             CALL Aij_plus_Bij(depsin_i,ddepsin,depsin_f)
             CALL Aij_plus_Bij(sigmaT_i,dsig,sigmaT_f)
             detec=0
             DO I=1,42
               drho_f(I)=drho_i(I)+ddrho(I)
               IF (drho_f(I).LT.0.0D0) THEN
                 activated(I)=0
                 detec=detec+1
               END IF
             END DO
             IF (detec.NE.0) THEN
               GOTO 10
             END IF
             DO I=1,3
               DO J=1,3
                 depsin_check(I,J)=0.0D0
               END DO
             END DO
             DO i_point=1,42
               IF (activated(i_point).EQ.1) THEN
                 CALL Dg_Dsigma(n_cosines(i_point,:),n_weight(i_point),
     & sigmaT_f,dg_dsig)
                 DO I=1,3
                   DO J=1,3
                     depsin_check(I,J)=depsin_check(I,J)
     & +drho_f(i_point)*dg_dsig(I,J)
                   END DO
                 END DO
               END IF
             END DO
             DO I=1,3
               DO J=1,3
                 R_f(I,J)=depsin_check(I,J)-depsin_f(I,J)
               END DO
             END DO
             Res=abs(R_f(1,1)+R_f(2,2)+R_f(3,3))/3
             fd_all=0.0D0
             CALL fdDP1(rho_o,drho_f,n_cosines,n_weight,sigmaT_f,
     & activated,fd_all)
             IF ((fd_all.LT.FTOL0).AND.(Res.LT.FTOL1)) THEN
               flag1=1
               GOTO 20
             ELSE
               DO I=1,3
                 DO J=1,3
                   sigmaT_i(I,J)=sigmaT_f(I,J)
                   depsin_i(I,J)=depsin_f(I,J)
                   R_i(I,J)=R_f(I,J)
                 END DO
               END DO
               DO i_point=1,42
                 drho_i(i_point)=drho_f(i_point)
               END DO
               INC=INC+1
             END IF
           END DO
         END DO       
C  Calculate increment of variables
20       DO I=1,3
           DO J=1,3
             sigmaT(I,J)=sigmaT_f(I,J)
             Depsin(I,J)=depsin_f(I,J)   
           END DO
         END DO
         DO i_point=1,42
           drho(i_point) = drho_f(i_point)
           rho(i_point) = rho_o(i_point) + drho_f(i_point)
         END DO
C
         CALL Aij_minus_Bij(sigmaT,sigma_o,dsig)
         CALL Aijkl_Bkl(matS0,dsig,Depsel)
         DO I=1,3
           DO J=1,3
             sigEner(I,J)=sigma_o(I,J)+theta*dsig(I,J)
             DepsE(I,J)=Deps(I,J)-Depsin(I,J)
           END DO
         END DO
         CALL Aij_plus_Bij(epsel_o,Depsel,epsel)
         CALL Aij_plus_Bij(epsE_o,DepsE,epsE)
         CALL Aij_plus_Bij(epsin_o,Depsin,epsin)
        
C - - - ENERGY CALCULATION
         CALL Aij_Bij(sigEner,Deps,energy(1))
         CALL Aij_Bij(sigEner,DepsE,energy(2))
         CALL Aij_Bij(sigEner,Depsel,energy(3))
         CALL Aij_Bij(sigEner,Depsin,energy(4))
         CALL fdDP(rho,n_cosines,n_weight,sigmaT,Yd,flag0,activated_1)
         DO i_point=1,42
           IF (activated(i_point).EQ.1) THEN
             energy(5)=energy(5)+Yd(i_point)*drho(i_point)
           END IF
         END DO
       END IF
C
C      write(7,*) 'OMEGA=',OMEGA
C	 write(7,*) 'SIGMA=',SIGMA
C	 write(7,*) 'epsin=',epsin
C	 write(7,*) 'DepsE=',DepsE 
C 
C --- UPDATING STATE VARIABLES
C
       CALL MAT2_MAT1(NTENS,epsel,epsel_v,TWO)
       CALL MAT2_MAT1(NTENS,epsE,epsE_v,TWO)
       CALL MAT2_MAT1(NTENS,epsin,epsin_v,TWO)
       CALL MAT2_MAT1(NTENS,sigmaT,STRESS,ONE)
C
C  STORE STATE VARIABLES
C
       STATEV(1:42)=rho(1:42)
       STATEV(43:48)=epsel_v(1:6)
       STATEV(49:54)=epsE_v(1:6)
       STATEV(55:60)=epsin_v(1:6)
       STATEV(61:65)=energy(1:5)
C
C  CREATE NEW JACOBIAN
C --- CONVERSION OF THE FORTH ORDER STIFFNESS TENSOR TO A SECOND
C     ORDER TENSOR
C
       CALL MAT4_MAT2(matDz,matDz_2,1)
C 	  write(7,*) 'DDSDDE=',matDz_2
C
C --- UPDATING THE JACOBIAN MATRIX FOR ACCELERATING CONVERGENCE
C
       DO I = 1 , NTENS
         DO J = 1 , NTENS
          DDSDDE (I , J) = matDz_2(I , J)
         ENDDO
       ENDDO
C
C
       RETURN
       END
C
C
C ====================================================================
C                              SUBROUTINES 
C ====================================================================
C
       SUBROUTINE matDO1(rho,n_cosines,n_weight,sigmaT,matDz,matS)
CC====================================================================
       INCLUDE 'ABA_PARAM.INC'
       DIMENSION rho(42),n_cosines(42,3),n_weight(42),
     & sigmaT(3,3),matDz(3,3,3,3),matDz_2(6,6),
     & matS(3,3,3,3),matS_2(6,6),NN(3,3,3,3),
     & Tri(3,3,3,3)
       INTEGER E(3,3)
       DOUBLE PRECISION rho,n_cosines,n_weight,sigmaT,
     & matDz,matDz_2,matS,matS_2,NN,Tri,E0,nu0,N_V0,a_0,
     & Kc,Ec,Alpha,Ko,Eo,c0,c1
C
       COMMON/parameters/E0,nu0,N_V0,a_0,Kc,Ec,Alpha,Ko,Eo,c0,c1
       DATA TPF,TWO /2.5D-1,2.0D0/
C    
       DO I=1,6
         DO J=1,6
           matS_2(I,J)=0.0D0
           matDz_2(I,J)=0.0D0
         END DO
       END DO
C
       DO I=1,3
         DO J=1,3
           E(I,J)=0.0D0
           DO K=1,3
             DO L=1,3
               NN(I,J,K,L)=0.0D0
               Tri(I,J,K,L)=0.0D0
               matS(I,J,K,L)=0.0D0
               matDz(I,J,K,L)=0.0D0
             END DO
           END DO 
         END DO
       END DO
C
      CALL DIDENTITY_2(E)
      b1=(1+nu0)/E0/2
      b2=nu0/E0
C 
       DO i_point=1,42
         DO i=1,3
           DO j=1,3
             DO k=1,3
               DO l=1,3
                 Call Ni_Aij_Nj(n_cosines(i_point, :), sigmaT, sigma_nn)
                 IF (sigma_nn.GE.0) THEN
                      NN(i,j,k,l) = NN(i,j,k,l)+ rho(i_point) 
     & *n_weight(i_point)*c0*n_cosines(i_point,i)*n_cosines(i_point,j)
     & *n_cosines(i_point,k)*n_cosines(i_point,l)
                 END IF
                 Tri(i,j,k,l) = Tri(i,j,k,l) + rho(i_point)
     & *n_weight(i_point)*c1*(TPF*(n_cosines(i_point,i)
     & *n_cosines(i_point,k)*E(j,l)+n_cosines(i_point,i)
     & *n_cosines(i_point,l)*E(j,k)+E(i,k) *n_cosines(i_point,j)
     & *n_cosines(i_point,l)+E(i,l)*n_cosines(i_point,j)
     & *n_cosines(i_point,k))-n_cosines(i_point,i)*n_cosines(i_point,j)
     & *n_cosines(i_point,k)*n_cosines(i_point,l))
                 END DO
             END DO
           END DO
         END DO
       END DO
C       
       DO i=1,3
         DO j=1,3
           DO k=1,3
             DO l=1,3
                matS(i,j,k,l) = b1*(E(i,k)*E(j,l)+E(i,l)*E(j,k))
     & -b2*E(i,j)*E(k,l)+NN(i,j,k,l)+Tri(i,j,k,l)
             END DO
          END DO
        END DO
      END DO
C
       CALL MAT4_MAT2(matS,matS_2,2)
       CALL INVERSE(matS_2,6,6,matDz_2)
       CALL MAT2_MAT4(matDz_2,matDz,1)
C	  
      RETURN
      END
C
       SUBROUTINE matDO2(rho,drho,n_cosines,n_weight,sigmaT,matDz)
CC ====================================================================
      INCLUDE 'ABA_PARAM.INC'
C
       DIMENSION rho(42),drho(42),n_cosines(42,3),n_weight(42),
     & sigmaT(3,3),matDz(3,3,3,3),matDz_2(6,6),
     & matS(3,3,3,3),matS_2(6,6),NN(3,3,3,3),
     & Tri(3,3,3,3)
       INTEGER E(3,3)
       DOUBLE PRECISION rho,drho,n_cosines,n_weight,sigmaT,
     & matDz,matDz_2,matS,matS_2,NN,Tri,E0,nu0,N_V0,a_0,Kc,
     & Ec,Alpha,Ko,Eo,c0,c1
C
      COMMON/parameters/E0,nu0,N_V0,a_0,Kc,Ec,Alpha,Ko,Eo,c0,c1
      DATA TPF,TWO /2.5D-1,2.0D0/
C
C      INITIALIZING TENSORS AND MATRIX
C    
       DO I=1,6
         DO J=1,6
           matS_2(I,J)=0.0D0
           matDz_2(I,J)=0.0D0
         END DO
       END DO
C
       DO I=1,3
         DO J=1,3
           E(I,J)=0.0D0
           DO K=1,3
             DO L=1,3
               NN(I,J,K,L)=0.0D0
               Tri(I,J,K,L)=0.0D0
               matS(I,J,K,L)=0.0D0
             END DO
           END DO 
         END DO
       END DO
C
      CALL DIDENTITY_2(E)
      b1=(1+nu0)/E0/2
      b2=nu0/E0
C 
      DO i_point=1,42
        DO i=1,3
          DO j=1,3
            DO k=1,3
              DO l=1,3
                Call Ni_Aij_Nj(n_cosines(i_point, 1:3),sigmaT,sigma_nn)
                IF (sigma_nn.GE.0) THEN
                  NN(i,j,k,l) = NN(i,j,k,l)+(rho(i_point)
     & +drho(i_point))*n_weight(i_point)*c0*n_cosines(i_point,i)
     & *n_cosines(i_point,j)*n_cosines(i_point,k)*n_cosines(i_point,l)
                 END IF
                 Tri(i,j,k,l) = Tri(i,j,k,l)+(rho(i_point)
     & +drho(i_point))*n_weight(i_point)*c1*(TPF*
     & (n_cosines(i_point,i)*n_cosines(i_point,k)*E(j,l) 
     & +n_cosines(i_point,i)*n_cosines(i_point,l)*E(j,k)+E(i,k)
     & *n_cosines(i_point,j)*n_cosines(i_point,l)+E(i,l)
     & *n_cosines(i_point,j)*n_cosines(i_point,k))-n_cosines(i_point,i)
     & *n_cosines(i_point,j)*n_cosines(i_point,k)*n_cosines(i_point,l))
                 END DO
              END DO
           END DO
         END DO
       END DO
C       
       DO i=1,3
         DO j=1,3
           DO k=1,3
             DO l=1,3
                matS(i,j,k,l) = b1*(E(i,k)*E(j,l)+E(i,l)*E(j,k))
     & -b2*E(i,j)*E(k,l)+NN(i,j,k,l)+Tri(i,j,k,l)
            END DO
          END DO
        END DO
      END DO
C
       CALL MAT4_MAT2(matS,matS_2,2)
       CALL INVERSE(matS_2,6,6,matDz_2)
       CALL MAT2_MAT4(matDz_2,matDz,1)
C	  
      RETURN
      END
C
       SUBROUTINE MAT1_MAT2(NTENS,VECTOR,TENSOR,FACT)
C ====================================================================
C 
C                MAT1_MAT2 : VECTOR TO TENSOR  
C
C ====================================================================
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION VECTOR(NTENS),TENSOR(3,3)
      DOUBLE PRECISION VECTOR,TENSOR
C
      DO I=1,3
        DO J=1,3
          TENSOR(I,J)=0.0D0
        END DO
      END DO
C
      TENSOR(1 , 1) = VECTOR( 1 )
      TENSOR(2 , 2) = VECTOR( 2 )
	  TENSOR(3 , 3) = VECTOR( 3 )
      TENSOR(1 , 2) = VECTOR( 4 )*FACT
      TENSOR(2 , 1) = VECTOR( 4 )*FACT
      TENSOR(1 , 3) = VECTOR( 5 )*FACT
      TENSOR(3 , 1) = VECTOR( 5 )*FACT
      TENSOR(2 , 3) = VECTOR( 6 )*FACT
      TENSOR(3 , 2) = VECTOR( 6 )*FACT
C
      RETURN
      END
C
       SUBROUTINE MAT2_MAT1(NTENS,TENSOR,VECTOR,FACT)
C
C ====================================================================
C
C =================== MAT2_MAT1: TENSOR TO VECTOR=====================
C
C ====================================================================
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION VECTOR(NTENS),TENSOR(3,3)
      DOUBLE PRECISION VECTOR,TENSOR
C
      DO I=1,NTENS
          VECTOR(I)=0.0D0
      END DO
C
      VECTOR( 1 ) = TENSOR(1 , 1)
      VECTOR( 2 ) = TENSOR(2 , 2)
      VECTOR( 3 ) = TENSOR(3 , 3)
      VECTOR( 4 ) = TENSOR(1 , 2)*FACT
      VECTOR( 5 ) = TENSOR(1 , 3)*FACT
      VECTOR( 6 ) = TENSOR(2 , 3)*FACT
C
      RETURN
      END
C
       SUBROUTINE MAT2_MAT4(DMATRIX,TENSOR,ICOE)
C======================================================================
C
C                             MAT2_MAT4
C
C======================================================================
      INCLUDE 'ABA_PARAM.INC'
C  
      DIMENSION DMATRIX(6,6),TENSOR(3,3,3,3)
      DOUBLE PRECISION DMATRIX,TENSOR
C
C     INITALIZATION
C
       DO I=1,3
         DO J=1,3
           DO K=1,3
             DO L=1,3 
               TENSOR(I,J,K,L)=0.0D0
             END DO
           END DO
         END DO
       END DO
C
      IF (ICOE.EQ.1) THEN
         COE1=1.
         COE2=1.
      ELSEIF(ICOE.EQ.2) THEN
         COE1=2.
         COE2=4.
      END IF
C      
      TENSOR(1,1,1,1) = DMATRIX(1,1)
      TENSOR(1,1,2,2) = DMATRIX(1,2)
      TENSOR(1,1,3,3) = DMATRIX(1,3)
      TENSOR(1,1,1,2) = DMATRIX(1,4)/COE1
      TENSOR(1,1,2,1) = DMATRIX(1,4)/COE1
      TENSOR(1,1,2,3) = DMATRIX(1,5)/COE1
      TENSOR(1,1,3,2) = DMATRIX(1,5)/COE1
      TENSOR(1,1,1,3) = DMATRIX(1,6)/COE1
      TENSOR(1,1,3,1) = DMATRIX(1,6)/COE1
C
      TENSOR(2,2,1,1) = DMATRIX(2,1)
      TENSOR(2,2,2,2) = DMATRIX(2,2)
      TENSOR(2,2,3,3) = DMATRIX(2,3)
      TENSOR(2,2,1,2) = DMATRIX(2,4)/COE1
      TENSOR(2,2,2,1) = DMATRIX(2,4)/COE1
      TENSOR(2,2,2,3) = DMATRIX(2,5)/COE1
      TENSOR(2,2,3,2) = DMATRIX(2,5)/COE1
      TENSOR(2,2,1,3) = DMATRIX(2,6)/COE1
      TENSOR(2,2,3,1) = DMATRIX(2,6)/COE1
C
      TENSOR(3,3,1,1) = DMATRIX(3,1)
      TENSOR(3,3,2,2) = DMATRIX(3,2)
      TENSOR(3,3,3,3) = DMATRIX(3,3)
      TENSOR(3,3,1,2) = DMATRIX(3,4)/COE1
      TENSOR(3,3,2,1) = DMATRIX(3,4)/COE1
      TENSOR(3,3,2,3) = DMATRIX(3,5)/COE1
      TENSOR(3,3,3,2) = DMATRIX(3,5)/COE1
      TENSOR(3,3,1,3) = DMATRIX(3,6)/COE1
      TENSOR(3,3,3,1) = DMATRIX(3,6)/COE1
C
      TENSOR(1,2,1,1) = DMATRIX(4,1)/COE1
      TENSOR(1,2,2,2) = DMATRIX(4,2)/COE1
      TENSOR(1,2,3,3) = DMATRIX(4,3)/COE1
      TENSOR(1,2,1,2) = DMATRIX(4,4)/COE2
      TENSOR(1,2,2,1) = DMATRIX(4,4)/COE2
      TENSOR(1,2,2,3) = DMATRIX(4,5)/COE2
      TENSOR(1,2,3,2) = DMATRIX(4,5)/COE2
      TENSOR(1,2,1,3) = DMATRIX(4,6)/COE2
      TENSOR(1,2,3,1) = DMATRIX(4,6)/COE2
C
      TENSOR(2,3,1,1) = DMATRIX(5,1)/COE1
      TENSOR(2,3,2,2) = DMATRIX(5,2)/COE1
      TENSOR(2,3,3,3) = DMATRIX(5,3)/COE1
      TENSOR(2,3,1,2) = DMATRIX(5,4)/COE2
      TENSOR(2,3,2,1) = DMATRIX(5,4)/COE2
      TENSOR(2,3,2,3) = DMATRIX(5,5)/COE2
      TENSOR(2,3,3,2) = DMATRIX(5,5)/COE2
      TENSOR(2,3,1,3) = DMATRIX(5,6)/COE2
      TENSOR(2,3,3,1) = DMATRIX(5,6)/COE2
C
      TENSOR(1,3,1,1) = DMATRIX(6,1)/COE1
      TENSOR(1,3,2,2) = DMATRIX(6,2)/COE1
      TENSOR(1,3,3,3) = DMATRIX(6,3)/COE1
      TENSOR(1,3,1,2) = DMATRIX(6,4)/COE2
      TENSOR(1,3,2,1) = DMATRIX(6,4)/COE2
      TENSOR(1,3,2,3) = DMATRIX(6,5)/COE2
      TENSOR(1,3,3,2) = DMATRIX(6,5)/COE2
      TENSOR(1,3,1,3) = DMATRIX(6,6)/COE2
      TENSOR(1,3,3,1) = DMATRIX(6,6)/COE2
C      
      TENSOR(2,1,1,1) = DMATRIX(4,1)/COE1
      TENSOR(2,1,2,2) = DMATRIX(4,2)/COE1
      TENSOR(2,1,3,3) = DMATRIX(4,3)/COE1
      TENSOR(2,1,1,2) = DMATRIX(4,4)/COE2
      TENSOR(2,1,2,1) = DMATRIX(4,4)/COE2
      TENSOR(2,1,2,3) = DMATRIX(4,5)/COE2
      TENSOR(2,1,3,2) = DMATRIX(4,5)/COE2
      TENSOR(2,1,1,3) = DMATRIX(4,6)/COE2
      TENSOR(2,1,3,1) = DMATRIX(4,6)/COE2
C
      TENSOR(3,2,1,1) = DMATRIX(5,1)/COE1
      TENSOR(3,2,2,2) = DMATRIX(5,2)/COE1
      TENSOR(3,2,3,3) = DMATRIX(5,3)/COE1
      TENSOR(3,2,1,2) = DMATRIX(5,4)/COE2
      TENSOR(3,2,2,1) = DMATRIX(5,4)/COE2
      TENSOR(3,2,2,3) = DMATRIX(5,5)/COE2
      TENSOR(3,2,3,2) = DMATRIX(5,5)/COE2
      TENSOR(3,2,1,3) = DMATRIX(5,6)/COE2
      TENSOR(3,2,3,1) = DMATRIX(5,6)/COE2
C
      TENSOR(3,1,1,1) = DMATRIX(6,1)/COE1
      TENSOR(3,1,2,2) = DMATRIX(6,2)/COE1
      TENSOR(3,1,3,3) = DMATRIX(6,3)/COE1
      TENSOR(3,1,1,2) = DMATRIX(6,4)/COE2
      TENSOR(3,1,2,1) = DMATRIX(6,4)/COE2
      TENSOR(3,1,2,3) = DMATRIX(6,5)/COE2
      TENSOR(3,1,3,2) = DMATRIX(6,5)/COE2
      TENSOR(3,1,1,3) = DMATRIX(6,6)/COE2
      TENSOR(3,1,3,1) = DMATRIX(6,6)/COE2
C
C
      RETURN
      END
C
       SUBROUTINE MAT4_MAT2(TENSOR,DMATRIX,ICOE)
C
C ====================================================================
C                        MAT4_MAT2                                                                  
C I        THIS PROGRAM TRANSFORMS THE FOURTH ORDER COMPLIANCE       I
C I        TENSOR TO A SECOND ORDER MATRIX                           I
C I                                                                  I
C ====================================================================
C
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION TENSOR(3,3,3,3),DMATRIX(6,6)
      DOUBLE PRECISION DMATRIX,TENSOR
C
      DATA ZERO,TWO /0.0D0,2.0D0/
C
C     D2 = THE SECOND ORDER STIFFNESS MATRIX
C
      DO I=1,6
        DO J=1,6
          DMATRIX(I,J)=0.0D0
        END DO
      END DO

      IF (ICOE.EQ.1) THEN
         COE1=1.
         COE2=1.
      ELSEIF(ICOE.EQ.2) THEN
         COE1=2.
         COE2=4.
      END IF
C
      DMATRIX(1,1)=TENSOR(1,1,1,1)
      DMATRIX(1,2)=TENSOR(1,1,2,2)
      DMATRIX(1,3)=TENSOR(1,1,3,3)
      DMATRIX(1,4)=TENSOR(1,1,1,2)*COE1
      DMATRIX(1,5)=TENSOR(1,1,2,3)*COE1
      DMATRIX(1,6)=TENSOR(1,1,1,3)*COE1
C
      DMATRIX(2,1)=TENSOR(2,2,1,1)
      DMATRIX(2,2)=TENSOR(2,2,2,2)
      DMATRIX(2,3)=TENSOR(2,2,3,3)
      DMATRIX(2,4)=TENSOR(2,2,1,2)*COE1
      DMATRIX(2,5)=TENSOR(2,2,2,3)*COE1
      DMATRIX(2,6)=TENSOR(2,2,1,3)*COE1
C
      DMATRIX(3,1)=TENSOR(3,3,1,1)
      DMATRIX(3,2)=TENSOR(3,3,2,2)
      DMATRIX(3,3)=TENSOR(3,3,3,3)
      DMATRIX(3,4)=TENSOR(3,3,1,2)*COE1
      DMATRIX(3,5)=TENSOR(3,3,2,3)*COE1
      DMATRIX(3,6)=TENSOR(3,3,1,3)*COE1
C 
C - — - here, engineering shear strain is used
C
      DMATRIX(4,1)=TENSOR(1,2,1,1)*COE1
      DMATRIX(4,2)=TENSOR(1,2,2,2)*COE1
      DMATRIX(4,3)=TENSOR(1,2,3,3)*COE1
      DMATRIX(4,4)=TENSOR(1,2,1,2)*COE2
      DMATRIX(4,5)=TENSOR(1,2,2,3)*COE2
      DMATRIX(4,6)=TENSOR(1,2,1,3)*COE2
C
      DMATRIX(5,1)=TENSOR(2,3,1,1)*COE1
      DMATRIX(5,2)=TENSOR(2,3,2,2)*COE1
      DMATRIX(5,3)=TENSOR(2,3,3,3)*COE1
      DMATRIX(5,4)=TENSOR(2,3,1,2)*COE2
      DMATRIX(5,5)=TENSOR(2,3,2,3)*COE2
      DMATRIX(5,6)=TENSOR(2,3,1,3)*COE2
C
      DMATRIX(6,1)=TENSOR(1,3,1,1)*COE1
      DMATRIX(6,2)=TENSOR(1,3,2,2)*COE1
      DMATRIX(6,3)=TENSOR(1,3,3,3)*COE1
      DMATRIX(6,4)=TENSOR(1,3,1,2)*COE2
      DMATRIX(6,5)=TENSOR(1,3,2,3)*COE2
      DMATRIX(6,6)=TENSOR(1,3,1,3)*COE2
C
C
      RETURN
      END
C
       SUBROUTINE INVERSE(A,N,NP,AINV)
C========================================================================
C
C    CALCULATE THE SECOND ORDER TENSOR A'S INVERSE, AINV
C    A^{-1} = AINV    
C    this subroutine inverses a (n x n) A matrix
C	 following a Gauss-Jordan elimination process
C
C========================================================================
       INCLUDE 'ABA_PARAM.INC'
       DIMENSION A(NP,NP),IPIV(NP),INDXR(NP),INDXC(NP),
     & A0(NP,NP),AINV(NP,NP)
       DOUBLE PRECISION A,IPIV,INDXR,INDXC,A0,AINV
C
      DO J=1,N
        IPIV(J)=0
      END DO
C
C     storage of the original A matrix
C
      DO I=1,N
        DO J=1,N
          A0(I,J)=A(I,J)
        END DO
      END DO
C
C	find a pivot (largest absolute value) among the rows of A that have not already been reduced
C
      DO I=1,N
        BIG=0.0D0
        DO J=1,N
          IF(IPIV(J).NE.1)THEN
            DO K=1,N
                IF(IPIV(K).EQ.0)THEN
                  IF(ABS(A(J,K)).GE.BIG)THEN
                    BIG=ABS(A(J,K))
                    IROW=J
                    ICOL=K
                    PIV=A(J,K)
                  END IF
                ELSEIF(IPIV(K).GT.1)THEN
                  write (7,*) 'Singular Matrix'
                END IF
            END DO
          END IF
        END DO
C
        IPIV(ICOL)=IPIV(ICOL)+1
        INDXR(I)=IROW
        INDXC(I)=ICOL
C	  
C     interchange the rows to put the pivot on the diagonal
C
        IF(IROW.NE.ICOL)THEN
          DO L=1,N
            DUM=A(IROW,L)
            A(IROW,L)=A(ICOL,L)
            A(ICOL,L)=DUM
          END DO
        END IF
C
C     reduction of the row of the pivot
C       
        IF(PIV.EQ.0) write (7,*) 'Singular Matrix2'
C       
        PIVIN=1./PIV          ! numerical stabilization
C
        A(ICOL,ICOL)=1.       ! numerical stabilization
        DO L=1,N
           A(ICOL,L)=A(ICOL,L)*PIVIN
        END DO
C
C     reduction of the column of the pivot
C
        DO LL=1,N
           IF(LL.NE.ICOL)THEN
             DUM=A(LL,ICOL)
             A(LL,ICOL)=0.    ! numerical stabilization
             DO L=1,N
                A(LL,L)=A(LL,L)-A(ICOL,L)*DUM
             END DO
           END IF
        END DO
      END DO
C
C     unscramble the columns to get A^{-1}
C		
      DO J=N,1,-1   ! reverse DO loop
        DO K=1,N
          DUM=A(K,INDXR(J))
          A(K,INDXR(J))=A(K,INDXC(J))
          A(K,INDXC(J))=DUM
        END DO
      END DO
C
C	restitution process of A and Ainv
C
      DO I=1,N
        DO J=1,N
          AINV(I,J)=A(I,J)
        END DO
      END DO
C
      DO I=1,N
        DO J=1,N
          A(I,J)=A0(I,J)
        END DO
      END DO
C     
      RETURN
      END
C
       SUBROUTINE DIDENTITY_2(DELTA)
C========================================================================
C                                                                       =
C                               DIDENTITY_2                             =
C                                                                       =
C========================================================================
C
      INTEGER DELTA(3,3)
      DATA ZERO,ONE /0,1/
C         
      DO I=1,3
        DO J=1,3
          DELTA(I,J)=ZERO
        END DO
		DELTA(I,I)=ONE
      END DO
C
      RETURN
      END
C
       SUBROUTINE fdDP(rho,n_cosines,n_weight,sigmaT,Yd,flag0,activated)
C=========================================================================
C       Damage Yield Function
C=========================================================================
       INCLUDE 'ABA_PARAM.INC'
C
       DIMENSION rho(42),n_cosines(42,3),n_weight(42),
     & sigmaT(3,3),fd(42),Yd(42)
       INTEGER activated(42),flag0
       DOUBLE PRECISION rho,n_cosines,n_weight,
     & sigmaT,fd,Yd
C       parameter(FTOL0=1.D-3)
C
      DO i_point=1,42
        fd(i_point)=0
        Yd(i_point)=0
        activated(i_point)=0
      END DO
      flag0=0
      DO  i_point=1,42
         call  YD_F(n_cosines(i_point,:), n_weight(i_point), sigmaT,
     & rho(i_point),fd(i_point),Yd(i_point))     
         IF (fd(i_point).GE.0) THEN
            flag0=flag0+1
            activated(i_point)=1
         END IF
      END DO
      RETURN
      END
C
       SUBROUTINE fdDP1(rho,drho,n_cosines,n_weight,sigmaT,
     & activated,fd_all)
C=========================================================================
       INCLUDE 'ABA_PARAM.INC'
C
       DIMENSION rho(42),drho(42),n_cosines(42,3),n_weight(42),
     & sigmaT(3,3)
       INTEGER activated(42),flag0
       DOUBLE PRECISION rho,drho,n_cosines,n_weight,sigmaT
C
      fd_al=0.0D0
      DO i_point=1,42
        IF (activated(i_point).EQ.1) THEN
          rho_new=rho(i_point)+drho(i_point)
          CALL  YD_F(n_cosines(i_point,:), n_weight(i_point), sigmaT,
     & rho_new,fd,Yd)
          fd_al=fd_al+fd**2
        END IF
      END DO
      fd_all=SQRT(fd_al)
C
      RETURN
      END
C
       SUBROUTINE YD_F(n,w,sigmaT,rho,fd,Yd)
C=========================================================================
C       Damage Yield single surface
C=========================================================================
       INCLUDE 'ABA_PARAM.INC'
	 DIMENSION n(3),sigmaT(3,3),Temp(3,3),M(3,3,3,3)
       INTEGER E(3,3)
       DOUBLE PRECISION n,sigmaT,Temp,M,E0,nu0,N_V0,a_0,Kc,
     & Ec,Alpha,Ko,Eo,c0,c1
       DATA TPF,ONE /2.5D-1,1.0D0/
       COMMON/parameters/E0,nu0,N_V0,a_0,Kc,Ec,Alpha,Ko,Eo,c0,c1
C
      DO I =1,3
        DO J=1,3
          E(I,J)=0
          Temp(I,J)=0
        END DO
      END DO
C
      CALL DIDENTITY_2(E)
      Trsigma=sigmaT(1,1)+sigmaT(2,2)+sigmaT(3,3)
      Call Ni_Aij_Nj(n,sigmaT,sigma_nn)
      IF (sigma_nn.GE.0) THEN
        DO i=1,3
          DO j=1,3
            DO k=1,3
              DO l=1,3
                M(i,j,k,l) = w*c0*n(i)*n(j)*n(k)*n(l)
     & +w*c1*((TPF*(n(i)*n(k)*E(j,l)+n(i)*n(l)*E(j,k)+E(i,k)*n(j)*n(l)
     & +E(i,l)*n(j)*n(k))-n(i)*n(j)*n(k)*n(l)))
              END DO
            END DO
          END DO
        END DO
        CALL Aijkl_Bij(M,sigmaT,Temp)
        CALL Aij_Bij(Temp,sigmaT,double_Yd)
        Yd= 0.5D0*double_Yd
    	  fd = Yd + alpha*Trsigma - Ko*(ONE+Eo*rho)
      ELSE
        DO i=1,3
          DO j=1,3
            DO k=1,3
              DO l=1,3
                M(i,j,k,l) = w*c1*((0.25*(n(i)*n(k)*E(j,l)
     & +n(i)*n(l)*E(j,k)+E(i,k)*n(j)*n(l)+E(i,l)*n(j)*n(k))
     & -n(i)*n(j)*n(k)*n(l)))
              END DO
            END DO
          END DO
        END DO
        CALL Aijkl_Bij(M,sigmaT,Temp)
        CALL Aij_Bij(Temp,sigmaT,double_Yd)
     	  Yd = 0.5*double_Yd
    	  fd = Yd + alpha*Trsigma - Kc*(1+Ec*rho)        
      END IF
C
      RETURN
      END
C
       SUBROUTINE Aijkl_Bkl(A,B,C)
C========================================================================
C                                                                       =
C                              Aijkl_Bkl                                =
C                                                                       =
C========================================================================
C
      DIMENSION A(3,3,3,3),B(3,3),C(3,3)
      DOUBLE PRECISION A,B,C
      DATA ZERO /0.0D0/
C
      DO I=1,3
        DO J=1,3
          C(I,J)=ZERO
          DO K=1,3
            DO L=1,3
              C(I,J)=C(I,J)+A(I,J,K,L)*B(K,L)
            END DO
          END DO
        END DO
      END DO
C
      RETURN
      END
C
       SUBROUTINE Aij_minus_Bij(A,B,C)
C====================================================================================
C                                                                      *
C                                  Aij_minus_Bij                        *
C                                                                      *
C====================================================================================
C
      DIMENSION A(3,3),B(3,3),C(3,3)
      DOUBLE PRECISION A,B,C
C
      DO I=1,3
        DO J=1,3
          C(I,J)=0
          C(I,J)=A(I,J)-B(I,J)
        END DO
      END DO
C
      RETURN
      END
C
       SUBROUTINE Aij_plus_Bij(A,B,C)
C=======================================================================
C                                                                      *
C                                  Aij_plus_Bij                        *
C                                                                      *
C=======================================================================
C
      DIMENSION A(3,3),B(3,3),C(3,3)
      DOUBLE PRECISION A,B,C
C
      DO I=1,3
        DO J=1,3
          C(I,J)=0
          C(I,J)=A(I,J)+B(I,J)
        END DO
      END DO
C
      RETURN
      END
C
       SUBROUTINE Aij_Bij(A,B,C)
C====================================================================================
C                                                                                   =
C                                   Aij_Bij                                         =
C                                                                                   =
C====================================================================================
       INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION A(3,3),B(3,3)
      DOUBLE PRECISION A,B,C
      DATA ZERO /0.0D0/
C
      C=ZERO
      DO I=1,3
        DO J=1,3
          C=C+A(I,J)*B(I,J)
        END DO
      END DO
C
      RETURN
      END
C
       SUBROUTINE A_Bij(A,B,C)
CC===========================================================
C
C                          A*Bij=Cij	  
C
C============================================================
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION B(3,3),C(3,3)
      DOUBLE PRECISION A,B,C
C
      DO I=1,3
        DO J=1,3
          C(I,J)=A*B(I,J)
        END DO
      END DO
C
      RETURN
      END
C
       SUBROUTINE Ni_Aij_Nj(N,A,C)
C
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION N(3),A(3,3)
      DOUBLE PRECISION A,N,C
C
      C=0
      DO I=1,3
        DO J=1,3
          C= C + N(I)*A(I,J)*N(J)
        END DO
      END DO
C
      RETURN
      END   
C
       SUBROUTINE Aijkl_Bij(A,B,C)
C=======================================================================
C                                                                      *
C                          Aijkl_Bij                                   *
C                                                                      *
C=======================================================================
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION A(3,3,3,3),B(3,3),C(3,3)
      DOUBLE PRECISION A,B,C
      DATA ZERO /0.0D0/
C
      DO K=1,3
        DO L=1,3
          C(K,L)=ZERO
          DO I=1,3
            DO J=1,3
              C(K,L)=C(K,L)+A(I,J,K,L)*B(I,J)
            END DO
          END DO
        END DO
      END DO
C
      RETURN
      END
C
       SUBROUTINE Aij_Bj(A,B,K,C)
C=======================================================================
C                                                                      *
C                           Aij_Bj                                    *
C                                                                      *
C=======================================================================
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION A(K,K),B(K),C(K)
      DOUBLE PRECISION A,B,C
C
      DO I=1,K
        C(I)=0.0D0
        DO J=1,K
          C(I) = C(I)+A(I,J)*B(J)
        END DO
      END DO
C
      RETURN
      END
C
       SUBROUTINE DF_DSIG_Drho(n,w,sigmaT,rho,drho,Df_Dsig,Df_Drho,fd)
C=======================================================================
C                                                                      *
C                           DF_DSIG_Drho                               *
C                                                                      *
C=======================================================================
       INCLUDE 'ABA_PARAM.INC'
C
       DIMENSION n(3),sigmaT(3,3),Temp(3,3),
     & Df_Dsig(3,3),M(3,3,3,3)
       INTEGER E(3,3)
       DOUBLE PRECISION n,sigmaT,Temp,Df_Dsig,M,E0,nu0,N_V0,
     & a_0,Kc,Ec,Alpha,Ko,Eo,c0,c1
       COMMON/parameters/E0,nu0,N_V0,a_0,Kc,Ec,Alpha,Ko,Eo,c0,c1
       DATA TPF,PFIVE /2.5D-1,0.5D0/
C
	  DO I=1,3
        DO J=1,3
          Df_Dsig(I,J)=0.0D0
          Temp(I,J)=0.0D0
          E(I,J)=0.0D0
          DO K=1,3
            DO L=1,3
              M(I,J,K,L)=0.0D0
            END DO
          END DO
        END DO
      END DO
C
      CALL DIDENTITY_2(E)
      Trsigma=sigmaT(1,1)+sigmaT(2,2)+sigmaT(3,3)
	  Call Ni_Aij_Nj(n,sigmaT,sigma_nn)
	  IF (sigma_nn.GE.0) THEN
        DO i=1,3
          DO j=1,3
            DO k=1,3
              DO l=1,3
                M(i,j,k,l) = w*c0*n(i)*n(j)*n(k)*n(l)
     & +w*c1*((TPF*(n(i)*n(k)*E(j,l)+n(i)*n(l)*E(j,k)+E(i,k)*n(j)*n(l)
     & +E(i,l)*n(j)*n(k))-n(i)*n(j)*n(k)*n(l)))
              END DO
            END DO
          END DO
        END DO
        CALL Aijkl_Bij(M,sigmaT,Temp)
        CALL Aij_Bij(Temp,sigmaT,Yd)
        fd = PFIVE*Yd+alpha*Trsigma-Ko*(1+Eo*(rho+drho))
        Df_Drho=-Ko*Eo
        DO i=1,3
          DO j=1,3
            IF (i.EQ.j) THEN
              Df_Dsig(i,j)=Temp(i,j)+alpha
            ELSE
              Df_Dsig(i,j)=Temp(i,j)
            END IF
          END DO
        END DO
	  ELSE
        DO i=1,3
          DO j=1,3
            DO k=1,3
              DO l=1,3
                M(i,j,k,l) = w*c1*((TPF*(n(i)*n(k)*E(j,l)
     & +n(i)*n(l)*E(j,k)+E(i,k)*n(j)*n(l)+E(i,l)*n(j)*n(k))
     & -n(i)*n(j)*n(k)*n(l)))
              END DO
            END DO
          END DO
        END DO
        CALL Aijkl_Bij(M,sigmaT,Temp)
        CALL Aij_Bij(Temp,sigmaT,Yd)
        fd = PFIVE*Yd+alpha*Trsigma-Kc*(1+Ec*(rho+drho))
        Df_Drho=-Kc*Ec
        DO i=1,3
          DO j=1,3
            IF (i.EQ.j) THEN
              Df_Dsig(i,j)=Temp(i,j)+alpha
            ELSE
              Df_Dsig(i,j)=Temp(i,j)
            END IF
          END DO
        END DO
	  END IF
C
      RETURN
      END
C
       SUBROUTINE Dg_Dsigma(n,w,sigmaT,dg_dsig)
C=======================================================================
C                                                                      *
C                           DF_DSIG_Drho                               *
C                                                                      *
C=======================================================================
       INCLUDE 'ABA_PARAM.INC'
C
       DIMENSION n(3),sigmaT(3,3),Dg_Dsig(3,3),M(3,3,3,3)
       INTEGER E(3,3)
       DOUBLE PRECISION n,sigmaT,Temp,Dg_Dsig,M,E0,nu0,N_V0,
     & a_0,Kc,Ec,Alpha,Ko,Eo,c0,c1
       COMMON/parameters/E0,nu0,N_V0,a_0,Kc,Ec,Alpha,Ko,Eo,c0,c1
       DATA TPF /2.5D-1/
C
C
	  DO I=1,3
        DO J=1,3
          Dg_Dsig(I,J)=0.0D0
          E(I,J)=0.0D0
          DO K=1,3
            DO L=1,3
              M(I,J,K,L)=0.0D0
            END DO
          END DO
        END DO
      END DO
C
      CALL DIDENTITY_2(E)
      Call Ni_Aij_Nj(n,sigmaT,sigma_nn)
      IF (sigma_nn.GE.0) THEN
          DO i=1,3
            DO j=1,3
              DO k=1,3
                DO l=1,3
                M(i,j,k,l) = w*c0*n(i)*n(j)*n(k)*n(l)
     & +w*c1*((TPF*(n(i)*n(k)*E(j,l)+n(i)*n(l)*E(j,k)+E(i,k)*n(j)*n(l)
     & +E(i,l)*n(j)*n(k))-n(i)*n(j)*n(k)*n(l)))
                END DO
              END DO
            END DO
          END DO
      ELSE
          DO i=1,3
            DO j=1,3
              DO k=1,3
                DO l=1,3
                M(i,j,k,l) = w*c1*((TPF*(n(i)*n(k)*E(j,l)
     & +n(i)*n(l)*E(j,k)+E(i,k)*n(j)*n(l)+E(i,l)*n(j)*n(k))
     & -n(i)*n(j)*n(k)*n(l)))
                END DO
              END DO
            END DO
          END DO
      END IF
      CALL Aijkl_Bij(M,sigmaT,dg_dsig)
C
      RETURN
      END
C
       SUBROUTINE Dg_Dsig_Sum_F(ddrho,activated,n_cosines,n_weight,
     & sigmaT,g_sig_sum)
C=======================================================================
C                                                                      *
C                          Dg_Dsig_Sum                                 *
C                                                                      *
C=======================================================================
       INCLUDE 'ABA_PARAM.INC'
C
       DIMENSION ddrho(42),n_cosines(42,3),n_weight(42),
     & sigmaT(3,3),g_sig_sum(3,3),NN(3,3,3,3)
       INTEGER E(3,3),activated(42)
       DOUBLE PRECISION ddrho,n_cosines,n_weight,sigmaT,g_sig_sum,NN,
     & E0,nu0,N_V0,a_0,Kc,Ec,Alpha,Ko,Eo,c0,c1
       COMMON/parameters/E0,nu0,N_V0,a_0,Kc,Ec,Alpha,Ko,Eo,c0,c1
       DATA TPF/2.5D-1/
C
       DO I=1,3
         DO J=1,3
           g_sig_sum(I,J)=0.0D0
           E(I,J)=0
           DO K=1,3
             DO L=1,3
               NN(I,J,K,L)=0.0D0
             END DO
           END DO 
         END DO
       END DO
C 
       CALL DIDENTITY_2(E)
       DO i_point=1,42
         IF (activated(i_point).EQ.1) THEN
           DO i=1,3
             DO j=1,3
               DO k=1,3
                 DO l=1,3
                   Call Ni_Aij_Nj(n_cosines(i_point, 1:3),sigmaT
     & ,sigma_nn)
                   IF  (sigma_nn.GE.0) THEN
                         NN(i,j,k,l) = NN(i,j,k,l)+ddrho(i_point)
     & *n_weight(i_point)*c0*n_cosines(i_point,i)*n_cosines(i_point,j)
     & *n_cosines(i_point,k)*n_cosines(i_point,l)
                   END IF
                   NN(i,j,k,l) = NN(i,j,k,l)+ddrho(i_point)
     & *n_weight(i_point)*c1*(TPF*(n_cosines(i_point,i)
     & *n_cosines(i_point,k)*E(j,l)+n_cosines(i_point,i)
     & *n_cosines(i_point,l)*E(j,k)+E(i,k)*n_cosines(i_point,j)
     & *n_cosines(i_point,l)+E(i,l)*n_cosines(i_point,j)
     & *n_cosines(i_point,k))-n_cosines(i_point,i)*n_cosines(i_point,j)
     & *n_cosines(i_point,k)*n_cosines(i_point,l))
                 END DO
               END DO
             END DO
           END DO
         END IF
       END DO
C
      CALL Aijkl_Bij(NN,sigmaT,g_sig_sum)
C
      RETURN
      END
C
       SUBROUTINE Inc_In_strain(ddrho,drho,n_cosines,n_weight,sigmaT,
     & dsig,activated,Inc_instrain)
C=====================================================================
C                                                                     *
C                          Inc_In_strain                              *
C                                                                     *
C=====================================================================
       INCLUDE 'ABA_PARAM.INC'
C
       DIMENSION ddrho(42),drho(42),n_cosines(42,3),n_weight(42),
     & sigmaT(3,3),dsig(3,3),Inc_instrain(3,3),
     & NN1(3,3,3,3),NN2(3,3,3,3),In_str1(3,3),In_str2(3,3)
       INTEGER E(3,3),activated(42)
       DOUBLE PRECISION ddrho,drho,n_cosines,n_weight,sigmaT,
     & dsig,Inc_instrain,NN1,NN2,In_str1,In_str2,E0,nu0,N_V0,
     & a_0,Kc,Ec,Alpha,Ko,Eo,c0,c1
       COMMON/parameters/E0,nu0,N_V0,a_0,Kc,Ec,Alpha,Ko,Eo,c0,c1
       DATA TPF/2.5D-1/
C
       DO I=1,3
         DO J=1,3
           E(I,J)=0.0D0
           In_str1(I,J)=0.0D0
           In_str2(I,J)=0.0D0
           Inc_instrain(I,J)=0.0D0
           DO K=1,3
             DO L=1,3
               NN1(I,J,K,L)=0.0D0
               NN2(I,J,K,L)=0.0D0
             END DO
           END DO 
         END DO
       END DO
       CALL DIDENTITY_2(E)
C
       DO i_point=1,42
         IF (activated(i_point).EQ.1) THEN
           DO i=1,3
             DO j=1,3
               DO k=1,3
                 DO l=1,3
                   Call Ni_Aij_Nj(n_cosines(i_point, 1:3),
     & sigmaT, sigma_nn)
                   IF  (sigma_nn.GE.0) THEN
                         NN1(i,j,k,l) = NN1(i,j,k,l)+drho(i_point)
     & *n_weight(i_point)*c0*n_cosines(i_point,i)*n_cosines(i_point,j)
     & *n_cosines(i_point,k)*n_cosines(i_point,l)
                         NN2(i,j,k,l) = NN2(i,j,k,l)+ddrho(i_point)
     & *n_weight(i_point)*c0*n_cosines(i_point,i)*n_cosines(i_point,j)
     & *n_cosines(i_point,k)*n_cosines(i_point,l)
                   END IF
                   NN1(i,j,k,l) = NN1(i,j,k,l)+drho(i_point)
     & *n_weight(i_point)*c1*(TPF*(n_cosines(i_point,i)
     & *n_cosines(i_point,k)*E(j,l)+n_cosines(i_point,i)
     & *n_cosines(i_point,l)*E(j,k)+E(i,k)*n_cosines(i_point,j)
     & *n_cosines(i_point,l)+E(i,l)*n_cosines(i_point,j)
     & *n_cosines(i_point,k))-n_cosines(i_point,i)*n_cosines(i_point,j)
     & *n_cosines(i_point,k)*n_cosines(i_point,l))
                   NN2(i,j,k,l) = NN2(i,j,k,l)+ddrho(i_point)
     & *n_weight(i_point)*c1*(TPF*(n_cosines(i_point,i)
     & *n_cosines(i_point,k)*E(j,l)+n_cosines(i_point,i)
     & *n_cosines(i_point,l)*E(j,k)+E(i,k)*n_cosines(i_point,j)
     & *n_cosines(i_point,l)+E(i,l)*n_cosines(i_point,j)
     & *n_cosines(i_point,k))-n_cosines(i_point,i)*n_cosines(i_point,j)
     & *n_cosines(i_point,k)*n_cosines(i_point,l))
                 END DO
               END DO
             END DO
           END DO
         END IF
       END DO  
C
      CALL Aijkl_Bij(NN1,dsig,In_str1)
      CALL Aijkl_Bij(NN2,sigmaT,In_str2)
      Call Aij_plus_Bij(In_str1,In_str2,Inc_instrain)
C
      RETURN
      END
C
       SUBROUTINE fd_lam(rho, drho, sigmaT, R, n_cosines, 
     & n_weight,ddepsin,ddrho,dsig,activated,NACT)
C=======================================================================
C                                                                      *
C                           DLAMBDA                                    *
C                                                                      *
C=======================================================================
       INCLUDE 'ABA_PARAM.INC'
C
       INTEGER activated(42),NACT
       DIMENSION rho(42),drho(42),sigmaT(3,3),n_cosines(42,3),
     & n_weight(42),R(3,3),ddepsin(3,3),ddrho(42),
     & dsig(3,3),Temp(3,3),Fv(NACT),Coef(NACT,NACT),Df_Dsig(3,3),
     & Dg_Dsig(3,3),C_Dg_Dsig(3,3),Coef_Inv(NACT,NACT),
     & ddrho_temp(NACT),g_sig_sum(3,3),Residual(3,3),
     & Inc_instrain(3,3),matDz(3,3,3,3)
       DOUBLE PRECISION rho,drho,n_cosines,n_weight,sigmaT,
     & R,ddepsin,ddrho,dsig,Temp,Fv,Coef,Df_Dsig,
     & Dg_Dsig,C_Dg_Dsig,Coef_Inv,ddrho_temp,g_sig_sum,Residual,
     & Inc_instrain,matDz,E0,nu0,N_V0,
     & a_0,Kc,Ec,Alpha,Ko,Eo,c0,c1
       COMMON/parameters/E0,nu0,N_V0,a_0,Kc,Ec,Alpha,Ko,Eo,c0,c1
       DATA TWO /2.0D0/
C
      DO I=1,3
         DO J=1,3
           Temp(I,J)=0.0D0
           Df_Dsig(I,J)=0.0D0
           Dg_Dsig(I,J)=0.0D0
           C_Dg_Dsig(I,J)=0.0D0
           g_sig_sum(I,J)=0.0D0
           Residual(I,J)=0.0D0
           dsig(I,J)=0.0D0
           Inc_instrain(I,J)=0.0D0
           ddepsin(I,J)=0.0D0
           DO K=1,3
             DO L=1,3
               matDz(I,J,K,L)=0.0D0
             END DO
           END DO 
         END DO
      END DO
      DO I=1,42
          ddrho(I)=0.0D0
      End DO
C
      DO I=1,NACT
        ddrho_temp(I)=0.0D0
        Fv(I)=0.0D0
        DO J=1,NACT
          Coef(I,J)=0.0D0
          Coef_Inv(I,J)=0.0D0
        END DO
      END DO
C
      CALL matDO2(rho,drho,n_cosines,n_weight,sigmaT,matDz)
      CALL Aijkl_Bkl(matDz,R,Temp)
      K=0
      Do I=1,42
        IF (activated(I).EQ.1) THEN
          CALL DF_DSIG_Drho(n_cosines(I,:), n_weight(I), sigmaT,
     & rho(I),drho(I),Df_Dsig,Df_Drho,fd)
          CALL Aij_Bij(Df_Dsig,Temp,Df_Dsig_S_sig)
          K=K+1
          Fv(K)=fd-Df_Dsig_S_sig
          L=0
          Do J=1,42
            IF (activated(J).EQ.1) THEN
              CALL Dg_Dsigma(n_cosines(J,:),n_weight(J),sigmaT,Dg_Dsig)
              CALL Aijkl_Bkl(matDz,Dg_Dsig,C_Dg_Dsig)
              CALL Aij_Bij(Df_Dsig,C_dg_dsig,Coef_temp)
              L=L+1
              Coef(K,L)=2*Coef_temp
              IF (K.EQ.L) THEN
                Coef(K,L)=Coef(K,L)-Df_Drho
              END IF 
            END IF
          END DO
        END IF
      END DO
      CALL INVERSE(Coef,K,K,Coef_Inv)  
      CALL Aij_Bj(Coef_Inv,Fv,K,ddrho_temp)
      K=0
      DO I=1,42
        IF (activated(I).EQ.1) THEN
          K=K+1
          ddrho(I)=ddrho_temp(K)
        END IF
      END DO
      CALL Dg_Dsig_Sum_F(ddrho,activated,n_cosines,n_weight,
     & sigmaT,g_sig_sum)
      DO I=1,3
        DO J=1,3
          Residual(I,J)=-R(I,J)-TWO*g_sig_sum(I,J)
        END DO
      END DO
      Call Aijkl_Bkl(matDz,Residual,dsig)
      CALL Inc_In_strain(ddrho,drho,n_cosines,n_weight,
     & sigmaT,dsig,activated,Inc_instrain)
      Call Aij_plus_Bij(R,Inc_instrain,ddepsin)
C      
C      WRITE(7,*),'Depsin=',Depsin
C 	   WRITE(7,*),'DOMEGA=',DOMEGA
C      WRITE(7,*) 'DSIG=',DSIG      
C
      RETURN
      END
C