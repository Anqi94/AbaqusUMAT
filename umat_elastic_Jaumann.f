C==========================================================================================
C ZERO-GRADED HYPOELASTIC CONSTITUTIVE MODEL BASED ON THE JAUMANN RATE OF THE CAUCHY STRESS 
C==========================================================================================
C PROPS(1) -> YOUNG'S MODULUS
C PROPS(2) -> POISSON'S RATIO
C
C Ee      = YOUNG'S MODULUS
C Nu      = POISSON'S RATIO
C Lambda  = LAME'S CONSTANT
C Mu      = LAME'S CONSTANT
C IDEN    = KRONECKER DELTA
C STRESS  = UMAT-SUBROUTINE STRESS ARRAY
C STRESSV = VECTOR VERSION OF "STRESS"
C STRESSM = MATRIX VERSION OF "STRESS"
C DSTRAN  = UMAT-SUBROUTINE INCREMENTAL STRAIN ARRAY
C DSTRANV = VECTOR VERSION OF "DSTRAN"
C DSTRANM = MATRIX VERSION OF "DSTRAN"
C C1      = FOURTH-ORDER ELASTICITY TENSOR
C DDSDDE  = UMAT-SUBROUTINE CONSISTENT JACOBIAN ARRAY
C C2      = FOURTH-ORDER TENSOR VERSION OF "DDSDDE"
C==========================================================================================
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),
     2 DDSDDT(NTENS),DRPLDE(NTENS),
     3 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     4 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3),
     5 JSTEP(4)
C
      INTEGER :: i,j,k,l    
      REAL*8 :: Ee,Nu,Lambda,Mu
      REAL*8, DIMENSION(NTENS) :: STRESSV,DSTRANV
      REAL*8, DIMENSION(3,3) :: IDEN,STRESSM,DSTRANM
      REAL*8, DIMENSION(3,3,3,3) :: C1,C2
      PARAMETER(ZERO=0.D0,ONE=1.D0,TWO=2.D0)
	  WRITE(*,*) ('Element Number NOEL=') 	
      WRITE(*,*) (noel) 
	  WRITE(*,*) ('Step number KSTEP=') 	
      WRITE(*,*) (kstep) 
	  WRITE(*,*) ('Increment number kinc=') 	
      WRITE(*,*) (kinc)
      WRITE(*,*) ('stress 0 =') 
      WRITE(*,*) (stress) 
      WRITE(*,*) ('strain=') 
      WRITE(*,*) (stran)  
      WRITE(*,*) ('dstrain=') 
      WRITE(*,*) (dstran)   
      WRITE(*,*) ('deformation gradient 0 =') 
      WRITE(*,*) (dfgrd0)   
      WRITE(*,*) ('deformation gradient 1 =') 
      WRITE(*,*) (dfgrd1)   
      WRITE(*,*) ('incremental rotation =') 
      WRITE(*,*) (drot) 
C
      Ee=PROPS(1)
      Nu=PROPS(2)
      Lambda=Nu*Ee/((ONE+Nu)*(ONE-TWO*Nu))
      Mu=Ee/(TWO*(ONE+Nu))
C
      DO i=1,3
      DO j=1,3
      IF (i.EQ.j) THEN
      IDEN(i,j)=ONE
      ELSE
      IDEN(i,j)=ZERO
      ENDIF 
      END DO
      END DO
C
      STRESSV=STRESS
      STRESSM(1,1)=STRESSV(1)
      STRESSM(2,2)=STRESSV(2)
      STRESSM(3,3)=STRESSV(3)
      STRESSM(1,2)=STRESSV(4)
      STRESSM(2,1)=STRESSV(4)
      STRESSM(1,3)=STRESSV(5)
      STRESSM(3,1)=STRESSV(5)
      STRESSM(2,3)=STRESSV(6)
      STRESSM(3,2)=STRESSV(6)
C
      DSTRANV=DSTRAN      
      DSTRANM(1,1)=DSTRANV(1)
      DSTRANM(2,2)=DSTRANV(2)
      DSTRANM(3,3)=DSTRANV(3)
      DSTRANM(1,2)=ONE/TWO*DSTRANV(4)
      DSTRANM(2,1)=ONE/TWO*DSTRANV(4)
      DSTRANM(1,3)=ONE/TWO*DSTRANV(5)
      DSTRANM(3,1)=ONE/TWO*DSTRANV(5)
      DSTRANM(2,3)=ONE/TWO*DSTRANV(6)
      DSTRANM(3,2)=ONE/TWO*DSTRANV(6)
C
      DO i=1,3
      DO j=1,3
      DO k=1,3
      DO l=1,3                      
      C1(i,j,k,l)=
     1 Lambda*IDEN(i,j)*IDEN(k,l)+
     2 Mu*IDEN(i,k)*IDEN(j,l)+
     3 Mu*IDEN(i,l)*IDEN(j,k)
      END DO
      END DO
      END DO
      END DO
C
      DO i=1,3
      DO j=1,3
      DO k=1,3
      DO l=1,3
      STRESSM(i,j)=STRESSM(i,j)+C1(i,j,k,l)*DSTRANM(k,l)
      END DO
      END DO
      END DO
      END DO
C
      DO i=1,3
      DO j=1,3
      DO k=1,3
      DO l=1,3
      C2(i,j,k,l)=C1(i,j,k,l)+STRESSM(i,j)*IDEN(k,l)
      END DO
      END DO
      END DO
      END DO
C
      STRESS(1)=STRESSM(1,1)
      STRESS(2)=STRESSM(2,2)
      STRESS(3)=STRESSM(3,3)
      STRESS(4)=STRESSM(1,2)
      STRESS(5)=STRESSM(1,3)
      STRESS(6)=STRESSM(2,3)
C
      DO i=1,3
      DDSDDE(i,1)=C2(i,i,1,1)
      DDSDDE(i,2)=C2(i,i,2,2)
      DDSDDE(i,3)=C2(i,i,3,3)
      DDSDDE(i,4)=C2(i,i,1,2)
      DDSDDE(i,5)=C2(i,i,1,3)
      DDSDDE(i,6)=C2(i,i,2,3)
      END DO
      DDSDDE(4,1)=C2(1,2,1,1)
      DDSDDE(4,2)=C2(1,2,2,2)
      DDSDDE(4,3)=C2(1,2,3,3)
      DDSDDE(4,4)=C2(1,2,1,2)
      DDSDDE(4,5)=C2(1,2,1,3)
      DDSDDE(4,6)=C2(1,2,2,3)
      DDSDDE(5,1)=C2(1,3,1,1)
      DDSDDE(5,2)=C2(1,3,2,2)
      DDSDDE(5,3)=C2(1,3,3,3)
      DDSDDE(5,4)=C2(1,3,1,2)
      DDSDDE(5,5)=C2(1,3,1,3)
      DDSDDE(5,6)=C2(1,3,2,3)
      DDSDDE(6,1)=C2(2,3,1,1)
      DDSDDE(6,2)=C2(2,3,2,2)
      DDSDDE(6,3)=C2(2,3,3,3)
      DDSDDE(6,4)=C2(2,3,1,2)
      DDSDDE(6,5)=C2(2,3,1,3)
      DDSDDE(6,6)=C2(2,3,2,3)
      WRITE(*,*)
      WRITE(*,*) (ddsdde)  
      WRITE(*,*)
      WRITE(*,*) (stress)
	  
      RETURN
      END
