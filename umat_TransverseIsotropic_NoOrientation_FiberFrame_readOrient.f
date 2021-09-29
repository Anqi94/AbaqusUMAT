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
      REAL*8 ::E11, E22, E33, v12, v13, v23, v21, v31, v32, G12, G13, G23, gamma, 
     1 D1111, D2222, D3333, D1122, D1133, D2233, D1212, D1313, D2323, 
     2 normf1, normf2, normf3
      real, dimension(3) :: e1,e2,e3,en1,en2,en3,f1,f2,f3,fn1,fn2,fn3,tmpfn2
	  dimension dstran_FF(ntens),
     1 stress0_FF(ntens),stress1_FF(ntens),stress1_NG(ntens),ddsdde_tmp(ntens,ntens)
	 
      REAL*8, DIMENSION(NTENS) :: STRESSV,DSTRANV,DSTRANV_rot,STRESSV_rot,STRESSV_rot_n
      REAL*8, DIMENSION(3,3) :: IDEN,STRESSM,DSTRANM,DSTRANM_rot,STRESSM_rot,STRESSM_rot_ini_tmp,
     1 STRESSM_rot_ini,STRESSM_rot_n,STRESSM_rot_n_rot,STRESSM_rot_tmp, Q,Q_tmp,Q_ini,PHI
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
	  do i = 1, 3
        do j = 1, 3
		  do k = 1, 3
		    do l = 1, 3
              C1(i,j,k,l)=0.0d0
			  C2(i,j,k,l)=0.0d0
		    end do
		  end do
        end do
      end do	  
	  do i = 1, 3
		do j = 1, 3
		  Q(i,j) = 0.0d0
		  Q_ini(i,j) = 0.0d0
		  DSTRANM_rot(i,j)= 0.0d0
		  STRESSM_rot_ini(i,j)= 0.0d0
		  STRESSM_rot_ini_tmp(i,j)= 0.0d0
		  STRESSM_rot(i,j)= 0.0d0
		  STRESSM_rot_n(i,j) = 0.0d0
		  STRESSM_rot_n_rot(i,j) = 0.0d0
		  STRESSM_rot_tmp(i,j) = 0.0d0
		end do
	  end do
	  DO i = 1,6
	    STRESSV_rot_n(i) = 0.0D0
	  end do
	  Q_ini(1,2) = -1.0d0
	  Q_ini(2,1) = 1.0d0
	  Q_ini(3,3) = 1.0d0

	  
      DO i=1,3
	    Q_tmp(i,1) = STATEV(i)
		Q_tmp(i,2) = STATEV(3+i)
		Q_tmp(i,3) = STATEV(6+i)
      END DO
!      WRITE(*,*) ('initial fiber direction Q') 
!	  WRITE(*,*) (Q_ini) 
!      WRITE(*,*) (Q_tmp) 
	  
	  do i = 1, 3
        do j = 1, 3
		  do n = 1, 3
            Q(i,j) = Q(i,j)+Q_tmp(i,n)*drot(n,j)
		  end do
        end do
      end do
C
      WRITE(*,*) ('Q_sum') 
      WRITE(*,*) (Q)
      DO i=1,3
	    STATEV(i) = Q(i,1)
		STATEV(3+i) = Q(i,2)
		STATEV(6+i) = Q(i,3)
      END DO
C---------------------------------------------------------------------------
C        Fiber direction initiation
C---------------------------------------------------------------------------     	 
      e1= (/1.0,0.0,0.0/)
      e2= (/0.0,1.0,0.0/)
      e3= (/0.0,0.0,1.0/)
C      f1= (/0.0,1.0,0.0/)
C      f2= (/-1.0,0.0,0.0/)
C      f3= (/0.0,0.0,1.0/) 
      f1= (/0.196117290802068,-0.980580444557435,0.0/)
      f2= (/0.980580444557435,0.196117290802068,0.0/)
      f3= (/0.0,0.0,1.0/) 
      do i = 1, ndi
        en1(i) = 0.0d0
        en2(i) = 0.0d0
        en3(i) = 0.0d0
        fn1(i) = 0.0d0
        fn2(i) = 0.0d0
        fn3(i) = 0.0d0
		tmpfn2(i) = 0.0d0
      end do

C---------------------------------------------------------------------------
C        Rotate initial base ei to eni
C---------------------------------------------------------------------------     
      do i = 1, ndi
        do j = 1, ndi
          en1(i) = en1(i) + Q(i,j) * e1(j)
          en2(i) = en2(i) + Q(i,j) * e2(j)
          en3(i) = en3(i) + Q(i,j) * e3(j)
        end do 
      end do
      WRITE(*,*) ('en1 = ')
      WRITE(*,*) (en1)
      WRITE(*,*) ('en2 = ')
      WRITE(*,*) (en2)		
      WRITE(*,*) ('en3 = ')
      WRITE(*,*) (en3)		  
C---------------------------------------------------------------------------
C        Rotate fiber base fi = F * ei
C---------------------------------------------------------------------------   

      do i = 1, ndi
        do j = 1, ndi
          fn1(i) = fn1(i) + DFGRD1(i,j) * f1(j)
          tmpfn2(i) = tmpfn2(i) + DFGRD1(i,j) * f2(j)
        end do 
      end do 
	  
	  normf1 = sqrt(fn1(1)**2.0 + fn1(2)**2.0 + fn1(3)**2.0 )
	  fn1(1) = fn1(1) / normf1
	  fn1(2) = fn1(2) / normf1
	  fn1(3) = fn1(3) / normf1
	  
	  
	  fn3(1) = fn1(2)*tmpfn2(3)-fn1(3)*tmpfn2(2)
	  fn3(2) = fn1(3)*tmpfn2(1)-fn1(1)*tmpfn2(3)
	  fn3(3) = fn1(1)*tmpfn2(2)-fn1(2)*tmpfn2(1)
	  
	  normf3 = sqrt(fn3(1)**2.0 + fn3(2)**2.0 + fn3(3)**2.0 )
	  fn3(1) = fn3(1) / normf3
	  fn3(2) = fn3(2) / normf3
	  fn3(3) = fn3(3) / normf3
	  
	  fn2(1) = fn1(2)*fn3(3)-fn1(3)*fn3(2)
	  fn2(2) = fn1(3)*fn3(1)-fn1(1)*fn3(3)
	  fn2(3) = fn1(1)*fn3(2)-fn1(2)*fn3(1)
	  
	  normf2 = sqrt(fn2(1)**2.0 + fn2(2)**2.0 + fn2(3)**2.0 )	
	  fn2(1) = -fn2(1) / normf2  
	  fn2(2) = -fn2(2) / normf2	  
	  fn2(3) = -fn2(3) / normf2
	  
	  
      WRITE(*,*) ('fn1 = ')
      WRITE(*,*) (fn1)
      WRITE(*,*) ('fn2 = ')
      WRITE(*,*) (fn2)		
      WRITE(*,*) ('fn3 = ')
      WRITE(*,*) (fn3)	
	  
!	  THETA(1,1) = fn1(1)*en1(1)+fn1(2)*en1(2)+fn1(3)*en1(3)
!	  THETA(1,2) = fn1(1)*en2(1)+fn1(2)*en2(2)+fn1(3)*en2(3)
!	  THETA(1,3) = fn1(1)*en3(1)+fn1(2)*en3(2)+fn1(3)*en3(3)
!	  THETA(2,1) = fn2(1)*en1(1)+fn2(2)*en1(2)+fn2(3)*en1(3)
!	  THETA(2,2) = fn2(1)*en2(1)+fn2(2)*en2(2)+fn2(3)*en2(3)
!	  THETA(2,3) = fn2(1)*en3(1)+fn2(2)*en3(2)+fn2(3)*en3(3)
!	  THETA(3,1) = fn3(1)*en1(1)+fn3(2)*en1(2)+fn3(3)*en1(3)
!	  THETA(3,2) = fn3(1)*en2(1)+fn3(2)*en2(2)+fn3(3)*en2(3)
!	  THETA(3,3) = fn3(1)*en3(1)+fn3(2)*en3(2)+fn3(3)*en3(3)	  
!     WRITE(*,*) ('THETA = ')
!      WRITE(*,*) (THETA)	

      PHI(1,1) = fn1(1)
	  PHI(2,1) = fn1(2)
	  PHI(3,1) = fn1(3)
      PHI(1,2) = fn2(1)
	  PHI(2,2) = fn2(2)
	  PHI(3,2) = fn2(3)
      PHI(1,3) = fn3(1)
	  PHI(2,3) = fn3(2)
	  PHI(3,3) = fn3(3)	  
      WRITE(*,*) ('PHI = ')
      WRITE(*,*) (PHI)	
C---------------------------------------------------------------------------
C        Rotate fiber base fi = F * ei
C---------------------------------------------------------------------------   

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
      do m = 1, 3
        do n = 1, 3
          do i = 1, 3
            do j = 1, 3
			  DSTRANM_rot(m,n) = DSTRANM_rot(m,n) + PHI(j,n)*DSTRANM(i,j)*PHI(i,m)
			  STRESSM_rot(m,n) = STRESSM_rot(m,n) + PHI(j,n)*STRESSM(i,j)*PHI(i,m)
            end do
          end do
        end do
      end do
	  
      WRITE(*,*) ('DSTRANM_rot') 
      WRITE(*,*) (DSTRANM_rot)
      WRITE(*,*) ('STRESSM_rot') 
      WRITE(*,*) (STRESSM_rot)	 
	  
      DSTRANV_rot(1)=DSTRANM_rot(1,1)
      DSTRANV_rot(2)=DSTRANM_rot(2,2)
      DSTRANV_rot(3)=DSTRANM_rot(3,3)
      DSTRANV_rot(4)=DSTRANM_rot(1,2)+DSTRANM_rot(2,1)
      DSTRANV_rot(5)=DSTRANM_rot(1,3)+DSTRANM_rot(3,1)
      DSTRANV_rot(6)=DSTRANM_rot(2,3)+DSTRANM_rot(3,2)

      STRESSV_rot_n(1)=STRESSM_rot(1,1)
      STRESSV_rot_n(2)=STRESSM_rot(2,2)
      STRESSV_rot_n(3)=STRESSM_rot(3,3)
      STRESSV_rot_n(4)=STRESSM_rot(1,2)
      STRESSV_rot_n(5)=STRESSM_rot(1,3)
      STRESSV_rot_n(6)=STRESSM_rot(2,3)
	  
C material properties
      E11 = props(1)      ! Young's modulus E11 - fiber longtitude direction
      E22 = props(2)      ! Young's modulus E22 - fiber transverse direction
      E33 = props(3)      ! Young's modulus E33 - fiber transverse direction
      v12 = props(4)      ! Poisson's ratio v12 - out of plane possion ratio 
      v13 = props(5)      ! Poisson's ratio v13 - out of plane possion ratio 
      v23 = props(6)      ! Poisson's ratio v23 - in plane possion ratio 
      G12 = props(7)      ! Shear Stiffness G12 - out of plane direction
      G13 = props(8)      ! Shear Stiffness G13 - out of plane direction
      G23 = props(9)      ! Shear Stiffness G23 - in plane direction


C Tangent matrix C - parameters
      v32 = v23 
      v31 = v13*E33/E11
      v21 = v31
      gamma = 1.0d0/(1-v12*v21-v23*v32-v13*v31-2*v21*v32*v13)


      D1111 = E11*(1.0d0-v23*v32)*gamma
      D2222 = E22*(1.0d0-v13*v31)*gamma 
      D3333 = E33*(1.0d0-v12*v21)*gamma 
      D1122 = E11*(v21+v31*v23)*gamma       ! = E22*(v12+v13*v23)*gamma
      D1133 = E11*(v31+v21*v32)*gamma       ! = E33*(v13+v12*v23)*gamma
      D2233 = E22*(v32+v12*v31)*gamma       ! = E33*(v23+v12*v13)*gamma
      D1212 = G12 
      D1313 = G13 
      D2323 = G23 
	  
C Stiffness matrix
      do i = 1, ntens
        do j = 1, ntens
          ddsdde(i,j) = 0.0d0
        end do
      end do

      ddsdde(1, 1) =  D1111 
      ddsdde(2, 2) =  D2222
      ddsdde(3, 3) =  D3333 

      ddsdde(1, 2) =  D1122 
      ddsdde(2, 1) =  D1122 
      ddsdde(1, 3) =  D1133 
      ddsdde(3, 1) =  D1133 
      ddsdde(2, 3) =  D2233 
      ddsdde(3, 2) =  D2233 

C Shear contribution
      ddsdde(4, 4) =  D1212 
      ddsdde(5, 5) =  D1313 
      ddsdde(6, 6) =  D2323 


! for debugging see the vumat version
      WRITE(*,*)
      WRITE(*,*) (ddsdde)   

C Stress increment evaluation
	  
      do i = 1, ntens
        do j = 1, ntens
          STRESSV_rot_n(i) = STRESSV_rot_n(i) + ddsdde(i,j) * DSTRANV_rot(j)
        end do 
      end do 
	  WRITE(*,*) (' updated stress in the fiber frame  ') 
      WRITE(*,*) (STRESSV_rot_n)
	  
      STRESSM_rot_n(1,1)=STRESSV_rot_n(1)
      STRESSM_rot_n(2,2)=STRESSV_rot_n(2)
      STRESSM_rot_n(3,3)=STRESSV_rot_n(3)
      STRESSM_rot_n(1,2)=STRESSV_rot_n(4)
      STRESSM_rot_n(2,1)=STRESSV_rot_n(4)
      STRESSM_rot_n(1,3)=STRESSV_rot_n(5)
      STRESSM_rot_n(3,1)=STRESSV_rot_n(5)
      STRESSM_rot_n(2,3)=STRESSV_rot_n(6)
      STRESSM_rot_n(3,2)=STRESSV_rot_n(6)
	  
      do m = 1, 3
        do n = 1, 3
          do i = 1, 3
            do j = 1, 3
			  STRESSM_rot_n_rot(m,n) = STRESSM_rot_n_rot(m,n) + PHI(n,j)*STRESSM_rot_n(i,j)*PHI(m,i)
            end do
          end do
        end do
      end do

      STRESS(1)=STRESSM_rot_n_rot(1,1)
      STRESS(2)=STRESSM_rot_n_rot(2,2)
      STRESS(3)=STRESSM_rot_n_rot(3,3)
      STRESS(4)=STRESSM_rot_n_rot(1,2)
      STRESS(5)=STRESSM_rot_n_rot(1,3)
      STRESS(6)=STRESSM_rot_n_rot(2,3)	
	  
	  WRITE(*,*) (' final stress ') 
      WRITE(*,*) (STRESS)	  
C
C
C
     
    
      C1(1,1,1,1)= ddsdde(1, 1);
      C1(2,2,2,2)= ddsdde(2, 2);
      C1(3,3,3,3)= ddsdde(3, 3);
      C1(1,1,2,2)= ddsdde(1, 2);
      C1(2,2,1,1)= ddsdde(2, 1);
      C1(1,1,3,3)= ddsdde(1, 3);
      C1(3,3,1,1)= ddsdde(3, 1);
      C1(2,2,3,3)= ddsdde(2, 3);
      C1(3,3,2,2)= ddsdde(3, 2);
    
      C1(1,2,2,1)= 2*ddsdde(4, 4);
      C1(1,3,3,1)= 2*ddsdde(5, 5);
      C1(2,3,3,2)= 2*ddsdde(6, 6);
      C1(2,1,1,2)= 2*ddsdde(4, 4);
      C1(3,1,1,3)= 2*ddsdde(5, 5);
      C1(3,2,2,3)= 2*ddsdde(6, 6);

      do i = 1, 3
        do j = 1, 3
          do k = 1, 3
            do l = 1, 3
			  do m = 1, 3
				do n = 1, 3
				  do o = 1, 3
					do p = 1, 3
					  C2(i,j,k,l) = C2(i,j,k,l) + PHI(i,m)*PHI(j,n)*PHI(k,o)*PHI(l,p)*C1(m,n,o,p)
					end do
				  end do
				end do
			  end do
            end do
          end do
        end do
      end do
      ddsdde(1, 1) =  C2(1,1,1,1) 
      ddsdde(2, 2) =  C2(2,2,2,2)
      ddsdde(3, 3) =  C2(3,3,3,3)

      ddsdde(1, 2) =  C2(1,1,2,2)
      ddsdde(2, 1) =  C2(2,2,1,1)
      ddsdde(1, 3) =  C2(1,1,3,3)
      ddsdde(3, 1) =  C2(3,3,1,1)
      ddsdde(2, 3) =  C2(2,2,3,3)
      ddsdde(3, 2) =  C2(3,3,2,2)

C Shear contribution
      ddsdde(4, 4) =  C2(1,2,2,1)/2.0d0 
      ddsdde(5, 5) =  C2(1,3,3,1)/2.0d0 
      ddsdde(6, 6) =  C2(2,3,3,2)/2.0d0 	  
      WRITE(*,*)
      WRITE(*,*) (ddsdde)  

	  
      RETURN
      END
