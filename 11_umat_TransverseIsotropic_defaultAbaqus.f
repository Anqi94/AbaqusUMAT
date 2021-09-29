c ======================================================================
c User Subroutine UMAT for Abaqus linear elastic material
c ======================================================================
      subroutine umat(stress,statev,ddsdde,sse,spd,scd,
     1 rpl,ddsddt,drplde,drpldt,
     2 stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname,
     3 ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
     4 celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc)
C
      include 'aba_param.inc'
C
      character*80 cmname
      dimension stress(ntens),statev(nstatv),
     1 ddsdde(ntens,ntens),
     2 ddsddt(ntens),drplde(ntens),
     3 stran(ntens),dstran(ntens),time(2),predef(1),dpred(1),
     4 props(nprops),coords(3),drot(3,3),dfgrd0(3,3),dfgrd1(3,3)
C
      integer i, j
      real E11, E22, E33, v12, v13, v23, v21, v31, v32, G12, G13, G23, gamma, 
     1 D1111, D2222, D3333, D1122, D1133, D2233, D1212, D1313, D2323
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
          stress(i) = stress(i) + ddsdde(i,j) * dstran(j)
        end do 
      end do 
C
      WRITE(*,*) ('stress 1 =') 
      WRITE(*,*) (stress) 
      return
      end