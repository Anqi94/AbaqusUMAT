	  function func(v) result(w)
		real, intent (in) :: v ! input
		real              :: w ! output
        WRITE(*,*) ('v =') 
        WRITE(*,*) (v)
		w = v**2 + v**3
!		w(1) = v(1)**2 + v(1)**3
!		w(2) = v(2)**2 + v(2)**3
!		w(3) = v(3)**2 + v(3)**3
		WRITE(*,*) ('w =') 
        WRITE(*,*) (w)
	  end function	  
	! ======================================================
	! Function returning an array: matrix * vector
	! ======================================================
	  function MatVecMult(A, v) result (w)

	   real, dimension(:,:), intent(in) :: A
	   real, dimension(:), intent(in)   :: v
	   real, dimension( SIZE(A,1) )     :: w   !! Define result using input param

	   integer :: i, j
	   integer :: N

	   N = size(v)

	   w = 0.0       !! clear whole vector
	   DO i = 1, N
		  w = w + v(i) * A( :, i )
	   END DO
	  end function
	! ======================================================
	! Function returning the inverse of the matrix: B =inv(A)
	! ======================================================
	  function inv3(G) result (B)

	   real*8, dimension(:,:), intent(in) :: G
	   real*8, dimension( 3,3 )     :: B   !! Define result using input param

      real*8,dimension(3,3) :: COG(3,3),ADJ(3,3)
      real :: det_G
	  integer :: i, j
	  
c
c	cofactors and determinant of g
      COG(1,1)=G(2,2)*G(3,3)-G(2,3)*G(3,2)
      COG(2,1)=G(1,3)*G(3,2)-G(1,2)*G(3,3)
      COG(3,1)=G(1,2)*G(2,3)-G(1,3)*G(2,2)
      COG(1,2)=G(2,3)*G(3,1)-G(2,1)*G(3,3)
      COG(2,2)=G(1,1)*G(3,3)-G(1,3)*G(3,1)
      COG(3,2)=G(1,3)*G(2,1)-G(1,1)*G(2,3)
      COG(1,3)=G(2,1)*G(3,2)-G(2,2)*G(3,1)
      COG(2,3)=G(1,2)*G(3,1)-G(1,1)*G(3,2)
      COG(3,3)=G(1,1)*G(2,2)-G(1,2)*G(2,1)
  
	  
      det_G=G(1,1)*COG(1,1)+G(2,1)*COG(2,1)+G(3,1)*COG(3,1)
C	  WRITE(*,*) ('determinant =') 
C      WRITE(*,*) (det_G)
c
C calculate Adjugate/Adjoint - transpose of COG
c
      ADJ = transpose(COG)
C	  WRITE(*,*) ('Adjoint =') 
C      WRITE(*,*) (ADJ)
C
C Multiply by 1/det
C
      B(1,1) = (1.0d0/det_G)*ADJ(1,1)
	  B(2,1) = (1.0d0/det_G)*ADJ(2,1)
	  B(3,1) = (1.0d0/det_G)*ADJ(3,1)
	  B(1,2) = (1.0d0/det_G)*ADJ(1,2)
	  B(2,2) = (1.0d0/det_G)*ADJ(2,2)
	  B(3,2) = (1.0d0/det_G)*ADJ(3,2)
	  B(1,3) = (1.0d0/det_G)*ADJ(1,3)
	  B(2,3) = (1.0d0/det_G)*ADJ(2,3)
	  B(3,3) = (1.0d0/det_G)*ADJ(3,3)
!      WRITE(*,*) ('(1.0d0/det_G)*ADJ(1,1)') 
!	  WRITE(*,*) ((1.0d0/det_G)*ADJ(1,1)) 
!	  WRITE(*,*) (B(1,1)) 
c	done
	  end function	  