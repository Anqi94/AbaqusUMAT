C     Subroutine to canculate the interse of a 3x3 matrix
C  
c	DATE:		2021/09/19  - by Anqi Li a.li@tue.nl
c********1*********2*********3*********4*********5*********6*********7**
c
c	inv_G = inv3(G)
c
c	calculate and return the inversed 3x3 matrix  G
c
c	INPUTS
c	G = 3x3 general matrix
c
c	OUTPUTS
c	inv_G = 3x3 inversed matrix of G
C
C	PRECISION:	single
C	COMMONS:	none
C	CALLS:		none
C	FUNCTIONS:	ABS, SQRT
C********1*********2*********3*********4*********5*********6*********7**
	  SUBROUTINE inv3(G,GINV)
c
c	declarations
      INCLUDE 'ABA_PARAM.INC'
	  
      dimension G(3,3),GINV(3,3),COG(3,3),ADJ(3,3)
      real det_G
	  
	  integer i, j
	  WRITE(*,*) ('inv_G =') 
      WRITE(*,*) (GINV)
	  WRITE(*,*) ('G =') 
      WRITE(*,*) (G)	  
	  
	  do i = 1, 3
        do j = 1, 3
          GINV(i,j) = 0.0d0
        end do
      end do	
	  WRITE(*,*) ('inv_G =') 
      WRITE(*,*) (GINV)
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
C      WRITE(*,*) ('matrix of cofactors =') 
C      WRITE(*,*) (COG) 	  
	  
      det_G=G(1,1)*COG(1,1)+G(2,1)*COG(2,1)+G(3,1)*COG(3,1)
C	  WRITE(*,*) ('determinant =') 
C      WRITE(*,*) (det_G)
c
C calculate Adjugate/Adjoint - transpose of COG
c
      ADJ = transpose(COG)
	  WRITE(*,*) ('Adjoint =') 
      WRITE(*,*) (ADJ)
C
C Multiply by 1/det
C
      GINV(1,1) = (1.0d0/det_G)*ADJ(1,1)
	  GINV(2,1) = (1.0d0/det_G)*ADJ(2,1)
	  GINV(3,1) = (1.0d0/det_G)*ADJ(3,1)
	  GINV(1,2) = (1.0d0/det_G)*ADJ(1,2)
	  GINV(2,2) = (1.0d0/det_G)*ADJ(2,2)
	  GINV(3,2) = (1.0d0/det_G)*ADJ(3,2)
	  GINV(1,3) = (1.0d0/det_G)*ADJ(1,3)
	  GINV(2,3) = (1.0d0/det_G)*ADJ(2,3)
	  GINV(3,3) = (1.0d0/det_G)*ADJ(3,3)
      WRITE(*,*) ('(1.0d0/det_G)*ADJ(1,1)') 
	  WRITE(*,*) ((1.0d0/det_G)*ADJ(1,1)) 
	  WRITE(*,*) (GINV(1,1)) 
c	done
      RETURN
      END