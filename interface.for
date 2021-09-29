	  interface
	  
		function func(v) result (w)
	     real, intent(in)   :: v
	     real               :: w   !! Define result using input param
		end function
		
		function MatVecMult(A, v) result (w)
	     real, dimension(:,:), intent(in) :: A
	     real, dimension(:), intent(in)   :: v
	     real, dimension( SIZE(A,1) )     :: w   !! Define result using input param
		end function
		
	    function inv3(G) result (B)
	     real*8, dimension(:,:), intent(in) :: G
	     real*8, dimension( 3,3 )       :: B   !! Define result using input param	
        end function
		
	  end interface 