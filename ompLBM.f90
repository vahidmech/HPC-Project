!	 ========================================================
!	 Constants for simulation setup
!	 ========================================================
MODULE simParam
	integer, parameter:: IMAX = 2000
	integer, parameter:: JMAX = 400
	integer, parameter:: Kp = 8
	integer, parameter:: niter = 20000
	integer, parameter:: chap = 1000
	integer, parameter:: nt = 16

	double precision, parameter:: H = 2.d0
	double precision, parameter:: L = 10.d0
	double precision, parameter:: nu= 4.5d-03
	double precision, parameter:: P0= 1.d0
	double precision, parameter:: Re= 20.d0
	double precision, parameter:: Uave= Re*nu/H
	double precision, parameter:: DX = 1.d0/dble(IMAX-1)*L
	double precision, parameter:: DT = DX
	
 	double precision :: BB(IMAX,JMAX,0:Kp)
	double precision :: Uexact(JMAX)
	double precision :: X(IMAX), Y(JMAX)
CONTAINS
!	integer, parameter:: nt = 16
! 	double precision :: BB(IMAX,JMAX,0:Kp)
	SUBROUTINE ExactSolution
		INTEGER :: I, J

		DO J=1,JMAX
			Y(J)=dble(J-1)/dble(JMAX-1)*H
			Uexact(J) = 1.5d0*Uave*( 1.d0-ABS(2.d0*Y(J)/H-1.d0)**2 )
		END DO
		DO I=1,IMAX
			X(I)=dble(I-1)/dble(IMAX-1)*L
		END DO
		
		Uexact(:) = 1.5d0*Uave*( 1.d0-ABS(2.d0*Y(:)/H-1.d0)**2 )
	
	END SUBROUTINE ExactSolution
END MODULE simParam

!	 ========================================================
!	 The main program, implementing a flow in a channel
!	 ========================================================

PROGRAM LBM

USE simParam
use omp_lib
IMPLICIT NONE

	double precision :: U(IMAX,JMAX), V(IMAX,JMAX), P(IMAX,JMAX)
	double precision :: DELTA1(50000)
	double precision :: F(IMAX,JMAX,0:Kp), Feq(IMAX,JMAX,0:Kp)
	double precision :: t1, t0, time
	integer          :: I, J, K, Iter
	call omp_set_num_threads(nt)

	PRINT*, DT, Uave
	CALL ExactSolution

!------- Initial Condition
	P = 1.d0
	U = Uave
	V = 0.d0
	CALL Fequi(F, U, V, P)

!------- Main Loop
	t0 = omp_get_wtime()
	DO Iter=1,niter
		!$omp parallel
		CALL Fequi(Feq, U, V, P)
		CALL StreamCollision(F, Feq)
		CALL BC(F, U, V, P)
		CALL MacroVar(F, U, V, P)
		!$omp end parallel
	
!		DELTA1(Iter) = MAXVAL(ABS( U(IMAX/2,:)-Uexact ))
!		IF ( MOD(Iter,chap)==0 ) PRINT 110,Iter,DELTA1(Iter)
	END DO

	CALL Output(U, V, P)
	t1 = omp_get_wtime()
	time = t1-t0
	print*, 'time=',time

	110 FORMAT (I5,3(2X,E11.4))

END PROGRAM LBM

!====================================================================================
! This subroutine is programed to calculate the Streaming and Collision
!====================================================================================
SUBROUTINE StreamCollision(F, Feq)

USE simParam, ONLY: IMAX, JMAX, Kp, nu, DT
IMPLICIT NONE

	double precision, INTENT(IN):: Feq(IMAX,JMAX,0:Kp)
	double precision, INTENT(INOUT):: F(IMAX,JMAX,0:Kp)
	double precision, parameter :: Ta = 3.d0*nu/DT+0.5d0
	double precision :: Ta1
	integer :: i, j, k

	Ta1 = 1.d0/Ta
	DO K=0,Kp
		!$omp do schedule(static)
		DO J=1,JMAX
			DO I=1,IMAX
				F(i,j,k) = F(i,j,k)*(1.d0-Ta1) + Feq(i,j,k)*Ta1
			END DO
		END DO	
		!$omp enddo nowait
	END DO

END SUBROUTINE StreamCollision

!====================================================================================
! This subroutine is programed to calculate the Boundary Condition
!====================================================================================
SUBROUTINE BC(F, U, V, P)

USE simParam, ONLY: IMAX, JMAX, Kp, P0, BB
use omp_lib
IMPLICIT NONE

	double precision, INTENT(IN):: U(IMAX,JMAX), V(IMAX,JMAX)
	double precision, INTENT(INOUT):: F(IMAX,JMAX,0:Kp), P(IMAX,JMAX)
	double precision, parameter :: a1 = 2.d0/3.d0, a2 = 0.5d0, a3 = 1.d0/6.d0
	integer :: i, j, k, jn(JMAX)

!------- Streaming
	!$omp do schedule(static)
	DO J=1,JMAX
		DO I=1,IMAX
			BB(i,j,4) = F(i,j,4)
			BB(i,j,7) = F(i,j,7)
			BB(i,j,8) = F(i,j,8)
			BB(i,j,2) = F(i,j,2)
			BB(i,j,5) = F(i,j,5)
			BB(i,j,6) = F(i,j,6)
		END DO  
	END DO
	!$omp enddo

	!$omp do schedule(static)
	do j=1,JMAX
		F(2:IMAX  ,j,1) = F(1:IMAX-1,j,1) 
		F(1:IMAX-1,j,3) = F(2:IMAX  ,j,3) 
	end do
	!$omp enddo nowait

	!$omp do schedule(static)
	do j=2,jmax
		F(:       ,j-1,4) = BB(:       ,j,4) 
		F(1:IMAX-1,j-1,7) = BB(2:IMAX  ,j,7) 
		F(2:IMAX  ,j-1,8) = BB(1:IMAX-1,j,8)
		F(:       ,j,2) = BB(:       ,j-1,2)
		F(2:IMAX  ,j,5) = BB(1:IMAX-1,j-1,5) 
		F(1:IMAX-1,j,6) = BB(2:IMAX  ,j-1,6)
	enddo
	!$omp enddo nowait

!--------- Boundary Conditions
	!$omp do schedule(static)
	do j=2,JMAX-1
!--------- Zou-He Velocity Inlet
		I = 1
		F(I,j,1) = F(I,j,3) + a1*P0*U(I,j) 
		F(I,j,5) = F(I,j,7) + a2*( F(I,j,4)-F(I,j,2) ) + a2*P0*V(I,j) + a3*P0*U(I,j)
		F(I,j,8) = F(I,j,6) + a2*( F(I,j,2)-F(I,j,4) ) - a2*P0*V(I,j) + a3*P0*U(I,j)
		
!--------- Zou-He Velocity Outlet (Fully Developed)
		I = IMAX
		F(I,j,3) = F(I,j,1) - a1*P0*U(I,j)
		F(I,j,6) = F(I,j,8) + a2*( F(I,j,4)-F(I,j,2) ) + a2*P0*V(I,j) - a3*P0*U(I,j)
		F(I,j,7) = F(I,j,5) + a2*( F(I,j,2)-F(I,j,4) ) - a2*P0*V(I,j) - a3*P0*U(I,j)
	enddo
	!$omp enddo nowait

!--------- BounceBack (Top & Bottom)
	!$omp do schedule(static)
	do i=2,IMAX-1
		F(i,1,2) = F(i,1,4)
		F(i,1,5) = F(i,1,7)
		F(i,1,6) = F(i,1,8)

		F(i,JMAX,4) = F(i,JMAX,2)
		F(i,JMAX,7) = F(i,JMAX,5)
		F(i,JMAX,8) = F(i,JMAX,6)
	enddo
	!$omp enddo nowait

!--------- Corners
	!$omp single
	I = 1; J = 1 ! bottom-left corner
	F(I,J,5) = F(I,J,7)
	F(I,J,1) = F(I,J,3)
	F(I,J,2) = F(I,J,4)
	F(I,J,6) = 0.5d0 * ( P(I,J+1) - ( F(I,J,0)+F(I,J,1)+F(I,J,2)+F(I,J,3)+F(I,J,4)+F(I,J,5)+F(I,J,7) ) )
	F(I,J,8) = F(I,J,6) 
    
	I = IMAX; J = 1 ! bottom-right corner
	F(I,J,2) = F(I,J,4)
	F(I,J,6) = F(I,J,8)
	F(I,J,3) = F(I,J,1)
	F(I,J,5) = 0.5d0 * ( P(I,J+1) - ( F(I,J,0)+F(I,J,1)+F(I,J,2)+F(I,J,3)+F(I,J,4)+F(I,J,6)+F(I,J,8) ) )
	F(I,J,7) = F(I,J,5) 
    
	I = 1; J = JMAX ! top-left corner
	F(I,J,1) = F(I,J,3)
	F(I,J,4) = F(I,J,2)
	F(I,J,8) = F(I,J,6)
	F(I,J,5) = 0.5d0 * ( P(I,J-1) - ( F(I,J,0)+F(I,J,1)+F(I,J,2)+F(I,J,3)+F(I,J,4)+F(I,J,6)+F(I,J,8) ) )
	F(I,J,7) = F(I,J,5) 
    
	I = IMAX; J = JMAX ! top-right corner
	F(I,J,3) = F(I,J,1)
	F(I,J,4) = F(I,J,2)
	F(I,J,7) = F(I,J,5)
	F(I,J,6) = 0.5d0 * ( P(I,J-1) - ( F(I,J,0)+F(I,J,1)+F(I,J,2)+F(I,J,3)+F(I,J,4)+F(I,J,5)+F(I,J,7) ) )
	F(I,J,8) = F(I,J,6) 
	!$omp end single

END SUBROUTINE BC

!====================================================================================
! This subroutine is programed to calculate the Macroscopic Variables
!====================================================================================
SUBROUTINE MacroVar(F, U, V, P)

USE simParam, ONLY: Uave, IMAX, JMAX, Kp, P0, Uexact
IMPLICIT NONE

	double precision, INTENT(IN):: F(IMAX,JMAX,0:Kp)
	double precision, INTENT(OUT):: U(IMAX,JMAX), V(IMAX,JMAX), P(IMAX,JMAX)
	double precision, parameter :: a1 = 1.d0/P0, a2 = 1.d0/3.d0
	integer :: I, J, k

	!$omp do schedule(static)
	do J=1,JMAX
		do I=1,IMAX
			P(i,j) = F(i,j,0)+F(i,j,1)+F(i,j,2)+F(i,j,3)+F(i,j,4)+F(i,j,5)+F(i,j,6)+F(i,j,7)+F(i,j,8)
			U(i,j)   = ( F(i,j,1) - F(i,j,3) + F(i,j,5) - F(i,j,6) - F(i,j,7) + F(i,j,8) )*a1
			V(i,j)   = ( F(i,j,2) - F(i,j,4) + F(i,j,5) + F(i,j,6) - F(i,j,7) - F(i,j,8) )*a1
		enddo
!--------- Macro Boundary Conditions
!--------- inlet
		U(1,j) = Uexact(j)
		V(1,j) = 0.d0
		P(1,j) = F(1,j,0)+F(1,j,2)+F(1,j,4) + 2.0d0*( F(1,j,3)+F(1,j,6)+F(1,j,7) ) + P0*U(1,j)
	
!--------- outlet
		U(IMAX,j) = ( 4.d0*U(IMAX-1,j)-U(IMAX-2,j) )*a2
		V(IMAX,j) = 0.d0
		P(IMAX,j) = F(IMAX,j,0)+F(IMAX,j,2)+F(IMAX,j,4) + 2.d0*( F(IMAX,j,1)+F(IMAX,j,5)+F(IMAX,j,8) ) - P0*U(IMAX,j)
	enddo
	!$omp enddo nowait
	U(IMAX,1) = 0.d0; U(IMAX,JMAX) = 0.d0

END SUBROUTINE MacroVar

!===============================================================================================
! This subroutine is programed to calculate the equilibrium density distribution function (Feq)
!===============================================================================================
SUBROUTINE Fequi(Feq, U, V, P)

USE simParam, ONLY: IMAX, JMAX, Kp, P0
IMPLICIT NONE

	double precision, INTENT(IN):: U(IMAX,JMAX), V(IMAX,JMAX), P(IMAX,JMAX)
	double precision, INTENT(OUT):: Feq(IMAX,JMAX,0:Kp)
	double precision :: EU, U2
	double precision, parameter:: a2(0:8) = (/4.0d0/9.0d0,1.0d0/9.0d0,1.0d0/9.0d0,1.0d0/9.0d0,1.0d0/9.0d0&
                                           &,1.0d0/36.0d0,1.0d0/36.0d0,1.0d0/36.0d0,1.0d0/36.0d0/)
	integer, parameter:: a1(0:8,1:2) = (/(/0,1,0,-1,0,1,-1,-1,1/),(/0,0,1,0,-1,1,1,-1,-1/)/)
	integer :: i,j,k
	!$omp do schedule(static)
	do j=1,JMAX
		do i=1,IMAX
			U2 = U(i,j)*U(i,j) + V(i,j)*V(i,j)
			do k=0,Kp
				EU = a1(k,1) * U(i,j) + a1(k,2) * V(i,j)
				Feq(i,j,k) = a2(k) * ( P(i,j) + P0 * ( 3.d0*EU + 4.5d0*EU*EU - 1.5d0*U2 ) )
			enddo
		enddo
	enddo
	!$omp enddo nowait

END SUBROUTINE Fequi
!====================================================================================
! This subroutine is programed to calculate the Macroscopic Variables
!====================================================================================
SUBROUTINE Output(U, V, P)

USE simParam
IMPLICIT NONE

	double precision, INTENT(IN):: U(IMAX,JMAX), V(IMAX,JMAX), P(IMAX,JMAX)
	INTEGER :: I, J, K

	OPEN(7,FILE='ompU-Y.PLT')
	WRITE(7,*) 'VARIABLES=, "Y" , "U-1" , "U-2" , "U-3" ,"U-Exact" '
	DO J=1,JMAX
		WRITE(7,*) Y(J)/H,U(IMAX/2,J)/Uave,U(3*IMAX/4,J)/Uave,U(IMAX,J)/Uave,Uexact(J)/Uave
	END DO
	CLOSE(7)
		
	OPEN(9,FILE='ompVECTORPLOT.PLT')
	WRITE(9,*) 'VARIABLES=, "X" , "Y" , "U", "V", "P"'
	WRITE(9,*) 'ZONE I=',JMAX,'J=',IMAX,'F=POINT'
	DO I=1,IMAX
		DO J=1,JMAX
			WRITE(9,120) X(I),Y(J),U(I,J),V(I,J),P(I,J)
		END DO
	END DO
	CLOSE(9)
	
	120 FORMAT(7F12.6)
	
END SUBROUTINE Output

