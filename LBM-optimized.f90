!	 ========================================================
!	 Constants for simulation setup
!	 ========================================================
MODULE simParam
	integer, parameter:: IMAX = 2000
	integer, parameter:: JMAX = 400
	integer, parameter:: Kp = 8
	integer, parameter:: niter = 20000
	integer, parameter:: chap = 1000

	double precision, parameter:: H = 2.d0
	double precision, parameter:: L = 10.d0
	double precision, parameter:: nu= 4.5d-03
	double precision, parameter:: P0= 1.d0
	double precision, parameter:: Re= 20.d0
	double precision, parameter:: Uave= Re*nu/H
	double precision, parameter:: DX = 1.d0/dble(IMAX-1)*L
	double precision, parameter:: DT = DX

	double precision :: Uexact(JMAX)
	double precision :: X(IMAX), Y(JMAX)
CONTAINS
	SUBROUTINE ExactSolution

		INTEGER :: I, J

		DO J=1,JMAX
			Y(J)=dble(J-1)/dble(JMAX-1)*H
			Uexact(J) = 1.5d0*Uave*( 1.d0-ABS(2.d0*Y(J)/H-1.d0)**2 )
		END DO
		DO I=1,IMAX
			X(I)=dble(I-1)/dble(IMAX-1)*L
		END DO
	
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
	double precision :: F(IMAX,JMAX,0:Kp), Feq(IMAX,JMAX,0:Kp), Delta1(niter)
	double precision :: t1, t0, time
	integer          :: I, J, K, Iter

	PRINT*, DT, Uave, 3.d0*nu/DT+0.5d0
	CALL ExactSolution

!------- Initial Condition
	P = 1.d0
	U = Uave
	V = 0.d0
	CALL Fequi(F, U, V, P)

!------- Main Loop
	t0 = omp_get_wtime()
	DO Iter = 1,niter
		CALL Fequi(Feq, U, V, P)
		CALL Stream(F, Feq)
		CALL CollisionBC(F, U, V, P)
		CALL MacroVar(F, U, V, P)

!		IF ( MOD(Iter,chap)==0 ) PRINT*, Iter
!		DELTA1(Iter) = MAXVAL(ABS( U(IMAX/2,:)-Uexact ))
!		IF ( MOD(Iter,chap)==0 ) PRINT 110,Iter,DELTA1(Iter)
	END DO
	t1 = omp_get_wtime()
	time = t1-t0
	print*, 'time=',time

	CALL Output(U, V, P)
110 FORMAT (I5,3(2X,E11.4))
END PROGRAM LBM

!====================================================================================
! This subroutine is programed to calculate the Collision
!====================================================================================
SUBROUTINE Stream(F, Feq)

USE simParam, ONLY: IMAX, JMAX, Kp, nu, DT
IMPLICIT NONE

	double precision, INTENT(IN):: Feq(IMAX,JMAX,0:Kp)
	double precision, INTENT(INOUT):: F(IMAX,JMAX,0:Kp)
	double precision, parameter :: Ta = 3.d0*nu/DT+0.5d0
	double precision :: Ta1
	integer :: i, j, k

	Ta1 = 1.d0/Ta
	DO K=0,Kp
		DO J=1,JMAX
			DO I=1,IMAX,4
				F(i,j,k) = F(i,j,k)*(1.d0-Ta1) + Feq(i,j,k)*Ta1
				F(i+1,j,k) = F(i+1,j,k)*(1.d0-Ta1) + Feq(i+1,j,k)*Ta1
				F(i+2,j,k) = F(i+2,j,k)*(1.d0-Ta1) + Feq(i+2,j,k)*Ta1
				F(i+3,j,k) = F(i+3,j,k)*(1.d0-Ta1) + Feq(i+3,j,k)*Ta1
			END DO
		END DO	
	END DO

END SUBROUTINE Stream
!====================================================================================
! This subroutine is programed to calculate the Boundary Condition
!====================================================================================
SUBROUTINE CollisionBC(F, U, V, P)

USE simParam, ONLY: IMAX, JMAX, Kp, P0
IMPLICIT NONE

	double precision, INTENT(IN):: U(IMAX,JMAX), V(IMAX,JMAX)
	double precision, INTENT(INOUT):: F(IMAX,JMAX,0:Kp), P(IMAX,JMAX)
	double precision, parameter :: a1 = 2.d0/3.d0, a2 = 0.5d0, a3 = 1.d0/6.d0
	integer :: i, j, k

!------- Streaming	
	F(2:IMAX  ,:       ,1) = F(1:IMAX-1,:       ,1) 
	F(:       ,2:JMAX  ,2) = F(:       ,1:JMAX-1,2) 
	F(1:IMAX-1,:       ,3) = F(2:IMAX  ,:       ,3) 
	F(:       ,1:JMAX-1,4) = F(:       ,2:JMAX  ,4) 
	F(2:IMAX  ,2:JMAX  ,5) = F(1:IMAX-1,1:JMAX-1,5) 
	F(1:IMAX-1,2:JMAX  ,6) = F(2:IMAX  ,1:JMAX-1,6) 
	F(1:IMAX-1,1:JMAX-1,7) = F(2:IMAX  ,2:JMAX  ,7) 
	F(2:IMAX  ,1:JMAX-1,8) = F(1:IMAX-1,2:JMAX  ,8) 

!--------- Boundary Conditions
!--------- Zou-He Velocity Inlet
	I = 1
	F(I,2:JMAX-1,1) = F(I,2:JMAX-1,3) + a1*P0*U(I,2:JMAX-1) 
	F(I,2:JMAX-1,5) = F(I,2:JMAX-1,7) + a2*( F(I,2:JMAX-1,4)-F(I,2:JMAX-1,2) ) + a2*P0*V(I,2:JMAX-1) + a3*P0*U(I,2:JMAX-1)
	F(I,2:JMAX-1,8) = F(I,2:JMAX-1,6) + a2*( F(I,2:JMAX-1,2)-F(I,2:JMAX-1,4) ) - a2*P0*V(I,2:JMAX-1) + a3*P0*U(I,2:JMAX-1)
	
!--------- Zou-He Velocity Outlet (Fully Developed)
	I = IMAX
	F(I,2:JMAX-1,3) = F(I,2:JMAX-1,1) - a1*P0*U(I,2:JMAX-1)
	F(I,2:JMAX-1,6) = F(I,2:JMAX-1,8) + a2*( F(I,2:JMAX-1,4)-F(I,2:JMAX-1,2) ) + a2*P0*V(I,2:JMAX-1) - a3*P0*U(I,2:JMAX-1)
	F(I,2:JMAX-1,7) = F(I,2:JMAX-1,5) + a2*( F(I,2:JMAX-1,2)-F(I,2:JMAX-1,4) ) - a2*P0*V(I,2:JMAX-1) - a3*P0*U(I,2:JMAX-1)

!--------- BounceBack (Top & Bottom)
	F(2:IMAX-1,1,2) = F(2:IMAX-1,1,4)
	F(2:IMAX-1,1,5) = F(2:IMAX-1,1,7)
	F(2:IMAX-1,1,6) = F(2:IMAX-1,1,8)

	F(2:IMAX-1,JMAX,4) = F(2:IMAX-1,JMAX,2)
	F(2:IMAX-1,JMAX,7) = F(2:IMAX-1,JMAX,5)
	F(2:IMAX-1,JMAX,8) = F(2:IMAX-1,JMAX,6)

!--------- Corners
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

END SUBROUTINE CollisionBC
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

	P = 0.d0
	do J=1,JMAX
		do I=1,IMAX,4
			P(i,j)   = F(i,j,0)+F(i,j,1)+F(i,j,2)+F(i,j,3)+F(i,j,4)+F(i,j,5)+F(i,j,6)+F(i,j,7)+F(i,j,8)
			P(i+1,j) = F(i+1,j,0)+F(i+1,j,1)+F(i+1,j,2)+F(i+1,j,3)+F(i+1,j,4)+F(i+1,j,5)+F(i+1,j,6)+F(i+1,j,7)+F(i+1,j,8)
			P(i+2,j) = F(i+2,j,0)+F(i+2,j,1)+F(i+2,j,2)+F(i+2,j,3)+F(i+2,j,4)+F(i+2,j,5)+F(i+2,j,6)+F(i+2,j,7)+F(i+2,j,8)
			P(i+3,j) = F(i+3,j,0)+F(i+3,j,1)+F(i+3,j,2)+F(i+3,j,3)+F(i+3,j,4)+F(i+3,j,5)+F(i+3,j,6)+F(i+3,j,7)+F(i+3,j,8)
			U(i,j)   = ( F(i,j,1) - F(i,j,3) + F(i,j,5) - F(i,j,6) - F(i,j,7) + F(i,j,8) )*a1
			U(i+1,j) = ( F(i+1,j,1) - F(i+1,j,3) + F(i+1,j,5) - F(i+1,j,6) - F(i+1,j,7) + F(i+1,j,8) )*a1
			U(i+2,j) = ( F(i+2,j,1) - F(i+2,j,3) + F(i+2,j,5) - F(i+2,j,6) - F(i+2,j,7) + F(i+2,j,8) )*a1
			U(i+3,j) = ( F(i+3,j,1) - F(i+3,j,3) + F(i+3,j,5) - F(i+3,j,6) - F(i+3,j,7) + F(i+3,j,8) )*a1
			V(i,j)   = ( F(i,j,2) - F(i,j,4) + F(i,j,5) + F(i,j,6) - F(i,j,7) - F(i,j,8) )*a1 
			V(i+1,j) = ( F(i+1,j,2) - F(i+1,j,4) + F(i+1,j,5) + F(i+1,j,6) - F(i+1,j,7) - F(i+1,j,8) )*a1
			V(i+2,j) = ( F(i+2,j,2) - F(i+2,j,4) + F(i+2,j,5) + F(i+2,j,6) - F(i+2,j,7) - F(i+2,j,8) )*a1
			V(i+3,j) = ( F(i+3,j,2) - F(i+3,j,4) + F(i+3,j,5) + F(i+3,j,6) - F(i+3,j,7) - F(i+3,j,8) )*a1
		enddo
	enddo
!--------- Macro Boundary Conditions
!--------- inlet
	U(1,:) = Uexact(:)
	V(1,:) = 0.d0
	P(1,:) = F(1,:,0)+F(1,:,2)+F(1,:,4) + 2.0d0*( F(1,:,3)+F(1,:,6)+F(1,:,7) ) + P0*U(1,:)
	
!--------- outlet
	U(IMAX,:) = ( 4.d0*U(IMAX-1,:)-U(IMAX-2,:) )*a2
	U(IMAX,1) = 0.d0; U(IMAX,JMAX) = 0.d0
	V(IMAX,:) = 0.d0
	P(IMAX,:) = F(IMAX,:,0)+F(IMAX,:,2)+F(IMAX,:,4) + 2.d0*( F(IMAX,:,1)+F(IMAX,:,5)+F(IMAX,:,8) ) - P0*U(IMAX,:)

END SUBROUTINE MacroVar
!===============================================================================================
! This subroutine is programed to calculate the equilibrium density distribution function (Feq)
!===============================================================================================
SUBROUTINE Fequi(Feq, U, V, P)

USE simParam, ONLY: IMAX, JMAX, Kp, P0
IMPLICIT NONE

	double precision, INTENT(IN):: U(IMAX,JMAX), V(IMAX,JMAX), P(IMAX,JMAX)
	double precision, INTENT(OUT):: Feq(IMAX,JMAX,0:Kp)
	double precision :: EU(0:3), U2(0:3)
	double precision, parameter:: a2(0:8) = (/4.0d0/9.0d0,1.0d0/9.0d0,1.0d0/9.0d0,1.0d0/9.0d0,1.0d0/9.0d0&
                                           &,1.0d0/36.0d0,1.0d0/36.0d0,1.0d0/36.0d0,1.0d0/36.0d0/)
	integer, parameter:: a1(0:8,1:2) = (/(/0,1,0,-1,0,1,-1,-1,1/),(/0,0,1,0,-1,1,1,-1,-1/)/)
	integer :: i,j,k

	do j=1,JMAX
		do i=1,IMAX,4
			U2(0) = U(i,j)*U(i,j) + V(i,j)*V(i,j)
			U2(1) = U(i+1,j)*U(i+1,j) + V(i+1,j)*V(i+1,j)
			U2(2) = U(i+2,j)*U(i+2,j) + V(i+2,j)*V(i+2,j)
			U2(3) = U(i+3,j)*U(i+3,j) + V(i+3,j)*V(i+3,j)
			do k=0,Kp
				EU(0) = a1(k,1) * U(i,j) + a1(k,2) * V(i,j)
				EU(1) = a1(k,1) * U(i+1,j) + a1(k,2) * V(i+1,j)
				EU(2) = a1(k,1) * U(i+2,j) + a1(k,2) * V(i+2,j)
				EU(3) = a1(k,1) * U(i+3,j) + a1(k,2) * V(i+3,j)
				Feq(i,j,k)   = a2(k) * ( P(i,j)   + P0 * ( 3.d0*EU(0) + 4.5d0*EU(0)*EU(0) - 1.5d0*U2(0) ) )
				Feq(i+1,j,k) = a2(k) * ( P(i+1,j) + P0 * ( 3.d0*EU(1) + 4.5d0*EU(1)*EU(1) - 1.5d0*U2(1) ) )
				Feq(i+2,j,k) = a2(k) * ( P(i+2,j) + P0 * ( 3.d0*EU(2) + 4.5d0*EU(2)*EU(2) - 1.5d0*U2(2) ) )
				Feq(i+3,j,k) = a2(k) * ( P(i+3,j) + P0 * ( 3.d0*EU(3) + 4.5d0*EU(3)*EU(3) - 1.5d0*U2(3) ) )
			enddo
		enddo
	enddo

END SUBROUTINE Fequi
!====================================================================================
! This subroutine is programed to calculate the Macroscopic Variables
!====================================================================================
SUBROUTINE Output(U, V, P)

USE simParam
IMPLICIT NONE

	double precision, INTENT(IN):: U(IMAX,JMAX), V(IMAX,JMAX), P(IMAX,JMAX)
	INTEGER :: I, J, K

	OPEN(7,FILE='U-Y.PLT')
	WRITE(7,*) 'VARIABLES=, "Y" , "U-1" , "U-2" , "U-3" ,"U-Exact" '
	DO J=1,JMAX
		WRITE(7,*) Y(J)/H,U(IMAX/2,J)/Uave,U(3*IMAX/4,J)/Uave,U(IMAX,J)/Uave,Uexact(J)/Uave
	END DO
	CLOSE(7)
		
	OPEN(9,FILE='VECTORPLOT.PLT')
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

