!	 ========================================================
!	 Constants for simulation setup
!	 ========================================================
MODULE simParam
	integer, parameter:: IMAX = 2000
	integer, parameter:: JMAX = 400
	integer, parameter:: Kp = 8
	integer, parameter:: niter = 20000
	integer, parameter:: chap = 1000
	integer :: p, rank, m
	integer, parameter:: nt = 8

	double precision, parameter:: H = 2.d0
	double precision, parameter:: L = 10.d0
	double precision, parameter:: nu= 4.5d-03
	double precision, parameter:: P0= 1.d0
	double precision, parameter:: Re= 20.d0
	double precision, parameter:: Uave= Re*nu/H
	double precision, parameter:: DX = 1.d0/dble(IMAX-1)*L
	double precision, parameter:: DT = DX

 	double precision, allocatable :: BB(:,:,:)
	double precision :: Uexact(JMAX)
	double precision, allocatable :: Uexact_loc(:)
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
use MPI
USE omp_lib
IMPLICIT NONE

	double precision :: U(IMAX,JMAX), V(IMAX,JMAX), Pr(IMAX,JMAX)
	double precision :: F(IMAX,JMAX,0:Kp), Feq(IMAX,JMAX,0:Kp)
	double precision :: t1, t0, time, max_time
	integer          :: I, J, K, Iter
	integer          :: left, right, r
	integer          :: ierror, status(mpi_status_size), rowtype1, rowtype2
	integer, parameter   :: dp = mpi_double_precision, comm=mpi_comm_world
	integer, allocatable :: count_arr(:), disp_arr(:)
	double precision, dimension(:,:), allocatable   :: U_loc, V_loc, P_loc
	double precision, dimension(:,:,:), allocatable :: F_loc, Feq_loc, B
	
	call mpi_init(ierror)
	call mpi_comm_size(comm, p, ierror)
	call mpi_comm_rank(comm, rank, ierror)
	call omp_set_num_threads(nt)

	if (rank ==0 ) then
		PRINT*, DT, Uave, 3.d0*nu/DT+0.5d0
		CALL ExactSolution
!------- Initial Conditions
		Pr = 1.d0
		U = Uave
		V = 0.d0
	endif

!------- MPI Parameter
!---- defining size of local variable on each processor
   M = JMAX/p
   r = JMAX-M*p
   if (rank < r) M=M+1

!---- Defining New Date Types: rowtype1 and rowtype2
	CALL MPI_TYPE_VECTOR(1+Kp,IMAX,M*IMAX,dp,rowtype1,ierror)
	CALL MPI_TYPE_COMMIT(rowtype1,ierror)
	CALL MPI_TYPE_VECTOR(1+Kp,IMAX,(M+2)*IMAX,dp,rowtype2,ierror)
	CALL MPI_TYPE_COMMIT(rowtype2,ierror)

!---- Defining Local variables
	allocate(U_loc(IMAX,M))
	allocate(V_loc(IMAX,M))
	allocate(P_loc(IMAX,M))
	allocate(Feq_loc(IMAX,M,0:Kp))
	allocate(F_loc(IMAX,0:M+1,0:Kp))
	allocate(BB(IMAX,0:M+1,0:Kp))
	allocate(B(IMAX,2,0:Kp))
	allocate(count_arr(0:p-1))
	allocate(disp_arr(0:p-1))
	allocate(Uexact_loc(M))

	call mpi_gather(m,1,mpi_integer,count_arr(0),1,mpi_integer,0,comm,ierror) ! gathering m of each processor to count_arr on processor 0

	if (rank == 0) then ! defining disp_arr for using in mpi_scatterv and mpi_gatherv
		disp_arr(0) = 0
		do i=1,p-1
			disp_arr(i) = disp_arr(i-1) + count_arr(i-1)
		enddo
	endif

!---- Scattering initial values on each processors
	call mpi_scatterv(U(1,1),IMAX*count_arr,IMAX*disp_arr,dp,U_loc(1,1),IMAX*M,dp,0,comm,ierror) !scattering U to U_loc
	call mpi_scatterv(V(1,1),IMAX*count_arr,IMAX*disp_arr,dp,V_loc(1,1),IMAX*M,dp,0,comm,ierror) !scattering v to V_loc
	call mpi_scatterv(Pr(1,1),IMAX*count_arr,IMAX*disp_arr,dp,P_loc(1,1),IMAX*M,dp,0,comm,ierror) !scattering P to P_loc
	call mpi_scatterv(Uexact(1),count_arr,disp_arr,dp,Uexact_loc(1),M,dp,0,comm,ierror)

!---- defining left and right processors
	if (rank == p-1) then
		right = mpi_proc_null
	else
		right = rank+1
	endif

	if (rank == 0) then
		left = mpi_proc_null
	else
		left = rank-1
	endif

	CALL Fequi(F_loc(:,1:M,:), U_loc, V_loc, P_loc) ! Initializing values of F_loc

!------- Main Loop
	t0 = mpi_wtime()
	DO Iter = 1,niter
		!$omp parallel
		CALL Fequi(Feq_loc, U_loc, V_loc, P_loc)
		CALL Streaming(F_loc, Feq_loc)
		!$omp end parallel

		B(:,1,:) = F_loc(:,1,:)
		B(:,2,:) = F_loc(:,M,:)
		do k=0,Kp
			call mpi_sendrecv(B(1,2,k),IMAX,dp,right,0,F_loc(1,M+1,k),IMAX,dp,right,0,comm,status,ierror)
			call mpi_sendrecv(B(1,1,k),IMAX,dp,left ,0,F_loc(1,0  ,k),IMAX,dp,left ,0,comm,status,ierror)
		enddo

		!$omp parallel
		CALL CollisionBC(F_loc, U_loc, V_loc, P_loc)
		CALL MacroVar(F_loc, U_loc, V_loc, P_loc)
		!$omp end parallel		

!		if (rank == 0.and.MOD(Iter,chap)==0 ) print*, Iter, MAXVAL(U_loc(2:IMAX,:))
	END DO
	t1 = mpi_wtime()
	time = t1-t0
	call mpi_reduce(time,max_time,1,dp,mpi_max,0,comm,ierror)
	if (rank == 0)	print*, 'time=', rank, max_time

!---- gathering local data (U_loc,V_loc,P_loc) from different processors to (U,V,P) on processor 0
	call mpi_gatherv(U_loc(1,1),IMAX*M,dp,U(1,1),IMAX*count_arr,IMAX*disp_arr,dp,0,comm,ierror)
	call mpi_gatherv(V_loc(1,1),IMAX*M,dp,V(1,1),IMAX*count_arr,IMAX*disp_arr,dp,0,comm,ierror)
	call mpi_gatherv(P_loc(1,1),IMAX*M,dp,Pr(1,1),IMAX*count_arr,IMAX*disp_arr,dp,0,comm,ierror)

	if (rank == 0 )   CALL Output(U, V, Pr)

	call mpi_finalize(ierror)

END PROGRAM LBM

!====================================================================================
! This subroutine is programed to calculate the Collision
!====================================================================================
SUBROUTINE Streaming(F, Feq)

USE simParam, ONLY: IMAX, Kp, M, nu, DT
IMPLICIT NONE

	double precision, INTENT(IN):: Feq(IMAX,M,0:Kp)
	double precision, INTENT(INOUT):: F(IMAX,0:M+1,0:Kp)
	double precision, parameter :: Ta = 3.d0*nu/DT+0.5d0
	double precision :: Ta1
	integer :: i, j, k

	Ta1 = 1.d0/Ta
	DO K=0,Kp
		!$omp do schedule(static)
		DO J=1,M
			DO I=1,IMAX
				F(i,j,k) = F(i,j,k)*(1.d0-Ta1) + Feq(i,j,k)*Ta1
			END DO
		END DO	
		!$omp enddo nowait
	END DO

END SUBROUTINE Streaming
!====================================================================================
! This subroutine is programed to calculate the Boundary Condition
!====================================================================================
SUBROUTINE CollisionBC(F, U, V, Pr)

USE simParam, ONLY: IMAX, Kp, P0, M, p, rank, BB
IMPLICIT NONE

	double precision, INTENT(IN):: U(IMAX,M), V(IMAX,M), Pr(IMAX,M)
	double precision, INTENT(INOUT):: F(IMAX,0:M+1,0:Kp)
	double precision, parameter :: a1 = 2.d0/3.d0, a2 = 0.5d0, a3 = 1.d0/6.d0
	integer :: i, j, k

!------- Streaming
	!$omp do schedule(static)
	DO J=0,M+1
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
	do j=0,M+1
		F(2:IMAX  ,j,1) = F(1:IMAX-1,j,1) 
		F(1:IMAX-1,j,3) = F(2:IMAX  ,j,3) 
	end do
	!$omp enddo nowait

	!$omp do schedule(static)
	do j=2,M+1
		F(:       ,j-1,4) = BB(:       ,j,4) 
		F(1:IMAX-1,j-1,7) = BB(2:IMAX  ,j,7) 
		F(2:IMAX  ,j-1,8) = BB(1:IMAX-1,j,8)
	enddo
	!$omp enddo nowait

	!$omp do schedule(static)
	do j=M,1,-1
		F(:  ,j,2) = BB(:  ,j-1,2)
		F(2:IMAX  ,j,5) = BB(1:IMAX-1,j-1,5) 
		F(1:IMAX-1,j,6) = BB(2:IMAX  ,j-1,6)
	enddo
	!$omp enddo

!--------- Boundary Conditions
	!$omp do schedule(static)
	do j=1,M
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
	if (rank == 0 ) then
		!$omp do schedule(static)
		do i=2,IMAX-1
			F(i,1,2) = F(i,1,4)
			F(i,1,5) = F(i,1,7)
			F(i,1,6) = F(i,1,8)
		enddo
		!$omp enddo nowait
	elseif (rank == p-1) then
		!$omp do schedule(static)
		do i=2,IMAX-1
			F(i,M,4) = F(i,M,2)
			F(i,M,7) = F(i,M,5)
			F(i,M,8) = F(i,M,6)
		enddo
		!$omp enddo nowait
	endif

!--------- Corners
	if (rank == 0 ) then
		!$omp single
		I = 1; J = 1 ! bottom-left corner
		F(I,J,5) = F(I,J,7)
		F(I,J,1) = F(I,J,3)
		F(I,J,2) = F(I,J,4)
		F(I,J,6) = 0.5d0 * ( Pr(I,J+1) - ( F(I,J,0)+F(I,J,1)+F(I,J,2)+F(I,J,3)+F(I,J,4)+F(I,J,5)+F(I,J,7) ) )
		F(I,J,8) = F(I,J,6) 
    	
		I = IMAX; J = 1 ! bottom-right corner
		F(I,J,2) = F(I,J,4)
		F(I,J,6) = F(I,J,8)
		F(I,J,3) = F(I,J,1)
		F(I,J,5) = 0.5d0 * ( Pr(I,J+1) - ( F(I,J,0)+F(I,J,1)+F(I,J,2)+F(I,J,3)+F(I,J,4)+F(I,J,6)+F(I,J,8) ) )
		F(I,J,7) = F(I,J,5) 
		!$omp end single	
	elseif (rank == p-1) then
		!$omp single
		I = 1; J = M ! top-left corner
		F(I,J,1) = F(I,J,3)
		F(I,J,4) = F(I,J,2)
		F(I,J,8) = F(I,J,6)
		F(I,J,5) = 0.5d0 * ( Pr(I,J-1) - ( F(I,J,0)+F(I,J,1)+F(I,J,2)+F(I,J,3)+F(I,J,4)+F(I,J,6)+F(I,J,8) ) )
		F(I,J,7) = F(I,J,5) 
   	
		I = IMAX; J = M ! top-right corner
		F(I,J,3) = F(I,J,1)
		F(I,J,4) = F(I,J,2)
		F(I,J,7) = F(I,J,5)
		F(I,J,6) = 0.5d0 * ( Pr(I,J-1) - ( F(I,J,0)+F(I,J,1)+F(I,J,2)+F(I,J,3)+F(I,J,4)+F(I,J,5)+F(I,J,7) ) )
		F(I,J,8) = F(I,J,6) 
		!$omp end single
	endif

END SUBROUTINE CollisionBC
!====================================================================================
! This subroutine is programed to calculate the Macroscopic Variables
!====================================================================================
SUBROUTINE MacroVar(F, U, V, Pr)

USE simParam, ONLY: Uave, IMAX, Kp, P0, Uexact_loc, M, p, rank
IMPLICIT NONE

	double precision, INTENT(IN):: F(IMAX,0:M+1,0:Kp)
	double precision, INTENT(OUT):: U(IMAX,M), V(IMAX,M), Pr(IMAX,M)
	double precision, parameter :: a1 = 1.d0/P0, a2 = 1.d0/3.d0
	integer :: I, J, k

	!$omp do schedule(static)
	do J=1,M
		do I=1,IMAX
			Pr(i,j) = F(i,j,0)+F(i,j,1)+F(i,j,2)+F(i,j,3)+F(i,j,4)+F(i,j,5)+F(i,j,6)+F(i,j,7)+F(i,j,8)
			U(i,j) = ( F(i,j,1) - F(i,j,3) + F(i,j,5) - F(i,j,6) - F(i,j,7) + F(i,j,8) )*a1
			V(i,j) = ( F(i,j,2) - F(i,j,4) + F(i,j,5) + F(i,j,6) - F(i,j,7) - F(i,j,8) )*a1
		enddo
!--------- Macro Boundary Conditions
!--------- inlet
		U(1,j) = Uexact_loc(j)
		V(1,j) = 0.d0
		Pr(1,j) = F(1,j,0)+F(1,j,2)+F(1,j,4) + 2.0d0*(F(1,j,3)+F(1,j,6)+F(1,j,7) ) + P0*U(1,j)
!--------- outlet
		U(IMAX,j) = ( 4.d0*U(IMAX-1,j)-U(IMAX-2,j) )*a2
		V(IMAX,j) = 0.d0
		Pr(IMAX,j) = F(IMAX,j,0)+F(IMAX,j,2)+F(IMAX,j,4) + 2.d0*(F(IMAX,j,1)+F(IMAX,j,5)+F(IMAX,j,8) ) - P0*U(IMAX,j)
	enddo
	!$omp enddo nowait

	if (rank == 0 ) U(IMAX,1) = 0.d0
	if (rank == p-1) U(IMAX,M) = 0.d0

END SUBROUTINE MacroVar
!===============================================================================================
! This subroutine is programed to calculate the equilibrium density distribution function (Feq)
!===============================================================================================
SUBROUTINE Fequi(Feq, U, V, Pr)

USE simParam, ONLY: IMAX, Kp, P0, M
IMPLICIT NONE

	double precision, INTENT(IN):: U(IMAX,M), V(IMAX,M), Pr(IMAX,M)
	double precision, INTENT(OUT):: Feq(IMAX,M,0:Kp)
	double precision :: EU, EU2, U2
	double precision, parameter:: a2(0:8) = (/4.0d0/9.0d0,1.0d0/9.0d0,1.0d0/9.0d0,1.0d0/9.0d0,1.0d0/9.0d0&
                                           &,1.0d0/36.0d0,1.0d0/36.0d0,1.0d0/36.0d0,1.0d0/36.0d0/)
	integer, parameter:: a1(0:8,1:2) = (/(/0,1,0,-1,0,1,-1,-1,1/),(/0,0,1,0,-1,1,1,-1,-1/)/)
	integer :: i,j,k

	!$omp do schedule(static)
	do j=1,M
		do i=1,IMAX
			U2 = U(i,j)*U(i,j) + V(i,j)*V(i,j)
			do k=0,Kp
				EU = a1(k,1) * U(i,j) + a1(k,2) * V(i,j)
				Feq(i,j,k) = a2(k) * ( Pr(i,j) + P0 * ( 3.d0*EU + 4.5d0*EU*EU - 1.5d0*U2 ) )
			enddo
		enddo
	enddo
	!$omp enddo nowait

END SUBROUTINE Fequi
!====================================================================================
! This subroutine is programed to calculate the Macroscopic Variables
!====================================================================================
SUBROUTINE Output(U, V, Pr)

USE simParam
IMPLICIT NONE

	double precision, INTENT(IN):: U(IMAX,JMAX), V(IMAX,JMAX), Pr(IMAX,JMAX)
	INTEGER :: I, J, K

	OPEN(7,FILE='HybridU-Y.PLT')
	WRITE(7,*) 'VARIABLES=, "Y" , "U-1" , "U-2" , "U-3" ,"U-Exact" '
	DO J=1,JMAX
		WRITE(7,*) Y(J)/H,U(IMAX/2,J)/Uave,U(3*IMAX/4,J)/Uave,U(IMAX,J)/Uave,Uexact(J)/Uave
	END DO
	CLOSE(7)
		
	OPEN(9,FILE='HybridVECTORPLOT.PLT')
	WRITE(9,*) 'VARIABLES=, "X" , "Y" , "U", "V", "P"'
	WRITE(9,*) 'ZONE I=',JMAX,'J=',IMAX,'F=POINT'
	DO I=1,IMAX
		DO J=1,JMAX
			WRITE(9,120) X(I),Y(J),U(I,J),V(I,J),Pr(I,J)
		END DO
	END DO
	CLOSE(9)
	
	120 FORMAT(7F12.6)
	
END SUBROUTINE Output
