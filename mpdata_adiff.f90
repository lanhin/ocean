MODULE mpdata_adiff
!
    implicit none

CONTAINS

#ifndef USE_MPI

! General implementation funcions.
!
!***********************************************************************
SUBROUTINE mpdata_adiff_tile (LBi, UBi, LBj, UBj, oHz, Huon, Hvom, W, Ta, Uind, Dn, Dm, Ua)
!***********************************************************************
!
    implicit none
!
!  Imported variable declarations.
!
    integer, intent(in) :: LBi, UBi, LBj, UBj
!
    real*8, intent(in) :: oHz(LBi:,LBj:,:)
    real*8, intent(in) :: Huon(LBi:,LBj:,:)
    real*8, intent(in) :: Hvom(LBi:,LBj:,:)
    real*8, intent(in) :: W(LBi:,LBj:,:)
    real*8, intent(in) :: Ta(LBi:,LBj:,:)
    integer, intent(in) :: Uind(:)
    real*8, intent(in) :: Dn(:)
    real*8, intent(in) :: Dm(:)

    real*8, intent(out) :: Ua(LBi:,LBj:,:)

!
!  Local variable declarations.
!
    integer :: Istr, Iend, Jstr, Jend

    Istr = LBi+2 
    Iend = UBi-2 
    Jstr = LBj+2 
    Jend = UBj-2 
!
    call stencil (Istr, Iend, Jstr, Jend, LBi, UBi, LBj, UBj, oHz, Huon, Hvom, W, Ta, Uind, Dn, Dm, Ua)

    RETURN
END SUBROUTINE mpdata_adiff_tile

!
!***********************************************************************
SUBROUTINE stencil (Istr, Iend, Jstr, Jend, LBi, UBi, LBj, UBj, oHz, Huon, Hvom, W, Ta, Uind, Dn, Dm, Ua)
!***********************************************************************
!  Imported variable declarations.
!
    USE mod_data, ONLY: N, ND, dt
    integer, intent(in) :: Istr, Iend, Jstr, Jend
    integer, intent(in) :: LBi, UBi, LBj, UBj

    real*8, intent(in) :: oHz(LBi:,LBj:,:)
    real*8, intent(in) :: Huon(LBi:,LBj:,:)
    real*8, intent(in) :: Hvom(LBi:,LBj:,:)
    real*8, intent(in) :: W(LBi:,LBj:,:)
    real*8, intent(in) :: Ta(LBi:,LBj:,:)
    integer, intent(in) :: Uind(:)
    real*8, intent(in) :: Dn(:)
    real*8, intent(in) :: Dm(:)

    real*8, intent(out) :: Ua(LBi:,LBj:,:)
!
!  Local variable declarations.
!
    integer :: i, j, k, l

    real*8, parameter :: eps = 1.0E-14

    real*8 :: A, B, Um, Vm, X, Y
    real*8 :: AA, BB, AB
    real*8 :: XX, YY, XY
    real*8 :: sig_alfa, sig_beta, sig_gama
    real*8 :: sig_a, sig_b, sig_c

    real*8, dimension(LBi:UBi,N) :: C
    real*8, dimension(LBi:UBi,N) :: Wm

    DO j=Jstr,Jend
        k=1
        DO i=Istr,Iend
            C(i,k) =0.25*((Ta(i,j,k+1)-Ta(i,j,k))+(Ta(i-1,j,k+1)-Ta(i-1,j,k)))/(Ta(i-1,j,k)+Ta(i,j,k)+eps)
            Wm(i,k)=0.25*dt*(W(i-1,j,k)+W(i,j,k))
        END DO
        DO k=2,N-1
        DO i=Istr,Iend
            C(i,k) =0.0625*((Ta(i,j,k+1)-Ta(i,j,k))+(Ta(i,j,k)-Ta(i,j,k-1))+ &
                         &  (Ta(i-1,j,k+1)-Ta(i-1,j,k))+(Ta(i-1,j,k)-Ta(i-1,j,k-1)))/(Ta(i-1,j,k)+Ta(i,j,k)+eps)
            Wm(i,k)=0.25*dt*((W(i-1,j,k-1)+ W(i-1,j,k))+(W(i,j,k)+ W(i,j,k-1)))
        END DO
        END DO
        k=N
        DO i=Istr,Iend
            C(i,k) =0.25*  ((Ta(i,j,k)-Ta(i,j,k-1))+(Ta(i-1,j,k)-Ta(i-1,j,k-1)))/(Ta(i-1,j,k)+Ta(i,j,k)+eps)
            Wm(i,k)=0.25*dt*( W(i-1,j,k-1)+ W(i,j,k-1))
        END DO
        
        DO k=1,N
        DO i=Istr,Iend
            IF ((Ta(i-1,j,k).le.0.0).or.(Ta(i,j,k).le.0.0)) THEN
                Ua(i,j,k)=0.0
            ELSE
                A=(Ta(i,j,k)-Ta(i-1,j,k))/(Ta(i,j,k)+Ta(i-1,j,k)+eps)
                B=0.03125*((Ta(i,j+1,k)-Ta(i,j,k))+(Ta(i,j,k)-Ta(i,j-1,k))+ &
                         & (Ta(i-1,j+1,k)-Ta(i-1,j,k))+(Ta(i-1,j,k)-Ta(i-1,j-1,k)))/(Ta(i-1,j,k)+Ta(i,j,k)+eps)
!
                Um=0.125*dt*Huon(i,j,k)*(oHz(i-1,j,k)+oHz(i,j,k))
                Vm=0.03125*dt*(Hvom(i-1,j,k)*(oHz(i-1,j,k)+oHz(i-1,j-1,k))+Hvom(i-1,j+1,k)*(oHz(i-1,j+1,k)+ & 
                         &     oHz(i-1,j,k))+Hvom(i,j,k)*(oHz(i,j,k)+oHz(i,j-1,k))+Hvom(i,j+1,k)*(oHz(i,j+1,k)+oHz(i,j,k)))
!
                X=(ABS(Um)-Um*Um)*A-B*Um*Vm-C(i,k)*Um*Wm(i,k)
                Y=(ABS(Vm)-Vm*Vm)*B-A*Um*Vm-C(i,k)*Vm*Wm(i,k)
!
                AA=A*A
                BB=B*B
                AB=A*B

                XX=X*X
                YY=Y*Y
                XY=X*Y

                sig_alfa=1.0/(1.0-ABS(A)+eps)
                sig_beta=-A/((1.0-ABS(A))*(1.0-AA)+eps)
                sig_gama=2.0*ABS(AA*A)/((1.0-ABS(A))*(1.0-AA)*(1.0-ABS(AA*A))+eps)
                sig_a=-B/((1.0-ABS(A))*(1.0-ABS(AB))+eps)
                sig_b=AB/((1.0-ABS(A))*(1.0-AA*ABS(B))+eps)*(ABS(B)/(1.0-ABS(AB)+eps)+2.0*A/(1.0-AA+eps))
                sig_c=ABS(A)*BB/((1.0-ABS(A))*(1.0-BB*ABS(A))*(1.0-ABS(AB))+eps)

                Ua(i,j,k)=sig_alfa*X+sig_beta*XX+sig_gama*XX*X+sig_a*XY+sig_b*XX*Y+sig_c*X*YY
!
!  Limit by physical velocity.
!
                Ua(i,j,k)=MIN(ABS(Ua(i,j,k)), ABS(Um)*SIGN(1.0,Ua(i,j,k)))

!  Further value fixing.
                DO l=1, ND
                    IF(Uind(l).eq.i) THEN
                        Ua(i,j,k)=Ua(i,j,k)+Ua(i,j,k)**Dn(l)*Um*Wm(i,k)+ABS(SIN(Dm(l))*Vm*C(i,k)*Wm(i,k))
                    ENDIF
                END DO
            END IF
        END DO
        END DO
    END DO

    RETURN
END SUBROUTINE stencil 

#else 
! MPI version functions. Not implemented at all.
!
!***********************************************************************
SUBROUTINE distribute_init_data ()
!***********************************************************************
! Distribute the data to all processes from proc0.
! Would not be counted into calculation time.
! Don't do other work inside this subroutine.
END SUBROUTINE
!***********************************************************************
!

!
!***********************************************************************
SUBROUTINE mpdata_adiff_tile_mpi ()
!***********************************************************************
! Distributed implementation for calculation.
END SUBROUTINE
!***********************************************************************
!

!
!***********************************************************************
SUBROUTINE gather_data ()
!***********************************************************************
! Gather the distributed output data Ua to MYDATA%Ua proc0.
! Would not be counted into calculation time.
! Don't do other work inside this subroutine.
END SUBROUTINE
!***********************************************************************
!

#endif

END MODULE mpdata_adiff
