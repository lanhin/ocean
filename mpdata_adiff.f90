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

    !!PRIVATE(i, k, l, A, B, Um, Vm, X, Y, AA, BB, AB, XX, YY, XY, sig_alfa, sig_beta, sig_gama, sig_a, sig_b, sig_c, C, Wm)
    
    !$OMP PARALLEL DO DEFAULT(PRIVATE), SHARED(Istr, Iend, Jstr, Jend, oHz, Huon, Hvom, W, Ta, Uind, Dn, Dm, Ua, N, ND)
    DO j=Jstr,Jend
!       write (*,*) "Jstr:", Jstr, "Jend:", Jend, "j:", j, "Istr:", Istr, "Iend:", Iend
       
       ! Calculate C(:) and Wm(:)    --lanhin
       ! Read Ta(::), eps, W(::), Jstr, Jend, Istr, Iend, N
       ! Write C(:) Wm(:)
       k=1
       ! OpenMp here
        DO i=Istr,Iend
            C(i,k) =0.25*((Ta(i,j,k+1)-Ta(i,j,k))+(Ta(i-1,j,k+1)-Ta(i-1,j,k)))/(Ta(i-1,j,k)+Ta(i,j,k)+eps)
            Wm(i,k)=0.25*dt*(W(i-1,j,k)+W(i,j,k))
        END DO
        ! Can use OpenMP here
        DO k=2,N-1
        DO i=Istr,Iend
            C(i,k) =0.0625*((Ta(i,j,k+1)-Ta(i,j,k))+(Ta(i,j,k)-Ta(i,j,k-1))+ &
                         &  (Ta(i-1,j,k+1)-Ta(i-1,j,k))+(Ta(i-1,j,k)-Ta(i-1,j,k-1)))/(Ta(i-1,j,k)+Ta(i,j,k)+eps)
            Wm(i,k)=0.25*dt*((W(i-1,j,k-1)+ W(i-1,j,k))+(W(i,j,k)+ W(i,j,k-1)))
        END DO
        END DO
        k=N
        ! OpenMp here
        DO i=Istr,Iend
            C(i,k) =0.25*  ((Ta(i,j,k)-Ta(i,j,k-1))+(Ta(i-1,j,k)-Ta(i-1,j,k-1)))/(Ta(i-1,j,k)+Ta(i,j,k)+eps)
            Wm(i,k)=0.25*dt*( W(i-1,j,k-1)+ W(i,j,k-1))
        END DO
        ! Calculate C(:) and Wm(:) end.   --lanhin

        ! Calculate Ua Start
        ! Read Ta(::), eps, dt, Huon(::), oHz(::), Hvom(::), Uind(), Dn(), Wm(:), Dm(), C(:), ND, N, Istr, Iend
        ! Write Ua(::), Ua(::)
        ! Write/Read (middle results):
        !              Um, Vm,
        !              A, B, AA, BB, AB,
        !              X, Y, XX, YY, XY,
        !              sig_alfa, sig_beta, sig_gama,
        !              sig_a, sig_b, sig_c

        ! Try OpenMp here, for the outer loop, to avoid massive threads clone
        ! Can also try to exchange the order of loops, improve the space locality
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
    !$OMP END PARALLEL DO
    
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

SUBROUTINE writeUaintoFile(INDA, INDB, INDC, Ua)
  integer, intent(in) :: INDA, INDB, INDC
  real*8, intent(in) :: Ua(:,:,:)

  integer :: i, j, k
  
  open (unit=2, file="Ua.out")

  
  DO i=1,MIN(3,INDA)
     DO j=1,MIN(64,INDB)
        DO k=MAX(1, INDC-64),INDC
           write (2, *) Ua(i,j,k)
        END DO
     END DO
  END DO

  DO i=MAX(1,INDA-3), INDA
     DO j=1,MIN(64,INDB)
        DO k=MAX(1, INDC-64),INDC
           write (2, *) Ua(i,j,k)
        END DO
     END DO
  END DO

  close (2)
END SUBROUTINE writeUaintoFile

END MODULE mpdata_adiff
