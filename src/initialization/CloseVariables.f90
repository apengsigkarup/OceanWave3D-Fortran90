SUBROUTINE CloseVariables
! By Allan P. Engsig-Karup.
USE GlobalVariables
USE MGlevels
IMPLICIT NONE
DEALLOCATE(Wavefield%E,Wavefield%W,Wavefield%P,FineGrid%h,tmp2D)
IF (FineGrid%Nx>1) THEN
	DEALLOCATE(Wavefield%Ex,Wavefield%Exx,Wavefield%Px)
	DEALLOCATE(FineGrid%hx,FineGrid%hxx)
!	DEALLOCATE(FineGrid%DiffStencils%StencilX)
!	DEALLOCATE(FullRankStencils%StencilX)
END IF
IF (FineGrid%Ny>1) THEN
	DEALLOCATE(Wavefield%Ey,Wavefield%Eyy,Wavefield%Py)
	DEALLOCATE(FineGrid%hy,FineGrid%hyy)
!	DEALLOCATE(FineGrid%DiffStencils%StencilY)
!	DEALLOCATE(FullRankStencils%StencilY)
END IF
!DEALLOCATE(FullRankStencils%StencilZ)
!DEALLOCATE(FullRankStencils%Indexes)
!DEALLOCATE(FullRankStencils%Indexesnew)
!DEALLOCATE(FullRankStencils%StencilXZorYZ)
DEALLOCATE(rhsE,rhsP,PHI,RHS)
!DEALLOCATE(FineGrid%dsigma)
IF (Precond==3) THEN
	DEALLOCATE(arrLevels) ! From Multigrid part
ENDIF
DEALLOCATE(workspace)
! Deallocate grid information (FIXME)
IF (relaxONOFF==1) THEN
	DEALLOCATE(RelaxZones)
ENDIF
! GD: New variable definition, SWENSE part
IF(swenseONOFF==1) THEN !We deallocate the incident wavefields values
    DEALLOCATE(Wavefield%E_I,Wavefield%Et_I,Wavefield%P_I_s,Wavefield%Pz_I_s,Wavefield%Pt_I_s,Wavefield%E_I_bp,Wavefield%Pz_I_bp)
    IF (FineGrid%Nx>1) THEN
    	DEALLOCATE(Wavefield%Ex_I, Wavefield%Exx_I, Wavefield%Px_I_s, Wavefield%Ex_I_bp,Wavefield%Px_I_bp)
    ENDIF
    IF (FineGrid%Ny>1) THEN
    	DEALLOCATE(Wavefield%Ey_I, Wavefield%Eyy_I, Wavefield%Py_I_s, Wavefield%Ey_I_bp, Wavefield%Py_I_bp)
	ENDIF
! Deallocate the surfacic source terms
! FIXME: Define even if non-SWENSE ?
! APEK: moved lines below inside loop as these are only allocated in SWENSE?
IF (FineGrid%Nx>1) THEN
    DEALLOCATE(Wavefield%SourceEx, Wavefield%SourcePx)
ENDIF
IF (FineGrid%Ny>1) THEN
    DEALLOCATE(Wavefield%SourceEy, Wavefield%SourcePy)
ENDIF
ENDIF
END SUBROUTINE CloseVariables
