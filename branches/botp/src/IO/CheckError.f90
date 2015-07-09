SUBROUTINE CheckError(STAT, err)
! By Allan P. Engsig-Karup.
IMPLICIT NONE
INTEGER :: STAT, err
IF (STAT/=0) THEN
   SELECT CASE (err)
      CASE (1)
	     WRITE (*,'(A)') 'Error: allocation of tmpArr failed.'
      CASE (2)
	     WRITE (*,'(A)') 'Error: memory allocation failed.'
	  CASE (3)
	     WRITE (*,'(A)') 'Error: could not allocate array for coarse-fine levels.'
	  CASE (4)
	  	 WRITE (*,'(A)') 'Error: Could not allocate variables in LUfactor.f90.'
	  CASE (5)
	     WRITE (*,'(A)') 'Error: Could not allocate SS(MAXS) for LU factorization.'
	  CASE (6)
	     WRITE (*,'(A)') 'Error: Could not allocate workspace memory for GMRES.'
	  CASE (7)
	     WRITE (*,'(A)') 'Error: Could not allocate necessary memory storage for SUBROUTINE PrepareFullOperatorStencils.'
	  CASE (8)
	     WRITE (*,'(A)') 'Error: Problem allocating sparse matrix arrays in SUBROUTINE GetSparseLinearLaplaceMatrix.'
      CASE (9) 
         WRITE (*,'(A)') 'Error: Problem allocating an index map for ghost points in SUBROUTINE MGPreProcess.'
      CASE (10)
         WRITE (*,'(A)') 'Error: PRoblem allocating FineGrid%dsigmanew in SUBROUTINE PreparePreconditioner.'      
      CASE (11)
         WRITE (*,'(A)') 'Error: PRoblem allocating FineGrid%dsigmanew in SUBROUTINE PreparePreconditioner.'
	  CASE DEFAULT
	     WRITE (*,'(A)') 'Error: unspecified error number (ERR) from calling routine.'
   END SELECT
   STOP ! FORCED STOP
ENDIF
END SUBROUTINE CheckError
