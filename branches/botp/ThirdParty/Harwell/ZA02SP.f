C COPYRIGHT (c) 1970 AEA Technology and Council for the Central
C                    Laboratory of the Research Councils

      REAL FUNCTION ZA02AS(X)
      REAL X

C  Returns the current CPU time in seconds when compiled
C  on a Generic Unix machine.  This version is appropriate
C  for Compaq/Dec, HP, Silicon Graphics, Sun and Intel(Linux)
C  machines, using default compilers.

      REAL RRTIME(2)
      REAL ETIME
      EXTERNAL ETIME

C  Note: DTIME can be used for time difference
      ZA02AS = ETIME(RRTIME)

      END
