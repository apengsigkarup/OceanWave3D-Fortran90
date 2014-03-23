      SUBROUTINE setDeltaTOCW(dtOF)
      ! Subroutine setting setting time increment dt in OceanWave3D
      !
      ! Input: dtOF:= Time increment from OpenFOAM. dt =
      ! min(min(OFnew,OFMAX),OCWdt)
      !
      ! Output: Void
      !
      ! Written by Bo Terp Paulsen, botp (botp@mek.dtu.dk)

      USE GlobalVariables
      IMPLICIT NONE
      REAL(KIND=long) :: dtOF

      ! Setting time increment in OceanWave3D
      !
      dt = dtOF

      END SUBROUTINE setDeltaTOCW

