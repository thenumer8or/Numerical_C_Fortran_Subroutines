C-----------------------------------------------------------------------
C  IMSL Name:  TDATE (Single precision version)
C
C  Computer:   PCDSMS/SINGLE
C
C  Revised:    November 5, 1985
C
C  Purpose:    Get today's date.
C
C  Usage:      CALL TDATE (IDAY, MONTH, IYEAR)
C
C  Arguments:
C     IDAY   - Day of the month.  (Output)
C              IDAY is between 1 and 31 inclusive.
C     MONTH  - Month of the year.  (Output)
C              MONTH is between 1 and 12 inclusive.
C     IYEAR  - AD year.  (Output)
C              For example, IYEAR = 1985.
C
C  GAMS:       R
C
C  Chapters:   STAT/LIBRARY Utilities
C              MATH/LIBRARY Utilities
C
C  Copyright:  1985 by IMSL, Inc.  All Rights Reserved.
C
C  Warranty:   IMSL warrants only that IMSL testing has been applied
C              to this code.  No other warranty, expressed or implied,
C              is applicable.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE TDATE (IDAY, MONTH, IYEAR)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER    IDAY, MONTH, IYEAR
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER    *2 DY, MH, YR
C                                  SPECIFICATIONS FOR SUBROUTINES
      EXTERNAL   GETDAT
C
      CALL GETDAT(YR, MH, DY)
      IDAY = DY
      MONTH = MH
      IYEAR = YR
C
      RETURN
      END
