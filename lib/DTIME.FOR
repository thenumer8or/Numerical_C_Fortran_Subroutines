C-----------------------------------------------------------------------
C  IMSL Name:  DTIME (Single precision version)
C
C  Computer:   PCDSMS/SINGLE
C
C  Revised:    November 5, 1985
C
C  Purpose:    Get time of day.
C
C  Usage:      CALL DTIME (IHOUR, MINUTE, ISEC)
C
C  Arguments:
C     IHOUR  - Hour of the day.  (Output)
C              IHOUR is between 0 and 23 inclusive.
C     MINUTE - Minute within the hour.  (Output)
C              MINUTE is between 0 and 59 inclusive.
C     ISEC   - Second within the minute.  (Output)
C              ISEC is between 0 and 59 inclusive.
C
C  GAMS:       R
C
C  Chapter:    Reference Material - Time utilities
C
C  Copyright:  1985 by IMSL, Inc.  All Rights Reserved.
C
C  Warranty:   IMSL warrants only that IMSL testing has been applied
C              to this code.  No other warranty, expressed or implied,
C              is applicable.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE DTIME (IHOUR, MINUTE, ISEC)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER    IHOUR, MINUTE, ISEC
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER    *2 HR, MN, SC, IH
C                                  SPECIFICATIONS FOR SUBROUTINES
      EXTERNAL   GETTIM
C
      CALL GETTIM(HR, MN, SC, IH)
      IHOUR = HR
      MINUTE = MN
      ISEC = SC
C
      RETURN
      END
