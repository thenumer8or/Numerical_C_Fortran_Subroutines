C-----------------------------------------------------------------------
C  IMSL Name:  CTIME (Single precision version)
C
C  Computer:   PCDSMS/SINGLE
C
C  Revised:    May 7, 1985
C
C  Purpose:    Return CPU time used in seconds.
C
C  Usage:      CTIME ()
C
C  Argument:
C     CTIME  - CPU time used (in seconds) since first call to CTIME.
C              (Output)
C
C  Remarks:
C  1. The first call to CTIME returns 0.0.
C
C  2. The accuracy of this routine depends on the hardware and the
C     operating system.  On some systems, identical runs can produce
C     timings differing by more than 10 percent.
C
C  GAMS:       S3
C
C  Chapter:    Reference Material - Time Utilities
C
C  Copyright:  1985 by IMSL, Inc.  All Rights Reserved.
C
C  Warranty:   IMSL warrants only that IMSL testing has been applied
C              to this code.  No other warranty, expressed or implied,
C              is applicable.
C
C-----------------------------------------------------------------------
C
      REAL FUNCTION CTIME ()
C                                  SPECIFICATIONS FOR SAVE VARIABLES
      REAL       BEGTIM
      SAVE       BEGTIM
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER    IER
      REAL       CPU
      INTEGER    *2 IHOUR, MINUTE, ISEC, IHSEC
C                                  SPECIFICATIONS FOR SUBROUTINES
      EXTERNAL   GETTIM
C
      DATA BEGTIM/-1.0E0/
C                                  Set defaults
      CTIME = 0.0E0
      IER = 0
C                                  Get CPU time in seconds
      CALL GETTIM(IHOUR, MINUTE, ISEC, IHSEC)
      CPU = (IHOUR*60+MINUTE)*60 + ISEC + IHSEC/100.0
      IF (CPU .LT. BEGTIM) CPU = CPU + 86400.0E0
C                                  Error are flagged by IER .NE. 0
      IF (IER .NE. 0) THEN
         CALL E1PSH ('CTIME ')
         CALL E1MES (5, 1, 'An error occured in obtaining the '//
     &               'runtime statistics.')
         CALL E1POP ('CTIME ')
C                                  Save initial CPU time
      ELSE IF (BEGTIM .LT. 0.0) THEN
         BEGTIM = CPU
C                                  Subtract initial CPU time
      ELSE
         CTIME = CPU - BEGTIM
      END IF
C
      RETURN
      END
