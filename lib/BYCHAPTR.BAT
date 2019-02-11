echo off
   if not exist MATHCORE.LIB goto note
   copy MATHCORE.LIB \LIB >nul
   del MATHCORE.LIB
   goto next
:note
   echo ÿ
   echo      You must install MATHCORE.LIB to access the IMSL MATH/LIBRARY!
:next
rem                      To include the BLAS in the chapter 9 library,
rem                      uncomment the next three lines.
rem   if exist ?MATH9.LIB errout LIB \LIB\MATH9 +DMATH9.LIB+SMATH9.LIB+BLAS.LIB; >nul
rem   if exist ?MATH9.LIB del ?MATH9.LIB
rem   if exist BLAS.LIB del BLAS.LIB
   if exist BLAS.LIB copy BLAS.LIB \LIB >nul
   if exist BLAS.LIB del BLAS.LIB
   if exist ?MATH10.LIB errout LIB \LIB\MATH10 +DMATH10.LIB+SMATH10.LIB; >nul
   if exist ?MATH10.LIB del ?MATH10.LIB
   if exist ?MATH9.LIB errout LIB \LIB\MATH9 +DMATH9.LIB+SMATH9.LIB; >nul
   if exist ?MATH9.LIB del ?MATH9.LIB
   if exist ?MATH8?.LIB errout LIB \LIB\MATH8/P:32 +DMATH8A.LIB+DMATH8B.LIB+SMATH8A.LIB+SMATH8B.LIB; >nul
   if exist ?MATH8?.LIB del ?MATH8?.LIB
   if exist ?MATH7.LIB errout LIB \LIB\MATH7 +DMATH7.LIB+SMATH7.LIB; >nul
   if exist ?MATH7.LIB del ?MATH7.LIB
   if exist ?MATH6.LIB errout LIB \LIB\MATH6 +DMATH6.LIB+SMATH6.LIB; >nul
   if exist ?MATH6.LIB del ?MATH6.LIB
   if exist ?MATH5?.LIB errout LIB \LIB\MATH5/P:32 +DMATH5A.LIB+DMATH5B.LIB+DMATH5C.LIB+SMATH5A.LIB+SMATH5B.LIB+SMATH5C.LIB; >nul
   if exist ?MATH5?.LIB del ?MATH5?.LIB
   if exist ?MATH4?.LIB errout LIB \LIB\MATH4 +DMATH4A.LIB+DMATH4B.LIB+SMATH4A.LIB+SMATH4B.LIB; >nul
   if exist ?MATH4?.LIB del ?MATH4?.LIB
   if exist ?MATH3?.LIB errout LIB \LIB\MATH3/P:32 +DMATH3A.LIB+DMATH3B.LIB+DMATH3C.LIB+SMATH3A.LIB+SMATH3B.LIB+SMATH3C.LIB; >nul
   if exist ?MATH3?.LIB del ?MATH3?.LIB
   if exist ?MATH2?.LIB errout LIB \LIB\MATH2 +DMATH2A.LIB+DMATH2B.LIB+SMATH2A.LIB+SMATH2B.LIB; >nul
   if exist ?MATH2?.LIB del ?MATH2?.LIB
   if exist ?MATH1?.LIB errout LIB \LIB\MATH1/P:32 +DMATH1A.LIB+DMATH1B.LIB+DMATH1C.LIB+SMATH1A.LIB+SMATH1B.LIB+SMATH1C.LIB; >nul
   if exist ?MATH1?.LIB del ?MATH1?.LIB
   echo ÿ
   echo      IMSL MATH/LIBRARY installation complete!
   exit
