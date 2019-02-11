echo off
   if not exist MATHCORE.LIB goto nocore
   echo \LIB\IMSLIB/P:128 >LIB.BLD
   echo Y >>LIB.BLD
   if exist DMATH1A.LIB echo +DMATH1A.LIB& >>LIB.BLD
   if exist DMATH1B.LIB echo +DMATH1B.LIB& >>LIB.BLD
   if exist DMATH1C.LIB echo +DMATH1C.LIB& >>LIB.BLD
   if exist SMATH1A.LIB echo +SMATH1A.LIB& >>LIB.BLD
   if exist SMATH1B.LIB echo +SMATH1B.LIB& >>LIB.BLD
   if exist SMATH1C.LIB echo +SMATH1C.LIB& >>LIB.BLD
   if exist DMATH2A.LIB echo +DMATH2A.LIB& >>LIB.BLD
   if exist DMATH2B.LIB echo +DMATH2B.LIB& >>LIB.BLD
   if exist SMATH2A.LIB echo +SMATH2A.LIB& >>LIB.BLD
   if exist SMATH2B.LIB echo +SMATH2B.LIB& >>LIB.BLD
   if exist DMATH3A.LIB echo +DMATH3A.LIB& >>LIB.BLD
   if exist DMATH3B.LIB echo +DMATH3B.LIB& >>LIB.BLD
   if exist DMATH3C.LIB echo +DMATH3C.LIB& >>LIB.BLD
   if exist SMATH3A.LIB echo +SMATH3A.LIB& >>LIB.BLD
   if exist SMATH3B.LIB echo +SMATH3B.LIB& >>LIB.BLD
   if exist SMATH3C.LIB echo +SMATH3C.LIB& >>LIB.BLD
   if exist DMATH4A.LIB echo +DMATH4A.LIB& >>LIB.BLD
   if exist DMATH4B.LIB echo +DMATH4B.LIB& >>LIB.BLD
   if exist SMATH4A.LIB echo +SMATH4A.LIB& >>LIB.BLD
   if exist SMATH4B.LIB echo +SMATH4B.LIB& >>LIB.BLD
   if exist DMATH5A.LIB echo +DMATH5A.LIB& >>LIB.BLD
   if exist DMATH5B.LIB echo +DMATH5B.LIB& >>LIB.BLD
   if exist DMATH5C.LIB echo +DMATH5C.LIB& >>LIB.BLD
   if exist SMATH5A.LIB echo +SMATH5A.LIB& >>LIB.BLD
   if exist SMATH5B.LIB echo +SMATH5B.LIB& >>LIB.BLD
   if exist SMATH5C.LIB echo +SMATH5C.LIB& >>LIB.BLD
   if exist DMATH6.LIB  echo +DMATH6.LIB&  >>LIB.BLD
   if exist SMATH6.LIB  echo +SMATH6.LIB&  >>LIB.BLD
   if exist DMATH7.LIB  echo +DMATH7.LIB&  >>LIB.BLD
   if exist SMATH7.LIB  echo +SMATH7.LIB&  >>LIB.BLD
   if exist DMATH8A.LIB echo +DMATH8A.LIB& >>LIB.BLD
   if exist DMATH8B.LIB echo +DMATH8B.LIB& >>LIB.BLD
   if exist SMATH8A.LIB echo +SMATH8A.LIB& >>LIB.BLD
   if exist SMATH8B.LIB echo +SMATH8B.LIB& >>LIB.BLD
   if exist DMATH9.LIB  echo +DMATH9.LIB&  >>LIB.BLD
   if exist SMATH9.LIB  echo +SMATH9.LIB&  >>LIB.BLD
   if exist DMATH10.LIB echo +DMATH10.LIB& >>LIB.BLD
   if exist SMATH10.LIB echo +SMATH10.LIB& >>LIB.BLD
   if exist BLAS.LIB    echo +BLAS.LIB&    >>LIB.BLD
   echo +MATHCORE.LIB >>LIB.BLD
   echo NUL >>LIB.BLD
   echo ÿ
   echo       IMSL MATH/LIBRARY build in progress ...
   errout LIB @LIB.BLD >nul
   echo ÿ
   echo       IMSL MATH/LIBRARY installation complete!
   del LIB.BLD >nul
   goto complete
:nocore
   echo ÿ
   echo       The library build cannot complete because
   echo ÿ     the core library has not been installed.
   echo ÿ
   echo       Please install the core library before building the combined library.
   goto complete
:complete
   exit
