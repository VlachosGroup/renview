@ECHO OFF

IF %ERRORLEVEL% EQU 0 (

	echo.
	echo Generating Visualization Files now!
	cd OUT.d
	COPY rpa_visualization.out ..\graphviz\bin\. >NUL
	COPY rpa_visualizationNormalized.out ..\graphviz\bin\. >NUL
	COPY rpa_visualizationMaxRateNormalized.out ..\graphviz\bin\. >NUL
	cd ..
	
	cd Species.d
	COPY *.out ..\graphviz\bin\. >NUL
	cd ..\graphviz\bin\
	FOR %%i IN (*.out) do (unflatten -c 10 %%i | dot -Tsvg %%i -o %%~ni.svg)
	FOR %%j IN (*.svg) do (COPY *.svg ..\..\Species.d\. >NUL)
	FOR %%k IN (*.svg) do (del %%~nk.out)
	FOR %%l IN (*.svg) do (del %%l)
	cd ..\..\
	
 IF %ERRORLEVEL% EQU 0 (
  echo.
  echo Build Succeeded.
 ) ELSE (
  echo.
  echo Build failed.
 )

) ELSE (
 echo.
 echo 'No Visualization files generated'
)
