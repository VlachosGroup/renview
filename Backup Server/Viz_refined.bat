@ECHO OFF

IF %ERRORLEVEL% EQU 0 (

	echo.
	echo Generating Refined Visualization File now!
	COPY RefinedVisualization.out graphviz\bin\. >NUL
	
	cd graphviz\bin\
	unflatten -c 10 RefinedVisualization.out | dot -Tsvg RefinedVisualization.out -o RefinedVisualization.svg
	COPY RefinedVisualization.svg ..\..\. >NUL
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
