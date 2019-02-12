#!/bin/bash

#  This script is designed to copy input and output files into a unique
#+ directory within the parent directory in order to facilitate 
#+ keeping track of individual job runs. This should allow a user to
#+ recreate a given run with a minimum amount of effort. The user may
#+ optionally specify a name for the output directory as a command-line
#+ argument, or the script will supply a default based on the current
#+ date and time.
#+ Initially written by Jonathan E. Sutton, Vlachos Research Group,
#+ University of Delaware, 2012/08/23.
#+ Last updated: 2013/06/13.

#  Define the path to the binary. For testing, this will be the local directory,
#+ but for production runs, this will be the location of the compiled binary
#+ on the server/cluster.
CKPATH='/home/chemkin/chemkin_repo_testing_new/ugupta/UploadFiles'

#printf 'Should previous results be discarded? Y/N: '
#read delete_output_str
#if [[ 'Y' == "$delete_output_str" || 'y' == "$delete_output_str" ]]; then
  delete_output=true
#else
#  delete_output=false
#fi

#make

#make clean_obj

#Convert all tabs to spaces in the input directory files
cd INP.d
for i in *; do
  sed -i 's/	/  /g' "$i" #First whitespace character is a literal tab
done
cd ..

#Create the output directory if needed
if [ ! -d "OUT.d" ]; then
  mkdir "OUT.d"
fi

#Create the species directory if needed
if [ ! -d "Species.d" ]; then
  mkdir "Species.d"
fi

#Remove old linking files, etc. and redirect error messages to /dev/null
#in order to suppress them.
if $delete_output; then
  rm INP.d/*link 2> /dev/null
  rm OUT.d/* 2> /dev/null
fi

#Prepare the model input
echo 'Preparing Chemkin binary linking files'
"$CKPATH/CHEMKIN.d/ckinterp.x"
"$CKPATH/CHEMKIN.d/skinterp.x"

echo 'Checking gas.inp'
grep -q "NO ERRORS" OUT.d/gas.out
if [ $? -ne 0 ]; then
  echo 'There is an error in gas.inp, exiting'
  exit
else
  echo 'No errors found in gas.inp'
fi

echo 'Checking surf.inp'
grep -q "NO ERRORS" OUT.d/surf.out
if [ $? -ne 0 ]; then
  echo 'There is an error in surf.inp, exiting'
  exit
else
  echo 'No errors found in surf.inp'
fi

echo -n 'Do you wish to run the model now? Y/N: '
read run_flag
#run_flag="y"
if [ 'Y' = $run_flag ] || [ 'y' = $run_flag ]; then
  echo 'Running model'
  "$CKPATH""/reactors.x"
  echo 'Done running model'

  #Optionally create the directory for the input and output.
  echo -n 'Do you wish to save the output to a new directory? Y/N: '
  read save_flag
#save_flag='n'
  if [ 'Y' = $save_flag ] || [ 'y' = $save_flag ]; then
    echo 'Creating a save directory'
    date_time=$(date +%Y.%m.%d_%H.%M.%S.%N)
    echo -n 'Enter a directory name or press return for default: '
    read dir_name
    if [ ${#dir_name} -eq 0 ]; then #Use the default value
      dir_name="run_""$date_time"
    fi
    while [ -e "$dir_name" ]; do
      echo 'Directory exists already, try again'
      date_time=$(date +%Y.%m.%d_%H.%M.%S.%N)
      echo -n 'Enter a directory name or press return for default: '
      read dir_name
      if [ ${#dir_name} -eq 0 ]; then #Use the default value
        dir_name="run_""$date_time"
      fi
    done
    mkdir $dir_name

    #Copy the input & output files.
    echo "Copying files to ${dir_name}"
    cp 'chemkin.sh' "${dir_name}/"
    cp -r *.d "${dir_name}/"
    echo 'Done'
  else
    echo 'Output will not be saved to a unique directory'
  fi

  echo 'Done running script'
else
  echo 'Model will not be run, exiting'
  exit
fi

echo -n 'Do you wish to visualize the reaction network? Y/N: '
read viz_flag
if [ 'Y' = $viz_flag ] || [ 'y' = $viz_flag ]; then
	#generate the svg file for visualization
	echo 'Generating Visualization Files now!'
	cd OUT.d
	cp rpa_visualizationBasis100.out ../graphviz/bin/
	cp rpa_visualizationNormalized.out ../graphviz/bin/
	cp rpa_visualizationMaxRateNormalized.out ../graphviz/bin/
	cd ..
	cp Species.d/* graphviz/bin/
	cd graphviz/bin/
	#unflatten -c 10 rpa_visualization.out | dot -Tsvg -o rpa_visualization.svg
	for i in *.out; do 
		unflatten -c 10 $i | dot -Tsvg -o ${i%.*}.svg
		rm $i 2> /dev/null
		cp ${i%.*}.svg ../../Species.d/
		rm ${i%.*}.svg 2> /dev/null
		rm $i 2> /dev/null
	done
	cd ../../Species.d/
	for i in *.out; do
		rm $i 2> /dev/null
	done
else
	echo 'No Visualization files generated'
	exit
fi



