if [ -e ./Species.d ]; then
	if [ -n "$CHEMKIN_OUT_RENAME" ]; then
		mv ./Species.d ./Species.d-$(date +%Y%m%d%H%M%S.%N)
		if [ $? -ne 0 ]; then
			echo "ERROR:  failed to rename $(pwd)/OUT.d"
			exit 1
		fi
	else
		rm -rf ./Species.d
		if [ $? -ne 0 ]; then
			echo "ERROR:  failed to remove $(pwd)/OUT.d"
			exit 1
		fi
	fi
fi
mkdir ./Species.d
if [ $? -ne 0 ]; then
	echo "ERROR:  failed to create $(pwd)/OUT.d"
	exit 1
fi

CHEMKIN_GRAPHVIZ_POSTPROCESS() {
	local input_file exec_output rc rc_out

	if [ ! -d OUT.d ]; then
		[ -z "$UD_QUIET_JOB_SETUP" ] && echo "ERROR: OUT.d directory is not present"
		return 1
	fi
	[ -z "$UD_QUIET_JOB_SETUP" ] && echo "-- Generating visualization files now"
	pushd OUT.d >/dev/null 2>&1
	rc_out=0
	for input_file in rpa_visualizationBasis100.out rpa_visualizationNormalized.out rpa_visualizationMaxRateNormalized.out; do
		if [ -f "$input_file" ]; then
			echo -n "   + processing $input_file : "
			exec_output=$( { unflatten -c 10 "$input_file" | dot -Tsvg -o "${input_file%.out}.svg" ; } 2>&1)
			rc=$?
			if [ $rc -ne 0 ]; then
				[ -z "$UD_QUIET_JOB_SETUP" ] && echo "failed"
				[ -z "$UD_QUIET_JOB_SETUP" ] && printf "%s\n" "$exec_output"
				rc_out=$rc
			else
    			[ -z "$UD_QUIET_JOB_SETUP" ] && echo 'ok'
    		fi
		fi
	done
	popd >/dev/null 2>&1
	return $rc_out
}

CHEMKIN_GRAPHVIZ_POSTPROCESS1() {
	echo 'Generating Visualization Files now!'	
	cd OUT.d
	unflatten -c 10 rpa_visualizationBasis100.out | dot -Tsvg -o rpa_visualizationBasis100.svg
	unflatten -c 10 rpa_visualizationNormalized.out | dot -Tsvg -o rpa_visualizationNormalized.svg
	unflatten -c 10 rpa_visualizationMaxRateNormalized.out | dot -Tsvg -o rpa_visualizationMaxRateNormalized.svg
	cp rpa_visualizationBasis100.svg ../Species.d/
	cp rpa_visualizationNormalized.svg ../Species.d/
	cp rpa_visualizationMaxRateNormalized.svg ../Species.d/
	cd ..
	cd Species.d
	for i in *.out; do 
		unflatten -c 10 $i | dot -Tsvg -o ${i%.*}.svg
		#rm $i 2> /dev/null
	done
	echo 'Done generating visualizations!'
}

