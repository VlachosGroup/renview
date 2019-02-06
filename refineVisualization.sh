#!/bin/bash

# This script is designed to refine the visualization master file generated 
#+ during simulation in CHEMKIN. Chemkin generates a visualization master file
#+ including edges and nodes corresponding to the Cutoffrate and EquilibriumTol 
#+ provided in this file.
#+ Author: Udit Gupta, Vlachos Research Group,
#+ University of Delaware, 2018/06/21.

Cutoffrate=5.0E-06
EquilibriumTol=0.05
NormalizationIndex=1
EquilLower=0.5-$EquilibriumTol
EquilUpper=0.5+$EquilibriumTol
NodeSep=0.25
RankSep=0.25
SplinesType="ortho" #Specify either ortho or spline
ActiveSpecies=" "
ActiveSpeciesNumber="0"
ReactionRate=" "
GasSpecies=" "
GasReactions=" "
NormalizationFactor=" "
TotalSpecies=" "
c=0
count=0
if [[ "$NormalizationIndex" =~ 1 ]]; then
	NormalizationFactor=" NetRateInitial "
elif [[ "$NormalizationIndex" =~ 2 ]]; then
	NormalizationFactor=" MaxRateNetwork "
fi
input="OUT.d/rpa_visualizationMaster2.out"
if [ ! -d "RefinedSpecies.d" ]; then
  mkdir "RefinedSpecies.d"
fi
echo "digraph G {" > RefinedVisualization.out
echo "graph [splines = "$SplinesType", bgcolor=lightgray, resolution=64, fontname=Arial, fontcolor=blue, fontsize=48];" >> RefinedVisualization.out
echo "node [shape=rectangle, style=filled, fontsize=35,width=0,height=0,fillcolor=cornsilk,shape=plaintext];" >> RefinedVisualization.out
echo "edge [fontsize=30];" >> RefinedVisualization.out
echo "label = \"Reaction Path Analysis\";" >> RefinedVisualization.out
echo "labelloc = \"t\";" >> RefinedVisualization.out
echo "center=1;" >> RefinedVisualization.out
echo "size=\"20,20\";" >> RefinedVisualization.out
echo "ranksep=\""$RankSep" equally\";" >> RefinedVisualization.out
echo "nodesep=\""$NodeSep" equally\";" >> RefinedVisualization.out
echo "rankdir=TB;" >> RefinedVisualization.out
echo "bgcolor=white;" >> RefinedVisualization.out



while IFS= read -r line
do 
	arr=($line)
	if [[ " ${arr[0]} " =~ "$NormalizationFactor" ]]; then
			ReactionRate="${arr[1]}"
	elif [[ " ${arr[0]} " =~ " GasSpecies " ]]; then
			GasSpecies="${arr[1]}"
	elif [[ " ${arr[0]} " =~ " GasReactions " ]]; then
			GasReactions="${arr[1]}"
	fi
	if [[ " ${arr[0]} " =~ " TotalSpecies " ]]; then
			TotalSpecies="${arr[1]}"
			for (( i=1; i<=TotalSpecies; i++ ))
			do 
				echo "digraph G {" > "RefinedSpecies.d/"$i".out"
				echo "graph [splines = "$SplinesType", bgcolor=lightgray, resolution=64, fontname=Arial, fontcolor=blue, fontsize=36];" >> "RefinedSpecies.d/"$i".out"
				echo "node [fontsize=12];" >> "RefinedSpecies.d/"$i".out"
				echo "edge [fontsize=24];" >> "RefinedSpecies.d/"$i".out"
				echo "center=1;" >> "RefinedSpecies.d/"$i".out"
				echo "size=\"10,10\";" >> "RefinedSpecies.d/"$i".out"
				echo "ranksep=\""$RankSep" equally\";" >> "RefinedSpecies.d/"$i".out"
				echo "nodesep=\""$NodeSep" equally\";" >> "RefinedSpecies.d/"$i".out"
				echo "rankdir=LR;" >> "RefinedSpecies.d/"$i".out"
				echo "bgcolor=white;" >> "RefinedSpecies.d/"$i".out"
			done
	fi
	if [[ " ${arr[0]} " =~ " Species " ]]; then
		if [[ ActiveSpecies != "${arr[1]}" ]]; then
			if [[ $count > 0 ]]; then
				if (( $(echo $(awk 'BEGIN{print ('$ActiveSpeciesNumber'<='$GasSpecies')?1:0}')) )); then 
					echo $ActiveSpecies"[shape=rectangle,style=filled,fontsize=35,width=0,height=0,fillcolor=cornsilk,URL=\""$ActiveSpeciesNumber".svg\",shape=plaintext];" >> RefinedVisualization.out
					echo $ActiveSpecies"[shape=rectangle,style=filled,fontsize=35,width=0,height=0,fillcolor=cornsilk,URL=\""$ActiveSpeciesNumber".svg\",shape=plaintext];" >> "RefinedSpecies.d/"$ActiveSpeciesNumber".out"
				else
					echo $ActiveSpecies"[shape=rectangle,style=filled,fontsize=35,width=0,height=0,fillcolor=azure,URL=\""$ActiveSpeciesNumber".svg\",shape=plaintext];" >> RefinedVisualization.out
					echo $ActiveSpecies"[shape=rectangle,style=filled,fontsize=35,width=0,height=0,fillcolor=azure,URL=\""$ActiveSpeciesNumber".svg\",shape=plaintext];" >> "RefinedSpecies.d/"$ActiveSpeciesNumber".out"
				fi
			fi
			count=0
			ActiveSpecies="${arr[1]}"
			ActiveSpeciesNumber="${arr[3]}"
			
		fi
	fi
	if [[ " ${arr[0]} " =~ " Rate " ]]; then
		if (( $(echo $(awk 'BEGIN{print ('${arr[1]}'>'$Cutoffrate')?1:0}')) )); then
			count=count+1
			RESULT1=$(awk -v dividend="${arr[1]}" -v divisor="$ReactionRate" 'BEGIN {printf "%.0f", ((dividend*100)/(divisor*10)+1); exit(0)}')
			RESULT2=$(awk -v dividend="${arr[1]}" -v divisor="$ReactionRate" 'BEGIN {printf "%.0f", ((dividend*100)/(divisor*120)+1); exit(0)}')
			RESULT3=$(awk -v dividend="${arr[1]}" -v divisor="$ReactionRate" 'BEGIN {printf "%.0f", ((dividend*100/divisor)); exit(0)}')
			RESULT4=$(awk -v dividend="${arr[11]}" -v divisor=20 'BEGIN {printf "%.0f", (-dividend/divisor)+1; exit(0)}')
			RESULT5=$(awk -v dividend="${arr[11]}" -v divisor=60 'BEGIN {printf "%.0f", (-dividend/divisor)+1; exit(0)}')
			RESULT6=$(awk -v dividend="${arr[11]}" -v divisor=1 'BEGIN {printf "%.2f", (dividend/divisor); exit(0)}')
			RESULT7=$(awk -v dividend="${arr[7]}" -v divisor=1 'BEGIN {printf "%.2f", (dividend/divisor); exit(0)}')
			RESULT8=$(awk -v dividend="${arr[17]}" -v divisor=1 'BEGIN {printf "%.4f", (dividend/divisor); exit(0)}')
			RESULT9=$(awk -v dividend="${arr[15]}" -v divisor=1 'BEGIN {printf "%.4f", (dividend/divisor); exit(0)}')
			RESULT10=$(awk -v dividend="${arr[13]}" -v divisor=1 'BEGIN {printf "%.0f", (dividend/divisor); exit(0)}')
			RESULT11=$(awk -v dividend="${arr[19]}" -v divisor=20 'BEGIN {printf "%.0f", (dividend/divisor)+1; exit(0)}')
			RESULT12=$(awk -v dividend="${arr[19]}" -v divisor=60 'BEGIN {printf "%.0f", (dividend/divisor)+1; exit(0)}')
			RESULT13=$(awk -v dividend="${arr[19]}" -v divisor=1 'BEGIN {printf "%.2f", (dividend/divisor); exit(0)}')
			
			if (( $(echo $(awk 'BEGIN{print ('${arr[7]}'>='$EquilLower')?1:0}')) )) && (( $(echo $(awk 'BEGIN{print ('${arr[7]}'<='$EquilUpper')?1:0}')) )); then
				echo "edge[dir=\"forward\",style=\"setlinewidth("$RESULT1")\",color=green,weight=2,arrowsize="$RESULT2",label=\"  "$RESULT3"%\"];" >> RefinedVisualization.out
				echo "edge[dir=\"forward\",style=\"setlinewidth("$RESULT4")\",color=green,weight=2,arrowsize="$RESULT5",label=\"  "${arr[9]}"   "$RESULT6" %  "$RESULT8" (mol/cm3)^"$RESULT10" "$RESULT7" \"];" >> "RefinedSpecies.d/"$ActiveSpeciesNumber".out"
				echo "edge[dir=\"forward\",style=\"setlinewidth("$RESULT11")\",color=green,weight=2,arrowsize="$RESULT12",label=\"  "${arr[9]}"   "$RESULT13" %  "$RESULT8" (mol/cm3)^"$RESULT10" "$RESULT7" \"];" >> "RefinedSpecies.d/"${arr[5]}".out"
			elif (( $(echo $(awk 'BEGIN{print ('${arr[9]}'<='$GasReactions')?1:0}')) )); then 
				echo "edge[dir=\"forward\",style=\"setlinewidth("$RESULT1")\",color=red,weight=2,arrowsize="$RESULT2",label=\"  "$RESULT3"%\"];" >> RefinedVisualization.out
				echo "edge[dir=\"forward\",style=\"setlinewidth("$RESULT4")\",color=red,weight=2,arrowsize="$RESULT5",label=\"  "${arr[9]}"   "$RESULT6" %  "$RESULT9" mol/s "$RESULT7" \"];" >> "RefinedSpecies.d/"$ActiveSpeciesNumber".out"
				echo "edge[dir=\"forward\",style=\"setlinewidth("$RESULT11")\",color=red,weight=2,arrowsize="$RESULT12",label=\"  "${arr[9]}"   "$RESULT13" %  "$RESULT9" mol/s "$RESULT7" \"];" >> "RefinedSpecies.d/"${arr[5]}".out"
			else 
				echo "edge[dir=\"forward\",style=\"setlinewidth("$RESULT1")\",color=black,weight=2,arrowsize="$RESULT2",label=\"  "$RESULT3"%\"];" >> RefinedVisualization.out
				echo "edge[dir=\"forward\",style=\"setlinewidth("$RESULT4")\",color=black,weight=2,arrowsize="$RESULT5",label=\"  "${arr[9]}"   "$RESULT6" %  "$RESULT9" mol/s "$RESULT7" \"];" >> "RefinedSpecies.d/"$ActiveSpeciesNumber".out"
				echo "edge[dir=\"forward\",style=\"setlinewidth("$RESULT11")\",color=black,weight=2,arrowsize="$RESULT12",label=\"  "${arr[9]}"   "$RESULT13" %  "$RESULT9" mol/s "$RESULT7" \"];" >> "RefinedSpecies.d/"${arr[5]}".out"
			fi
			echo "$ActiveSpecies->${arr[3]}" >> RefinedVisualization.out
			echo "$ActiveSpecies->${arr[3]}" >> "RefinedSpecies.d/"$ActiveSpeciesNumber".out"
			echo "$ActiveSpecies->${arr[3]}" >> "RefinedSpecies.d/"${arr[5]}".out"
			if (( $(echo $(awk 'BEGIN{print ('${arr[5]}'<='$GasSpecies')?1:0}')) )); then
				echo ${arr[3]}"[shape=rectangle,style=filled,fontsize=35,width=0,height=0,fillcolor=cornsilk,URL=\""${arr[5]}".svg\",shape=plaintext];" >> RefinedVisualization.out
				echo ${arr[3]}"[shape=rectangle,style=filled,fontsize=35,width=0,height=0,fillcolor=cornsilk,URL=\""${arr[5]}".svg\",shape=plaintext];" >> "RefinedSpecies.d/"$ActiveSpeciesNumber".out"
				echo ${arr[3]}"[shape=rectangle,style=filled,fontsize=35,width=0,height=0,fillcolor=cornsilk,URL=\""${arr[5]}".svg\",shape=plaintext];" >> "RefinedSpecies.d/"${arr[5]}".out"
				echo $ActiveSpecies"[shape=rectangle,style=filled,fontsize=35,width=0,height=0,fillcolor=cornsilk,URL=\""$ActiveSpeciesNumber".svg\",shape=plaintext];" >> "RefinedSpecies.d/"${arr[5]}".out"
			else
				echo ${arr[3]}"[shape=rectangle,style=filled,fontsize=35,width=0,height=0,fillcolor=azure,URL=\""${arr[5]}".svg\",shape=plaintext];" >> RefinedVisualization.out
				echo ${arr[3]}"[shape=rectangle,style=filled,fontsize=35,width=0,height=0,fillcolor=azure,URL=\""${arr[5]}".svg\",shape=plaintext];" >> "RefinedSpecies.d/"$ActiveSpeciesNumber".out"
				echo ${arr[3]}"[shape=rectangle,style=filled,fontsize=35,width=0,height=0,fillcolor=azure,URL=\""${arr[5]}".svg\",shape=plaintext];" >> "RefinedSpecies.d/"${arr[5]}".out"
				echo $ActiveSpecies"[shape=rectangle,style=filled,fontsize=35,width=0,height=0,fillcolor=azure,URL=\""$ActiveSpeciesNumber".svg\",shape=plaintext];" >> "RefinedSpecies.d/"${arr[5]}".out"
			fi
		fi
	fi
done < "$input"
echo "}" >> RefinedVisualization.out
for (( i=1; i<=TotalSpecies; i++ ))
do 
	echo "}" >> "RefinedSpecies.d/"$i".out"
done

echo "Finished generating visualization file!"
cp RefinedVisualization.out graphviz/bin/
cp RefinedSpecies.d/* graphviz/bin/
cd graphviz/bin/
#unflatten -c 10 RefinedVisualization.out | dot -Tsvg -o RefinedVisualization.svg
for i in *.out; do 
	unflatten -c 10 $i | dot -Tsvg -o ${i%.*}.svg
	rm $i 2> /dev/null
	cp ${i%.*}.svg ../../RefinedSpecies.d/
	rm ${i%.*}.svg 2> /dev/null
	rm $i 2> /dev/null
done
#cp RefinedVisualization.svg ../../.
#rm RefinedVisualization.out 2> /dev/null
#rm RefinedVisualization.svg 2> /dev/null
cd ../../RefinedSpecies.d/
for i in *.out; do
	rm $i 2> /dev/null
done
echo "Created Visualization file RefinedVisualization.svg!"