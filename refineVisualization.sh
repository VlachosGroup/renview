#!/bin/bash

# This script is designed to refine the visualization master file generated 
#+ during simulation in CHEMKIN. Chemkin generates a visualization master file
#+ including edges 
#+ Author: Udit Gupta, Vlachos Research Group,
#+ University of Delaware, 2018/06/21.

Cutoffrate=1.0E-08
EquilibriumTol=0.1
EquilLower=0.5-$EquilibriumTol
EquilUpper=0.5+$EquilibriumTol
ActiveSpecies=" "
ActiveSpeciesNumber=" "
MaxReactionRate=" "
GasSpecies=" "
GasReactions=" "
c=0
count=0
input="OUT.d/rpa_visualizationMaster.out"
echo "strict digraph G {" > RefinedVisualization.out
echo "graph [bgcolor=lightgray, resolution=64, fontname=Arial, fontcolor=blue, fontsize=48];" >> RefinedVisualization.out
echo "node [shape=rectangle, style=filled, fontsize=35,width=0,height=0,fillcolor=cornsilk,shape=plaintext];" >> RefinedVisualization.out
echo "edge [fontsize=30];" >> RefinedVisualization.out
echo "label = \"Reaction Path Analysis\";" >> RefinedVisualization.out
echo "labelloc = \"t\";" >> RefinedVisualization.out
echo "center=1;" >> RefinedVisualization.out
echo "size=\"20,20\";" >> RefinedVisualization.out
echo "ranksep=\"1.0 equally\";" >> RefinedVisualization.out
echo "nodesep=\"1.0 equally\";" >> RefinedVisualization.out
echo "rankdir=TB;" >> RefinedVisualization.out
echo "bgcolor=white;" >> RefinedVisualization.out

while IFS= read -r line
do 
	arr=($line)
	if [[ " ${arr[0]} " =~ " MaxRate " ]]; then
			MaxReactionRate="${arr[1]}"
	elif [[ " ${arr[0]} " =~ " GasSpecies " ]]; then
			GasSpecies="${arr[1]}"
	elif [[ " ${arr[0]} " =~ " GasReactions " ]]; then
			GasReactions="${arr[1]}"
	fi
	if [[ " ${arr[0]} " =~ " Species " ]]; then
		if [[ ActiveSpecies != "${arr[1]}" ]]; then
			if [[ $count > 0 ]]; then
				if (( $(echo $(awk 'BEGIN{print ('$ActiveSpeciesNumber'<='$GasSpecies')?1:0}')) )); then 
					echo $ActiveSpecies"[shape=rectangle,style=filled,fontsize=35,width=0,height=0,fillcolor=cornsilk,URL=\""$ActiveSpeciesNumber".svg\",shape=plaintext];" >> RefinedVisualization.out
				else
					echo $ActiveSpecies"[shape=rectangle,style=filled,fontsize=35,width=0,height=0,fillcolor=azure,URL=\""$ActiveSpeciesNumber".svg\",shape=plaintext];" >> RefinedVisualization.out
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
			if (( $(echo $(awk 'BEGIN{print ('${arr[7]}'>='$EquilLower')?1:0}')) )) && (( $(echo $(awk 'BEGIN{print ('${arr[7]}'<='$EquilUpper')?1:0}')) )); then
				RESULT1=$(awk -v dividend="${arr[1]}" -v divisor="$MaxReactionRate" 'BEGIN {printf "%.0f", (dividend*100)/(divisor*5); exit(0)}')
				RESULT2=$(awk -v dividend="${arr[1]}" -v divisor="$MaxReactionRate" 'BEGIN {printf "%.0f", ((dividend*100)/(divisor*30)+1); exit(0)}')
				RESULT3=$(awk -v dividend="${arr[1]}" -v divisor="$MaxReactionRate" 'BEGIN {printf "%.0f", dividend*100/divisor; exit(0)}')
				echo "edge[dir=\"forward\",style=\"setlinewidth("$RESULT1")\",color=green,weight=2,arrowsize="$RESULT2",label=\"  "$RESULT3"%\"];" >> RefinedVisualization.out
			elif (( $(echo $(awk 'BEGIN{print ('${arr[9]}'<='$GasReactions')?1:0}')) )); then 
				RESULT1=$(awk -v dividend="${arr[1]}" -v divisor="$MaxReactionRate" 'BEGIN {printf "%.0f", (dividend*100)/(divisor*5); exit(0)}')
				RESULT2=$(awk -v dividend="${arr[1]}" -v divisor="$MaxReactionRate" 'BEGIN {printf "%.0f", ((dividend*100)/(divisor*30)+1); exit(0)}')
				RESULT3=$(awk -v dividend="${arr[1]}" -v divisor="$MaxReactionRate" 'BEGIN {printf "%.0f", dividend*100/divisor; exit(0)}')
				echo "edge[dir=\"forward\",style=\"setlinewidth("$RESULT1")\",color=red,weight=2,arrowsize="$RESULT2",label=\"  "$RESULT3"%\"];" >> RefinedVisualization.out
			else 
				RESULT1=$(awk -v dividend="${arr[1]}" -v divisor="$MaxReactionRate" 'BEGIN {printf "%.0f", (dividend*100)/(divisor*5); exit(0)}')
				RESULT2=$(awk -v dividend="${arr[1]}" -v divisor="$MaxReactionRate" 'BEGIN {printf "%.0f", ((dividend*100)/(divisor*30)+1); exit(0)}')
				RESULT3=$(awk -v dividend="${arr[1]}" -v divisor="$MaxReactionRate" 'BEGIN {printf "%.0f", dividend*100/divisor; exit(0)}')
				echo "edge[dir=\"forward\",style=\"setlinewidth("$RESULT1")\",color=black,weight=2,arrowsize="$RESULT2",label=\"  "$RESULT3"%\"];" >> RefinedVisualization.out
			fi
			echo "$ActiveSpecies->${arr[3]}" >> RefinedVisualization.out
			if (( $(echo $(awk 'BEGIN{print ('${arr[5]}'<='$GasSpecies')?1:0}')) )); then
				echo ${arr[3]}"[shape=rectangle,style=filled,fontsize=35,width=0,height=0,fillcolor=cornsilk,URL=\""${arr[5]}".svg\",shape=plaintext];" >> RefinedVisualization.out
			else
				echo ${arr[3]}"[shape=rectangle,style=filled,fontsize=35,width=0,height=0,fillcolor=azure,URL=\""${arr[5]}".svg\",shape=plaintext];" >> RefinedVisualization.out
			fi
		fi
	fi
done < "$input"

echo "}" >> RefinedVisualization.out
echo "Finished generating visualization file!"
cp RefinedVisualization.out graphviz/bin/
cd graphviz/bin/
unflatten -c 10 RefinedVisualization.out | dot -Tsvg -o RefinedVisualization.svg
cp RefinedVisualization.svg ../../.
rm RefinedVisualization.out 2> /dev/null
rm RefinedVisualization.svg 2> /dev/null
echo "Created Visualization file RefinedVisualization.svg!"