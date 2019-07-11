# Visualization Tool
The **Re**action **N**etwork **View**er (ReNView) generates a graphic representation of the reaction fluxes within the system essential for identifying dominant reaction pathways and mechanism reduction.

# Developer
<ul>
  <li> Udit Gupta (<ugupta@udel.edu>) <\li>
<\ul>

# Dependencies
<ul>
  <li> Python3 <\li>
  <li> Numpy: <\li>
  <li> Pandas: <\li>
  <li> Graphviz <\li>
<\ul>

Files for using the Visualization tool
1) species_comp.out - A species composition file specifying species name, phase, and the elemental composition of the molecule. In case of heterogeneous systems, surface coverages can also be provided for node coloring.

2) reaction_rates.out - A reactions file specifying the forward, reverse, net rate, partial equilibrium index, and reaction string. The reaction string should contain species names as mentioned in the species_comp.out file.



Features in the visualization tool:
1) 
