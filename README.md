Reaction Network Viewer (ReNView)
=================================

The **Re**\ action **N**\ etwork **View**\ er (ReNView) generates a graphic representation of the reaction fluxes within the system essential for identifying dominant reaction pathways and mechanism reduction.

Documentation
-------------

See our <documentation page> for examples, and equations used.

Developer
---------
Udit Gupta (<ugupta@udel.edu>)

Dependencies
------------

- Python3
- <Numpy> : Used for vector and matrix operations
- Pandas:
- Graphviz:

Files for using the Visualization tool
1) species_comp.out - A species composition file specifying species name, phase, and the elemental composition of the molecule. In case of heterogeneous systems, surface coverages can also be provided for node coloring.

2) reaction_rates.out - A reactions file specifying the forward, reverse, net rate, partial equilibrium index, and reaction string. The reaction string should contain species names as mentioned in the species_comp.out file.



Features in the visualization tool:
1) 
