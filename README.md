Reaction Network Viewer (ReNView)
=================================

The **Re**action **N**etwork **View**er (ReNView) generates a graphic representation of the reaction fluxes within the system essential for identifying dominant reaction pathways and mechanism reduction.

Documentation
-------------

See our [documentation page](https://github.com/VlachosGroup/ReNView/wiki/Documentation) for examples, and equations used.

Developer
---------
Udit Gupta (uditgupta0912@gmail.com)

Dependencies
------------

- Python3
- [Numpy](http://www.numpy.org/) : Used for vector and matrix operations
- [Pandas](https://pandas.pydata.org/): Used to import data from input files and process headers
- [Graphviz](https://www.graphviz.org/): Used to generate visualizations from text files

Getting Started
---------------
1. Install using pip::

  >  pip install ReNView
  
Files for using the Visualization tool
--------------------------------------
1) species_comp.out - A species composition file specifying species name, phase, and the elemental composition of the molecule. In case of heterogeneous systems, surface coverages can also be provided for node coloring.

2) reaction_rates.out - A reactions file specifying the forward, reverse, net rate, partial equilibrium index, and reaction string. The reaction string should contain species names as mentioned in the species_comp.out file.

Features in the visualization tool:
-----------------------------------
- Network Visualization
- Species Visualization
- Legend generation

License
-------

This project is licensed under the GNU LGPL License - see the [LICENSE.md](https://github.com/VlachosGroup/ReNView/blob/master/LICENSE.md) file for details

Publications
------------
- Gupta, U.; Vlachos, D. G. Reaction Network Viewer (ReNView): An Open-Source Framework for Reaction Path Visualization of Chemical Reaction Systems. SoftwareX 2020, 11, 100442. https://doi.org/10.1016/J.SOFTX.2020.100442.

Contributing
------------

If you have a suggestion or find a bug, please post to our Issues page with 
the enhancement or bug tag respectively.

Finally, if you would like to add to the body of code, please:

- fork the development branch
- make the desired changes
- write the appropriate unit tests
- submit a pull request.

Questions
---------

If you are having issues, please post to our Issues page with the 
help wanted or question tag. We will do our best to assist.

Funding
-------

This material is based upon work supported by the Department of Energy's Office 
of Energy Efficient and Renewable Energy's Advanced Manufacturing Office under 
Award Number DE-EE0007888-9.5.

Special Thanks
--------------

-  Dr. Jeffrey Frey (capsule compatibility)
-  Gerhard Wittreich (testing)
-  Hilal Ezgi Toraman (testing)
-  Jonathan Lym (testing)
-  Jaynell Keely (Logo design)
