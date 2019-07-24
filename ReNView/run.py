# -*- coding: utf-8 -*-

import os
import pydot
import numpy as np
import pandas as pd
import graphviz
import subprocess
from graphviz import ENGINES
from graphviz import Digraph
from subprocess import check_call

#from Reactions.reactions_import import *
from GraphGenerator.GraphGenerator import *

input_species_file('data/example_ammonia/species_comp.out')
input_reactions_file('data/example_ammonia/reaction_rates.out')
input_initial_reactant('NH3')
input_reaction_cutoffrate(1.0E-09)
input_elements_desired(['N', 'H'])
input_normalization(2)
input_output_directory('results/example_ammonia/')
generate_visualizations()

erase_data()