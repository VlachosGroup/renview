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

from renview import GraphGenerator as gg

def run(species_file, reactions_file, initial_reactant, reaction_cutoffrate,
        elements_desired, normalization, output_directory):
    gg.input_species_file(species_file)
    gg.input_reactions_file(reactions_file)
    gg.input_initial_reactant(initial_reactant)
    gg.input_reaction_cutoffrate(reaction_cutoffrate)
    gg.input_elements_desired(elements_desired)
    gg.input_normalization(normalization)
    gg.input_output_directory(output_directory)
    gg.generate_visualizations()
    gg.erase_data()

if __name__ == "__main__":
    try:
        os.chdir(os.path.dirname(__file__))
    except:
        pass
    default_kwargs = {
        'species_file': r'./data/example_ammonia/species_comp.out',
        'reactions_file': r'./data/example_ammonia/reaction_rates.out',
        'initial_reactant': 'NH3',
        'reaction_cutoffrate': 1.0E-09,
        'elements_desired': ['N', 'H'],
        'normalization': 2,
        'output_directory': r'./results/example_ammonia/',
        }
    run(**default_kwargs)
