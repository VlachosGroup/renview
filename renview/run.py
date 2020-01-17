# -*- coding: utf-8 -*-

import os

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

def run_for_omkm(omkm_directory, initial_reactant, elements_desired,
                 reaction_cutoffrate=1.e-9, normalization=2):
    default_kwargs = {
        'species_file': os.path.join(omkm_directory, r'species.out'),
        'reactions_file': os.path.join(omkm_directory, r'rates_ss.out'),
        'initial_reactant': initial_reactant,
        'reaction_cutoffrate': reaction_cutoffrate,
        'elements_desired': elements_desired,
        'normalization': normalization,
        'output_directory': omkm_directory,
        }
    run(**default_kwargs)


if __name__ == "__main__":
    # Change the folder to the same directory as the script
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
