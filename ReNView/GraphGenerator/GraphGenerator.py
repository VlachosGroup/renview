# -*- coding: utf-8 -*-
"""
Demonstration script for reading in reaction rates and reactions to generate
visualizations
Function definitions in ReNView
"""

import os
import numpy as np
import pandas as pd
from graphviz import *
from subprocess import check_call

#from Legend.legend import generate_legend

#Global variables specified here that are used for storing data
species_filename = ''
reactions_filename = ''
output_directory_name = ''
openmkm_species_list = []
openmkm_species_phase = []
openmkm_species_surf_cov = []
openmkm_species_composition = []
openmkm_reaction_strings = []
openmkm_reaction_reactants = []
openmkm_reaction_products = []
openmkm_reaction_fwdrate = []
openmkm_reaction_revrate = []
openmkm_reaction_netrate = []
openmkm_reaction_pei = []
Elemental_List = []
Normalization_requested = 0
openmkm_species_consider = []
stoich_matrix = np.array([])
rpa_local = np.array([])
Elements_Available = [] #gets populated based on the system under study
is_surface_cov_def = False
header_count = 0

# Inputs from OpenMKM or any other chemical kinetics software
InitialReactant = ''
Reaction_Rate_Cutoff = ''
Equilibrium_Tolerance = 0.05
Equilibrium_Upper = 0.5 + Equilibrium_Tolerance
Equilibrium_Lower = 0.5 - Equilibrium_Tolerance
Rank_Sep = 0.25
Node_Sep = 0.25

def generate_legend(fname, cov_def):
    '''
    This function generates the legend.
    
        fname : Legend filename
        cov_def : Surface coverages specification either True or False
    '''
    f = open(fname, "w+")
    f.write('digraph test {\n')
    f.write('graph [ratio=fill];\n')
    f.write('node [label=\"\\N\", fontsize=15, shape=plaintext];\n')
    f.write('graph [bb="0,0,352,154"];\n')
    
    f.write('legend [width=2,label=<\n')
    f.write('    <TABLE ALIGN="LEFT">\n')
    f.write('<TR>\n')
    f.write('      <TD COLSPAN="2"><B>Legend</B></TD>\n')
    f.write('     </TR>\n')
    if cov_def:
        f.write('<TR>\n')
        f.write('            <TD bgcolor="coral"></TD>\n')
        f.write('            <TD >0.85+</TD>\n')
        f.write('</TR>\n')
        f.write('<TR>\n')
        f.write('            <TD bgcolor="orangered"></TD>\n')
        f.write('            <TD >0.7-0.85</TD>\n')
        f.write('</TR>\n')
        f.write('<TR>\n')
        f.write('            <TD bgcolor="red"></TD>\n')
        f.write('            <TD >0.55-0.7</TD>\n')
        f.write('</TR>\n')
        f.write('<TR>\n')
        f.write('            <TD bgcolor="darkorange"></TD>\n')
        f.write('            <TD >0.4-0.55</TD>\n')
        f.write('</TR>\n')
        f.write('<TR>\n')
        f.write('            <TD bgcolor="orange"></TD>\n')
        f.write('            <TD >0.2-0.4</TD>\n')
        f.write('</TR>\n')
        f.write('<TR>\n')
        f.write('            <TD bgcolor="gold"></TD>\n')
        f.write('            <TD >0.1-0.2</TD>\n')
        f.write('</TR>\n')
        f.write('<TR>\n')
        f.write('            <TD bgcolor="lawngreen"></TD>\n')
        f.write('            <TD >0.05-0.1</TD>\n')
        f.write('</TR>\n')
        f.write('<TR>\n')
        f.write('            <TD bgcolor="green"></TD>\n')
        f.write('            <TD >1E-5-0.05</TD>\n')
        f.write('</TR>\n')
        f.write('<TR>\n')
        f.write('            <TD bgcolor="dodgerblue"></TD>\n')
        f.write('            <TD >1E-10-1E-5</TD>\n')
        f.write('</TR>\n')
        f.write('<TR>\n')
        f.write('            <TD bgcolor="purple"></TD>\n')
        f.write('            <TD >0.0-1E-10</TD>\n')
        f.write('</TR>\n')
    else:
        f.write('<TR>\n')
        f.write('            <TD bgcolor="cornsilk"></TD>\n')
        f.write('            <TD >Gas species</TD>\n')
        f.write('</TR>\n') 
        f.write('<TR>\n')
        f.write('            <TD bgcolor="azure"></TD>\n')
        f.write('            <TD >Surface species</TD>\n')
        f.write('</TR>\n')
        
    f.write('        \n')
    f.write('    </TABLE>\n')
    f.write('> ];\n')
    f.write('}\n')    
    f.close()

def input_species_file(fname):
    '''
    This function is used to specify the path to the species_comp file. The species_comp.out file contains the information regarding the species in     the reaction network. The file should contain the following headers - Species_name, Phase, Surf_cov, Element symbols. An example has been           listed in the previous section. 
    '''
    global species_filename
    species_filename = fname

def input_reactions_file(fname):
    '''
    This function is used to specify the path to the reaction_rates file. The reaction_rates.out contains the information regarding the reactions       in the reaction network. The file should contain the following headers - forward rate (Fwd_Rate), reverse rate (Rev_Rate), net rate (Net_Rate),     partial equilibrium index (PEI), and reaction string (Reaction_string). An example has been listed in the previous section.
    '''
    global reactions_filename
    reactions_filename = fname
    
def input_initial_reactant(reac_string):
    '''
    This function is used to specify the inlet reactant used for the simulation of the microkinetic model. This should be same as the name             specified in the species list.
    '''
    global InitialReactant
    InitialReactant = reac_string

def input_reaction_cutoffrate(cutoffrate):
    '''
    This function is used to refine the edges in the reaction network visualization and helps in removing edges that have a net reaction rate lower     than the specified value.
    '''
    global Reaction_Rate_Cutoff
    Reaction_Rate_Cutoff = cutoffrate

def input_elements_desired(list):
    '''
    This function is used to refine the nodes in the reaction network visualization and helps in removing nodes that do not contain any of the          specified elements in the list.
    '''
    global Elemental_List
    Elemental_List = list

def input_normalization(norm):
    '''
    This function specifies the desired normalization desired by the user.
    
        1 - Normalization using Maximum reaction rate in the network
        2 - Normalization using Net rate of the inlet reactant
        3 - Visualization using local consumption of species
    '''
    global Normalization_requested
    Normalization_requested = norm
    
def input_output_directory(fname):
    '''
    This specifies the output directory desired where all the visualizations generated are stored.
    '''
    global output_directory_name
    output_directory_name = fname

def erase_data():
    '''
    This function is used to reset all the datastructures in the run if using multiple case studies at a time.
    '''
    global species_filename, reactions_filename, output_directory_name
    global openmkm_species_list, openmkm_species_phase, openmkm_species_surf_cov, openmkm_species_composition
    global openmkm_reaction_strings, openmkm_reaction_reactants, openmkm_reaction_products, openmkm_reaction_fwdrate
    global openmkm_reaction_revrate, openmkm_reaction_netrate, openmkm_reaction_pei
    global Elemental_List, Normalization_requested, openmkm_species_consider, Elements_Available
    global stoich_matrix, rpa_local
    global is_surface_cov_def, header_count
    species_filename = ''
    reactions_filename = ''
    output_directory_name = ''
    openmkm_species_list = []
    openmkm_species_phase = []
    openmkm_species_surf_cov = []
    openmkm_species_composition = []
    openmkm_reaction_strings = []
    openmkm_reaction_reactants = []
    openmkm_reaction_products = []
    openmkm_reaction_fwdrate = []
    openmkm_reaction_revrate = []
    openmkm_reaction_netrate = []
    openmkm_reaction_pei = []
    Elemental_List = []
    Normalization_requested = 0
    openmkm_species_consider = []
    stoich_matrix = np.array([])
    rpa_local = np.array([])
    Elements_Available = []
    is_surface_cov_def = False
    header_count = 0
    
def get_color(surf_cov, phase):
    '''
    This function is used to get the color for node coloring. 
    
        surf_cov : Surface coverage of a species lying between 0 - 1
        phase : Species phase either Gas or Surface
    '''
    global is_surface_cov_def
    if is_surface_cov_def:
        if phase == 'Gas':
            color = 'none'
        else:
            if float(surf_cov) >= 0.85:
                color = 'coral'
            elif 0.7 <= float(surf_cov) < 0.85:
                color = 'orangered'
            elif 0.55 <= float(surf_cov) < 0.7:
                color = 'red'
            elif 0.4 <= float(surf_cov) < 0.55:
                color = 'darkorange'
            elif 0.2 <= float(surf_cov) < 0.4:
                color = 'orange'
            elif 0.1 <= float(surf_cov) < 0.2:
                color = 'gold'
            elif 0.01 <= float(surf_cov) < 0.1:
                color = 'lawngreen'
            elif 1E-5 <= float(surf_cov) < 0.01:
                color = 'springgreen'
            elif 1E-10 <= float(surf_cov) < 1E-5:
                color = 'dodgerblue'
            else:
                color = 'purple'
    else:
        if phase == 'Gas':
            color = 'cornsilk'
        else:
            color = 'azure'
    
    return color

def normalization_type(norm):
    '''
    This function returns the filename for network visualization based on the 
    normalization specified by the user
    
        norm : Integer value either 1, 2, or 3
    '''
    if norm == 1:
        return 'Normalization_MaxReactionRate.txt'
    elif norm == 2:
        return 'Normalization_NetReactionRate.txt'
    elif norm == 3:
        return 'Normalization_LocalConsumption.txt'
    else:
        print('No such normalization exists')
        
def print_header(fname, Rank_Sep, Node_Sep):
    '''
    This function generates the code for the header of the network visualization file
    
        fname : Visualization filename
        Rank_Sep : Vertical distance between nodes (if rankdir = TB)
                    Horizontal distance between nodes (if rankdir = LR)
        Node_Sep : Horizontal distance between nodes (if rankdir = TB)
                    Vertical distance between nodes (if rankdir = LR)          
    '''
    f = open(fname, "w+")
    f.write('digraph G {\n')
    f.write('splines = true;\n')
    f.write('graph [bgcolor=lightgray, resolution=64, fontname=Arial, fontcolor=blue, fontsize=36];\n')
    f.write('node [fontsize=12];\n')
    f.write('edge [fontsize=30];\n')
    f.write('label = "Reaction Path Analysis";\n')
    f.write('labelloc = "t";\n')
    f.write('center=1;\n')
    f.write('size="10,10";\n')
    f.write('ranksep="' + str(Rank_Sep) + ' equally";\n')
    f.write('nodesep="' + str(Node_Sep) + ' equally";\n')
    f.write('rankdir=TB;\n')
    f.write('bgcolor=white;\n')
    f.close()
    
def print_header_species(fname, Rank_Sep, Node_Sep):
    '''
    This function generates the code for the header of the species visualization file
    
        fname : Visualization filename
        Rank_Sep : Vertical distance between nodes (if rankdir = TB)
                    Horizontal distance between nodes (if rankdir = LR)
        Node_Sep : Horizontal distance between nodes (if rankdir = TB)
                    Vertical distance between nodes (if rankdir = LR)          
    '''
    f = open(fname, "w+")
    f.write('digraph G {\n')
    f.write('splines = true;\n')
    f.write('graph [bgcolor=lightgray, resolution=64, fontname=Arial, fontcolor=blue, fontsize=36];\n')
    f.write('node [fontsize=12];\n')
    f.write('edge [fontsize=30];\n')
    f.write('label = "Reaction Path Analysis";\n')
    f.write('labelloc = "t";\n')
    f.write('center=1;\n')
    f.write('size="10,10";\n')
    f.write('ranksep="' + str(Rank_Sep) + ' equally";\n')
    f.write('nodesep="' + str(Node_Sep) + ' equally";\n')
    f.write('rankdir=LR;\n')
    f.write('bgcolor=white;\n')
    f.close()
    
def normalization(fname,normalization_rate,normalization_type,stm,rpa_l):
    '''
    This function generates the normalizations used for specifying linewidth, arrowsize, 
    and edge thickness in the visualization
    
        fname : Visualization filename
        normalization_rate : Reaction rate used for normalization based on user's input
        normalization_type : Integer value either 1, 2, or 3
        stm : matrix having element stm(i,j) as stoichiometric 
                coefficients of species i in reaction j
        rpa_l : matrix having element rpa_l(i,j) as % consumption/production 
                of species i in reaction j
                - element value is negative if species i is net consumed in the reaction j
                - element value is positive if species i is net produced in the reaction j
    '''
    global openmkm_species_consider
    stoich_matrix = stm
    rpa_local = rpa_l
    f = open(fname, "a+")
    for i in range(len(openmkm_species_list)):
        if openmkm_species_consider[i] == True:
            count = 0 #Keeping a count of how many reactions with reaction rate more than cutoff are present
            for j in range(len(openmkm_reaction_strings)):
                #First checking if a reaction exists above the cutoff rate
                if float(stoich_matrix[i][j] != 0.0):
                    if abs(float(openmkm_reaction_netrate[j])) >= float(Reaction_Rate_Cutoff):
                        count = count + 1
                #So if the species is being consumed in the reaction
                if float(rpa_local[i][j]) <= 0.0 and float(stoich_matrix[i][j]) != 0.0 and float(stoich_matrix[i][j]) < 0.0:
                    for prod in openmkm_reaction_products[j]:
                        #check if the product needs to be considered for visualization
                        if openmkm_species_consider[int(openmkm_species_list.index(prod))] == True:
                            if abs(float(openmkm_reaction_netrate[j])) >= float(Reaction_Rate_Cutoff):
                                if float(openmkm_reaction_pei[j]) >= Equilibrium_Lower and float(openmkm_reaction_pei[j]) <= Equilibrium_Upper:
                                    if normalization_type == 1 or normalization_type == 2:
                                        linewidth = abs(int(((float(openmkm_reaction_netrate[j]))*100)/(normalization_rate*20))) + int(1)
                                        arrowsize = abs(int(((float(openmkm_reaction_netrate[j]))*100)/(normalization_rate*120))) + int(1)
                                        edgelabel = abs(int(((float(openmkm_reaction_netrate[j]))*100)/(normalization_rate)))
                                    elif normalization_type == 3:
                                        linewidth = abs(int(rpa_local[i][j]/20)) + int(1)
                                        arrowsize = abs(int(rpa_local[i][j]/60)) + int(1)
                                        edgelabel = abs(int(rpa_local[i][j]))
                                    else:
                                        linewidth = 0
                                        arrowsize = 0
                                        edgelabel = 0
                                    if float(rpa_local[i][j]) < 0.0:
                                        f.write('edge[dir="forward",style="setlinewidth(' + str(linewidth) + ')",color=green,weight=2,arrowsize=' + str(arrowsize) + ',label="   ' + str(edgelabel) + '%"];\n')
                                    else:
                                        f.write('edge[dir="none",style="setlinewidth(' + str(linewidth) + ')",color=green,weight=2,arrowsize=' + str(arrowsize) + ',label="   ' + str(edgelabel) + '%"];\n')
                                    f.write('"' + str(openmkm_species_list[i]) + '"->"' + str(prod) + '"\n')
                                else:
                                    if normalization_type == 1 or normalization_type == 2:
                                        linewidth = abs(int(((float(openmkm_reaction_netrate[j]))*100)/(normalization_rate*20))) + int(1)
                                        arrowsize = abs(int(((float(openmkm_reaction_netrate[j]))*100)/(normalization_rate*120))) + int(1)
                                        edgelabel = abs(int(((float(openmkm_reaction_netrate[j]))*100)/(normalization_rate)))
                                    elif normalization_type == 3:
                                        linewidth = abs(int(rpa_local[i][j]/20)) + int(1)
                                        arrowsize = abs(int(rpa_local[i][j]/60)) + int(1)
                                        edgelabel = abs(int(rpa_local[i][j]))
                                    else:
                                        linewidth = 0
                                        arrowsize = 0
                                        edgelabel = 0
                                    if float(rpa_local[i][j]) < 0.0:
                                        f.write('edge[dir="forward",style="setlinewidth(' + str(linewidth) + ')",color=black,weight=2,arrowsize=' + str(arrowsize) + ',label="   ' + str(edgelabel) + '%"];\n')
                                    else:
                                        f.write('edge[dir="none",style="setlinewidth(' + str(linewidth) + ')",color=black,weight=2,arrowsize=' + str(arrowsize) + ',label="   ' + str(edgelabel) + '%"];\n')
                                    f.write('"' + str(openmkm_species_list[i]) + '"->"' + str(prod) + '"\n')
                if float(rpa_local[i][j]) >= 0.0 and float(stoich_matrix[i][j]) != 0.0 and float(stoich_matrix[i][j]) < 0.0:
                    for prod in openmkm_reaction_products[j]:
                        if openmkm_species_consider[int(openmkm_species_list.index(prod))] == True:
                            if abs(float(openmkm_reaction_netrate[j])) >= float(Reaction_Rate_Cutoff):
                                if float(openmkm_reaction_pei[j]) >= Equilibrium_Lower and float(openmkm_reaction_pei[j]) <= Equilibrium_Upper:
                                    if normalization_type == 1 or normalization_type == 2:
                                        linewidth = abs(int(((float(openmkm_reaction_netrate[j]))*100)/(normalization_rate*20))) + int(1)
                                        arrowsize = abs(int(((float(openmkm_reaction_netrate[j]))*100)/(normalization_rate*120))) + int(1)
                                        edgelabel = abs(int(((float(openmkm_reaction_netrate[j]))*100)/(normalization_rate)))
                                    elif normalization_type == 3:
                                        linewidth = abs(int(rpa_local[i][j]/20)) + int(1)
                                        arrowsize = abs(int(rpa_local[i][j]/60)) + int(1)
                                        edgelabel = abs(int(rpa_local[i][j]))
                                    else:
                                        linewidth = 0
                                        arrowsize = 0
                                        edgelabel = 0
                                    if float(rpa_local[i][j]) > 0.0:
                                        f.write('edge[dir="forward",style="setlinewidth(' + str(linewidth) + ')",color=green,weight=2,arrowsize=' + str(arrowsize) + ',label="   ' + str(edgelabel) + '%"];\n')
                                    else:
                                        f.write('edge[dir="none",style="setlinewidth(' + str(linewidth) + ')",color=green,weight=2,arrowsize=' + str(arrowsize) + ',label="   ' + str(edgelabel) + '%"];\n')
                                    f.write('"' + str(prod) + '"->"' + str(openmkm_species_list[i]) + '"\n')
                                else:
                                    if normalization_type == 1 or normalization_type == 2:
                                        linewidth = abs(int(((float(openmkm_reaction_netrate[j]))*100)/(normalization_rate*20))) + int(1)
                                        arrowsize = abs(int(((float(openmkm_reaction_netrate[j]))*100)/(normalization_rate*120))) + int(1)
                                        edgelabel = abs(int(((float(openmkm_reaction_netrate[j]))*100)/(normalization_rate)))
                                    elif normalization_type == 3:
                                        linewidth = abs(int(rpa_local[i][j]/20)) + int(1)
                                        arrowsize = abs(int(rpa_local[i][j]/60)) + int(1)
                                        edgelabel = abs(int(rpa_local[i][j]))
                                    else:
                                        linewidth = 0
                                        arrowsize = 0
                                        edgelabel = 0
                                    if float(rpa_local[i][j]) > 0.0:
                                        f.write('edge[dir="forward",style="setlinewidth(' + str(linewidth) + ')",color=black,weight=2,arrowsize=' + str(arrowsize) + ',label="   ' + str(edgelabel) + '%"];\n')
                                    else:
                                        f.write('edge[dir="none",style="setlinewidth(' + str(linewidth) + ')",color=black,weight=2,arrowsize=' + str(arrowsize) + ',label="   ' + str(edgelabel) + '%"];\n')
                                    f.write('"' + str(prod) + '"->"' + str(openmkm_species_list[i]) + '"\n')
                    
            if count > 0:
                color1 = get_color(openmkm_species_surf_cov[i],openmkm_species_phase[i])
                f.write('"' + str(openmkm_species_list[i]) + '"[shape=rectangle,style=filled,fontsize=35,width=0,height=0,fillcolor=' + str(color1) + ',URL="' + str("Species\\") + str(i) + '.svg",shape=plaintext];\n')


def generate_visualizations():
    '''
    This executes the primary function generating visualizations for the reaction network.
    '''
    #Implementing pandas to read in the header and initialize the element fields
    #delimiter = '\t'
    #pandas puts in the delimited \t between spaces and 
    #so it needs to be removed when reading in the file
    #Populating all elements provided in input species file
    global rpa_local
    global is_surface_cov_def, header_count
    input_species = pd.read_csv(species_filename, header = 0, delimiter = '\t')
    for col in input_species.columns:
        if col == 'Species_name':
            header_count = header_count + 1
        if col == 'Phase':
            header_count = header_count + 1
        if col == 'Surf_cov':
            is_surface_cov_def = True
            header_count = header_count + 1
        if col != 'Species_name' and col != 'Phase' and col != 'Surf_cov':
            Elements_Available.append(col) 
    #Populating species datastructures
    with open(species_filename, 'r') as f_ptr:
        next(f_ptr)
        for line in f_ptr:
            line = line.replace('\n', '')
            words = line.split()
            openmkm_species_list.append(words[0])
            openmkm_species_phase.append(words[1])
            openmkm_species_surf_cov.append(words[2])
            composition = {}
            for j in range(int(header_count),len(words)):
                if int(words[j]) > 0:
                    composition[Elements_Available[j-int(header_count)]] = words[j]
            openmkm_species_composition.append(composition)
    
    TotalSpecies = len(openmkm_species_list)        
    #Populating reactions datastructures 
    with open(reactions_filename, 'r') as f_ptr:
        for i in range(1):
            f_ptr.readline()
        for line in f_ptr:
            line = line.replace('\n', '')
            words = line.split()
            openmkm_reaction_fwdrate.append(words[0]) # fwd rate
            openmkm_reaction_revrate.append(words[1]) # rev rate
            openmkm_reaction_netrate.append(words[2]) # net rate
            openmkm_reaction_pei.append(words[3]) # pei
            
            reaction_string = []
            reaction_str = ''
            reactants = []
            products = []
            stoich_list = np.zeros(TotalSpecies)
            for j in range(4,len(words)):
                reaction_string.append(words[j])
                reaction_str = reaction_str + str(' ') + str(words[j])
            openmkm_reaction_strings.append(reaction_str)
            arrow_index = reaction_string.index('<=>')
            for j in range(0,arrow_index):
                reactants.append(reaction_string[j])
            for k in range(arrow_index+1, len(reaction_string)):
                products.append(reaction_string[k])
            reactant_without_coeff = []
            react1 = []
            temp_r = []
            for cell_r in reactants:
                if cell_r == '+':
                    react1.append(temp_r)
                    if len(temp_r) > 1:
                        reac_stoich = temp_r[0]
                        reac_string = temp_r[1]
                    else:
                        reac_stoich = 1
                        reac_string = temp_r[0]
                    species_index = openmkm_species_list.index(reac_string)
                    stoich_list[species_index] = -int(reac_stoich)
                    reactant_without_coeff.append(reac_string)
                    temp_r = []
                else:
                    temp_r.append(cell_r)
                if cell_r == reactants[-1]:
                    react1.append(temp_r)
                    if len(temp_r) > 1:
                        reac_stoich = temp_r[0]
                        reac_string = temp_r[1]
                    else:
                        reac_stoich = 1
                        reac_string = temp_r[0]
                    species_index = openmkm_species_list.index(reac_string)
                    stoich_list[species_index] = -int(reac_stoich)
                    reactant_without_coeff.append(reac_string)
            openmkm_reaction_reactants.append(reactant_without_coeff)
            
            product_without_coeff = []
            prod1 = []
            temp_p = []
            for cell_p in products:
                if cell_p == '+':
                    prod1.append(temp_p)
                    if len(temp_p) > 1:
                        prod_stoich = temp_p[0]
                        prod_string = temp_p[1]
                    else:
                        prod_stoich = 1
                        prod_string = temp_p[0]
                    species_index = openmkm_species_list.index(prod_string)
                    stoich_list[species_index] = int(prod_stoich)
                    product_without_coeff.append(prod_string)
                    temp_p = []
                else:
                    temp_p.append(cell_p)
                if cell_p == products[-1]:
                    prod1.append(temp_p)
                    if len(temp_p) > 1:
                        prod_stoich = temp_p[0]
                        prod_string = temp_p[1]
                    else:
                        prod_stoich = 1
                        prod_string = temp_p[0]
                    species_index = openmkm_species_list.index(prod_string)
                    stoich_list[species_index] = int(prod_stoich)
                    product_without_coeff.append(prod_string)
            openmkm_reaction_products.append(product_without_coeff)
    
            if len(openmkm_reaction_netrate) == 1:
                stoich_matrix = stoich_list
            else:
                stoich_matrix = np.vstack( (stoich_matrix, stoich_list) )
            
    stoich_matrix = np.transpose(stoich_matrix)        
    
    #General output specifying number of reactions and species in system
    TotalReactions = len(openmkm_reaction_netrate)
    
    for i in range(len(openmkm_species_list)):
        consider = False
        species_dict = openmkm_species_composition[i]
        #Checking if elements in user specified list is present in species
        for elem in Elemental_List:
            for elem1 in species_dict.keys():
                if elem == elem1:
                    consider = True
        openmkm_species_consider.append(consider)
    
    check_output = os.path.exists(output_directory_name)        
    if check_output == False:
        os.makedirs(output_directory_name)
    
    species_output_dir = str(output_directory_name) + str('Species/')
    #Check if species folder is present
    check_species = os.path.exists(species_output_dir)
    if check_species == False:
        os.mkdir(species_output_dir)
    
    filename = normalization_type(Normalization_requested)
    pre, ext = os.path.splitext(filename)
    filename = str(output_directory_name) + str(filename)
    filename_svg = str(output_directory_name) + str(pre) + str('.svg')
    print_header(filename, Rank_Sep, Node_Sep)
    f = open(filename, "a+")
    
    MaxReactionRate = max(abs(float(n)) for n in openmkm_reaction_netrate)
    reaction_deln = np.sum(stoich_matrix, axis=0)
    initial_reactant_index = openmkm_species_list.index(InitialReactant)
    
    NetRate = 0.0
    for i in range(len(stoich_matrix[initial_reactant_index])):
        NetRate = NetRate + float(stoich_matrix[initial_reactant_index][i])*float(openmkm_reaction_netrate[i])

    rpa_local = np.zeros( (TotalSpecies, TotalReactions) )
    rpa_consumption = []
    rpa_production = []
    
    #Calculate the net rate of formation and consumption for each species    
    for i in range(len(openmkm_species_list)):
        NetConsumptionRate = 0.0
        NetProductionRate = 0.0
        for j in range(len(openmkm_reaction_strings)):
            if float(stoich_matrix[i][j]) < 0 and float(openmkm_reaction_netrate[j]) > 0: # => consumption
                NetConsumptionRate = NetConsumptionRate + float(openmkm_reaction_netrate[j])*float(stoich_matrix[i][j])
            elif float(stoich_matrix[i][j]) < 0 and float(openmkm_reaction_netrate[j]) < 0: # => production
                NetProductionRate = NetProductionRate + float(openmkm_reaction_netrate[j])*float(stoich_matrix[i][j])
            elif float(stoich_matrix[i][j]) > 0 and float(openmkm_reaction_netrate[j]) > 0: # => production
                NetProductionRate = NetProductionRate + float(openmkm_reaction_netrate[j])*float(stoich_matrix[i][j])
            elif float(stoich_matrix[i][j]) > 0 and float(openmkm_reaction_netrate[j]) < 0: # => consumption
                NetConsumptionRate = NetConsumptionRate + float(openmkm_reaction_netrate[j])*float(stoich_matrix[i][j])
            else:
                #print('Either species does not participate in reaction or reaction rate is zero')
                pass
        rpa_consumption.append(NetConsumptionRate)
        rpa_production.append(NetProductionRate)
        
    #Next generate rpa_local for each species 
    for i in range(len(openmkm_species_list)):
        for j in range(len(openmkm_reaction_strings)):
            if float(stoich_matrix[i][j]) < 0 and float(openmkm_reaction_netrate[j]) > 0: # => consumption
                rpa_local[i][j] = float(openmkm_reaction_netrate[j])*float(stoich_matrix[i][j])*float(100)/float(abs(rpa_consumption[i]))
            elif float(stoich_matrix[i][j]) < 0 and float(openmkm_reaction_netrate[j]) < 0: # => production
                rpa_local[i][j] = float(openmkm_reaction_netrate[j])*float(stoich_matrix[i][j])*float(100)/float(abs(rpa_production[i]))
            elif float(stoich_matrix[i][j]) > 0 and float(openmkm_reaction_netrate[j]) > 0: # => production
                rpa_local[i][j] = float(openmkm_reaction_netrate[j])*float(stoich_matrix[i][j])*float(100)/float(abs(rpa_production[i]))
            elif float(stoich_matrix[i][j]) > 0 and float(openmkm_reaction_netrate[j]) < 0: # => consumption
                rpa_local[i][j] = float(openmkm_reaction_netrate[j])*float(stoich_matrix[i][j])*float(100)/float(abs(rpa_consumption[i]))
            else:
                #print('Either species does not participate in reaction or reaction rate is zero')
                pass
    if Normalization_requested == 1:
        normalization(filename, MaxReactionRate, Normalization_requested,stoich_matrix,rpa_local)
    elif Normalization_requested == 2:
        normalization(filename, abs(NetRate), Normalization_requested,stoich_matrix,rpa_local)
    elif Normalization_requested == 3:
        normalization(filename, abs(NetRate), Normalization_requested,stoich_matrix,rpa_local)
    else:
        f.write('This type of normalization does not exist\n')
    
    f.write('}')
    f.close()

    #Generate network visualization    
    check_call(['dot','-Tsvg',f.name,'-o',filename_svg])
    
    #Generate species visualizations in output/Species folder
    for i in range(len(openmkm_species_list)):
        speciescounter = int(i)
        fname = str(species_output_dir) + str(speciescounter) + str('.txt')
        fname_svg = str(species_output_dir) + str(speciescounter) + str('.svg')
        print_header_species(fname, Rank_Sep, Node_Sep)
        f = open(fname, "a+")
        color = get_color(openmkm_species_surf_cov[i],openmkm_species_phase[i])
        f.write('"' + str(openmkm_species_list[i]) + '"[shape=rectangle,style=filled,fontsize=35,width=0,height=0,fillcolor=' + str(color) + ',URL="' + str(i) + '.svg",shape=plaintext];\n')
        for j in range(len(openmkm_reaction_strings)):
            if float(rpa_local[i][j]) <= 0.0 and float(stoich_matrix[i][j]) != 0.0 and float(stoich_matrix[i][j]) < 0.0:
                for prod in openmkm_reaction_products[j]:
                    if openmkm_species_consider[int(openmkm_species_list.index(prod))] == True:
                        if abs(float(openmkm_reaction_netrate[j])) >= float(Reaction_Rate_Cutoff):
                            color1 = get_color(openmkm_species_surf_cov[int(openmkm_species_list.index(prod))],openmkm_species_phase[int(openmkm_species_list.index(prod))])
                            f.write('"' + str(openmkm_species_list[int(openmkm_species_list.index(prod))]) + '"[shape=rectangle,style=filled,fontsize=35,width=0,height=0,fillcolor=' + str(color1) + ',URL="' + str(int(openmkm_species_list.index(prod))) + '.svg",shape=plaintext];\n')
                            if float(openmkm_reaction_pei[j]) >= Equilibrium_Lower and float(openmkm_reaction_pei[j]) <= Equilibrium_Upper:
                                linewidth = abs(int(rpa_local[i][j]/20)) + int(1)
                                arrowsize = abs(int(rpa_local[i][j]/60)) + int(1)
                                edgelabel = abs(int(rpa_local[i][j]))
                                if float(rpa_local[i][j]) < 0.0:
                                    f.write('edge[dir="forward",style="setlinewidth(' + str(linewidth) + ')",color=green,weight=2,arrowsize=' + str(arrowsize) + ',label="   ' + str(j) + '   ' + str(edgelabel) + '%   ' + str(abs(float(openmkm_reaction_netrate[j]))) + ' mol/s    ' + str(openmkm_reaction_pei[j]) + '"];\n')
                                else:
                                    f.write('edge[dir="none",style="setlinewidth(' + str(linewidth) + ')",color=green,weight=2,arrowsize=' + str(arrowsize) + ',label="   ' + str(j) + '   ' + str(edgelabel) + '%   ' + str(abs(float(openmkm_reaction_netrate[j]))) + ' mol/s    ' + str(openmkm_reaction_pei[j]) + '"];\n')
                                f.write('"' + str(openmkm_species_list[i]) + '"->"' + str(prod) + '"\n')
                            else:
                                linewidth = abs(int(rpa_local[i][j]/20)) + int(1)
                                arrowsize = abs(int(rpa_local[i][j]/60)) + int(1)
                                edgelabel = abs(int(rpa_local[i][j]))
                                if float(rpa_local[i][j]) < 0.0:
                                    f.write('edge[dir="forward",style="setlinewidth(' + str(linewidth) + ')",color=black,weight=2,arrowsize=' + str(arrowsize) + ',label="   ' + str(j) + '   ' + str(edgelabel) + '%   ' + str(abs(float(openmkm_reaction_netrate[j]))) + ' mol/s    ' + str(openmkm_reaction_pei[j]) + '"];\n')
                                else:
                                    f.write('edge[dir="none",style="setlinewidth(' + str(linewidth) + ')",color=black,weight=2,arrowsize=' + str(arrowsize) + ',label="   ' + str(j) + '   ' + str(edgelabel) + '%   ' + str(abs(float(openmkm_reaction_netrate[j]))) + ' mol/s    ' + str(openmkm_reaction_pei[j]) + '"];\n')
                                f.write('"' + str(openmkm_species_list[i]) + '"->"' + str(prod) + '"\n')
            if float(rpa_local[i][j]) <= 0.0 and float(stoich_matrix[i][j]) != 0.0 and float(stoich_matrix[i][j]) > 0.0:
                for reac in openmkm_reaction_reactants[j]:
                    if openmkm_species_consider[int(openmkm_species_list.index(reac))] == True:
                        if abs(float(openmkm_reaction_netrate[j])) >= float(Reaction_Rate_Cutoff):
                            color1 = get_color(openmkm_species_surf_cov[int(openmkm_species_list.index(reac))],openmkm_species_phase[int(openmkm_species_list.index(reac))])
                            f.write('"' + str(openmkm_species_list[int(openmkm_species_list.index(reac))]) + '"[shape=rectangle,style=filled,fontsize=35,width=0,height=0,fillcolor=' + str(color1) + ',URL="' + str(int(openmkm_species_list.index(reac))) + '.svg",shape=plaintext];\n')
                            if float(openmkm_reaction_pei[j]) >= Equilibrium_Lower and float(openmkm_reaction_pei[j]) <= Equilibrium_Upper:
                                linewidth = abs(int(rpa_local[i][j]/20)) + int(1)
                                arrowsize = abs(int(rpa_local[i][j]/60)) + int(1)
                                edgelabel = abs(int(rpa_local[i][j]))
                                if float(rpa_local[i][j]) < 0.0:
                                    f.write('edge[dir="forward",style="setlinewidth(' + str(linewidth) + ')",color=green,weight=2,arrowsize=' + str(arrowsize) + ',label="   ' + str(j) + '   ' + str(edgelabel) + '%   ' + str(abs(float(openmkm_reaction_netrate[j]))) + ' mol/s    ' + str(openmkm_reaction_pei[j]) + '"];\n')
                                else:
                                    f.write('edge[dir="none",style="setlinewidth(' + str(linewidth) + ')",color=green,weight=2,arrowsize=' + str(arrowsize) + ',label="   ' + str(j) + '   ' + str(edgelabel) + '%   ' + str(abs(float(openmkm_reaction_netrate[j]))) + ' mol/s    ' + str(openmkm_reaction_pei[j]) + '"];\n')
                                f.write('"' + str(openmkm_species_list[i]) + '"->"' + str(reac) + '"\n')
                            else:
                                linewidth = abs(int(rpa_local[i][j]/20)) + int(1)
                                arrowsize = abs(int(rpa_local[i][j]/60)) + int(1)
                                edgelabel = abs(int(rpa_local[i][j]))
                                if float(rpa_local[i][j]) < 0.0:
                                    f.write('edge[dir="forward",style="setlinewidth(' + str(linewidth) + ')",color=black,weight=2,arrowsize=' + str(arrowsize) + ',label="   ' + str(j) + '   ' + str(edgelabel) + '%   ' + str(abs(float(openmkm_reaction_netrate[j]))) + ' mol/s    ' + str(openmkm_reaction_pei[j]) + '"];\n')
                                else:
                                    f.write('edge[dir="none",style="setlinewidth(' + str(linewidth) + ')",color=black,weight=2,arrowsize=' + str(arrowsize) + ',label="   ' + str(j) + '   ' + str(edgelabel) + '%   ' + str(abs(float(openmkm_reaction_netrate[j]))) + ' mol/s    ' + str(openmkm_reaction_pei[j]) + '"];\n')
                                f.write('"' + str(openmkm_species_list[i]) + '"->"' + str(reac) + '"\n')
            if float(rpa_local[i][j]) >= 0.0 and float(stoich_matrix[i][j]) != 0.0 and float(stoich_matrix[i][j]) < 0.0:
                for prod in openmkm_reaction_products[j]:
                    if openmkm_species_consider[int(openmkm_species_list.index(prod))] == True:
                        if abs(float(openmkm_reaction_netrate[j])) >= float(Reaction_Rate_Cutoff):
                            color1 = get_color(openmkm_species_surf_cov[int(openmkm_species_list.index(prod))],openmkm_species_phase[int(openmkm_species_list.index(prod))])
                            f.write('"' + str(openmkm_species_list[int(openmkm_species_list.index(prod))]) + '"[shape=rectangle,style=filled,fontsize=35,width=0,height=0,fillcolor=' + str(color1) + ',URL="' + str(int(openmkm_species_list.index(prod))) + '.svg",shape=plaintext];\n')
                            if float(openmkm_reaction_pei[j]) >= Equilibrium_Lower and float(openmkm_reaction_pei[j]) <= Equilibrium_Upper:
                                linewidth = abs(int(rpa_local[i][j]/20)) + int(1)
                                arrowsize = abs(int(rpa_local[i][j]/60)) + int(1)
                                edgelabel = abs(int(rpa_local[i][j]))
                                if float(rpa_local[i][j]) > 0.0:
                                    f.write('edge[dir="forward",style="setlinewidth(' + str(linewidth) + ')",color=green,weight=2,arrowsize=' + str(arrowsize) + ',label="   ' + str(j) + '   ' + str(edgelabel) + '%   ' + str(abs(float(openmkm_reaction_netrate[j]))) + ' mol/s    ' + str(openmkm_reaction_pei[j]) + '"];\n')
                                else:
                                    f.write('edge[dir="none",style="setlinewidth(' + str(linewidth) + ')",color=green,weight=2,arrowsize=' + str(arrowsize) + ',label="   ' + str(j) + '   ' + str(edgelabel) + '%   ' + str(abs(float(openmkm_reaction_netrate[j]))) + ' mol/s    ' + str(openmkm_reaction_pei[j]) + '"];\n')
                                f.write('"' + str(prod) + '"->"' + str(openmkm_species_list[i]) + '"\n')
                            else:
                                linewidth = abs(int(rpa_local[i][j]/20)) + int(1)
                                arrowsize = abs(int(rpa_local[i][j]/60)) + int(1)
                                edgelabel = abs(int(rpa_local[i][j]))
                                if float(rpa_local[i][j]) > 0.0:
                                    f.write('edge[dir="forward",style="setlinewidth(' + str(linewidth) + ')",color=black,weight=2,arrowsize=' + str(arrowsize) + ',label="   ' + str(j) + '   ' + str(edgelabel) + '%   ' + str(abs(float(openmkm_reaction_netrate[j]))) + ' mol/s    ' + str(openmkm_reaction_pei[j]) + '"];\n')
                                else:
                                    f.write('edge[dir="none",style="setlinewidth(' + str(linewidth) + ')",color=black,weight=2,arrowsize=' + str(arrowsize) + ',label="   ' + str(j) + '   ' + str(edgelabel) + '%   ' + str(abs(float(openmkm_reaction_netrate[j]))) + ' mol/s    ' + str(openmkm_reaction_pei[j]) + '"];\n')
                                f.write('"' + str(prod) + '"->"' + str(openmkm_species_list[i]) + '"\n')
            if float(rpa_local[i][j]) >= 0.0 and float(stoich_matrix[i][j]) != 0.0 and float(stoich_matrix[i][j]) > 0.0:
                for reac in openmkm_reaction_reactants[j]:
                    if openmkm_species_consider[int(openmkm_species_list.index(reac))] == True:
                        if abs(float(openmkm_reaction_netrate[j])) >= float(Reaction_Rate_Cutoff):
                            color1 = get_color(openmkm_species_surf_cov[int(openmkm_species_list.index(reac))],openmkm_species_phase[int(openmkm_species_list.index(reac))])
                            f.write('"' + str(openmkm_species_list[int(openmkm_species_list.index(reac))]) + '"[shape=rectangle,style=filled,fontsize=35,width=0,height=0,fillcolor=' + str(color1) + ',URL="' + str(int(openmkm_species_list.index(reac))) + '.svg",shape=plaintext];\n')
                            if float(openmkm_reaction_pei[j]) >= Equilibrium_Lower and float(openmkm_reaction_pei[j]) <= Equilibrium_Upper:
                                linewidth = abs(int(rpa_local[i][j]/20)) + int(1)
                                arrowsize = abs(int(rpa_local[i][j]/60)) + int(1)
                                edgelabel = abs(int(rpa_local[i][j]))
                                if float(rpa_local[i][j]) > 0.0:
                                    f.write('edge[dir="forward",style="setlinewidth(' + str(linewidth) + ')",color=green,weight=2,arrowsize=' + str(arrowsize) + ',label="   ' + str(j) + '   ' + str(edgelabel) + '%   ' + str(abs(float(openmkm_reaction_netrate[j]))) + ' mol/s    ' + str(openmkm_reaction_pei[j]) + '"];\n')
                                else:
                                    f.write('edge[dir="none",style="setlinewidth(' + str(linewidth) + ')",color=green,weight=2,arrowsize=' + str(arrowsize) + ',label="   ' + str(j) + '   ' + str(edgelabel) + '%   ' + str(abs(float(openmkm_reaction_netrate[j]))) + ' mol/s    ' + str(openmkm_reaction_pei[j]) + '"];\n')
                                f.write('"' + str(reac) + '"->"' + str(openmkm_species_list[i]) + '"\n')
                            else:
                                linewidth = abs(int(rpa_local[i][j]/20)) + int(1)
                                arrowsize = abs(int(rpa_local[i][j]/60)) + int(1)
                                edgelabel = abs(int(rpa_local[i][j]))
                                if float(rpa_local[i][j]) > 0.0:
                                    f.write('edge[dir="forward",style="setlinewidth(' + str(linewidth) + ')",color=black,weight=2,arrowsize=' + str(arrowsize) + ',label="   ' + str(j) + '   ' + str(edgelabel) + '%   ' + str(abs(float(openmkm_reaction_netrate[j]))) + ' mol/s    ' + str(openmkm_reaction_pei[j]) + '"];\n')
                                else:
                                    f.write('edge[dir="none",style="setlinewidth(' + str(linewidth) + ')",color=black,weight=2,arrowsize=' + str(arrowsize) + ',label="   ' + str(j) + '   ' + str(edgelabel) + '%   ' + str(abs(float(openmkm_reaction_netrate[j]))) + ' mol/s    ' + str(openmkm_reaction_pei[j]) + '"];\n')
                                f.write('"' + str(reac) + '"->"' + str(openmkm_species_list[i]) + '"\n')
        
        f.write('}')
        f.close()
        check_call(['dot','-Tsvg',f.name,'-o',fname_svg])
    
    #Generate legend for the visualizations
    legend_output_file = str(output_directory_name) + str('legend_cov.out')
    legend_output_file_svg = str(output_directory_name) + str('legend.svg')
#    generate_legend(legend_output_file, is_surface_cov_def)
    generate_legend(legend_output_file, is_surface_cov_def)
    #generate_legend(legend_output_file, True)
    check_call(['dot','-Tsvg',legend_output_file,'-o',legend_output_file_svg])