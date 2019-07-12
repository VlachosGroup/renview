# -*- coding: utf-8 -*-
"""
Created on Fri May 10 15:13:35 2019

@author: ugupta
"""

import os
os.environ["PATH"] += os.pathsep + 'C:/Program Files (x86)/Graphviz2.38/bin/'
import pydot
import numpy as np
import pandas as pd
from graphviz import *
from subprocess import check_call

from Reactions.reactions_import import *
from Legend.legend import *

# Inputs from user include species_comp, reactions file
# along with cutoff rate, elements desired, rank sep and node sep if needed

'''
Demonstration script for reading in reaction rates and reactions to generate
visualizations
'''

def get_color(surf_cov, phase):
#    print(surf_cov, phase)
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
            color = 'green'
        elif 1E-10 <= float(surf_cov) < 1E-5:
            color = 'dodgerblue'
        else:
            color = 'purple'
    
    return color

def normalization_type(norm):
    if norm == 1:
#        print('Normalization requested is Max. Reaction Rate')
        return 'Normalization_MaxReactionRate.txt'
    elif norm == 2:
#        print('Normalization requested is Net Reaction Rate')
        return 'Normalization_NetReactionRate.txt'
    elif norm == 3:
#        print('Normalization requested is Local consumption')
        return 'Normalization_LocalConsumption.txt'
    else:
        print('No such normalization exists')
        
def print_header(fname, Rank_Sep, Node_Sep):
#    print('Inside print header')
#    print(fname)
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
#    print('Inside print header')
#    print(fname)
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
    
def normalization(fname,normalization_rate):
#    print('IN here',openmkm_species_list)
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
#                        print(openmkm_species_list[i],prod)
                        #check if the product needs to be considered for visualization
                        if openmkm_species_consider[int(openmkm_species_list.index(prod))] == True:
                            if abs(float(openmkm_reaction_netrate[j])) >= float(Reaction_Rate_Cutoff):
                                if float(openmkm_reaction_pei[j]) >= Equilibrium_Lower and float(openmkm_reaction_pei[j]) <= Equilibrium_Upper:
                                    linewidth = abs(int(((float(openmkm_reaction_netrate[j]))*100)/(normalization_rate*20))) + int(1)
                                    arrowsize = abs(int(((float(openmkm_reaction_netrate[j]))*100)/(normalization_rate*120))) + int(1)
                                    edgelabel = abs(int(((float(openmkm_reaction_netrate[j]))*100)/(normalization_rate)))
                                    if float(rpa_local[i][j]) < 0.0:
                                        f.write('edge[dir="forward",style="setlinewidth(' + str(linewidth) + ')",color=green,weight=2,arrowsize=' + str(arrowsize) + ',label="   ' + str(edgelabel) + '%"];\n')
                                    else:
                                        f.write('edge[dir="none",style="setlinewidth(' + str(linewidth) + ')",color=green,weight=2,arrowsize=' + str(arrowsize) + ',label="   ' + str(edgelabel) + '%"];\n')
                                    f.write('"' + str(openmkm_species_list[i]) + '"->"' + str(prod) + '"\n')
                                else:
                                    linewidth = abs(int(((float(openmkm_reaction_netrate[j]))*100)/(normalization_rate*20))) + int(1)
                                    arrowsize = abs(int(((float(openmkm_reaction_netrate[j]))*100)/(normalization_rate*120))) + int(1)
                                    edgelabel = abs(int(((float(openmkm_reaction_netrate[j]))*100)/(normalization_rate)))
                                    if float(rpa_local[i][j]) < 0.0:
                                        f.write('edge[dir="forward",style="setlinewidth(' + str(linewidth) + ')",color=black,weight=2,arrowsize=' + str(arrowsize) + ',label="   ' + str(edgelabel) + '%"];\n')
                                    else:
                                        f.write('edge[dir="none",style="setlinewidth(' + str(linewidth) + ')",color=black,weight=2,arrowsize=' + str(arrowsize) + ',label="   ' + str(edgelabel) + '%"];\n')
                                    f.write('"' + str(openmkm_species_list[i]) + '"->"' + str(prod) + '"\n')
                if float(rpa_local[i][j]) >= 0.0 and float(stoich_matrix[i][j]) != 0.0 and float(stoich_matrix[i][j]) < 0.0:
                    for prod in openmkm_reaction_products[j]:
#                        print(openmkm_species_list[i], prod, 'new')
                        if openmkm_species_consider[int(openmkm_species_list.index(prod))] == True:
                            if abs(float(openmkm_reaction_netrate[j])) >= float(Reaction_Rate_Cutoff):
                                if float(openmkm_reaction_pei[j]) >= Equilibrium_Lower and float(openmkm_reaction_pei[j]) <= Equilibrium_Upper:
                                    linewidth = abs(int(((float(openmkm_reaction_netrate[j]))*100)/(normalization_rate*20))) + int(1)
                                    arrowsize = abs(int(((float(openmkm_reaction_netrate[j]))*100)/(normalization_rate*120))) + int(1)
                                    edgelabel = abs(int(((float(openmkm_reaction_netrate[j]))*100)/(normalization_rate)))
                                    if float(rpa_local[i][j]) > 0.0:
                                        f.write('edge[dir="forward",style="setlinewidth(' + str(linewidth) + ')",color=green,weight=2,arrowsize=' + str(arrowsize) + ',label="   ' + str(edgelabel) + '%"];\n')
                                    else:
                                        f.write('edge[dir="none",style="setlinewidth(' + str(linewidth) + ')",color=green,weight=2,arrowsize=' + str(arrowsize) + ',label="   ' + str(edgelabel) + '%"];\n')
                                    f.write('"' + str(prod) + '"->"' + str(openmkm_species_list[i]) + '"\n')
                                else:
                                    linewidth = abs(int(((float(openmkm_reaction_netrate[j]))*100)/(normalization_rate*20))) + int(1)
                                    arrowsize = abs(int(((float(openmkm_reaction_netrate[j]))*100)/(normalization_rate*120))) + int(1)
                                    edgelabel = abs(int(((float(openmkm_reaction_netrate[j]))*100)/(normalization_rate)))
                                    if float(rpa_local[i][j]) > 0.0:
                                        f.write('edge[dir="forward",style="setlinewidth(' + str(linewidth) + ')",color=black,weight=2,arrowsize=' + str(arrowsize) + ',label="   ' + str(edgelabel) + '%"];\n')
                                    else:
                                        f.write('edge[dir="none",style="setlinewidth(' + str(linewidth) + ')",color=black,weight=2,arrowsize=' + str(arrowsize) + ',label="   ' + str(edgelabel) + '%"];\n')
                                    f.write('"' + str(prod) + '"->"' + str(openmkm_species_list[i]) + '"\n')
                    
            if count > 0:
                color1 = get_color(openmkm_species_surf_cov[i],openmkm_species_phase[i])
                f.write('"' + str(openmkm_species_list[i]) + '"[shape=rectangle,style=filled,fontsize=35,width=0,height=0,fillcolor=' + str(color1) + ',URL="' + str("Species\\") + str(i) + '.svg",shape=plaintext];\n')

def normalization_local(fname):
    print('now doing local consumption visualization')
#    print('IN here',openmkm_species_list)
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
#                        print(openmkm_species_list[i],prod)
                        #check if the product needs to be considered for visualization
                        if openmkm_species_consider[int(openmkm_species_list.index(prod))] == True:
                            if abs(float(openmkm_reaction_netrate[j])) >= float(Reaction_Rate_Cutoff):
                                if float(openmkm_reaction_pei[j]) >= Equilibrium_Lower and float(openmkm_reaction_pei[j]) <= Equilibrium_Upper:
                                    linewidth = abs(int(rpa_local[i][j]/20)) + int(1)
                                    arrowsize = abs(int(rpa_local[i][j]/60)) + int(1)
                                    edgelabel = abs(int(rpa_local[i][j]))
                                    if float(rpa_local[i][j]) < 0.0:
                                        f.write('edge[dir="forward",style="setlinewidth(' + str(linewidth) + ')",color=green,weight=2,arrowsize=' + str(arrowsize) + ',label="   ' + str(edgelabel) + '%"];\n')
                                    else:
                                        f.write('edge[dir="none",style="setlinewidth(' + str(linewidth) + ')",color=green,weight=2,arrowsize=' + str(arrowsize) + ',label="   ' + str(edgelabel) + '%"];\n')
                                    f.write('"' + str(openmkm_species_list[i]) + '"->"' + str(prod) + '"\n')
                                else:
                                    linewidth = abs(int(rpa_local[i][j]/20)) + int(1)
                                    arrowsize = abs(int(rpa_local[i][j]/60)) + int(1)
                                    edgelabel = abs(int(rpa_local[i][j]))
                                    if float(rpa_local[i][j]) < 0.0:
                                        f.write('edge[dir="forward",style="setlinewidth(' + str(linewidth) + ')",color=black,weight=2,arrowsize=' + str(arrowsize) + ',label="   ' + str(edgelabel) + '%"];\n')
                                    else:
                                        f.write('edge[dir="none",style="setlinewidth(' + str(linewidth) + ')",color=black,weight=2,arrowsize=' + str(arrowsize) + ',label="   ' + str(edgelabel) + '%"];\n')
                                    f.write('"' + str(openmkm_species_list[i]) + '"->"' + str(prod) + '"\n')
                if float(rpa_local[i][j]) >= 0.0 and float(stoich_matrix[i][j]) != 0.0 and float(stoich_matrix[i][j]) < 0.0:
                    for prod in openmkm_reaction_products[j]:
#                        print(openmkm_species_list[i], prod, 'new')
                        if openmkm_species_consider[int(openmkm_species_list.index(prod))] == True:
                            if abs(float(openmkm_reaction_netrate[j])) >= float(Reaction_Rate_Cutoff):
                                if float(openmkm_reaction_pei[j]) >= Equilibrium_Lower and float(openmkm_reaction_pei[j]) <= Equilibrium_Upper:
                                    linewidth = abs(int(rpa_local[i][j]/20)) + int(1)
                                    arrowsize = abs(int(rpa_local[i][j]/60)) + int(1)
                                    edgelabel = abs(int(rpa_local[i][j]))
                                    if float(rpa_local[i][j]) > 0.0:
                                        f.write('edge[dir="forward",style="setlinewidth(' + str(linewidth) + ')",color=green,weight=2,arrowsize=' + str(arrowsize) + ',label="   ' + str(edgelabel) + '%"];\n')
                                    else:
                                        f.write('edge[dir="none",style="setlinewidth(' + str(linewidth) + ')",color=green,weight=2,arrowsize=' + str(arrowsize) + ',label="   ' + str(edgelabel) + '%"];\n')
                                    f.write('"' + str(prod) + '"->"' + str(openmkm_species_list[i]) + '"\n')
                                else:
                                    linewidth = abs(int(rpa_local[i][j]/20)) + int(1)
                                    arrowsize = abs(int(rpa_local[i][j]/60)) + int(1)
                                    edgelabel = abs(int(rpa_local[i][j]))
                                    if float(rpa_local[i][j]) > 0.0:
                                        f.write('edge[dir="forward",style="setlinewidth(' + str(linewidth) + ')",color=black,weight=2,arrowsize=' + str(arrowsize) + ',label="   ' + str(edgelabel) + '%"];\n')
                                    else:
                                        f.write('edge[dir="none",style="setlinewidth(' + str(linewidth) + ')",color=black,weight=2,arrowsize=' + str(arrowsize) + ',label="   ' + str(edgelabel) + '%"];\n')
                                    f.write('"' + str(prod) + '"->"' + str(openmkm_species_list[i]) + '"\n')
                    
            if count > 0:
                #If surface coverages provided change fillcolor based on the legend
                color1 = get_color(openmkm_species_surf_cov[i],openmkm_species_phase[i])
                f.write('"' + str(openmkm_species_list[i]) + '"[shape=rectangle,style=filled,fontsize=35,width=0,height=0,fillcolor=' + str(color1) + ',URL="' + str("Species\\") + str(i) + '.svg",shape=plaintext];\n')


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


#Inputs from OpenMKM or any other chemical kinetics software
InitialReactant = 'NH3'
Equilibrium_Tolerance = 0.05
Reaction_Rate_Cutoff = 1.0E-09
Elemental_List = ['N', 'H']
Elements_Available = []
Rank_Sep = 0.25
Node_Sep = 0.25
Normalization_requested = 1

Equilibrium_Upper = 0.5 + Equilibrium_Tolerance
Equilibrium_Lower = 0.5 - Equilibrium_Tolerance


#Implementing pandas to read in the header and initialize the element fields
#delimiter = '\t'
#pandas puts in the delimited \t between spaces and 
#so it needs to be removed when reading in the file
#Populating all elements provided in input species file
input_species = pd.read_csv('..\\example\\species_comp.out', header = 0, delimiter = '\t')
for col in input_species.columns:
    if col != 'SpeciesName' and col != 'Phase':
        Elements_Available.append(col) 
#print(Elements_Available)
#Populating species datastructures
with open('..\\example\\species_comp.out', 'r') as f_ptr:
    next(f_ptr)
    for line in f_ptr:
        line = line.replace('\n', '')
        words = line.split()
        openmkm_species_list.append(words[0])
        openmkm_species_phase.append(words[1])
        openmkm_species_surf_cov.append(words[2])
        composition = {}
        for j in range(3,len(words)):
            if int(words[j]) > 0:
                composition[Elements_Available[j-3]] = words[j]
        openmkm_species_composition.append(composition)
        
#print('Species stored')
#Populating reactions datastructures 
#print(openmkm_species_surf_cov)
        
#reactions_read('rates_ss.out')
with open('..\\example\\rates_ss.out', 'r') as f_ptr:
    for i in range(2):
        f_ptr.readline()
    for line in f_ptr:
        line = line.replace('\n', '')
        words = line.split()
        openmkm_reaction_fwdrate.append(words[0])
        openmkm_reaction_revrate.append(words[1])
        openmkm_reaction_netrate.append(words[2])
        openmkm_reaction_pei.append(words[3])
        openmkm_reaction_strings.append(words[4])

openmkm_species_consider = []
for i in range(len(openmkm_species_list)):
#    print(openmkm_species_list[i])
    consider = False
    species_dict = openmkm_species_composition[i]
#    print(species_dict)
    #Checking if elements in user specified list is present in species
    for elem in Elemental_List:
        for elem1 in species_dict.keys():
            if elem == elem1:
                consider = True
    openmkm_species_consider.append(consider)

filename = normalization_type(Normalization_requested)
#print(filename)
pre, ext = os.path.splitext(filename)
#print(pre,ext)
filename_svg = str(pre) + str('.svg')
#print(filename_svg)
print_header(filename, Rank_Sep, Node_Sep)
f = open(filename, "a+")

#General output specifying number of reactions and species in system
TotalReactions = len(openmkm_reaction_netrate)
TotalSpecies = len(openmkm_species_list)
#print('The networks includes', TotalReactions, ' reactions and ',TotalSpecies, ' species!')

#print('Checking if initial reactant is present: ', openmkm_species_list.count(InitialReactant))
MaxReactionRate = max(abs(float(n)) for n in openmkm_reaction_netrate)
#print('Maximum reaction rate is: ',MaxReactionRate)

stoich_matrix = np.zeros( (TotalSpecies, TotalReactions) )

total_corelations = 0
openmkm_reactions_list = []
with open('..\\example\\reactions.out','r') as f_ptr:
    for line in f_ptr:
        line = line.replace('\n', '')
        reac_prod = line.split('<=>')
        reactants = reac_prod[0]
        products = reac_prod[1]
        reaction_index = openmkm_reaction_strings.index(line)

        reactant = reactants.split('+')
        reactant_without_coeff = []
        for reac in reactant:
            reac_split = reac.split('|')
            if len(reac_split) > 1:
                reac_stoich = reac_split[0]
                reac_string = reac_split[1]
            else:
                reac_stoich = 1
                reac_string = reac_split[0]
            species_index = openmkm_species_list.index(reac_string)
            stoich_matrix[species_index][reaction_index] = -int(reac_stoich)
            reactant_without_coeff.append(reac_string)
            total_corelations = total_corelations + 1
        openmkm_reaction_reactants.append(reactant_without_coeff)
        
        product = products.split('+')
        product_without_coeff = []
        for prod in product:
            prod_split = prod.split('|')
            if len(prod_split) > 1:
                prod_stoich = prod_split[0]
                prod_string = prod_split[1]
            else:
                prod_stoich = 1
                prod_string = prod_split[0]
            species_index = openmkm_species_list.index(prod_string)
            stoich_matrix[species_index][reaction_index] = int(prod_stoich)
            product_without_coeff.append(prod_string)
            total_corelations = total_corelations + 1
        openmkm_reaction_products.append(product_without_coeff)
            
#print('total number of rows in rpa_local matrix ',total_corelations)

reaction_deln = np.sum(stoich_matrix, axis=0)
#print('Listing all deln of reactions: ',reaction_deln)

initial_reactant_index = openmkm_species_list.index(InitialReactant)
#print('Index of initial reactant: ',initial_reactant_index)

NetRate = 0.0
#print(stoich_matrix[initial_reactant_index])
for i in range(len(stoich_matrix[initial_reactant_index])):
    NetRate = NetRate + float(stoich_matrix[initial_reactant_index][i])*float(openmkm_reaction_netrate[i])
    
#print('NetRate of starting species: ',abs(NetRate),' species: ',InitialReactant)

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
#            print('Either species does not participate in reaction or reaction rate is zero')
            pass
    rpa_consumption.append(NetConsumptionRate)
    rpa_production.append(NetProductionRate)
    
#print(rpa_consumption)
#print(rpa_production)

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
#            print('Either species does not participate in reaction or reaction rate is zero')
            pass
#print(rpa_local)    
if Normalization_requested == 1:
    normalization(filename, MaxReactionRate)
elif Normalization_requested == 2:
    normalization(filename, abs(NetRate))
elif Normalization_requested == 3:
    normalization_local(filename)
else:
    f.write('This type of normalization does not exist\n')




f.write('}')
f.close()

#Species Visualizations need to be added

check_call(['dot','-Tsvg',f.name,'-o',filename_svg])

generate_legend('legend_cov.out')
check_call(['dot','-Tsvg','legend_cov.out','-o','legend.svg'])

#check_call(['dot','-Tpng','RefinedVisualization.dot','-o','rpa_visualizationMaxRateNormalized1.png'])
#check_call(['dot','-Tsvg','RefinedVisualization.dot','-o','rpa_visualizationMaxRateNormalized1.svg'])
#dot -Ttiff rpa_visualizationMaxRateNormalized.out -o rpa_visualizationMaxRateNormalized1.svg

#(graph,) = pydot.graph_from_dot_file('rpa_visualizationMaxRateNormalized.dot')
#graph.write_png('rpa_visualizationMaxRateNormalized1.png')
        

#Check if species folder is present
check = os.path.exists("Species")
#print(os.path.exists("Species"))
if check == False:
    os.mkdir('Species')
#print(os.path.exists("Species"))

#Generate species visualizations
for i in range(len(openmkm_species_list)):
    speciescounter = int(i)
#    print(speciescounter)
    fname = str('Species/') + str(speciescounter) + str('.txt')
    fname_svg = str('Species/') + str(speciescounter) + str('.svg')
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
