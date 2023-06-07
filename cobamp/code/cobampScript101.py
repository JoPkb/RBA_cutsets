__author__ = 'Basile Sugranes', 'Martial Scavino', 'Ravy Leon Foun Lin', 'Jérémie Muller'
__version__ = '1.0.1'

#### GESTION DES ARGUMENTS AVEC ARGPARSER ####
import argparse

parser = argparse.ArgumentParser(description="Recherche des cutsets minimaux à l'aide de la librairie cobamp")
output = parser.add_argument_group('output')
output.add_argument('-o', '--output', type=str, nargs= 1, required= True, help="Définit le chemin du fichier de sortie, le nom d\'un dossier se termine par un '/'. Exemple: ../Data/")
dataset = parser.add_argument_group('dataset')
dataset.add_argument('-i', '--input', type= str, nargs= '+', required= True, help="Définit le fichier d'entrée au format sbml|xml")
dataset.add_argument('-m', '--metabolites', type= str, nargs= '?', default= False, help="Définit le fichier de métabolites à utiliser en tant que réaction de biomasse")
parameters = parser.add_argument_group('parameters')
parameters.add_argument('-I', '--iterative', action="store_true", help= "Définit ITERATIVE comme le type d'algorithme à utiliser (au lieu de POPULATE)")
parameters.add_argument('-s', '--stop_criteria', type= int, nargs= '?', default= 2, help= "Critère d'arrêt du programme (taille des cutsets pour POPULATE ou nombre de cutsets pour ITERATIVE")
parameters.add_argument('-n', '--thread', type= int, default= 1, nargs= '?', help= "Nombre de coeur à utiliser pour le programme")
parser.add_argument('-v', '--version', action='version', version=__version__)
args = parser.parse_args()

### IMPORTATION DES LIBRAIRIES ###

import cobra
import cbmpy
import time
import pickle
import os, sys
from threading import Thread, Event
from cobamp.wrappers import KShortestMCSEnumeratorWrapper
from cobra.io import validate_sbml_model

#### IMPORTATION DU MODELE ####

if len(args.input) == 1:
    print('\n\n-----Importation et conversion du modèle par CBMPY-----\n\n')
else: print('\n\n-----Importation et conversion des modèles par CBMPY-----\n\n')

#### Parsing du nom du modèle d'entrée pour la création des fichiers de sortie:
model_file_name = args.input.split("/")[-1].split(".xml")[0]

for fichier in args.input:

    converted_model, errors = validate_sbml_model(args.input)
    if len(errors["COBRA_ERROR"]) > 0:
        print(errors["COBRA_ERROR"])
    

    #### Ajout de la réaction de biomasse si le fichier texte nécessaire est rempli/présent ####

    if args.metabolites :                                                           # Vérification de la présence du fichier contenant les métabolites à ajouter. 
                                                                                    # S'il est absent, le modèle conservera sa biomasse originale.
        biomass_name = 'artificial_biomass'
        biomass_components_manual = []
        with open(args.metabolites, 'r') as input_metabolites :
            metabolites_list = input_metabolites.readlines()
        
            for metabolite_to_add in metabolites_list :
                biomass_components_manual.append(metabolite_to_add.strip('\n'))

        reaction_biom = cobra.Reaction("artificial_biomass")
        ### Ajout des bornes pour définir la réaction comme irréversible
        reaction_biom.lower_bound = 0.0
        reaction_biom.upper_bound = 1000.0
        temp001x= cobra.Metabolite(
            'temp001x',
            formula='',
            name='biomass_',
            compartment='C_c')
        #Ajout de la réaction de sortie de la biomasse :
        output_biom = cobra.Reaction('EX_temp001x')
        output_biom.lower_bound = 0.0
        output_biom.upper_bound = 1000.0
        # Ajout du métabolite qui sera sorti par la réaction de sortie de biomasse.
        output_biom.add_metabolites({temp001x :-1.0})
        converted_model.add_reaction(output_biom)
        converted_model.add_reaction(reaction_biom)
        converted_model.reactions.artificial_biomass.add_metabolites({temp001x : 1.0})

        for metabolite in biomass_components_manual :
            converted_model.reactions.artificial_biomass.add_metabolites({converted_model.metabolites.get_by_id(str(metabolite)) : -1.0}) # Ajout des métabolites un par un.
                                                                                                                                          # Il faut bien sûr que les id des métabolites correspondent
                                                                                                                                          # à des métabolites déjà présents dans le modèle.
    else :
        biomass_name = "biomass_components"
    #### COBAMP ####    

    print('\n\n-----Looking for cutsets-----\n\n')
    yield_space = {
        # Not used here, so we kept it empty
            }

    flux_space = {
       # flux id     : (lower_bound, upper_bound)
        biomass_name: (1.0, None)
    }
    # Récupération du type d'algorithme à utiliser (Si l'argument commence par p alors on prend POPULATE)    
    if args.iterative:
        algorithme = 'ITERATIVE'
    else: algorithme = 'POPULATE'
    print(f'-----with ALGORITHM_TYPE_{algorithme}-----')
    # Récupération des solutions en fonction du type d'algorithme
    if not args.iterative:
        ksefm = KShortestMCSEnumeratorWrapper(
            model=converted_model,
            target_flux_space_dict=flux_space,
            target_yield_space_dict=yield_space,
            algorithm_type=KShortestMCSEnumeratorWrapper.ALGORITHM_TYPE_POPULATE,
            n_threads=args.thread,
            stop_criteria=args.stop_criteria,
            solver='CPLEX') # or 'GUROBI', 'GLPK'

    else:
        ksefm = KShortestMCSEnumeratorWrapper(
            model=converted_model,
            target_flux_space_dict=flux_space,
            target_yield_space_dict=yield_space,
            algorithm_type=KShortestMCSEnumeratorWrapper.ALGORITHM_TYPE_ITERATIVE,
            n_threads=args.thread,
            stop_criteria=args.stop_criteria,
            solver='CPLEX') # or 'GUROBI', 'GLPK'


    enumerator = ksefm.get_enumerator()
    t1 = time.time()

    solution = []
    count=0

    print('-----Computing-----')

    for sol in enumerator :
        solution.append(sol)
        t2 = time.time()
        if args.iterative:
            count+=1
        else: count+=len(sol)
        print(f'Search in progress, cutset found : {count} | Time passed since cutset search launch : {t2-t1}s')

    print('-----Creating the pickle-----')
    
    f=open(f'{args.output[0]}{model_file_name}_mcs.pkl' ,'wb')
    pickle.dump(solution, f)
    f.close()

    t2 = time.time()
    print(f'\ncutsets search complete, {model_file_name} took :, {t2-t1}\nFound {len(solution)} MCS.')
