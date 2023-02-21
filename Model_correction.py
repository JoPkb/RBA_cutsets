import cobra
from datetime import datetime
import requests
import json
from tqdm import tqdm

### Function to copy-paste genes from one model to another :

def copy_genes(origin_model, target_model) :
    
    for target_reaction in target_model.reactions :
        try :
            target_reaction.gene_reaction_rule = origin_model.reactions.get_by_id(target_reaction.id).gene_reaction_rule
            print(f"\nDEBUG -- Added gene_reaction_rule to reaction {target_reaction.id}")
        except KeyError :
            print(f"\nERROR -- KeyError, The reaction from target model {target_reaction.id} might not have a counterpart in origin_model.")
    return target_model

def fix_formulas(model) :

    for m in model.metabolites :
       
        m.formula = ""
    
    return model


### Gene id conversion :
def get_ids(model) :
    url="http://rest.ensembl.org/xrefs/id/"
    #print("\nCreating working copy of the model...")
    #converted = model.copy()
    #print("\nModel copied.")
    error_messages = []
    for gene in tqdm(model.genes) :
        gene_id = gene.id
    
        # Checking if the gene id is an ENSEMBL id :
        if len(gene_id) == 15 :

            #print(f"\nDEBUG : getting ncbi gene id for {gene_id}")
            response = requests.get(url+gene_id+"?external_db=EntrezGene")

            try :

                # Parsing of the ncbi gene ID from the text response
                new_id = response.text.split("primary_id: ")[1].split("\n")[0]
                gene.annotation["ncbigene"] = new_id
                #print(f"\nAdded ncbi gene id {new_id} to annotations of gene {gene_id}")

            except IndexError :
                    current_time = datetime.now().strftime("%H:%M:%S")
                    error_messages.append(f"\n[{current_time}] ERROR - No ncbi gene id found for gene {gene_id}.")
                    
                    pass
        else :
            current_time = datetime.now().strftime("%H:%M:%S")
            error_messages.append(f"\n[{current_time}]ERROR - Gene id {gene_id} is not an ENSEMBL id.")
            
    for message in error_messages :
        print(message)
    """
    for reaction in converted.reactions :
        print(f"\nChanging gene_reaction_rule for reaction {reaction.id}")
        reaction_genes = list(reaction.genes)
        string = ""
        for gene_name in reaction_genes :
            if reaction_genes.index(gene_name) != len(reaction_genes)-1 :
                string += f"{gene_name} or "
            else :
                string += f"{gene_name}"
        reaction.gene_reaction_rule = string
        #reaction.gene_name_reaction_rule = string"""

    #return converted
"""
### Function to add biomass_reaction :
def add_biomass_reaction(model, biomass_metabolites) :

    #biomass_metabolites = './Hep-G2/Hep-G2.xml-5a91f1955e8a9c8a5c5d2ff4a1737c21/utils_biomass_metabolites.txt'
    biomass_components_manual = []
    with open(biomass_metabolites, 'r') as input_metabolites :
        metabolites_list = input_metabolites.readlines()
        
        for metabolite_to_add in metabolites_list :
            biomass_components_manual.append(metabolite_to_add.strip('\n'))

    reaction = cobra.Reaction('artificial_biomass')
    reaction.lower_bound = 0.0
    reaction.upper_bound = 1000.0
    temp001x= cobra.Metabolite(
        'temp001x',
        formula='',
        name='biomass_',
        compartment='C_c')
    #Ajout de la sortie de la biomasse
    output_biom = cobra.Reaction('EX_temp001x')
    output_biom.lower_bound = 0.0
    output_biom.upper_bound = float("inf")

    output_biom.add_metabolites({temp001x :-1.0})
    model.add_reactions([output_biom])
    model.add_reactions([reaction])
    model.reactions.artificial_biomass.add_metabolites({temp001x : 1.0})

    for metabolite in biomass_components_manual :
        model.reactions.artificial_biomass.add_metabolites({model.metabolites.get_by_id(str(metabolite)) : -1.0})
    return model




def remove_gene_dupes(json_in, json_out) :
    with open(json_in) as input_file :
        data = json.load(input_file)

    unique_ids = []
    for g in data["genes"] :
        id = g["id"]
        if id not in unique_ids :
            unique_ids.append(id)
        else : 
            print(f"{id} already found somewhere else.")
    data["genes"] = [{"id" : id, "name" : ""} for id in unique_ids]

    with open(json_out, "w") as output_file :
        json.dump(data, output_file, indent="")
"""