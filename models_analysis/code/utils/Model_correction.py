import cobra
from datetime import datetime
import requests
import json
from tqdm import tqdm
from cobra.core.model import Model
from cobra.core.reaction import Reaction

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
def get_ids(model, external_db, annotation_key):
    url="http://rest.ensembl.org/xrefs/id/"
    #print("\nCreating working copy of the model...")
    #converted = model.copy()
    #print("\nModel copied.")
    error_messages = []
    if len(model.genes) == 0 :
        print(f"\nERROR -- No genes in model {model.name}")
        return 0
    for gene in tqdm(model.genes) :
        gene_id = gene.id
        #print(f"Checking out gene {gene_id}")
        # Checking if the gene id is an ENSEMBL id :
        if len(gene_id) == 15 :

            #print(f"\nDEBUG : getting ncbi gene id for {gene_id}")
            response = requests.get(url+gene_id+"?external_db="+external_db)

            try :

                # Parsing of the ncbi gene ID from the text response
                new_id = response.text.split("primary_id: ")[1].split("\n")[0]
                gene.annotation[annotation_key] = new_id
                #print(f"\nAdded ncbi gene id {new_id} to annotations of gene {gene_id}")

            except IndexError :
                    current_time = datetime.now().strftime("%H:%M:%S")
                    error_messages.append(f"\n[{current_time}] ERROR - No {annotation_key} gene id found for gene {gene_id}.")
                    
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


def get_subsystem(source_model, target_model) :
    for source_reaction in source_model.reactions :
        source_reaction_id = source_reaction.id

        if len(source_reaction.subsystem) != 0 :
            try :
                print(f"\nGetting subsystem information from source model for reaction {source_reaction_id}")
                target_model.reactions.get_by_id(source_reaction_id).notes["SUBSYSTEM"] = source_reaction.notes["SUBSYSTEM"]
                target_model.reactions.get_by_id(source_reaction_id).subsystem = source_reaction.notes["SUBSYSTEM"]
            except KeyError :
                print(f"\nKeyError for reaction {source_reaction_id} from source model.")
                pass
        else :
            print(f"\nNo subsystem found for reactions {source_reaction_id}")



def get_names_from_genes(target_model) :
    url = "http://rest.ensembl.org/lookup/id/"
    archive_url = "http://rest.ensembl.org/archive/id/"
    for reaction in tqdm(target_model.reactions) :
        try :
            gene_ids = [gene.id for gene in reaction.genes]
        except :
            reaction_name = "Null"
            reaction.name = reaction_name
            continue

        names = []
        for gene_id in gene_ids :
            response = requests.get(url+str(gene_id))
            try :
                reaction_name = response.text.split("display_name: ")[1].split("\n")[0]
                #print(f"{gene_id} : {reaction_name}")
            except IndexError :
                
                archive_response = requests.get(archive_url+gene_id)
                try :
                    new_gene_id = archive_response.text.split("stable_id: ")[1].split("\n")[0]
                    #print(f"\nGot an error with gene id {gene_id}. retrying with updated gene_id : {new_gene_id}", flush=True)
                
                    
                    new_response = requests.get(url+new_gene_id)
                
                    reaction_name = new_response.text.split("display_name: ")[1].split("\n")[0]
                except IndexError :
                    reaction_name = "Null"
            if len(reaction_name) != 0 :
                names.append(reaction_name)
        if len(names) != 0 :
            #print(f"\nFound Enzyme name {reaction_name} for reaction {reaction.id}")
            reaction.name = "/".join(names)
        else :
            continue


def get_exchanges_reactions(target_model) :
    boundary_reactions = target_model.boundary
    reactions_to_add = []
    for b_r in boundary_reactions :
        for metab, sto in b_r.metabolites.items() :
            # Creating boundary metabolite
            lb = b_r.lower_bound
            ub = b_r.upper_bound
            new_metab_x = cobra.core.metabolite.Metabolite(metab.id[:-1]+"x", name = metab.name, compartment = "C_x")
            if sto > 0.0 :
                b_r.add_metabolites({new_metab_x: -1.0}) # Intake

            else :
                b_r.add_metabolites({new_metab_x : 1.0}) # Secretion
            new_reaction_x = cobra.core.reaction.Reaction(id = "EX_" + new_metab_x.id, name = "EX_" + new_metab_x.id, subsystem= "Exchange reactions", lower_bound= lb, upper_bound = ub)
            new_reaction_x.add_metabolites({new_metab_x : 1.0})
            reactions_to_add.append(new_reaction_x)
    target_model.add_reactions(reactions_to_add)


def check_biomass_metabolites(target_model, metabolite_id = "_") :
    print(f"\nChecking if the metabolite [cofactors and vitamins] is included in the biomass_components reaction of the model.")

    biomass_reaction = target_model.reactions.biomass_components

    for m in biomass_reaction.metabolites :
        if "cofactors and vitamins" in m.name or m.id == metabolite_id:
            id = m.id
            biomass_reaction.add_metabolites({target_model.metabolites.get_by_id(id) : 0.0}, combine=False)
            print(f"\nFound unwanted metabolite in biomass_components reactions : {id}. Removed it.")

###- MEDIUM DEFINITION AND APPLICATION FUNCTIONS -###


def export_medium(model:Model, outfile:str):
    """
    Function exporting the cobra metabolic model's medium
    in a tabular text file with semicolon separators.
    Used to manually specify which intakes the model takes.
    Type 1 at the end of a line to include that reaction, 
    leave empty to leave it out.
    """
    str = ""
    for exchange_reaction in model.medium:
        reaction = model.reactions.get_by_id(exchange_reaction)
        metab = [m for m in reaction.metabolites.keys()][0]
        str += f"\n{reaction.id};{metab.name};"
    
    with open(outfile, "w") as out_medium :
        out_medium.write(str)


def parse_medium_definition(file_path:str, model:Model):
    """
    Function parsing the tabular text file defining the medium
    for the cobra metabolic model. Each line with a "1" at the end 
    will get included in the "reactions to keep" list, the others will
    not.
    Returns a tuple of two lists : exchange reactions to remove, exchange reactions to keep.
    """
    to_remove = []
    to_keep = []
    model.objective = "biomass_components"
    with open(file_path, "r") as medium_data:
        medium = medium_data.read()

    for line in medium.split("\n")[1:]:
        reaction_id = line.split(";")[0]
        reaction_presence = line.split(";")[2]

        if reaction_presence == "1":

            reaction_to_keep = model.reactions.get_by_id(reaction_id)
            m = [m.name for m in reaction_to_keep.metabolites][0]
            # print(f"\n{reaction_id} -- {m} is present in medium")
            to_keep.append(reaction_id)

        else:

            reaction_to_remove = model.reactions.get_by_id(reaction_id)
            m = [m.name for m in reaction_to_remove.metabolites][0]
            # print(f"\nReaction {reaction_id} -- {m} not in medium.")
            to_remove.append(reaction_id)

    return (to_remove, to_keep)


def change_reaction_bounds(model: Model, reactions: dict[str:tuple[float, float]]):
    """
    Function used to cut one or more reactions in a cobra metabolic model.
    Takes as input a cobra model object and a dictionary containing reaction
    IDs and their desired new bounds.
    """
    for reaction_id, bounds in reactions.items():
        reaction_to_change = model.reactions.get_by_id(reaction_id)
        reaction_to_change.bounds = bounds


def apply_medium_safe(model:Model, to_remove:list[str], objective_function:str, threshold:float):
    """
    Function used to apply the medium constraint on the model
    given as input, from the list of reactions to remove given as input.
    The list should contain strings of reactions IDs from the model.
    The function removes the reactions one by one, while checking that
    the objective function of the model stays solvable. If the objective_value
    decreases by more than the set threshold given as input, then the reaction 
    will be restored.
    The function returns a tuple of the modified model and 
    a list of all the reaction IDs that have been restored.
    """
    #Creating a copy of the model
    modified_model = model.copy()

    modified_model.objective = objective_function
    objective_value = modified_model.optimize().objective_value
    #For each reaction, we cut it, check the biomass production, then, if it does not
    #respect our criteria, it is restored.
    reactions_restored = []

    for reaction_id in to_remove:
        rdict = {reaction_id : (0.0,1000.0)}
        change_reaction_bounds(modified_model, rdict)
        obj_val_tmp = modified_model.optimize().objective_value
        rdict_restore = {reaction_id : (-1000.0,1000.0)} 

        if obj_val_tmp == 0 :
            print(f"\nWhen cutting reaction {reaction_id}, biomass fell to 0.")
            change_reaction_bounds(modified_model, rdict_restore)
            reactions_restored.append(reaction_id)

        elif objective_value - obj_val_tmp > threshold:
            print(f"\nWhen cutting reaction {reaction_id}, biomass was reduced by {objective_value-obj_val_tmp}")
            change_reaction_bounds(modified_model, rdict_restore)
            reactions_restored.append(reaction_id)

        else:
            print(f"Minimal consequences when cutting reaction {reaction_id}, biomass changed by {objective_value-obj_val_tmp}")
    print(f"\nMedium definition finished. Removed {len(to_remove)-len(reactions_restored)} reactions out of {len(to_remove)}")
    return modified_model, reactions_restored


def get_reac_id_from_genes(target_model, source_model):
    
    pass



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