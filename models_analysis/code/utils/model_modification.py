# Creates a dictionary of reactions associated to overexpressed genes.
def get_gene_associated_reaction(genes_IDs:list[str], model:Model) :
    """
    Function used to get all reactions associated to the given list of gene IDs and only to gene IDs contained in the given list
    gene_IDs:list[str]
    model:cobra.core.model.Model
    """
    # Definition of local variables to return
    gene_associated_reactions = {}
    n_errors = 0
    
    # Iterating over all model reactions. For each reaction, the function checks if their associated gene(s) are included in the list of genes.
    for reaction in model.reactions:
        associated = True
        try:
            model_reaction_genes = [g.annotation['ncbigene'] for g in reaction.genes]
        except KeyError:
            n_errors += 1
        
        if len(model_reaction_genes) > 0:
            for gene in model_reaction_genes:
                if gene not in genes_IDs:
                    associated = False
        else:
            associated = False

        # the boolean associated will only stay True if all the genes associated to a reaction are found in the given gene list.        
        if associated == True:
            gene_associated_reactions[reaction.id] = model_reaction_genes
    
    return gene_associated_reactions, n_errors

# Safely removing down-regulated genes-associated reactions from model :
def remove_and_check(model:Model, reactions:list[str], threshold = float(0.2), objective = "biomass_components"):

    reactions_removed = 0
    total_reactions = len(reactions)
    model.objective = objective

    for reaction in reactions:
        reaction_obj = model.reactions.get_by_id(reaction)
        model.remove_reactions([reaction_obj])
        new_val = model.optimize().objective_value
        print(f"obj_val : {new_val}; {type(new_val)}")
        if new_val <= threshold:
            model.add_reactions([reaction_obj])
        else:
            reactions_removed += 1
    print(f"{reactions_removed}/{total_reactions} reactions removed.")

# Reads medium definition file and applies the bounds to reactions of the model.
def apply_medium(model:Model, medium_file_path:str):
    new_medium = {}
    with open(medium_file_path, "r") as medium_buffer:
        medium_data = medium_buffer.read()
    
    for line in medium_data.split("\n"):
        print(line)
        reaction_id = line.split(";")[0]
        reaction_activity = line.split(";")[2]
        if reaction_activity == "1":
            new_medium[reaction_id] = abs(float(line.split(";")[3]))
        else:
            new_medium[reaction_id] = float(0.0)
    
    model.medium = new_medium


#### CONSTRAINT CREATION TAKING ONLY REACTIONS THAT ARE EXCLUSIVELY ASSOCIATED WITH OVEREXPRESSED GENES ####
#Writes lp4 constraint files by matching overexpressed gene-associated reactions to their compressed counterpart.
def write_constraint_file_overexpressed(reactions_list, reaction_subset_path, outfile_path):
    # Takes as input a list of reactions that are associated to overexpressed genes
    # The function then selects compressed reactions which only contain these reactions.

    with open(reaction_subset_path, "r") as reaction_subset:
        data = reaction_subset.read()
        r_subsets = eval(data)

    # downregulated_reactions_list = [reaction_id for reaction_id in reactions_list]
    compressed_gene_associated_reactions = []
    for compressed_reaction_name, compressed_reaction_data in r_subsets.items():
        associated = True
        #print(f"{compressed_reaction_name} :\n")
        for reaction in compressed_reaction_data["reacs"]:
            if "rev" in reaction :
                reaction = reaction[:-4]
            if reaction not in reactions_list :
                #print(f"{reaction} not in list")
                associated = False # If one of the compressed reaction subsets is not linked to an overexpressed reaction, the compressed reaction is not selected.

        if associated == True :
            print(f"{compressed_reaction_name} is gene-associated.")
            compressed_gene_associated_reactions.append(compressed_reaction_name)

    
    prefix = ":- not support"
    suffix = ".\n"

    for compressed_reaction_id in compressed_gene_associated_reactions:
        reaction_number = str(compressed_gene_associated_reactions.index(compressed_reaction_id))
        constr_out_str = ""
        constr_out_str += f'{prefix}("mcs_{compressed_reaction_id}"){suffix}'
        constr_out_str += f'target("mcs_{compressed_reaction_id}").'
        with open(outfile_path+reaction_number+"_posconstr.lp4", "w") as constr_out:
            constr_out.write(constr_out_str)


#### CONSTRAINT CREATION TAKING INTO ACCOUNT EVERY REACTION SUBSET THAT CONTAIN AT LEAST ONE REACTION LINKED TO AN OVEREXPRESSED GENE ####
def write_constraint_file(reactions_list, reaction_subset_path, outfile_path):
    # Takes as input a list of reactions that are associated to overexpressed genes
    # The function then selects compressed reactions which only contain these reactions.

    with open(reaction_subset_path, "r") as reaction_subset:
        data = reaction_subset.read()
        r_subsets = eval(data)

    # downregulated_reactions_list = [reaction_id for reaction_id in reactions_list]
    compressed_gene_associated_reactions = []
    for compressed_reaction_name, compressed_reaction_data in r_subsets.items():
        associated = False
        #print(f"{compressed_reaction_name} :\n")
        for reaction in compressed_reaction_data["reacs"]:
            if "rev" in reaction :
                reaction = reaction[:-4]
            if reaction in reactions_list :
                #print(f"{reaction} not in list")
                associated = True # If one of the compressed reaction subsets is not linked to an overexpressed reaction, the compressed reaction is not selected.

        if associated == True :
            print(f"{compressed_reaction_name} is gene-associated.")
            compressed_gene_associated_reactions.append(compressed_reaction_name)

    
    prefix = ":- not support"
    suffix = ".\n"

    for compressed_reaction_id in compressed_gene_associated_reactions:
        reaction_number = str(compressed_gene_associated_reactions.index(compressed_reaction_id))
        constr_out_str = ""
        constr_out_str += f'{prefix}("mcs_{compressed_reaction_id}"){suffix}'
        constr_out_str += f'target("mcs_{compressed_reaction_id}").'
        with open(outfile_path+reaction_number+"_posconstr.lp4", "w") as constr_out:
            constr_out.write(constr_out_str)
