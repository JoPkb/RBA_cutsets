import os
import pickle
from collections import defaultdict
from cobra.core.model import Model



def get_mcs_counters(mcs_directory, output_sorting) :
    """
    Getting counter dictionaries of MCS results.
    returns 3 dictionaries of MCS counts : 
    overall_dMCS_sizes : counter of the length of MCS over all positive constraints runs
        Note : this dictionary effectively contains the length distribution of compressed mcs.
    per_file_dMCS_sizes : counter of the length of MCS for each individual positive constraint run
    per_file_dMCS : Dictionary linking the name of the positive constraint and the list of MCS associated to it.
    """
    #mcs_directory="../../../HCC_compressed_2/mcs_posconstr_analysis_2/mcs_decomp/"

    mcs_file_list = os.listdir(mcs_directory)
    overall_MCS_sizes = {}
    per_file_dMCS_sizes= {}
    per_file_dMCS = {}
    per_file_sat = {}
    per_constraint_sat = {}
    with open(output_sorting, "r") as sat_file:
        sat_lines = sat_file.readlines()
        for line in sat_lines:
            mcs_filename = line.split(": ")[-1].strip("\n")
            sat = line.split(" file :")[0]
            per_file_sat[mcs_filename] = sat
        

    for mcs_file in mcs_file_list:
        print(f"Reading {mcs_file}")
        with open(f"{mcs_directory}{mcs_file}", "rb") as mcs_in:
            mcs_list = pickle.load(mcs_in)
            if len(mcs_list) > 0:
                dMCS_sizes = {}
                for mcs in mcs_list:
                    if len(mcs) > 0 :
                        dMCS_sizes[len(mcs[0])] = dMCS_sizes.get(len(mcs[0]), 0 ) + 1
                        overall_MCS_sizes[len(mcs[0])] = overall_MCS_sizes.get(len(mcs[0]), 0) + 1

                
                mcs_log_file = mcs_file.split("_decomp")[0] + ".txt"
                with open("../../../HCC_compressed_2/mcs_posconstr_analysis_3/mcs_res/"+mcs_log_file, "r") as mcs_log:
                    posconstr_line = mcs_log.readlines()[0]
                    posconstr = posconstr_line.split('support("')[-1].split('")')[0]
                try:
                    per_constraint_sat[posconstr] = per_file_sat[mcs_log_file]
                except KeyError:
                    per_constraint_sat[posconstr] = "NA"
                per_file_dMCS_sizes[posconstr] = dMCS_sizes
                per_file_dMCS[posconstr] = mcs_list
    return overall_MCS_sizes, per_file_dMCS_sizes, per_file_dMCS, per_file_sat, per_constraint_sat

def get_gene_associated_mcs(genes_IDs:list[str], model:Model, mcs_list) :
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
    
    gene_associated_reactions_list = [elem for elem in gene_associated_reactions.keys()]
    # MCS sorting
    wanted_mcs = []
    posconstr_found = []
    posconstr2mcs = defaultdict(list)
    #mcs2posconstr = defaultdict(list)
    posconstr2mcs_sizes = {}
    for mcs in mcs_list:
        match=False
        for reaction in mcs :
            if reaction in gene_associated_reactions_list:
                match=True
                posconstr_found.append(reaction)
                #mcs2posconstr[mcs].append(reaction)
                posconstr2mcs[reaction].append(mcs)
        if match:
            wanted_mcs.append(mcs)

    return wanted_mcs, posconstr2mcs

def parse_cnapy_res(res_file) :
    cnapy_mcs_reactions = []
    with open(res_file, "r") as cnapy_res_buffer:
        cnapy_mcs = eval(cnapy_res_buffer.read())

    for mcs_dict in cnapy_mcs:
        mcs = []
        for key in mcs_dict.keys():
            mcs.append(key)

        cnapy_mcs_reactions.append(mcs)
    return cnapy_mcs_reactions

