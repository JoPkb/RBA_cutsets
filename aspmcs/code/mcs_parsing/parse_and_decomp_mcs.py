from utils import parse_output, decompress_cutsets, sol_after_ko
import cobra
import pandas as pd
from itertools import combinations, combinations_with_replacement, permutations, chain, product
from tqdm import tqdm
import sys
import pickle

def itertools_product(listofcutsets, model):
    i = 1
    good_cutsets = []
    for cutset_orID in tqdm(listofcutsets):
        good_cutsets_subs = []
        #print(f"Cutset : {cutset_orID}\n")
        cuts_combinations = list(product(*cutset_orID)) # Décompression des réactions retrouvées dans les cutsets
        for comb in cuts_combinations:
            with model as m:
                retval = sol_after_ko(comb, m, obj=objective)
                # for reac in cs:
                #     m.reactions.get_by_id(reac.strip("_rev"))
                #opt = m.optimize()
            #print("<",opt.objective_value,">", cs )
            if retval < 1e-6:
                good_cutsets_subs.append(comb)
            else:
                print(f"{comb} is not a cutset.")

            if len(comb) > 1 :
                for nm1 in combinations(comb, len(comb)-1):
                    retval = sol_after_ko(nm1, m)
                    if retval < 1e-6:
                        print('[N-1]',"Objective value = ", retval,"For n-1 subset\
                           of cutset >>>", nm1, 'not MIN cutset')
        good_cutsets.apend(good_cutsets_subs)
        i+=1
    return good_cutsets

if __name__ == "__main__":
    model_path = sys.argv[1]
    cutsets_path = sys.argv[2]
    reaction_subsets = sys.argv[3]
    output_path = sys.argv[4]
    model = cobra.io.read_sbml_model(model_path)
    cutsets_list = parse_output(cutsets_path)
    with open(reaction_subsets, 'r') as buffer_subsets:
        reaction_subsets_dict = eval(buffer_subsets.read())
    listofcutsets_original_id = []
    for reac_list in cutsets_list:
        #Pour chaque cutset:

        cutset_original_IDs = []
        #print(f"\n\nFor cutset {reac_list}:\n\n")
        for comp_reac_id in reac_list:
            #Pour chaque réaction compressée du cutset
            #On veut récupérer le subset de réactions pré-compression.
            reac_subset = reaction_subsets_dict[comp_reac_id[4:]]["reacs"]
            cutset_original_IDs.append(reac_subset)
        listofcutsets_original_id.append(cutset_original_IDs)
    
    good_cutsets = itertools_product(listofcutsets_original_id, model)
    with open(output_path, "wb") as out:
        pickle.dump(good_cutsets,out)