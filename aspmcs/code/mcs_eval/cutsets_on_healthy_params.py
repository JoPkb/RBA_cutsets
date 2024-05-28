import cobra
import pickle
import sys
import os

# File formatted by Black extension.


def test_cutsets_comb(
    decompressed_cs_list, model, objective, bounds_modifications=None
):
    """
    Function used to test the effect of MCS on a model's metabolic functions.
    INPUT:
    - decompressed_cs_list: list of decompressed Minimal Cut Sets.
    - model: Cobra metabolic model
    - objective: string corresponding to the ID of a reaction from the cobra metabolic model.
    - bounds_modifications: Optional, if some additional constraints on reaction bounds
    need to be added in order to test a metabolic pathway's functionality, they are to be given here.

    OUTPUT:
    - invalid_cutsets_comb_index: list containing the index of MCS from the original decompressed_cs_list
    that prevent the realization of one of the tested metabolic function.
    - valid_mcs : list of list of tuples (following mcs input file format) of valid MCS combinations
    """

    i = 0
    invalid_cutsets_comb_index = []
    invalid_cutsets_index = {}
    #fluxes_per_cs = []
    not_in_healthy = []
    valid_mcs = []
    for cutset_reaction_combinations in decompressed_cs_list:
        j = 0
        cs_fluxes = {}
        valid_combinations = []
        # For each combination of the cutset:
        for cutset_combination in cutset_reaction_combinations:
            # print("reactions_id: " + str(reactions_id))
            # reactions_bounds = []
            
            with model as model_buffer:
                # Apply bounds modifications
                if bounds_modifications:
                    for rid, bounds in bounds_modifications:
                        model_buffer.reactions.get_by_id(rid).bounds = bounds

                # Cut reactions from MCS
                for rid in cutset_combination:
                    # print("rid: "+str(rid))
                    if rid.endswith("rev"):
                        rid = rid[:-4]
                    try:
                        # reactions_bounds.append(model_buffer.reactions.get_by_id(rid).bounds)
                        model_buffer.reactions.get_by_id(rid).bounds = (0.0, 0.0)

                    except KeyError:
                        not_in_healthy.append(rid)

                # Test that the target pathway works
                model_buffer.objective = objective
                sol = model_buffer.optimize()
                if sol.objective_value > 1e-6:
                    valid_combinations.append(cutset_combination)
                else:
                    print(f"\n[{decompressed_cs_list.index(cutset_reaction_combinations)}/{len(decompressed_cs_list)}] - Cutset {cutset_combination} prevents {objective}", flush=True)
                    invalid_cutsets_comb_index.append(
                        cutset_reaction_combinations.index(cutset_combination)
                    )
        if len(valid_combinations) > 0 :
            print(f"[{decompressed_cs_list.index(cutset_reaction_combinations)}/{len(decompressed_cs_list)}] - Found {len(valid_combinations)} valid MCS combinations. : {valid_combinations}", flush=True)
            valid_mcs.append(valid_combinations)
            
                # for reaction in model_buffer.reactions:
                #     cs_fluxes[reaction] = reaction.flux
                # fluxes_per_cs.append(cs_fluxes)
                # for reaction in model.reactions:
                #     cs_fluxes[reaction] = reaction.fluxes

                # restoring bounds

        if len(invalid_cutsets_index) > 0:
            invalid_cutsets_index[
                decompressed_cs_list.index(cutset_reaction_combinations)
            ] = invalid_cutsets_comb_index
        j += 1

    i += 1
    print(
        f"Reactions from cancer model absent from healthy_model :\n {set(not_in_healthy)}\n"
    )
    return invalid_cutsets_comb_index, valid_mcs


if __name__ == "__main__":
    #This scripts runs a check on the viability of each metabolic function in the healthy cell after MCS reactions KOs.
    #For each metabolic function, the results are stored as an entry of the dictionary objective2result. 
    #If MCS are contained in multiple files, they will be read one by one and one result dictionary will be created for each.
    #The result dictionary is written to a pickle at the end.

    decompressed_cs_list_path = sys.argv[1]
    model_path = sys.argv[2]
    metabolic_functions_path = sys.argv[3]
    output_folder = sys.argv[4]
    print(len(sys.argv))
    print(f"decomp_mcs : {decompressed_cs_list_path}\nmodel_path : {model_path}\nmetabolic_functions_path : {metabolic_functions_path}, \noutput_folder : {output_folder}", flush=True)
    if len(sys.argv) > 5 :
        print(f"\n5 : {sys.argv[5]}", flush=True)
        multiple = True
    else :
        print(f"\n5 : None", flush=True)
        multiple = False

    print(f"\nReading model...", flush=True)
    healthy = cobra.io.read_sbml_model(model_path)
    h_c = healthy.copy()
    print("Finished reading model...", flush=True)

    if multiple : 
        print("Starting MCS eval from multiple files", flush=True)
        # If we have multiple mcs files to deal with, the input path should be that of a folder.
        for mcs_file in os.listdir(decompressed_cs_list_path) :
            file_path = decompressed_cs_list_path+mcs_file
            file_number = mcs_file.split("decomp_mcs")[0][7:]

            with open(file_path, "rb") as input_mcs:
                decompressed_cs_list = pickle.load(input_mcs)

            print(f"\nReading metabolic functions...", flush=True)
            with open("metabolic_functions.txt", "r") as metab_func:
                metabolic_functions = eval(metab_func.read())

            objective2result = {}
            i = 1
            for objective_reaction, bounds_modification in metabolic_functions.items():
                print(f"\nTesting Objective n°{i}/{len(metabolic_functions)}", flush=True)
                invalid_cs_index, valid_mcs = test_cutsets_comb(
                    decompressed_cs_list, h_c, objective_reaction, bounds_modification
                )
                objective2result[objective_reaction] = (invalid_cs_index, valid_mcs)
                i += 1

            # In case the output folder str ends with a / :
            if output_folder.endswith("/"):
                output_folder = output_folder[:-1]
            with open(f"{output_folder}/mcs_res_{file_number}.p", "wb") as out_results:
                pickle.dump(objective2result, out_results)
            

    else :
        print("Starting MCS eval from a single file", flush=True)
        print(f"\nReading Cutsets...", flush=True)
        with open(decompressed_cs_list_path, "rb") as input_mcs:
            decompressed_cs_list = pickle.load(input_mcs)

        print(f"\nReading metabolic functions...", flush=True)
        with open("metabolic_functions.txt", "r") as metab_func:
            metabolic_functions = eval(metab_func.read())
        objective2result = {}
        i = 1
        for objective_reaction, bounds_modification in metabolic_functions.items():
            print(f"\nTesting Objective n°{i}/{len(metabolic_functions)}", flush=True)
            invalid_cs_index, valid_mcs = test_cutsets_comb(
                decompressed_cs_list, h_c, objective_reaction, bounds_modification
            )
            objective2result[objective_reaction] = (invalid_cs_index, valid_mcs)
            i += 1

        # In case the output folder str ends with a / :
        if output_folder.endswith("/"):
            output_folder = output_folder[:-1]
        with open(f"{output_folder}/mcs_res.p", "wb") as out_results:
            pickle.dump(invalid_cs_index, out_results)

