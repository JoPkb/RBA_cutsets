import cobra
import pickle

def test_cutsets(decompressed_cs_list, model, objective, bounds_modifications = None):
    i = 0
    invalid_cutsets_comb_index = []
    invalid_cutsets_index = {}
    fluxes_per_cs = []
    not_in_healthy = []
    # Before anything, we copy the model, and apply the bounds changes if needed
    model_working_copy = model.copy()
    if bounds_modifications:
        for reaction, bounds in bounds_modifications:
            model_working_copy.reactions.get_by_id(reaction).bounds = bounds

    for cutset_reaction_combinations in decompressed_cs_list:
        j = 0
        
        cs_fluxes = {}
        #Take the first cutset combination :
        cutset_combination = cutset_reaction_combinations[0]
  
        #reactions_bounds = []
        #modified_reactions_bounds = []
        with model as model_buffer:
                        
            # Apply bounds modifications
            if bounds_modifications:
                for rid, bounds in bounds_modifications:
                    model_buffer.reactions.get_by_id(rid).bounds = bounds

            for rid in cutset_combination:
                
                #print("rid: "+str(rid))
                if rid.endswith("rev"): 
                    rid = rid[:-4]
                try:
                    #reactions_bounds.append(model_buffer.reactions.get_by_id(rid).bounds)
                    model_buffer.reactions.get_by_id(rid).bounds = (0.0,0.0)
                    
                except KeyError:
                    not_in_healthy.append(rid)
            

            model_buffer.objective = objective
            sol = model_buffer.optimize()

            if sol.objective_value > 1e-6:
                #print(f"\nCutset {reactions_id} has no impact on neoglucogenesis")
                pass
            else:
                print(f"\nCutset {cutset_combination} prevents {objective}")
                invalid_cutsets_comb_index.append(cutset_reaction_combinations.index(cutset_combination))
        
            for reaction in model_buffer.reactions:
                cs_fluxes[reaction] = reaction.flux
            fluxes_per_cs.append(cs_fluxes)
        # #restoring bounds
        # for rid, rbounds in zip(cutset_combination, reactions_bounds):
        #     if rid.endswith("rev"):
        #         rid = rid[:-4]
        #     try:
        #         model.reactions.get_by_id(rid).bounds = rbounds
        #     except KeyError:
        #         pass
                
        if len(invalid_cutsets_comb_index) > 0:
            invalid_cutsets_index[decompressed_cs_list.index(cutset_reaction_combinations)] = invalid_cutsets_comb_index
        j+=1

    i+=1
    print(f"\nReactions in cancer model absent from healthy model : \n{set(not_in_healthy)}")
    return invalid_cutsets_comb_index, fluxes_per_cs


if __name__ == "__main__":
    print(f"\nReading Cutsets...")
    decompressed_cs_list_path = "../results/vanilla_run/res_final/merged_decompressed_cutsets.p"
    with open(decompressed_cs_list_path, "rb") as input_mcs:
        decompressed_cs_list = pickle.load(input_mcs)
    print(f"\nReading model...")
    healthy=cobra.io.read_sbml_model("../../models_storage/iHep_updated_v011.xml")

    h_c = healthy.copy()
    print(f"\nReading metabolic functions...")
    with open("metabolic_functions.txt", "r") as metab_func:
        metabolic_functions = eval(metab_func.read())

    objective2result = {}
    i=1
    for objective_reaction, bounds_modification in metabolic_functions.items():
        print(f"\nTesting Objective n°{i}/{len(metabolic_functions)}")
        invalid_cs_index, flux2cutset = test_cutsets(decompressed_cs_list, h_c, objective_reaction, bounds_modification)
        objective2result[objective_reaction] = (invalid_cs_index, flux2cutset)
        i+=1

    with open("../results/invalid_cs_index.p", "wb") as out_invalid_cs:
        pickle.dump(invalid_cs_index, out_invalid_cs)
    with open("../results/flux2cutsets.p", "wb") as out_flux2cutsets:
        pickle.dump(flux2cutset, out_flux2cutsets)