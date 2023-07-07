import cobra
import pickle
import sys
def test_cutsets_comb(decompressed_cs_list, model, objective, bounds_modifications=None):
    i = 0
    invalid_cutsets_comb_index = []
    invalid_cutsets_index = {}
    fluxes_per_cs = []
    not_in_healthy = []
    for cutset_reaction_combinations in decompressed_cs_list:
        j = 0
        cs_fluxes = {}
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
                        model_buffer.reactions.get_by_id(
                            rid).bounds = (0.0, 0.0)

                    except KeyError:
                        not_in_healthy.append(rid)

                # Test that the target pathway works
                model_buffer.objective = objective
                sol = model_buffer.optimize()
                if sol.objective_value > 1e-6:
                    # print(f"\nCutset {reactions_id} has no impact on neoglucogenesis")
                    pass
                else:
                    print(
                        f"\nCutset {cutset_combination} prevents {objective}")
                    invalid_cutsets_comb_index.append(
                        cutset_reaction_combinations.index(cutset_combination))

                for reaction in model_buffer.reactions:
                    cs_fluxes[reaction] = reaction.flux
                fluxes_per_cs.append(cs_fluxes)
                # for reaction in model.reactions:
                #     cs_fluxes[reaction] = reaction.fluxes

                # restoring bounds

        if len(invalid_cutsets_index) > 0:
            invalid_cutsets_index[decompressed_cs_list.index(
                cutset_reaction_combinations)] = invalid_cutsets_comb_index
        j += 1

    i += 1
    print(f"Reactions from cancer model absent from healthy_model :\n {set(not_in_healthy)}\n")
    return invalid_cutsets_comb_index, fluxes_per_cs



if __name__ == "__main__":
    output_1 = sys.argv[1]
    output_2 = sys.argv[2]
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
        invalid_cs_index, flux2cutset = test_cutsets_comb(decompressed_cs_list, h_c, objective_reaction, bounds_modification)
        objective2result[objective_reaction] = (invalid_cs_index, flux2cutset)
        i+=1

    with open(output_1, "wb") as out_invalid_cs:
        pickle.dump(invalid_cs_index, out_invalid_cs)
    with open(output_2, "wb") as out_flux2cutsets:
        pickle.dump(flux2cutset, out_flux2cutsets)
