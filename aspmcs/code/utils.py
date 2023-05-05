from itertools import combinations
import warnings
def knock_out(model, lr):
    model = model.copy()
    for r in lr:
        lb = False
        r = r[:] # remove quotes and ending bracket
        if r.startswith('mcs_'):
            r = r[4:]
        if r.startswith('R_'):
            r = r[2:]
        if r.endswith('_rev'):
            r = r[:-4]
            lb = True
        gr = model.reactions.get_by_id(r)
        if lb:
            gr.lower_bound = 0
        else:
            gr.upper_bound = 0
    return model

def sol_after_ko(lr, model, obj, check_reaction=None):
    try:
        model_c = knock_out(model, lr)
        model_c.objective = obj
        opt = model_c.optimize()
        if check_reaction:
            print(f"Reaction {check_reaction} : flux = {model_c.reactions.get_by_id(check_reaction).flux}")
        return opt.objective_value if opt.status != 'infeasible' else 0.0
    except KeyError as e:
        warnings.warn(str(e))
        return 0.0

def nm1_combinations(orig_lr, model, saves):
    all_lr = {}
    for lr in itertools.combinations(orig_lr, len(orig_lr)-1):
        fs = frozenset(lr)
        if fs not in saves:
            saves[fs] = sol_after_ko(lr, model)
        all_lr[fs] = saves[fs]
        if all_lr[fs] == 0.0:
            saves[frozenset(orig_lr)] = 0.0
            return 0.0, all_lr, fs
    orig = sol_after_ko(orig_lr, model)
    saves[frozenset(orig_lr)] = orig
    return orig, all_lr, None


def parse_output(path):
    with open(path) as output_mcs_buffer:
        output_mcs_data = output_mcs_buffer.readlines()

    cutset_list = []
    for line in output_mcs_data:
        # if "Answer" in line:
        #     print(line)
        #     answer_line_found = True
        # elif answer_line_found:
        #     print(line)
        #     answer_line_found = False
        if "cutset" in line:
            reacs = line.strip("\n").split("cutset(\"")[1:]

            cutset_list.append([rid.split("\")")[0] for rid in reacs])
    return cutset_list


def decompress_cutsets(subsets_file_path, comp_cutsets_list):
    with open(subsets_file_path, 'r') as buffer_subsets:
        reaction_subsets_dict = eval(buffer_subsets.read())

    decomp_cutsets_list = []
    for reac_list in comp_cutsets_list:
        #Itère sur chaque cutset

        cs_decomp = []
        cs_decomp_comb = []
        for reac in reac_list:
            reac_id = reac[4:]
            original_reac_subset = reaction_subsets_dict[reac_id]["reacs"]
            cs_decomp.append(original_reac_subset)
            cs_decomp_comb.append([c for c in combinations(cs_decomp,len(reac_list))])
        decomp_cutsets_list.append(cs_decomp_comb)
            
    return decomp_cutsets_list
