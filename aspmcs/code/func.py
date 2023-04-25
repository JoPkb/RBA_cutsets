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

def sol_after_ko(lr, model):
    try:
        model_c = knock_out(model, lr)
        opt = model_c.optimize()
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
    with open("../results/output_subset_local_13.txt") as output_mcs_buffer:
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
