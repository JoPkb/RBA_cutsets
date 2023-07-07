import numpy, cobra
import efmtool_link.efmtool_intern as efmtool_intern
import efmtool_link.efmtool_extern as efmtool_extern
import efmtool_link.efmtool4cobra as efmtool4cobra
from pathlib import Path
import pickle


def flux_variability_analysis(model: cobra.Model, loopless=False, fraction_of_optimum=0.0,
                              processes=None, results_cache_dir: Path=None, fva_hash=None, print_func=print):
    # all bounds in the model must be finite because the COBRApy FVA treats unbounded results as errors
    if results_cache_dir is not None:
        fva_hash.update(pickle.dumps((loopless, fraction_of_optimum, model.tolerance))) # integrate solver tolerances?
        if fraction_of_optimum > 0:
            fva_hash.update(pickle.dumps(model.reactions.list_attr("objective_coefficient")))
            fva_hash.update(model.objective_direction.encode())
        file_path = results_cache_dir / (model.id+"_FVA_"+fva_hash.hexdigest())
        fva_result = None
        if Path.exists(file_path):
            try:
                fva_result = pandas.read_pickle(file_path)
                print_func("Loaded FVA result from", str(file_path))
            except:
                print_func("Loading FVA result from", str(file_path), "failed, running FVA.")
        else:
            print_func("No cached result available, running FVA...")
        if fva_result is None:
            fva_result = cobra.flux_analysis.flux_variability_analysis(model, reaction_list=None, loopless=loopless,
                                                             fraction_of_optimum=fraction_of_optimum,
                                                             pfba_factor=None, processes=processes)
            try:
                fva_result.to_pickle(file_path)
                print_func("Saved FVA result to ", str(file_path))
            except:
                print_func("Failed to write FVA result to ", str(file_path))
        return fva_result
    else:
        return cobra.flux_analysis.flux_variability_analysis(model, reaction_list=None, loopless=loopless,
                                                             fraction_of_optimum=fraction_of_optimum,
                                                             pfba_factor=None, processes=processes)


def fva_model(compressed_model):
    fva_tolerance=1e-9
    with model as fva: # can be skipped when a compressed model is available
        # when include_model_bounds=False modify bounds so that only reversibilites are used?
        # fva.solver = 'glpk_exact' # too slow for large models
        fva.tolerance = fva_tolerance
        fva.objective = model.problem.Objective(0.0)
        if fva.problem.__name__ == 'optlang.glpk_interface':
            # should emulate setting an optimality tolerance (which GLPK simplex does not have)
            fva.solver.configuration._smcp.meth = GLP_DUAL
            fva.solver.configuration._smcp.tol_dj = fva_tolerance
        elif fva.problem.__name__ == 'optlang.coinor_cbc_interface':
            fva.solver.problem.opt_tol = fva_tolerance
        fva_res = flux_variability_analysis(fva, fraction_of_optimum=0.0, processes=1, 
            results_cache_dir=None, fva_hash=None)
    for i in range(fva_res.values.shape[0]): # assumes the FVA results are ordered same as the model reactions
        if abs(fva_res.values[i, 0]) > fva_tolerance: # resolve with glpk_exact?
            compressed_model.reactions[i].lower_bound = fva_res.values[i, 0]
        else:
            compressed_model.reactions[i].lower_bound = 0
        if abs(fva_res.values[i, 1]) > fva_tolerance: # resolve with glpk_exact?
            compressed_model.reactions[i].upper_bound = fva_res.values[i, 1]
        else:
            compressed_model.reactions[i].upper_bound = 0
    return compressed_model




# all parameters here for simpler utilisation
model_file = 'models_storage/Hep-G2_v2.xml' # model file name
dir_name = './' # parent directory of the directory where mparser_cli is
module_name = './' # directory where mparser_cli is
generated_files_path = '../compressed/' # path where the generated files go
model_name = 'new' # prefix name for all generated files (do not put a path there)

model = cobra.io.read_sbml_model(model_file)
original_model = model.copy()

# network compression (currently combination of reaction subsets only)
# IMPORTANT: the model is modified by this function, if you want to keep the full model copy it first
model = fva_model(model)
subT = efmtool4cobra.compress_model_sympy(model) # subT is a matrix for conversion of flux vectors between the full and compressed model
rd = cobra.util.array.create_stoichiometric_matrix(model, array_type='lil')
# model compression makes sure that irreversible reactions always point in the forward direction
rev_rd = [int(r.reversibility) for r in model.reactions]
len(model.reactions)