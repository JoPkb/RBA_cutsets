# commande de lancement bash 
# python run_straindesign.py ../../aspefm_mcs_checker/data/model.xml rsub_93 output_path.txt > output_cnapy.txt 2> output_cnapy_err.txt
import straindesign as sd
import cobra
import sys
from cobra.io import validate_sbml_model
import logging
import pickle
model_path = sys.argv[1]
objective = sys.argv[2]
output_path = sys.argv[3]
#print(sys.argv)

#model, errors = validate_sbml_model(model_path)
#model_copy = model.copy()

#model_copy.objective = objective

model = cobra.io.read_sbml_model(model_path)

#module_suppress = sd.SDModule(model_copy,sd.names.SUPPRESS,constraints='objective <= 1e-3')
module_suppress = sd.SDModule(model, sd.names.SUPPRESS, constraints=f'{objective}>=0.001')
module_protect = sd.SDModule(model, sd.names.PROTECT, constraints='')

logging.basicConfig(level=logging.INFO)
# Compute strain designs
sols = sd.compute_strain_designs(model,
                                 sd_modules = [module_suppress, module_protect],
                                 time_limit = 7200,
                                 max_cost = 5,
                                 solution_approach = sd.names.ANY)
# Print solutions
print(f"Solutions found, expanded to {len(sols.reaction_sd)} solutions in the uncompressed netork.")
print(f"Example knockout set: {[s for s in sols.reaction_sd[0]]}")

with open(output_path, "wb") as p_out :
    pickle.dump(sols, p_out)
print(f"Wrote solutions object to {output_path}.\nDONE.")
