import straindesign as sd
import cobra
import sys
from cobra.io import validate_sbml_model
import logging
import pickle
model_path = sys.argv[1]
objective = sys.argv[2]
output_path = sys.argv[3]
if sys.argv[4]:
    time_limit = sys.argv[4]
else:
    time_limit = 61200

model, errors = validate_sbml_model(model_path)
model_copy = model.copy()

model_copy.objective = objective

module_suppress = sd.SDModule(model_copy,sd.names.SUPPRESS,constraints='objective <= 1e-3')

logging.basicConfig(level=logging.INFO)
# Compute strain designs
sols = sd.compute_strain_designs(model_copy,
                                 time_limit = 61200,
                                 sd_modules = [module_suppress],
                                 max_cost = -5,
                                 solution_approach = sd.names.ANY)
# Print solutions
print(f"Solutions found, expanded to {len(sols.reaction_sd)} solutions in the uncompressed netork.")
print(f"Example knockout set: {[s for s in sols.reaction_sd[0]]}")

with open(output_path, "wb") as p_out :
    pickle.dump(sols, p_out)
print(f"Wrote solutions object to {output_path}.\nDONE.")
