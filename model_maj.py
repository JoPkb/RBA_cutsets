import cobra
from cobra.io import validate_sbml_model
import Model_correction as mc
import cbmpy
import sys

def maj(source_model, target_model, dir) :

    print("\nLoading model to update :\n")
    target_model_file_name = target_model.split(".xml")[0]

    target_model_temp = cbmpy.CBRead.readSBML2FBA(dir + target_model)
    cbmpy.CBWrite.writeSBML3FBCV2(target_model_temp, dir+target_model_file_name+"_cbmpy.xml")
    del(target_model_temp)
    
    target_model_cbmpy, errors = validate_sbml_model(dir+target_model_file_name+"_cbmpy.xml")
    print(errors["COBRA_ERROR"])
    # ^ Maybe not needed : can just use target_model_temp ? 

    print("\nLoading model model : \n")
    source_model_cobra, errors_source = validate_sbml_model(dir+source_model)

    print("\nCopying genes from model model to target model : \n")
    target_model_cbmpy = mc.copy_genes(source_model_cobra, target_model_cbmpy)


    target_model_cbmpy.reactions.biomass_components.add_metabolites({"m01602c" : 0.0}, combine=False)

    print("\nGetting ENSEMBL and ncbi genes IDs to the annotations of the target model : \n")
    for gene in target_model_cbmpy.genes :
        if len(gene.id) == 15 :
            gene.annotation["ensembl"] = str(gene.id)
    
    mc.get_ids(target_model_cbmpy, "EntrezGene", "ncbigene")

    print("\nCleaning up the metabolites list in the target model : \n")
    to_remove = []
    for metabolite in target_model_cbmpy.metabolites :
        if "m" not in metabolite.id :
            to_remove.append(metabolite)
    target_model_cbmpy.remove_metabolites(to_remove)
    
    print(f"\nNumber of metabolites in target model : {len(target_model_cbmpy.metabolites)}")
    print("\nSetting the reactions to finite bounds in the target model : \n")
    for reaction in target_model_cbmpy.reactions :
        ub = 1000.0 if reaction.upper_bound == float("inf") or reaction.upper_bound == 1000.0 else 0.0
        lb = -1000.0 if reaction.lower_bound == float("-inf") or reaction.lower_bound == -1000.0 else 0.0
        #print(f"DEBUG --- reaction {reaction.id}\noriginal bounds : ({reaction.lower_bound},{reaction.upper_bound})\nnew bounds : ({lb},{ub})\n\t########\n")
        reaction.bounds = (lb,ub)
    print("\n\nDONE\n\n")
    return target_model_cbmpy
if __name__ == "__main__" :
    converted_target_model = maj(sys.argv[1], sys.argv[2], sys.argv[3])

    cobra.io.write_sbml_model(converted_target_model, sys.argv[3].split(".xml")[0] + "_updated.xml")