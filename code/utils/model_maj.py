import cobra
from cobra.io import validate_sbml_model
import utils.Model_correction as mc
import cbmpy
import sys
import importlib
import os

importlib.reload(mc)
def maj(dir,source_model=None, target_model=None, bounds_check =True, genes_id_copy = True, alt_gene_ids = True, metab_id_check = True, bounds_value_check = True, subsystem_copy = True, add_exchanges = True, get_names = True, check_biom_reac = True) :

    print("\nLoading model to update...\n")
    target_model_file_name = target_model.split(".xml")[0]
    

    # If all bounds are the same, load the model with 
    if bounds_check :
        print("\nLoading with cbmpy")
        if not os.path.exists(dir+target_model_file_name+"_cbmpy.xml") :
            target_model_temp = cbmpy.CBRead.readSBML2FBA(dir + target_model)
            cbmpy.CBWrite.writeSBML3FBCV2(target_model_temp, dir+target_model_file_name+"_cbmpy.xml")
            del(target_model_temp)
            working_target_model, errors = validate_sbml_model(dir+target_model_file_name+"_cbmpy.xml")
        else :
            working_target_model, errors = validate_sbml_model(dir+target_model_file_name+"_cbmpy.xml")

        print("Cobra model loading errors :\n")
        print(errors["COBRA_ERROR"])
    # ^ Maybe not needed : can just use target_model_temp ? 
    else :

        working_target_model, errors = validate_sbml_model(dir+target_model)
        print("Cobra model loading errors :\n")
        print(errors["COBRA_ERROR"])
    
    if genes_id_copy or subsystem_copy:
        print("\nLoading source model for genes or subsystem copy : \n")
        source_model_cobra, errors_source = validate_sbml_model(dir+source_model)
    else :
        pass

    if genes_id_copy :
        # TEST # if no genes in model.genes : (else, dont do it)
        print("\nCopying genes from model model to target model : \n")
        working_target_model = mc.copy_genes(source_model_cobra, working_target_model)
    else :
        pass

    
    if alt_gene_ids :
        # TEST # if gene.annotation exists : (else : dont)
        # THEN TEST # if gene.annotation["ensembl"] and/or gene.annotation["ncbigene"] are present
        print("\nGetting ENSEMBL and ncbi genes IDs to the annotations of the target model : \n")
        for gene in working_target_model.genes :
            if len(gene.id) == 15 and "ENSG" in gene.id[:4] :
                gene.annotation["ensembl"] = str(gene.id)
        
        mc.get_ids(working_target_model, "EntrezGene", "ncbigene")
    else :
        pass

    if metab_id_check :
        print("\nCleaning up the metabolites list in the target model : \n")
        to_remove = []
        for metabolite in working_target_model.metabolites :
            if "m" not in metabolite.id :
                to_remove.append(metabolite)
        working_target_model.remove_metabolites(to_remove)
    else :
        pass
    
    if bounds_value_check :
        print(f"\nNumber of metabolites in target model : {len(working_target_model.metabolites)}")
        print("\nSetting the reactions to finite bounds in the target model : \n")

        for reaction in working_target_model.reactions :
            ub = 1000.0 if reaction.upper_bound == float("inf") or reaction.upper_bound == 1000.0 else 0.0
            lb = -1000.0 if reaction.lower_bound == float("-inf") or reaction.lower_bound == -1000.0 else 0.0
            #print(f"DEBUG --- reaction {reaction.id}\noriginal bounds : ({reaction.lower_bound},{reaction.upper_bound})\nnew bounds : ({lb},{ub})\n\t########\n")
            reaction.bounds = (lb,ub)
    else :
        pass

    if subsystem_copy :
        print(f"\nGetting subsystem information from source model.")
        mc.get_subsystem(source_model_cobra, working_target_model)
    else :
        pass

    print("\n\nDONE\n\n")


    if add_exchanges :
        print(f"\nCopying exchange reactions from source model.")
        mc.get_exchanges_reactions(target_model=working_target_model)

    else :
        pass
    
    if get_names :
        print(f"\nGetting reactions names from genes ID.")
        mc.get_names_from_genes(working_target_model)
    else :
        pass

    if check_biom_reac :
        print("\nChecking the presence of cofactors and vitamins in the biomass reaction of the target model.")
        mc.check_biomass_metabolites(target_model)
        if "m01602c" in [m.id for m in working_target_model.reactions.biomass_components.metabolites] :
            working_target_model.reactions.biomass_components.add_metabolites({"m01602c" : 0.0}, combine=False)
        else:
            pass
    else :
        pass
    return working_target_model


    

if __name__ == "__main__" :
    converted_target_model = maj(sys.argv[1], sys.argv[2], sys.argv[3])

    cobra.io.write_sbml_model(converted_target_model, sys.argv[3].split(".xml")[0] + "_updated.xml")




"""
MISC

    # TEST # if all bounds are te same : same_bounds = True, else False
    same_bounds = True
    for i in range(len(target_model.reactions)) :
        if i == 0 :
            bound_to_compare = target_model.reactions[i].bounds
        else :
            if bound_to_compare != target_model.reactions[i].bounds :
                same_bounds = False


"""