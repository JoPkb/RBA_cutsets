import cobra
import cbmpy
import requests

### Gene id conversion :
def convert_ids(model) :
    url="http://rest.ensembl.org/xrefs/id/"
    converted = model.copy()

    for gene in converted.genes :
        gene_id = gene.id

        # Checking if the gene id is an ENSEMBL id :
        if len(gene_id) == 15 :

            print(f"\nDEBUG : getting uniprot gene id for {gene_id}")
            response = requests.get(url+gene_id+"?external_db=Uniprot_gn")

            try :

                # Parsing of the uniprot ID from the text response
                uniprot_id = response.text.split("primary_id: ")[1].split("\n")[0]
                gene.id = uniprot_id
                print(f"\nGene ensembl id {gene_id} changed to uniprot id : {uniprot_id}.")

            except IndexError :
                    
                    print(f"\nNo uniprot id found for gene {gene_id}.")
                    pass
        else :
            
            print(f"\nGene already has uniprot id : {gene_id}")

    return converted

### Function to add biomass_reaction :
def add_biomass_reaction(model, biomass_metabolites) :

    biomass_metabolites = './Hep-G2/Hep-G2.xml-5a91f1955e8a9c8a5c5d2ff4a1737c21/utils_biomass_metabolites.txt'
    biomass_components_manual = []
    with open(biomass_metabolites, 'r') as input_metabolites :
        metabolites_list = input_metabolites.readlines()
        
        for metabolite_to_add in metabolites_list :
            biomass_components_manual.append(metabolite_to_add.strip('\n'))

    reaction = cobra.Reaction('artificial_biomass')
    reaction.lower_bound = 0.0
    reaction.upper_bound = 1000.0
    temp001x= cobra.Metabolite(
        'temp001x',
        formula='',
        name='biomass_',
        compartment='C_c')
    #Ajout de la sortie de la biomasse
    output_biom = cobra.Reaction('EX_temp001x')
    output_biom.lower_bound = 0.0
    output_biom.upper_bound = 1000.0

    output_biom.add_metabolites({temp001x :-1.0})
    model.add_reactions([output_biom])
    model.add_reactions([reaction])
    model.reactions.artificial_biomass.add_metabolites({temp001x : 1.0})

    for metabolite in biomass_components_manual :
        model.reactions.artificial_biomass.add_metabolites({model.metabolites.get_by_id(str(metabolite)) : -1.0})
    return model

### Function to copy-paste genes from one model to another :

def copy_genes(origin_model, target_model) :
    
    for target_reaction in target_model.reactions :
        try :
            target_reaction.genes = origin_model.reactions.get_by_id(target_reaction.id).genes
        except KeyError :
            print(f"ERROR -- KeyError, The reaction form target model {target_reaction.id} might not have a counterpart in origin_model.")
    return target_model