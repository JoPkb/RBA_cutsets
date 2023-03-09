import pandas as pd
from contextlib import redirect_stdout
import sys

def print_exchanges(optimized_model, filter, ) : 
    intakes = []
    secretions = []
    neutrals = []
    if filter == "non_null" :
        flux_comparison = 0.0
    elif filter == "all" :
        flux_comparison = 200000.0
    for reaction in optimized_model.boundary :
        if reaction.flux != flux_comparison :
            

            # Jolification 
            spaces = 12
            spaces_str = ""
            for i in str(round(reaction.flux)) :
                spaces -= 1
            if "EX_t" in reaction.id :
                spaces -= 1
            for j in range(spaces) :
                spaces_str += " "
            # Fin de la jolification

            ml = [metab for metab in reaction.metabolites]
            m = [metab for metab in reaction.metabolites][0]    
            if reaction.flux > 0.0 :
                secretions.append(f"{reaction.id} : {round(reaction.flux)}{spaces_str}ub : {reaction.upper_bound}\t---\t\tmetabolites : \t id : {m.id} --- metabolite name : {m.name} ; id : {m.id}")
            elif reaction.flux < 0.0 :
                intakes.append(f"{reaction.id} : {round(reaction.flux)}{spaces_str}ub : {reaction.upper_bound}\t---\t\tmetabolites : \t id : {m.id} --- metabolite name : {m.name} ; id : {m.id}")
            else :
                neutrals.append(f"{reaction.id} : {round(reaction.flux)}{spaces_str}ub : {reaction.upper_bound}\t---\t\tmetabolites : \t id : {m.id} --- metabolite name : {m.name} ; id : {m.id}")
                
    intakes.sort(key=lambda f : float(f.split(": ")[1].split("ub")[0]))
    secretions.sort(key=lambda f : float(f.split(": ")[1].split("ub")[0]), reverse=True)
    print("\n##########\nINTAKES :\n")
    for i in intakes :
        print(i)
    print("\n##########\nSECRETIONS :\n")
    for s in secretions :
        print(s)
    print("\n##########\nNEUTRALS : \n")
    for n in neutrals :
        print(n)


def print_reactions(reaction, flux = 0.0, v=True):
    '''
    Affiche les reactions d'un modele de façon plus lisible
    '''
    liste_reactif = [i for i in reaction.reactants]
    liste_produit = [i for i in reaction.products]
    string = ""
    string2 = ""
    if not "EX_" in reaction.id :
        if flux >= float('0'):
                fleche = "-->"
        elif flux <= float('0') :
                fleche = "<--"
        else:
                fleche = "<=>"
    else :
        if flux >= float('0') :
            fleche = "<--"
        elif flux <= float('0') :
             fleche = "-->"
        else :
             fleche = "<=>"
    
    for i in liste_reactif:
        if i != liste_reactif[-1]:
            string += str(float(abs(reaction.metabolites[i])))
            string += " "
            string += str(i.name)
            string += " + "
            
            string2 += str(float(abs(reaction.metabolites[i])))
            string2 += " "
            string2 += str(i.id)
            string2 += " + "
            
        else:
            string += str(float(abs(reaction.metabolites[i])))
            string += " "
            string += str(i.name)
            string += " "
            
            string2 += str(float(abs(reaction.metabolites[i])))
            string2 += " "
            string2 += str(i.id)
            string2 += " "
            
    string += fleche
    string += " "
    
    string2 += fleche
    string2 += " "
    
    
    for i in liste_produit:
        if i != liste_produit[-1]:
            string += str(float(abs(reaction.metabolites[i])))
            string += " "
            string += str(i.name)
            string += " + "
            
            string2 += str(float(abs(reaction.metabolites[i])))
            string2 += " "
            string2 += str(i.id)
            string2 += " + "
            
        else:
            string += str(float(abs(reaction.metabolites[i])))
            string += " "
            string += str(i.name)
            
            string2 += str(float(abs(reaction.metabolites[i])))
            string2 += " "
            string2 += str(i.id)
    if v :    
        print(string)
        print("")
        print(string2)
    else :
        return(string, string2)
    





def parcours_test(reaction, flux_dict, v = False) :
    
    #Metabolites list, used to determine which metabolites have already been visited.
    m_l = []
    
    for m in reaction.metabolites :
        print(f"\nChecking out metabolite {m.id}") 
        
        if not m in m_l :
            print(f"\nNew metabolite, adding {m.id} to list of visited metabolites.")
            m_l.append(m)

            for r in m.reactions :
                print(f"\nChecking out reaction {r.id}")

                if r.id not in flux_dict.keys() :

                    # This part checks if there is a flux for the given reaction, and if it is an exchange reaction.
                    # If there is a flux and it is not an exchange reaction, it is added to the flux dict.
                    # If there is a flux and it is an exchange reaction, it is only added to the flux dict and the recursion starts back at the next reaction.
                    # If it was an exchange reaction involving C_x or C_s, it behaves normally.
                    if r.flux != 0.0 and len(r.compartments) == 1:
                        print(f"\nNew reaction, adding {r.id} to flux dict.")
                        flux_dict[r.id] = r.flux
                        flux_dict = parcours(r, flux_dict)

                    elif r.flux != 0.0 and len(r.compartments) >1 :
                        print(f"\nNew exchange reaction, adding {r.id} : (exchange {r.compartments}) to flux dict.")
                        flux_dict[r.id] = r.flux
                        if "C_x" in r.compartments or "C_s" in r.compartments :
                            flux_dict = parcours(r, flux_dict) 
                    else :
                        print(f"\nERROR -- flux == 0 for {r.id}")
                        
                else :
                    print(f"\nERROR -- id in dict for {r.id}")
                    continue
        else :
            print(f"\nERROR -- metabolite {m.id} already visited")
            break
        
            
    return flux_dict




       
def parcours(reaction, flux_dict, max_iterations = 10000, i=0,v = True) :
    #Metabolites list, used to determine which metabolites have already been visited.
    m_l = []
    
    """
    if not v :
        orig_stdout = sys.stdout
        f = open('out.txt', 'w')
        sys.stdout = f"""
    for m in reaction.metabolites :
        print(f"\nChecking out metabolite {m.id}") 
        
        if not m in m_l and not i >= max_iterations:
            print(f"\nNew metabolite, adding {m.id} to list of visited metabolites.")
            m_l.append(m)

            for r in m.reactions :
                print(f"\nChecking out reaction {r.id}")

                if r.id not in flux_dict.keys() :

                    # This part checks if there is a flux for the given reaction, and if it is an exchange reaction.
                    # If there is a flux and it is not an exchange reaction, it is added to the flux dict.
                    # If there is a flux and it is an exchange reaction, it is only added to the flux dict and the recursion starts back at the next reaction.
                    # If it was an exchange reaction involving C_x or C_s, it behaves normally.
                    if r.flux != 0.0 :
                        print(f"\nNew reaction, adding {r.id} to flux dict.")
                        flux_dict[r.id] = r.flux
                        flux_dict = parcours(r, flux_dict, max_iterations, i=i)
                        i +=1
                    else :
                        print(f"\nERROR -- flux == 0 for {r.id}")
                        
                else :
                    print(f"\nERROR -- id in dict for {r.id}")
                    continue
        else :
            print(f"\nERROR -- metabolite {m.id} already visited")
            continue
    
    """if not v :
        f.close()
        sys.stdout = orig_stdout"""
            
    return flux_dict

def run_parcours(reaction, model, max_iterations=10000) :
    flux_dict = {}

    f = parcours(reaction, flux_dict, max_iterations)
    for reaction, flux in f.items() :
        print_reactions(model.reactions.get_by_id(reaction), flux)
        print(f"FLUX : {flux} --- ID : {model.reactions.get_by_id(reaction).id} --- compartment : {model.reactions.get_by_id(reaction).compartments} \n\n---\n\n")

def build_reaction_df(optimized_model) :

    reactions_list = [r for r in optimized_model.reactions]

    compartments_reactions_dict = {}

    
       
    for compartment in optimized_model.compartments : 
        compartments_reactions_dict[str(compartment)] = {
            "flux" : [abs(r.flux) for r in reactions_list if str(compartment) in r.compartments],\
            "subSystem" : [r.subsystem for r in reactions_list if str(compartment) in r.compartments],\
            "id" : [r.id for r in reactions_list if str(compartment) in r.compartments],\
            "name" : [r.name for r in reactions_list]
        } 
    
    return (pd.DataFrame(compartments_reactions_dict["C_c"]), \
            pd.DataFrame(compartments_reactions_dict["C_r"]), \
            pd.DataFrame(compartments_reactions_dict["C_s"]), \
            pd.DataFrame(compartments_reactions_dict["C_m"]), \
            pd.DataFrame(compartments_reactions_dict["C_p"]), \
            pd.DataFrame(compartments_reactions_dict["C_x"]), \
            pd.DataFrame(compartments_reactions_dict["C_l"]), \
            pd.DataFrame(compartments_reactions_dict["C_g"]), \
            pd.DataFrame(compartments_reactions_dict["C_n"]))
    


    

    