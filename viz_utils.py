
def print_exchanges(optimized_model, filter) : 
    intakes = []
    secretions = []
    neutrals = []
    if filter == "non_null" :
        flux_comparison = 0.0
    elif filter == "all" :
        flux_comparison = 200000.0
    for reaction in optimized_model.reactions :
        if reaction.flux != flux_comparison :
            if "EX_" in reaction.id :

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



def parcours(reaction, flux_dict, v = False) :
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
                    if r.flux != 0.0 :
                        print(f"\nNew reaction, adding {r.id} to flux dict.")
                        flux_dict[r.id] = r.flux
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

def run_parcours(reaction) :
    flux_dict = {}

    return parcours(reaction, flux_dict)