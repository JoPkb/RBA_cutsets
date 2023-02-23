
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