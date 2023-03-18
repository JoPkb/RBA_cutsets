import pandas as pd
from contextlib import redirect_stdout
from tqdm import tqdm
import sys
import plotly.express as px
from itertools import cycle
import matplotlib.pyplot as plt
import cobra

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
    for m in tqdm(reaction.metabolites) :
        #print(f"\nChecking out metabolite {m.id}") 
        
        if not m in m_l and not i >= max_iterations:
            #print(f"\nNew metabolite, adding {m.id} to list of visited metabolites.")
            m_l.append(m)

            for r in m.reactions :
                #print(f"\nChecking out reaction {r.id}")

                if r.id not in flux_dict.keys() :

                    # This part checks if there is a flux for the given reaction, and if it is an exchange reaction.
                    # If there is a flux and it is not an exchange reaction, it is added to the flux dict.
                    # If there is a flux and it is an exchange reaction, it is only added to the flux dict and the recursion starts back at the next reaction.
                    # If it was an exchange reaction involving C_x or C_s, it behaves normally.
                    if r.flux != 0.0 :
                        #print(f"\nNew reaction, adding {r.id} to flux dict.")
                        flux_dict[r.id] = r.flux
                        flux_dict = parcours(r, flux_dict, max_iterations, i=i)
                        i +=1
                    else :
                        #print(f"\nERROR -- flux == 0 for {r.id}")
                        pass
                        
                else :
                    #print(f"\nERROR -- id in dict for {r.id}")
                    continue
        else :
            #print(f"\nERROR -- metabolite {m.id} already visited")
            continue
    
    """if not v :
        f.close()
        sys.stdout = orig_stdout"""
            
    return flux_dict

def run_parcours(reaction : cobra.core.reaction.Reaction, model : cobra.core.model.Model, max_iterations=10000) :
    flux_dict = {}

    f = parcours(reaction, flux_dict, max_iterations)
    for reaction, flux in f.items() :
        print_reactions(model.reactions.get_by_id(reaction), flux)
        print(f"FLUX : {flux} --- ID : {model.reactions.get_by_id(reaction).id} --- compartment : {model.reactions.get_by_id(reaction).compartments} \n\n---\n\n")

def build_reaction_df(optimized_model, by_compartment = True) :

    reactions_list = [r for r in optimized_model.reactions]
    if by_compartment :
        compartments_reactions_dict = {}

        
        
        for compartment in optimized_model.compartments : 
            compartments_reactions_dict[str(compartment)] = {
                "flux" : [abs(r.flux) for r in reactions_list if str(compartment) in r.compartments],\
                "subSystem" : [r.subsystem for r in reactions_list if str(compartment) in r.compartments],\
                "id" : [r.id for r in reactions_list if str(compartment) in r.compartments],\
                "name" : [r.name for r in reactions_list if str(compartment) in r.compartments],\
                "compartment" : [compartment for r in reactions_list if str(compartment) in r.compartments],
                "direction" : [str(r.flux)[0] for r in reactions_list if str(compartment) in r.compartments],\
                "reactants" : [" + ".join([m.name for m in r.reactants]) for r in reactions_list if str(subsystem) in r.subsystem],\
                "products" : [" + ".join([m.name for m in r.products]) for r in reactions_list if str(subsystem) in r.subsystem]
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
    else :

        subsystem_reactions_dict = {}
        dataframes_to_return = []
        for subsystem in optimized_model.groups :
            subsystem_reactions_dict[str(subsystem)] = {
                "flux" : [abs(r.flux) for r in reactions_list if str(subsystem) in r.subsystem],\
                "subSystem" : [r.subsystem for r in reactions_list if str(subsystem) in r.subsystem],\
                "id" : [r.id for r in reactions_list if str(subsystem) in r.subsystem],\
                "name" : [r.name for r in reactions_list if str(subsystem) in r.subsystem],\
                "compartment" : [str([comp for comp in r.compartments][0]) for r in reactions_list if str(subsystem) in r.subsystem],\
                "direction" : [str(r.flux)[0] for r in reactions_list if str(subsystem) in r.subsystem],\
                "reactants" : [" + ".join([m.name for m in r.reactants]) for r in reactions_list if str(subsystem) in r.subsystem],\
                "products" : [" + ".join([m.name for m in r.products]) for r in reactions_list if str(subsystem) in r.subsystem]
            }

        return subsystem_reactions_dict
    



def get_subsystem_fluxes(dfs, multiple = True) :
    subsystem_fluxes = {}

    if multiple : 
        for df in dfs :
            for line in df.iterrows() :
                #print(line[1]["flux"])
                if "Transport" not in line[1]["subSystem"] and "Exchange" not in line[1]["subSystem"] and len(line[1]["subSystem"]) > 1:
                    try :
                        subsystem_fluxes[line[1]["subSystem"]] += abs(line[1]["flux"])
                    except KeyError :
                        subsystem_fluxes[line[1]["subSystem"]] = abs(line[1]["flux"])
            df_bar = pd.DataFrame(subsystem_fluxes,index=["Fluxes"]).T
        return df_bar
    else :
        for line in dfs.iterrows() :
            #print(line[1]["flux"])
            if "Transport" not in line[1]["subSystem"] and "Exchange" not in line[1]["subSystem"] and len(line[1]["subSystem"]) > 1:
                try :
                    subsystem_fluxes[line[1]["subSystem"]] += abs(line[1]["flux"])
                except KeyError :
                    subsystem_fluxes[line[1]["subSystem"]] = abs(line[1]["flux"])
        df_bar = pd.DataFrame(subsystem_fluxes,index=["Fluxes"]).T
        return df_bar
    



### DATA VISUALISATION ###

def plot_treemap(df, model, title, path=['subSystem', 'id'], flux_filter=0.0, color_by = "subsystem") :
    ### Building colormap :

    
    if color_by == "subsystem" :
        # cols = ['#E48F72','#FC6955','#7E7DCD','#BC7196','#86CE00','#E3EE9E','#22FFA7','#FF0092','#C9FBE5','#B68E00','#00B5F7','#6E899C',
        # '#D626FF','#AF0038','#0D2A63','#6C4516','#DA60CA','#1616A7','#620042','#A777F1','#862A16','#778AAE','#6C7C32','#B2828D',
        # '#FC0080','#00A08B','#511CFB','#EB663B','#750D86','#B68100','#222A2A','#FB0D0D','#1CA71C','#E15F99',
        # '#2E91E5','#DC587D','#EEA6FB','#479B55','#FF9616','#F6F926','#0DF9FF','#FE00CE','#FED4C4','#6A76FC','#00FE35','#FD3216',
        # '#2E91E5','#E15F99','#1CA71C','#FB0D0D','#222A2A','#B68100','#750D86','#EB663B','#511CFB','#00A08B','#FB00D1',
        # '#FC0080','#B2828D','#6C7C32','#778AAE','#862A16','#A777F1','#620042','#1616A7','#DA60CA','#6C4516','#0D2A63','#AF0038']
        cols = ["#4D455D", "#E96479", "#F5E9CF", "#7DB9B6", "#539165", "#820000", "#FFB100"]
        # Used to plot fluxes compartment by compartment --> coloring by the next highest hierarchical category = subsystems.
    
        
        subsystems = set()
        for r in model.reactions :
            if len(r.subsystem) >0 and not "Transport" in r.subsystem and not "Exchange" in r.subsystem :
                subsystems.add(r.subsystem)

        cmap = {}
        for subsystem, color in zip(subsystems, cycle(set(cols))) :
            cmap[subsystem] = color
        ###

        df = df.loc[(df["subSystem"] != "Transport, mitochondrial")& (df["subSystem"] != "Transport, extracellular" ) & (df["subSystem"] != "Exchange reactions")& (df["subSystem"] != "Transport, peroxisomal" )]
        fig = px.treemap(df.loc[(df["flux"] >= flux_filter) & (df["name"] != "Null")].to_dict() , path=path, 
                    values='flux', hover_name= "name" ,color='subSystem',hover_data = ["direction", "reactants", "products"], color_discrete_map=cmap)

        fig.update_layout(title_text=title, font_size=12)

    elif color_by == "compartment" :
        cols = ["#4D455D", "#E96479", "#F5E9CF", "#7DB9B6", "#539165", "#820000", "#FFB100"]
        # Used to plot fluxes subsystem by subsystem --> coloring by the next highest hierarchical category = compartments.

        compartments = set()
        for r in model.reactions :
            if len(r.compartments) > 0 and not "C_r" in r.compartments and not "C_l" in r.compartments and not "C_g" in r.compartments :
                for comp in r.compartments :
                    compartments.add(comp)
        cmap =  {}
        for compartment, color in zip(compartments, cycle(set(cols))) :
            cmap[compartment] = color
        ###

        df = df.loc[(df["subSystem"] != "Transport, mitochondrial")& (df["subSystem"] != "Transport, extracellular" ) & (df["subSystem"] != "Exchange reactions")& (df["subSystem"] != "Transport, peroxisomal" )]
        fig = px.treemap(df.loc[(df["flux"] >= flux_filter) & (df["name"] != "Null")].to_dict() , path=path, 
                    values='flux', hover_name= "name" ,color='compartment', hover_data = ["direction", "reactants", "products"], color_discrete_map=cmap)
    return fig 


def compartment_fluxes_barplots(model_1, model_2) :
    barplots = {}
    for compartments_iHep, compartments_G2 in zip(build_reaction_df(model_1), build_reaction_df(model_2)):

        color_Ekeley = 'salmon'
        color_Porter = 'lightblue'
        subS_HepG2 = get_subsystem_fluxes(compartments_G2, multiple=False)
        subS_iHep = get_subsystem_fluxes(compartments_iHep, multiple=False)
        df_both = pd.concat([subS_iHep, subS_HepG2],axis=1)
        sum_HepG2 = subS_HepG2.sum(0)
        sum_iHep = subS_iHep.sum(0)
        df_both.columns = ["iHep", "HepG2"]

        df_both = df_both.loc[(df_both["iHep"] != 0.0)& (df_both["HepG2"] != 0.0)]
        
        compartment_G2 = compartments_G2["compartment"][0]
        compartment_iHep = compartments_iHep["compartment"][0]

        

        if not compartment_iHep == compartment_G2 :
            print("Erreur compartements")
            break
        if not "C_r" in compartment_iHep and not "C_s" in compartment_iHep and not "C_x" in compartment_iHep and not "C_g" in compartment_iHep and not "C_l" in compartment_iHep :

            labels = list(df_both.index)
            normalized_fluxes_iHep = [float(val/sum_iHep) for val in df_both.loc[:,"iHep"].values]
            normalized_fluxes_G2 = [float(val/sum_HepG2) for val in df_both.loc[:,"HepG2"].values]

            xmax = max(max(normalized_fluxes_G2), max(normalized_fluxes_iHep))
            xmin = min(min(normalized_fluxes_G2), min(normalized_fluxes_iHep))
            
            fig, axes = plt.subplots(figsize=(5,7), ncols=2, sharey=True)
            fig.tight_layout()


            axes[0].barh(labels, normalized_fluxes_iHep, align='center', color=color_Ekeley, zorder=10)
            axes[0].set_title(f"iHep : {compartment_iHep}")
            axes[1].barh(labels, normalized_fluxes_G2 , align='center', color=color_Porter, zorder=10)
            axes[1].set_title(f"Hep_G2 : {compartment_G2}")
            axes[0].set_xlabel("Fraction des comptages totaux")
            axes[1].set_xlabel("Fraction des comptages totaux")
            axes[0].invert_yaxis() # labels read top-to-bottom
            axes[0].invert_xaxis() # mirror data for both duildings
            plt.xlim = (xmin, xmax)
            plt.close()
            barplots[compartment_G2] = (fig, axes)
    return barplots




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
        if flux > 0.0 :
                fleche = "-->"
        elif flux < 0.0 :
                fleche = "<--"
        else:
                fleche = "<=>"
    else :
        if flux > 0.0 :
            fleche = "<--"
        elif flux < 0.0 :
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
    