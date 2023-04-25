# coding: utf-8
txtfile='new_reactionSubsets.txt'
with open(txtfile, 'r') as f:
    dicto = eval(f.read())
    
import cobra
cobra.solver
cobra.util.solvers
txtfile='new_reactionSubsets.txt'
with open(txtfile, 'r') as f:
    dicto = eval(f.read())
    
dicto
l = ['rsub_459']
cobra.io.read_sbml_model('HepG2_medium.xml')
model = cobra.io.read_sbml_model('HepG2_medium.xml')
l
cs = dicto[l[0]]['reacs']]
cs = dicto[l[0]]['reacs']
cs
csets = []
cs
l2 = ['rsub_1289']
cs2 = dicto[l2[0]]['reacs']
cs2
csets
csets.append(cs); csets.append(cs2)
l3 = ['rsub_459', 'rsub_1289', 'rsub_1200']
tr = dicto[l[2]]['reacs']
tr = dicto[l3[2]]['reacs']
tr
cs3 = ['HMR_0432', 'HMR_5432', 'HMR_3875']
csets.append(cs3)
csets
from itertools import combinations
model
for cs in csets:
    with model as m:
        for reac in cs:
            m.reactions.get_by_id(reac).knock_out()
        opt = model.optimize()
        
for cs in csets:
    with model as m:
        for reac in cs:
            m.reactions.get_by_id(reac).knock_out()
        opt = model.optimize()
    print(opt, cs)
    if len(cs) > 1:
        for nm1 in combinations(cs, len(cs)-1):
            with model as m:
                for reac in cs:
                    m.reactions.get_by_id(reac).knock_out()
                opt = model.optimize()
            print('N-1', opt, cs)
            
for cs in csets:
    with model as m:
        for reac in cs:
            m.reactions.get_by_id(reac).knock_out()
        opt = model.optimize()
    print(opt, cs)
    if len(cs) > 1:
        for nm1 in combinations(cs, len(cs)-1):
            with model as m:
                for reac in nm1:
                    m.reactions.get_by_id(reac).knock_out()
                opt = model.optimize()
            print('N-1', opt, nm1)
            
for cs in csets:
    with model as m:
        for reac in cs:
            m.reactions.get_by_id(reac).knock_out()
        opt = model.optimize()
    print(opt, cs)
    if len(cs) > 1:
        for nm1 in combinations(cs, len(cs)-1):
            with model as m:
                for reac in nm1:
                    m.reactions.get_by_id(reac).knock_out()
                opt = model.optimize()
            print('N-1', opt, nm1, 'not MIN cutset' if opt == 0.0)
for cs in csets:
    with model as m:
        for reac in cs:
            m.reactions.get_by_id(reac).knock_out()
        opt = model.optimize()
    print(opt, cs)
    if len(cs) > 1:
        for nm1 in combinations(cs, len(cs)-1):
            with model as m:
                for reac in nm1:
                    m.reactions.get_by_id(reac).knock_out()
                opt = model.optimize()
            print('N-1', opt, nm1, 'not MIN cutset' if opt == 0.0 else '')
            
for cs in csets:
    with model as m:
        for reac in cs:
            m.reactions.get_by_id(reac).knock_out()
        opt = model.optimize()
    print(opt, cs)
    if len(cs) > 1:
        for nm1 in combinations(cs, len(cs)-1):
            with model as m:
                for reac in nm1:
                    m.reactions.get_by_id(reac).knock_out()
                opt = model.optimize()
            print('N-1', opt, nm1, 'not MIN cutset' if opt < 1e-6 else '')
            
for cs in csets:
    with model as m:
        for reac in cs:
            m.reactions.get_by_id(reac).knock_out()
        opt = model.optimize()
    print(opt, cs)
    if len(cs) > 1:
        for nm1 in combinations(cs, len(cs)-1):
            with model as m:
                for reac in nm1:
                    m.reactions.get_by_id(reac).knock_out()
                opt = model.optimize()
            print('N-1', opt, nm1, 'not MIN cutset' if opt.objective_value < 1e-6 else '')
            
cs2 = dicto[l2[0]]['reacs']
cs2
l2[0]
r.knock_out => 'rev' 'backwards' r.lower_bound = 0,  '(pas rev: default)', r.upper_bound = 0
for cs in csets:
    with model as m:
        for reac in cs:
            m.reactions.get_by_id(reac).knock_out()
        opt = m.optimize()
    print(opt, cs)
    if len(cs) > 1:
        for nm1 in combinations(cs, len(cs)-1):
            with model as m:
                for reac in nm1:
                    m.reactions.get_by_id(reac).knock_out()
                opt = m.optimize()
            print('N-1', opt, nm1, 'not MIN cutset' if opt.objective_value < 1e-6 else '')
            
model.reactions.get_by_id('HMR_5432').bounds
model.reactions.get_by_id('HMR_3875').bounds
csets.append(['HMR_0432', 'HMR_5432'])
for cs in csets:
    with model as m:
        for reac in cs:
            m.reactions.get_by_id(reac).knock_out()
        opt = m.optimize()
    print(opt, cs)
    if len(cs) > 1:
        for nm1 in combinations(cs, len(cs)-1):
            with model as m:
                for reac in nm1:
                    m.reactions.get_by_id(reac).knock_out()
                opt = m.optimize()
            print('N-1', opt, nm1, 'not MIN cutset' if opt.objective_value < 1e-6 else '')
            
cobra.objective
model.objective
model.objective_direction
print(model.objective)
dicto['rsub_1696']['reacs']
dicto['rsub_1159']['reacs']
csets.append(['HMR_0351'], ['HMR_5856'])
csets.append(['HMR_0351']); csets.append(['HMR_5856'])
csets
for cs in csets:
    with model as m:
        for reac in cs:
            m.reactions.get_by_id(reac).knock_out()
        opt = m.optimize()
    print(opt, cs)
    if len(cs) > 1:
        for nm1 in combinations(cs, len(cs)-1):
            with model as m:
                for reac in nm1:
                    m.reactions.get_by_id(reac).knock_out()
                opt = m.optimize()
            print('N-1', opt, nm1, 'not MIN cutset' if opt.objective_value < 1e-6 else '')
            
model.reactions.get_by_id('biomass_components').metabolites
