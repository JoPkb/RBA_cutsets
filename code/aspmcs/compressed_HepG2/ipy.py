# coding: utf-8
with open(txtfile, 'r') as f:
    dicto = eval(f.read())
    
txtfile='new_reactionSubsets.txt'
with open(txtfile, 'r') as f:
    dicto = eval(f.read())
    
dicto
for rsub, v in dicto.items():
    print(rsub, v['reacs'])
    
get_ipython().run_line_magic('log', '')
get_ipython().run_line_magic('save', '')
get_ipython().run_line_magic('save', 'ipy.py')
python parse_output.py --split-revs --csv output_mcs.csv
import pandas
csv = pandas.read_csv('output_mcs.csv')
python parse_output.py output_mcs.txt --split-revs --csv output_mcs.csv
for i, row in csv.iterrows
['mcs_rsub_1', 'mcs_rsub_2']
l = ['mcs_rsub_1', 'mcs_rsub_2']
l = list(map(lambda x: x[4:], l))
l
for i, row in csv.iterrows():
    row.loc[row != 0]
    
for i, row in csv.iterrows():
    row.loc[row != 0].index
    
l = ['mcs_rsub_1', 'mcs_rsub_2', 'mcs_M01xx', 'mcs_rsub_2616_tgt']
l = list(filter(lambda x: 'tgt' not in x and 'rsub' in x, l))
l
l = list(map(lambda x: x[4:], l))
l
decomp_reacs = []
for rs in l: # l = support de solution
    decomp_reacs.extend(dicto[rs]['reacs'])
    
decomp_reacs
decomp_cutsets = []
for rs in l: # l = support de solution
    cs = []
    for r in (dicto[rs]['reacs']):
        cs.append(r)
        
cs
drs1 = dicto['rsub_1']['reacs']
drs2 = dicto['rsub_2']['reacs']
from itertools import combinations
combinations(drs1, drs2, 2)
list(combinations([drs1, drs2], 2))
drs1
drs2
list(combinations((drs1, drs2), 2))
list(combinations((drs1), 2))
list(combinations(zip(drs1, drs2), 2))
l = ['rsub_459']
import cobra
