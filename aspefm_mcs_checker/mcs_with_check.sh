
aspmcs=(clingo clingoLP.py mcs\[LP\].lp4) # clingoLP.py with extensions enabled
aspmcs="${aspmcs[@]}"

mcschecker=mcs_checker.py

#### MODEL DATA

lp4Dual=data/model.mcs.lp4
xml=data/model.xml

#### ACTUAL CONSTRS

cstr=data/model.mcs.constr.lp4
target=data/model.mcs.target.lp4
tspLimit=data/model.mcs.transporterLimit.lp4
treatment=data/model.mcs.treatment.lp4

### PARAMS including epsilon

params=(-c nstrict=0 -c accuracy=10 --heuristic Domain --enum-mode domRec --stats=2 -c epsilon="(1,3)" --verbose=3)
params="${params[@]}"

### WRITE CONSTRAINTS

echo ":- 4 {transporter(R) : cutset(R)}." > $tspLimit # nb tsp max: 3
echo ":- 6 {cutset(R) : reaction(R)}." > $cstr # taille max: 5
echo ':- not support("mcs_rsub_93_tgt").' > $target
echo ':- not support("mcs_rsub_844").' > $treatment
echo 'target("mcs_rsub_844").' >> $treatment # specify reaction as wanted for mcs_checker

### RUN CODE

$aspmcs $mcschecker $lp4Dual $params $tspLimit $treatment $target $cstr -n 0 -c nb=100000 -c mcscheckfile=\"$xml\"  > output_mcs_statin.txt  

# --time-limit=129600
# (sortie err) 2> output_errors.mcs.txt





