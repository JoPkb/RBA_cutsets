#!/bin/bash
aspmcs=(clingo clingoLP.py mcs\[LP\].lp4) # clingoLP.py with extensions enabled
aspmcs="${aspmcs[@]}"

mcschecker=mcs_checker.py

#### MODEL DATA

lp4Dual=../HCC_compressed_2/new_reduced_mcs.lp4
xml=../HCC_compressed_2/new_reduced.xml

#### ACTUAL CONSTRS
treatments_folder=$1
cstr=data/model.mcs.constr.lp4
target=data/model.mcs.target.lp4
tspLimit=data/model.mcs.transporterLimit.lp4
#treatment=data/model.mcs.treatment.lp4

### PARAMS including epsilon

params=(-c nstrict=0 -c accuracy=10 --heuristic Domain --enum-mode domRec --stats=2 -c epsilon="(1,3)" --verbose=3)
params="${params[@]}"

### WRITE CONSTRAINTS

echo ":- 4 {transporter(R) : cutset(R)}." > $tspLimit # nb tsp max: 3
echo ":- 6 {cutset(R) : reaction(R)}." > $cstr # taille max: 5
echo ':- not support("mcs_rsub_36").' > $target
#echo ':- not support("mcs_rsub_844").' > $treatment
#echo 'target("mcs_rsub_844").' >> $treatment # specify reaction as wanted for mcs_checker

### RUN CODE
i=0
for treatment_file in $(ls $treatments_folder)
do
    i=$((i+1))
    treatment="$treatments_folder"/"$treatment_file"
    echo $i
    $aspmcs $mcschecker $lp4Dual $params $target $cstr $treatment -c nb=100000 -n 0 -c mcscheckfile=\"$xml\" --time-limit=43200  > res/output_${i}.txt &  

done

# IFS="_"
# read -ra name_arr <<< $filename
# #treatment_file=${treatment_folder}${treatment}
# echo $name_arr
# $aspmcs $mcschecker $lp4Dual $params $target $cstr $treatment -c nb=100000 -n 0 -c mcscheckfile=\"$xml\"  > data/posconstr/res/output_${name_arr}.txt &  

#$aspmcs $mcschecker $lp4Dual $params $tspLimit $treatment $target $cstr -n 0 -c mcscheckfile=\"$xml\"  > output_mcs_statin.txt  

# --time-limit=129600 # 1 jour et demi

# (sortie err) 2> output_errors.mcs.txt # inutile à moins que tu cherches à débugger avec -c trace=1 (déconseillé)





