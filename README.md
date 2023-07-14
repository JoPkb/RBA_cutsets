# RBA_cutsets


### Exécution d'ASPEFM avec contraintes positives et négatives et checker de minimalité :

chemin : aspefm_mcs_checker/mcs_with_check.sh \
Ce script permet une exécution plus ergonomique de la recherche de MCS. Il faut spécifier les contraintes à appliquer, qui sont ensuite stockées dans des fichier textes et donnés en argument à clingo ensuite.
A exécuter au sein d'un environnement clingo !

Dans aspefm_mcs_checker/data/ se trouvent les fichiers contenant le modèle au format lp4, ainsi que les contraintes.

### Traitement des résultats :

chemin aspmcs/code/*
* **compression** : Notebook exécutant la fonction de compression d'EFMtool. A utiliser avant la recherche de MCS dans un nouveau modèle.
* **mcs_parsing** : Parsing des résultats d'ASPEFM : parse_and_decomp_mcs.py permet de parser les logs pour récupérer les MCS, puis de les décompresser. Le script prend quatre arguments en entrée :
  - 1 : Chemin vers le modèle (SBML)
  - 2 : Chemin vers le fichier contenant les logs d'ASPEFM (TXT)
  - 3 : Chemin vers le fichier de reaction subsets permettant la décompression des MCS. Il se trouve dans models_storage/lp_models/{model_name}_reactionSubsets.txt
  - 4 : Chemin spécifiant l'emplacement et le nom du fichier de sortie. Le fichier de sortie contient les MCS parsés et décompressés. Ils se trouvent sous forme d'une liste python, stockée dans un fichier pickle.
* **mcs_eval** : Contient trois versions du même script. Selon les besoins, ils peuvent chacuns être utilisés et fonctionne globalement de la même façon. Le plus généraliste est *cutsets_on_healthy_params.py*. Ce script permet de tester l'effet des MCS sur les fonctions métaboliques du modèle sain. Pour son fonctionnement, il faut donner comme argument au script les éléments suivants :
  - Le chemin vers le fichier pickle contenant la liste de MCS décompressée.
  - Le chemin vers le modèle métabolique au format SBML sur lequel l'effet des MCS doit être vérifié.
  - Le chemin vers un fichier indiquant quelles fonctions métaboliques du modèle sont à tester, ainsi que les contraintes additionnelles sur les bornes de réactions à ajouter pour permettre le test.
  - Le chemin vers le dossier d'output où deux fichier pickle seront créés par le script.
  
  Le fichier spécifiant les fonctions métaboliques à vérifier doit avoir la forme suivante: 
```
{"biomass_components" : None, \
"HMR_9728" : [("HMR_4381", (0.0,0.0))]}
```
Il sera interprété comme un dictionnaire python par le script. Les clés sont les ID des réactions du modèle à tester. Les valeurs sont les contraintes additionnelles à ajouter pour tester la fonction métabolique. Par exemple, la réaction HMR_9728 n'est testable qu'en coupant la réaction HMR_4381.

Le fichier *aspmcs/code/mcs_eval/metabolic_functions.txt* contient les ID de réactions et les contraintes additionnelles nécessaires pour tester la fonctionalité des voies métaboliques décrites dans mon rapport de stage.

Le script retourne deux variables, stockées dans deux fichiers Pickle. Le premier (*{path}/invalid_cs.p*) contient les indices des MCS non valides dans la liste originale. Le second (*{path}/flux2cutsets.p*) contient, pour chaque MCS, un dictionnaire stockant pour chaque réaction, son flux lorsque le MCS est appliqué au modèle.
