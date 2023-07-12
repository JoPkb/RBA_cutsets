# RBA_cutsets


### Exécution d'ASPEFM avec contraintes positives et négatives et checker de minimalité :

chemin : aspefm_mcs_checker/mcs_with_check.sh \
Ce script permet une exécution plus ergonomique de la recherche de MCS. Il faut spécifier les contraintes à appliquer, qui sont ensuite stockées dans des fichier textes et donnés en argument à clingo ensuite.
A exécuter au sein d'un environnement clingo !

Dans aspefm_mcs_checker/data/ se trouvent les fichiers contenant le modèle au format lp4, ainsi que les contraintes.

### Traitement des résultats :

chemin aspmcs/code/*
* compression : Notebook exécutant la fonction de compression d'EFMtool. A utiliser avant la recherche de MCS dans un nouveau modèle.
* mcs_parsing : Parsing des résultats d'ASPEFM : parse_and_decomp_mcs.py permet de parser les logs pour récupérer les MCS, puis de les décompresser. Le script prend quatre arguments en entrée :
  - 1 : Chemin vers le modèle (SBML)
  - 2 : Chemin vers le fichier contenant les logs d'ASPEFM (TXT)
  - 3 : Chemin vers le fichier de reaction subsets permettant la décompression des MCS. Il se trouve dans models_storage/lp_models/{model_name}_reactionSubsets.txt
  - 4 : Chemin spécifiant l'emplacement et le nom du fichier de sortie. Le fichier de sortie contient les MCS parsés et décompressés. Ils se trouvent sous forme d'une liste python, stockée dans un fichier pickle.
* mcs_eval : 
