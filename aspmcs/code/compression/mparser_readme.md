==== Example usage of `mparser_cli.py` : ====

`python mparser_cli.py METATOOL ./data/sbml/toy/toy_model_biofuel.txt ASP ./data/sbml/toy/toy_model_biofuel.lp4`

`python mparser_cli.py SBML ./data/sbml/toy/toy_model_biofuel.xml ASP ./data/sbml/toy/toy_model_biofuel.lp4`

`python mparser_cli.py METATOOL ./data/sbml/toy/toy_model_biofuel.txt --to-dual-mcs ASP ./data/sbml/toy/toy_model_biofuel.mcs.lp4 --target-reactions R17`

==== Example usage of `mparser_gui.py` : ====

`python mparser_gui.py`

Needs installation of Python module `pysimplegui`.

==== Additional configuration ====

When reading SBML files and converting them to other formats, make sure to check the code in file `reading_parameters.py`.
