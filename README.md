# Initialisation
1. Clone this repository: `git clone git@github.com:lavakin/bioinformatics_tools.git`
2. move to the repository `cd bioinformatics_tools`
3. set up a virtual environment `python3 -m venv bioinf`
4. activate the environment `. bioinf/bin/activate`
5. build and install package `./setup.py build && ./setup.py install`

After you are done, just type `deactivate` to terminate the virtual environment

# Use
## Distance 
The package implements two distance measures. 

### Editing distance
The function returns the distance along with all best scoring alignments
```
>>> from bioinf.distance import *
>>> editing_distance('parek','uprk')
[['parek', 'upr-k', '3'], ['-parek', 'up-r-k', '3']]
```
