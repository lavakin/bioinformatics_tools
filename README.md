# Initialization
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
The editing distance function returns the distance along with all best scoring alignments
```python
>>> from bioinf.distance import *
>>> editing_distance('parek','uprk')
[['parek', 'upr-k', '3'], ['-parek', 'up-r-k', '3']]
```
### Hamming distance
```python
>>> from bioinf.distance import *
>>> hamming_distance('brok','brak')
1
```
The hamming distance is not defined for sequence with different lengths. Thus an error occurs in such case.
```python
>>> hamming_distance('brok','braky')
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
  File ".../bioinf/distance.py", line 40, in hamming_distance
    raise SequencesNotTheSameLength()
bioinf.distance.SequencesNotTheSameLength: Sequences does not have the same length
```
## Processing fasta files
Import file:
```python
>>> from bioinf.fasta import *
>>> fasta = Fasta('./bioinf/fasta_files/1R8Y.fasta')
```
Get ids:
```python
>>> fasta.get_ids()
['pdb|1R8Y|H', 'P01013']
```
Get sequence by id or index:
```python
>>> seq = fasta.get_sequence('P01013')
>>> seq = fasta.sequences[1]
```
Working with sequences:
```python
>>> len(seq)
232
>>> seq[1:10]
'IKDLLVSSS'
>>> seq.description
'GENE X PROTEIN (OVALBUMIN-RELATED)'
>>> seq.id
'P01013'
```
## Processing multiple sequence alignment
Parse MSA
```python
>>> from bioinf.MSA import *
>>> msa = MSA('bioinf/clust')
```
Get sequence and column:
```python
>>> msa.get_sequence('UniRef90_UPI000')
SeqRecord(seq=Seq('--------------------------------------------------MEPM...DSE'), id='UniRef90_UPI000', name='<unknown name>', description='UniRef90_UPI000', dbxrefs=[])
>>> msa.get_sequence(2)
SeqRecord(seq=Seq('--------------------------------------------------MEPM...DSE'), id='UniRef90_UPI000', name='<unknown name>', description='UniRef90_UPI000',dbxrefs=[])
>>> msa.get_column(50)
'MMMMMMMMMMMMMMMMMMM-MMMMMMMMMMMMMMMMMMMMMMMMMMM-MMMMMMM'
```
Sum of Pairs:
```python
>>> blosum = Align.substitution_matrices.load('BLOSUM62')
>>> msa.get_sum_of_pairs(blosum)
1806364.0
>>> msa.get_sum_of_pairs_column(74,blosum)
15741.0
```
Scoring matrix as an argument.

## Conservation determination from multiple aligned sequences

Get conservation score:
```python
>>> msa.get_conserv_for_seq(50, blosum)
68118.0
```
Get N best positions along with the respective scores:
```python
>>> msa.get_N_best_for_sequence(10,5,blosum)
[(283, 432.0), (261, 432.0), (224, 475.0), (370, 486.0), (368, 486.0), (332, 486.0), (328, 486.0), (259, 486.0), (218, 486.0), (74, 583.0)]
>>> msa.get_N_best(10,blosum)
[(283, 11880.0), (261, 11880.0), (224, 12771.0), (370, 13365.0), (368, 13365.0), (332, 13365.0), (328, 13365.0), (259, 13365.0), (218, 13365.0), (74, 15741.0)]
```
Get conservation score for a sequence
```python
>>> msa.get_conserv_for_seq(50,blosum)
68118.0
```
## Processing PDB files
Get substructures:
```python
>>> from bioinf.pdb_parser import *
>>> pdb = PDB('./pdbs/1b0b.pdb')
>>> pdb.structure[0]
<Model id=0>
>>> pdb.structure[0]["A"]
<Chain id=A>
>>> pdb.structure[0]["A"][100]
<Residue GLY het=  resseq=100 icode= >
>>> pdb.structure[0]["A"][100]["CA"]
<Atom CA>
```
Get information about the structure:
```python
>>> pdb.get_info()
{'number of models:': 1, 'number of chains': 1, 'number of residues': 346, 'number of atoms': 1289}
>>> pdb.get_diameter()
52.42697045264779
```
Get nearest neighbors:
```python
>>> pdb.get_nearest_residues(pdb.structure[0]['A'][111]['CA'],3)
[<Residue VAL het=  resseq=112 icode= >, <Residue PHE het=  resseq=110 icode= >, <Residue LYS het=  resseq=111 icode= >]
>>> pdb.get_nearest_atoms(pdb.structure[0]['A'][111]['CA'],2)
[<Atom C>, <Atom N>, <Atom CA>, <Atom CB>]
```
## Computing structure-related properties
Ratio of surface and buried amino acids
```python
>>> pdb.get_buried_ratio()
{'buried_ratio': 0.241, 'exposed_ratio': 0.759}
>>> pdb.get_buried_exposed_by_amk()
{'buried': Counter({'ALA': 6, 'LEU': 5, 'PHE': 4, 'TRP': 3, 'MET': 3, 'LYS': 2, 'VAL': 2, 'GLY': 2, 'ASN': 1, 'SER': 1, 'PRO': 1, 'CYS': 1, 'GLU': 1, 'THR': 1, 'ILE': 1}), 'exposed': Counter({'ALA': 22, 'GLY': 14, 'SER': 9, 'LYS': 
>>> pdb.get_buried_exposed_ratio_by_amk()
{'LYS': 0.182, 'ARG': 0, 'ILE': 0.5, 'VAL': 0.286, 'THR': 0.2, 'HIS': 0, 'TRP': 0.75, 'ASP': 0, 'GLU': 0.2, 'MET': 0.5, 'PRO': 0.333, 'SER': 0.1, 'PHE': 0.4, 'ASN': 0.167, 'ALA': 0.214, 'TYR': 0, 'CYS': 1, 'LEU': 0.556, 'GLY': 0.125, 'GLN': 0}
9, 'ASP': 8, 'PHE': 6, 'GLN': 5, 'VAL': 5, 'ASN': 5, 'LEU': 4, 'THR': 4, 'GLU': 4, 'MET': 3, 'HIS': 2, 'PRO': 2, 'ARG': 2, 'TRP': 1, 'ILE': 1, 'TYR': 1})}
```
Get histogram:
```python
>>> import matplotlib.pyplot as plt
>>> plt = pdb.get_histogram()
>>> plt.show()
```
Get ratio of polar aminoacids in exposed and buried:
```python
>>> pdb.get_polarity()
{'buried ratio': 0.20588235294117646, 'exposed ratio': 0.45794392523364486}
```
 ## A2a and caffeine receptor and hemoglobin (1b0b) differences
 The caffeine receptor has 35% of its aminoacids buried. From the buried aminoacids, 32% are polar and from the exposed 41% are polar. The number of atoms is 2410.
 
 The hemoglobin has 24% of its aminoacids buried. From the buried aminoacids, 21% are polar and from the exposed 46% are polar. The number of atoms is 1289.
 
 The hemoglobin molecule is thus smaller and less compact, whereas the caffeine receptor is larger and more compact. This is why caffeine rceptor has larger ratio of buried aminoacids. This is also why polar aminoacids are more often buried by caffeine receptor then by hemoglobin.
 




