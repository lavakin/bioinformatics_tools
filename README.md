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
>>> msa.get_sum_of_pairs(blosum)
1806364.0
>>> msa.get_sum_of_pairs_column(74,blosum)
15741.0
```






