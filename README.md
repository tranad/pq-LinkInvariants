# pq-LinkInvariants
Compute arbitrary link invariants up to three strands for <a href="http://www.codecogs.com/eqnedit.php?latex=\mathcal&space;Z(\textrm{Vec}_G^{\omega^u})" target="_blank"><img src="http://latex.codecogs.com/gif.latex?\mathcal&space;Z(\textrm{Vec}_G^{\omega^u})" title="\mathcal Z(\textrm{Vec}_G^{\omega^u})" /></a> where <a href="https://www.codecogs.com/eqnedit.php?latex=G=\mathbb&space;Z_{11}&space;\rtimes&space;\mathbb&space;Z_5" target="_blank"><img src="https://latex.codecogs.com/gif.latex?G=\mathbb&space;Z_{11}&space;\rtimes&space;\mathbb&space;Z_5" title="G=\mathbb Z_{11} \rtimes \mathbb Z_5" /></a>.

[(https://arxiv.org/abs/1805.05736)](https://arxiv.org/abs/1805.05736). 

Mathematica packages containing data on <a href="https://www.codecogs.com/eqnedit.php?latex=G=\mathbb&space;Z_{11}&space;\rtimes&space;\mathbb&space;Z_5" target="_blank"><img src="https://latex.codecogs.com/gif.latex?G=\mathbb&space;Z_{11}&space;\rtimes&space;\mathbb&space;Z_5" title="G=\mathbb Z_{11} \rtimes \mathbb Z_5" /></a> and <a href="http://www.codecogs.com/eqnedit.php?latex=\mathcal&space;Z(\textrm{Vec}_G^{\omega^u})" target="_blank"><img src="http://latex.codecogs.com/gif.latex?\mathcal&space;Z(\textrm{Vec}_G^{\omega^u})" title="\mathcal Z(\textrm{Vec}_G^{\omega^u})" /></a> are included, as well as a sample notebook that computes the S and T matrices and checks that it satisfies the modular relations.

```python
SetDirectory[NotebookDirectory[]];
<< "Z11xZ5-data.m"
<< "DZ11xZ5.m"
DistributeDefinitions[ComputeLinkInvariant];

u = 1;
braidword = {1};
numStrands = 2;

dT1 = ComputeLinkInvariant[u, {braidword, numStrands}];
```

The above computes a twisted unknot for the u=1 category. To compute the braid closure <a href="http://www.codecogs.com/eqnedit.php?latex=\widehat{\sigma_2^{-3}&space;\sigma_1&space;\sigma_2^{2}}" target="_blank"><img src="http://latex.codecogs.com/gif.latex?\widehat{\sigma_2^{-3}&space;\sigma_1&space;\sigma_2^{2}}" title="\widehat{\sigma_2^{-3} \sigma_1 \sigma_2^{2}}" /></a> on three strands for u=4 do:
```python
ComputeLinkInvariant[4, {{2,2,1,-2,-2,-2}, 3}]
```
