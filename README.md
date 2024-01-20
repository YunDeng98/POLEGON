![Logo](POLEGON.png)
# POLEGON
POLEGON stands for **P**rior-**O**blivious **L**ength **E**stimation in **G**enealogies with **O**riented **N**etwork. POLEGON works with inferred Ancestral Recombination Graph (ARG) to re-calibrate the branch length, **without the usage of any prior**. The inferred ARGs should be in tskit format, with mutations mapped to branches. It is also important that the genealogies in the ARG should be **linked**, in that the adjacent trees should differ relatively small. After the branch length has been inferred, the population size history can be subsequently inferred with the calibrated branch length.

The details of the algorithms can be found at: [add some link], which is also the citation source. 

# Input and output
POLEGON takes .trees files with tskit tree sequence format (of course there need to be mutations in it!). Its output is also in tree sequence format, for both posterior samples and posterior average. 

# Basic usage
Fixated on the topology, POLEGON can generate you the posterior samples of the ARG and the posterior average of them. 

The basic commands is:

```
polegon_master -m mutation_rate -input original_tree_sequence -output updated_tree_sequence -num_samples N -thin K -scaling_rep L
```

The following details to these arguments can be displayed if you simply type `polegon`

|flag|required?|details|  
|-------------------|-----|---|  
|**-input**|required|the prefix of the tree sequence file|
|**-output**|required|the prefix of the re-sampled tree sequence file|
|**-m**|conditionally required|per base pair per generation mutation rate|
|**-map**|conditionally required|mutation rate map for the region|
|**-num_samples**|optional|the number of posterior ARG samples. Default: 100|
|**-thinning**|optional|the number of thinning iterations in MCMC. Default: 10|
|**-scaling_rep**|optional|the number of rescaling steps after MCMC. Default: 5|

If you want to use a mutation map, rather than a constant mutation rate along the genome, the mutation map file should be formatted as follows:

```
0 1.2e-8
100000 2e-8
200000 1e-8
```

this means that the mutation rate between 0-100kb is 1.2e-8, and between 100-200kb is 2e-8. The coordinates must start from 0 and the last coordinate must be larger than (or equal to) the sequence length in the tree sequence file, so that mutation map is fully defined. 

# Suggestions from the developers
1. **TL;DR: Only SINGER, fast-SINGER and tsinfer work with POLEGON.**
Unfortunately, POLEGON doesn't work with every current ARG inference methods. First it needs to output in **tree sequence format**, so methods like ARGweaver wouldn't work here; Second of course mutations have to be **mapped** to the branches, so ARG-Needle wouldn't work here; Thirdly the tree sequence needs to have some **spatial regularity** (it should be!), that is, adjacent trees differ relatively small, so Relate's output is not ideal (you should use Relate's native branch length estimator).
2. 
