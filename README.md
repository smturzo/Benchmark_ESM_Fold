# Benchmark\_ESM\_Fold

Benchmark ESMFold predictions by RMSD and TMScore!  

## Example Usage

1. On terminal type the following:

   python benchmark\_esmfold\.py  -s /path\_to\_/native\_structure\.pdb -o /path\_to\_/store\_output/

2. Result is printed on screen as:
   ./path\_to\_/native\_structure.pdb seq\_length rmsd tmscore plddt

3. Output prediction is stored as: 
   native\_structure\_esm\_pred\.pdb

## Reference

```bibtex
@article{lin2022language,
  title={Language models of protein sequences at the scale of evolution enable accurate structure prediction},
  author={Lin, Zeming and Akin, Halil and Rao, Roshan and Hie, Brian and Zhu, Zhongkai and Lu, Wenting and dos Santos Costa, Allan and Fazel-Zarandi, Maryam and Sercu, Tom and Candido, Sal and others},
  journal={bioRxiv},
  year={2022},
  publisher={Cold Spring Harbor Laboratory}
}
@article{https://doi.org/10.1002/prot.20264,
author = {Zhang, Yang and Skolnick, Jeffrey},
title = {Scoring function for automated assessment of protein structure template quality},
journal = {Proteins: Structure, Function, and Bioinformatics},
volume = {57},
number = {4},
pages = {702-710},
doi = {https://doi.org/10.1002/prot.20264},
url = {https://onlinelibrary.wiley.com/doi/abs/10.1002/prot.20264},
eprint = {https://onlinelibrary.wiley.com/doi/pdf/10.1002/prot.20264},
abstract = {Abstract We have developed a new scoring function, the template modeling score (TM-score), to assess the quality of protein structure templates and predicted full-length models by extending the approaches used in Global Distance Test (GDT)1 and MaxSub.2 First, a protein size-dependent scale is exploited to eliminate the inherent protein size dependence of the previous scores and appropriately account for random protein structure pairs. Second, rather than setting specific distance cutoffs and calculating only the fractions with errors below the cutoff, all residue pairs in alignment/modeling are evaluated in the proposed score. For comparison of various scoring functions, we have constructed a large-scale benchmark set of structure templates for 1489 small to medium size proteins using the threading program PROSPECTOR\_3 and built the full-length models using MODELLER and TASSER. The TM-score of the initial threading alignments, compared to the GDT and MaxSub scoring functions, shows a much stronger correlation to the quality of the final full-length models. The TM-score is further exploited as an assessment of all ‘new fold’ targets in the recent CASP5 experiment and shows a close coincidence with the results of human-expert visual assessment. These data suggest that the TM-score is a useful complement to the fully automated assessment of protein structure predictions. The executable program of TM-score is freely downloadable at http://bioinformatics.buffalo.edu/TM-score. Proteins 2004. © 2004 Wiley-Liss, Inc.},
year = {2004}
}

```
