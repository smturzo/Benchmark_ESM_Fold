# Benchmark\_ESM\_Fold

Benchmark ESMFold predictions by RMSD and TMScore!  
Currently this script is limited to predictions with sequence length < 400 because of the [API](https://esmatlas.com/about#api) limitations.

TODO: Work around the 400 seq lenght limitation.  
TODO: Add argparse for pdb and path to save
TODO: Add documentation for usage
## Example Usage

1. On terminal type the following:


2. Result is printed on screen!


## Reference

```bibtex
@article{lin2022language,
  title={Language models of protein sequences at the scale of evolution enable accurate structure prediction},
  author={Lin, Zeming and Akin, Halil and Rao, Roshan and Hie, Brian and Zhu, Zhongkai and Lu, Wenting and dos Santos Costa, Allan and Fazel-Zarandi, Maryam and Sercu, Tom and Candido, Sal and others},
  journal={bioRxiv},
  year={2022},
  publisher={Cold Spring Harbor Laboratory}
}

```
