# `baktfold` Method

## ProstT5 3Di Inference

* `baktfold` begins by predicting the Foldseek 3Di tokens for every **hypothetical** protein input using the [ProstT5](https://github.com/mheinzinger/ProstT5) protein language model
* Alternatively, this step is skipped if you specify pre-computed protein structures in the .pdb or .cif formats using the `--structures` flag along with specifying the directory containing the structures with `--structure_dir`
* You can also annotate all proteins using `-a`

##  Foldseek Structural Comparison

* `baktfold` then creates a [Foldseek](https://github.com/steineggerlab/foldseek) database combining the AA and 3Di representations of each protein, and compares this to each of the four `baktfold` databases with Foldseek (detailed below)
* Alternatively, you can specify protein structures that you have pre-computed for your phage(s) instead of using ProstT5 with the parameter `--structure_dir`
* You can also specify a custom Foldseek database with `--custom-db`

## `baktfold` databases

1. Swiss-Prot (n=590,183)
2. AFDB Clusters (v6) (n=3,085,778)
3. PDB (n=294,848)
4. CATH (n=195,223)

## Downstream Annotation Processing

* The top hit for each database lower than the E-value threshold (default E=1e-03) is then kept for each of the four searches
* For each protein the following logic is conducted to select the annotation from the top hit(s) from the 4 database searches:
* The top hit for each database lower than the E-value threshold (default
E=1e-03) is then kept for each of the four searches. 
* If Baktfold finds hits for multiple databases, the overall functional label will be taken:

1. From the user’s custom database (if specified)
2. Swiss-Prot
3. AFDB if no Swiss-Prot hits were found,
4. PDB if neither Swiss-Prot nor AFDB hits were found and 
5. finally CATH if no other databases were hit. 

In practice, 4. and 5. above are extremely rare.

All top hit accessions are available in `baktfold.inference.tsv` while the detailed foldseek statistics are available for each database in the respective `_tophit.tsv` files in Baktfold’s output.




