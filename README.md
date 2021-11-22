# AFIG v0.1.0
Using BUSCO Lineage data to train HMM model and annotate fungi genome.  
By using taxonkit to identify the lineage of your fungi genome. This script will automatically find a busco dataset and map conserved ancestral proteins to this genome. SNAP will be used to train a species-specific HMM to annotate your genome. 
This script will also configure a bash script for further antiSMASH analysis.


### Requirements

```bash
Maker v3.01
snap
e2e
A BUSCO database (One may need to build it you self because of license issues.)
```


### Install

```bash
pip install e2e
git clone https://github.com/JinyuanSun/AFIG.git
cd AFIG && export PATH="$(pwd):$PATH"
```
Please set `path/to/snap/` and  `path/to/perl/` for the `run_snap()` function at line 30 in `genome_annotation_pipeline.py`.

### Usage
```bash
usage: genome_annotation_pipe.py [-h] [-n NAME] [-g GENOME_FILE_NAME]
                                 [-p THREAD_NUM]

Automate, RNA-seq free, fungi genome annotation pipeline.

optional arguments:
  -h, --help            show this help message and exit
  -n NAME, --name NAME
  -g GENOME_FILE_NAME, --genome_file_name GENOME_FILE_NAME
  -p THREAD_NUM, --thread_num THREAD_NUM

demo: python3 genome_annotation_pipe.py -n Apiotrichum_mycotoxinovorans -g Am_genome.fasta -p 16
```

### Change Log
**2021-11-22:**  
Replace `taxonkit` with `ETE Toolkit`.  
[ETE](http://etetoolkit.org/) is a Python programming toolkit that assists in the
 automated manipulation, analysis and visualization of phylogenetic trees.
 
### To Do List
1. Build a smaller fungal specific lineages dataset derived from BUSCO.
2. Find or make a better software to replace Maker. It is too heavy.

**Feel free to contact me or contribute to this open-source project to make it better!**
