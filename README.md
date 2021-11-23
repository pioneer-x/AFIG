# AFIG v0.1.0
Using BUSCO Lineage data to train HMM model and annotate fungi genome.  
By using e2e to identify the lineage of your fungi genome. This script will automatically find a busco dataset 
and map conserved ancestral proteins to this genome. SNAP will be used to train a species-specific HMM to annotate your genome. 
This pipeline will also configure a bash script for further antiSMASH analysis.


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
cd AFIG && tar vxzf fungi_busco_sub.tar.gz && export PATH="$(pwd):$PATH"
```


### Usage
```bash
usage: genome_annotation_pipe.py [-h] [-n NAME] [-g GENOME_FILE_NAME] [-e ENGINE] [-p THREAD_NUM]

Automate, RNA-seq free, fungi genome annotation pipeline.

optional arguments:
  -h, --help            show this help message and exit
  -n NAME, --name NAME
  -g GENOME_FILE_NAME, --genome_file_name GENOME_FILE_NAME
  -e ENGINE, --engine ENGINE
                        use maker or exonerate
  -p THREAD_NUM, --thread_num THREAD_NUM

demo: genome_annotation_pipe.py -n Apiotrichum_mycotoxinovorans -g Am_genome.fasta -p 16 -e maker
```

### Change Log
**2021-11-23:**  
Built a sub-database of released official [BUSCO databse](https://busco-data.ezlab.org/v4/data/lineages/). Only fungi are
included.

**2021-11-22:**  
Replace `taxonkit` with `ETE Toolkit`.  
[ETE](http://etetoolkit.org/) is a Python programming toolkit that assists in the
 automated manipulation, analysis and visualization of phylogenetic trees.
 
### To Do List
1. Build a smaller fungal specific lineages dataset derived from BUSCO.
2. Find or make a better software to replace Maker. It is too heavy.

### For further functional annotation and visulation, refer to [CMPP](https://github.com/JinyuanSun/CMPP).

**Feel free to contact me or contribute to this open-source project to make it better!**
