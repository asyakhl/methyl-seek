
#### Setting up the working environment

Before running the pipeline, certain packages are required to be installed within a custom conda environment.

```
module load python
source /data/$USER/conda/etc/profile.d/conda.sh
conda create --name meth
conda activate meth
mamba install -yc bioconda bwameth methyldackel
conda deactivate meth
```

#### Setting up the working files

Alter the config.yaml so the rawdata_dir is the absolute path of the directory containing all your fastqs.
Alter the result_dir so it is the absolute path of the working directory containing your snakemake pipeline, where results will be stored.
Alter samples in config.yaml to be the absolute path of your samples.txt. Check this is correct. The samples file should have the following tab-delimited format:

```
sample  group comp
S1  GA  GAvsGB
S2  GA  GAvsGB
S3  GB  GAvsGB
S4  GB  GAvsGB
S5  GC  GAvsGC
S6  GC  GAvsGC
S1  GA  GAvsGC
S2  GA  GAvsGC
```

Where GA, GB & GC are the groups these 6 samples belong to.


###### Need combine multiple Snakefiles into one while still keeping four steps different

The pipeline is divided into 4 steps:

 * bismark - which performs quality control steps and maps the data using Bismark, before extracting CpG profiles using MethylDackel.
 * bwa - which performs quality control steps and maps the data using BWA-Meth instead of Bismark, before extracting CpG profiles using MethylDackel.
 * dmr - which uses the previously generated CpG profiles to identify differentially methylated regions between groups.
 * dcv - which uses the previously generated CpG to deconvolute the data and identify which tissues samples belong to based on methylation profiles

