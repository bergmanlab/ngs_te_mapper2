# **ngs_te_mapper**: <sub><sup>A program to identify transposable element insertions using next generation sequencing data</sup></sub>

# Table of Contents
* [Introduction](#intro)
* [Installation](#install)
* [Usage](#run)
* [Output](#output)
* [Getting help](#help)
* [License](#license)

# <a name="intro"></a> Introduction
ngs_te_mapper is an python re-implementation of the method for detecting transposable element (TE) insertions from next-generation sequencing (NGS) data published in [Linheiro and Bergman (2012) PLoS ONE 7(2): e30008](http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0030008). The original method only detected non-reference (aka _de novo_) TE insertions, however it has been extended to identify TE insertions also found in the reference genome. Additionally, the current implementation uses [BWA](http://bio-bwa.sourceforge.net/) as a short read mapping engine instead of BLAT, as was used originally in Linheiro and Bergman (2012).

Non-reference TE insertions are detected using a two-stage process that relies on the presence of target site duplications (TSDs) in the region flanking the TE insertion. An overview of the two-stage mapping procedure is shown below, and is taken from [Figure 1](https://journals.plos.org/plosone/article/figure/image?size=large&id=10.1371/journal.pone.0030008.g001) of Linheiro and Bergman (2012). 
<p align="center">
<img src="https://raw.githubusercontent.com/bergmanlab/ngs_te_mapper/master/img/journal.pone.0030008.g001.jpg?raw=true"/>
</p>

In the first stage, raw reads from a whole genome shotgun sequence are used to query against a library of reference TE sequences. 'Junction reads' that span the start/end of TE and genomic flanking sequences are retained. Such reads are often referred as 'split-reads', although in reality these reads are not split in the resequenced genome.

In the second stage, unmodified junction reads on each side of TE identified in the first stage are separately aligned to the repeatmasked reference genome. Genome-wide coverage profiles are computed and regions with enriched coverage are annotated in bed format that represents 5' and 3' clusters of junction reads. The region of overlap between 5' and 3' clusters of junction reads defines the location of the target site duplication (TSD) for non-reference TE insertions. The orientation of the TE is determined from the relative orientation of alignments of the junction reads to the reference genome and TE library.

Non-reference TE insertion sites are annotated as the span of TSD on zero-based, half-open coordinates and orientation is assigned in the strand field, following the framework described in [Bergman (2012) Mob Genet Elements. 2:51-54](http://www.landesbioscience.com/journals/mge/article/19479/). 

ngs_te_mapper also includes a third stage mapping procedure to estimate the allele frequency of the insertion. In this stage the raw reads are used to query against the candidate-masked reference genome (all the non-candidate insertion regions are hard-masked). For each non-reference TE insertion candidate, 'Junction reads' covering 5' and 3' side of insertion were counted as number of soft-clipped reads overlapping 10bp window on 5' and 3' side of TSD, respectively. 'Reference reads' were counted as number of non-soft-clipped reads spanning the TSD with at least 3bp extension on both side. The allele frequency was estimated as max(5' Junction reads, 3' Junction reads)/'Reference reads'.

Reference TE insertions are detected using a similar strategy to non-reference insertions, independently of any reference TE annotation. The first stage in detecting reference TE insertions is identical to the first stage of detecting non-reference TE insertions described above. The second stage in identifying reference TE insertions involves alignment of the renamed, but otherwise unmodified, junction reads to the reference genome. Alignments of the complete junction read (i.e. non-TE and TE components) are clustered to identify the two ends of the reference TE insertion. The orientation of the reference TE is then determined from the relative orientation of alignments of the junction reads to the reference genome and TE library.

# <a name="install"></a> Installation
## Use Conda to install software dependencies
ngs_te_mapper is written in python3 and is designed to run on linux operating system. Installation of software dependencies for ngs_te_mapper is automated by Conda, thus a working installation of Conda is required to install ngs_te_mapper. Conda can be installed via the Miniconda installer.

### Installing Miniconda (Python 3.X)
```
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O $HOME//miniconda.sh
bash ~/miniconda.sh -b -p $HOME/miniconda # silent mode
echo "export PATH=\$PATH:\$HOME/miniconda/bin" >> $HOME/.bashrc # add to .bashrc
source $HOME/.bashrc
conda init
```
- `conda init` requires you to close and open a new terminal before it take effect

### Update Conda
```
conda update conda
```

## Install software dependencies
After installing and updating Conda, you can now use conda to install dependencies and create running environment for ngs_te_mapper.
### Clone ngs_te_mapper Repository
```
git clone git@github.com:bergmanlab/ngs_te_mapper.git
cd ngs_te_mapper
```
### Create ngs_te_mapper Conda Environment
```
conda env create -f envs/ngs_te_mapper.yml
```
- This installs all the software dependencies needed to run ngs_te_mapper into the ngs_te_mapper Conda environment.

### Activate ngs_te_mapper Conda Environment
```
conda activate ngs_te_mapper
```
- This adds the dependencies installed in the ngs_te_mapper conda environment to the environment PATH so that they can be used by the ngs_te_mapper scripts.
- This environment must always be activated prior to running any of the ngs_te_mapper scripts
- NOTE: Sometimes activating conda environments does not work via conda activate myenv when run through a script submitted to a queueing system, this can be fixed by activating the environment in the script as shown below
```
CONDA_BASE=$(conda info --base)
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate ngs_te_mapper
```
- For more on Conda: see the [Conda User Guide](https://docs.conda.io/projects/conda/en/latest/index.html).

## Running ngs_te_mapper on test dataset
- A test dataset is provided in the `test/` directory, you can test whether your ngs_te_mapper installation is successful by running ngs_te_mapper on this dataset, which should take less than one minute to finish on a single thread machine.
```
conda activate ngs_te_mapper
cd test
python3 ../sourceCode/ngs_te_mapper.py -o test_output -f reads.fastq -r ref_1kb.fasta -l library.fasta
```

# <a name="run"></a> Usage
## Command line help page
```
usage: ngs_te_mapper.py [-h] -f READS -l LIBRARY -r REFERENCE [-n REGION]
                        [-w WINDOW] [--min_mapq MIN_MAPQ] [--min_af MIN_AF]
                        [--tsd_max TSD_MAX] [--gap_max GAP_MAX] [-m MAPPER]
                        [-t THREAD] [-o OUT] [-p PREFIX] [-k]

Script to detect non-reference TEs from single end short read data

required arguments:
  -f READS, --reads READS
                        raw reads in fastq or fastq.gz format, separated by
                        comma
  -l LIBRARY, --library LIBRARY
                        TE concensus sequence
  -r REFERENCE, --reference REFERENCE
                        reference genome

optional arguments:
  -h, --help            show this help message and exit
  -n REGION, --region REGION
                        region to filter
  -w WINDOW, --window WINDOW
                        merge window for identifying TE clusters (default =
                        10)
  --min_mapq MIN_MAPQ   minimum mapping quality of alignment (default = 20)
  --min_af MIN_AF       minimum allele frequency (default = 0.1)
  --tsd_max TSD_MAX     maximum TSD size (default = 25)
  --gap_max GAP_MAX     maximum gap size (default = 5)
  -m MAPPER, --mapper MAPPER
                        read alignment program (default = 'bwa')
  -t THREAD, --thread THREAD
                        thread (default = 1)
  -o OUT, --out OUT     output dir (default = '.')
  -p PREFIX, --prefix PREFIX
                        output prefix
  -k, --keep_files      If provided then all intermediate files will be kept
                        (default: remove intermediate files)
```

# <a name="output"></a> Output
ngs_te_mapper outputs reference and non-referece TE insertion predictions in BED format (0-based).
- `<sample>.nonref.bed`: non-reference TE insertion annotation predicted by ngs_te_mapper pipeline in BED format (0-based).
- `<sample>.ref.bed`: reference TE insertion annotation predicted by ngs_te_mapper pipeline in BED format (0-based).

## TE insertion annotation in bed format
ngs_te_mapper generates standard BED file `<sample>.nonref.bed` and `<sample>.ref.bed` that have detailed information for each reference and non-reference TE insertion.

Column | Description
-- | --
chromosome | The chromosome name where the TE insertion occurred
position | Starting breakpoint position of the TE insertions.
end | Ending breakpoint position of the TE insertions.
info | Includes TE family, TSD, Allele Frequency, 3' support, 5' support and reference reads. Separated by '\|'.
score | '.'
strand | Strand that TE insertion occurs

## Log file output by ngs_te_mapper
For each ngs_te_mapper run, a log file called `<sample>.log` is generated that records all the major steps in the program and error messages.
`<sample>.log`: log file of ngs_te_mapper run.

# <a name="help"></a> Getting help
Please use the [Github Issue page](https://github.com/bergmanlab/ngs_te_mapper/issues) if you have questions.

# <a name="license"></a> License
Copyright (c) 2020 Shunhua Han, Raquel S. Linheiro, Michael G. Nelson, Casey M. Bergman

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.