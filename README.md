# Zhou Lab @ AGIS Pan_TE 
Both for linear or graph genome TE detection and annotation

<img width="800" alt="The Pan_TE workflow" src="https://github.com/unavailable-2374/Pan_TE/tree/main/image/PGTA.jpg?raw=true" >

## Requirements

### Tools
- [Look4LTRs ](https://github.com/BioinformaticsToolsmith/Look4LTRs)
- [RMblast]
- [RepeatScout]
- [RECON]
- [DeepTE]
- [Minigraph-Cactus]
- [TRF-finder]
- [CD-Hit]
- [MAFFT]
- [RepeatClassfier]
- [blast]
- [Pfam-scan]

## Software Installation 

### 1.Download the latest Pipeline:

  git clone https://github.com/unavailable-2374/Genome-Wide-annotation-pipeline.git

### 2.Install

    cd Pan_TE
    export PATH=/PATH/TO/bin >> ~/.bashrc
    mamba env create -f pgta.yml
    conda activate pgta
## Usage

        Usage:
          perl Pan_TE.pl [options]
        For example:
