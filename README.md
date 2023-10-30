# Zhou Lab @ AGIS Pan_TE 
Both for linear or graph genome TE detection and annotation

<img width="800" alt="The Pan_TE workflow" src="https://github.com/unavailable-2374/Pan_TE/tree/main/image/PGTA.jpg?raw=true" >

## Requirements

### Tools
- [Look4LTRs ](https://github.com/BioinformaticsToolsmith/Look4LTRs)
- [RMblast ] ()
- [RepeatScout ] ()
- [RECON ] ()
- [DeepTE ] ()
- [Minigraph-Cactus ] ()
- [TRF-finder ] ()
- [CD-Hit ] ()
- [MAFFT ] ()
- [RepeatClassfier ] ()
- [blast ] ()
- [hmmscan ] ()

## Software Installation 

### 1.Download the latest Pipeline:

  git clone https://github.com/unavailable-2374/Pan_TE.git

### 2.Install

    cd Pan_TE
    export PATH=/PATH/TO/bin >> ~/.bashrc
    mamba env create -f pgta.yml
    conda activate pgta
## Usage

  Usage:
    perl Pan_TE.pl [options]
  For example:
    Pan_TE.pl --genome gemome.fa  --cpu 40 --out demo --model P --model_dir path/to/model --hmmscan path/to/hmmscam
    Parameters:
        [General]
        --ref <string>     Required
    genome file in fasta format.

    [other]
        --list <string> default:NA
        path file for genome .
    
        --out <string>    default: .
        the work dir.
    
        -M <int>    
        memory limit (in MB) for the program, default 0; 0 for unlimitted;
    
        --model <string>
        P or M or F or O. P:Plants, M:Metazoans, F:Fungi, and O: Others.
    
        --model_dir <string>
        Provide model_dir that could be downloaded from website (optional requirements). 
    
        --hmmscan <int>
        path to hmmscan
    
        --cpu <int>    default: 4
        the number of threads, preferably in multiples of 4.
    
        --sensitive defaut:normal
        Sensitivity represents a parameter when merging genomes; 
        the more sensitive is more cautious for merging, 
        the longer the procedure takes and the more TEs may be obtained.
        super,hight,normal,low
    
        --help|-h Display this help info
    
    Version: 1.0.0
    USAGE
