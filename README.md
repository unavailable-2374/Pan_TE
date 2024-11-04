# Pan_TE 
Design for linear or graph genome TE detection and classification 

## Software Installation 

### 1.Download the latest Pipeline:

    git clone https://github.com/unavailable-2374/Pan_TE.git

### 2.Install:

    cd Pan_TE
    chmod 750 bin/*
    export PATH=/YOUR/PATH/TO/bin //* like: export PATH=/public/home/soft/Pan_TE/bin:$PATH
    mamba env create -f env/pgta.yml
    ln -s /PATH/TO/miniconda3/envs/PGTA/bin/x86_64-conda-linux-gnu-g++ /PATH/TO/miniconda3/envs/PGTA/bin/g++
    conda activate PGTA

### [Look4LTRs](https://github.com/BioinformaticsToolsmith/Look4LTRs) intallation section:
    git clone https://github.com/BioinformaticsToolsmith/Look4LTRs.git
    cd Look4LTRs
    mkdir bin
    cd bin
    cmake ..
    vim ../src/Util.h add "#include <cstdint>" on the top
    vim ../src/KmerHistogram.h add "#include <cstdint>" on the top
    make look4ltrs
    export PATH=/YOUR/PATH/TO/Look4LTRs/bin:$PATH

### [ClassifyTE](https://github.com/manisa/ClassifyTE/tree/master) intallation section:
    git clone https://github.com/manisa/ClassifyTE.git
go to [this link](https://drive.google.com/file/d/1CuDciG0Ru5zRBhffjQmgJdqSMQB89mfh/view?usp=sharing)
    
- Click on **ClassifyTE_Models.zip**. This will automatically download models built with TE datasets.
- Unzip and copy all the models from "ClassifyTE_Models" directory into the folder **model** inside the root folder **ClassifyTE**
 
### Usage:
    Parameters:
   [General]
    --genome <string>         Required. Genome file in FASTA format.
    --ClassifyTE_dir <string> Provide path to ClassifyTE for classification.
    --RM_dir <string>         Provide path to rmblastn for classification.

[Other]
    --list <string>           Default: NA. Path file for GFA, see gfa.list for format.
    --out <string>            Default: current directory. The work directory.
    -M <int>                  Memory limit (in MB), default: 0 (unlimited).
    --flag <string>           PGGB or Minigraph-Cactus.
    --cpu <int>               Default: 4. Number of threads, preferably in multiples of 4.
    --fragment_size <int>     Default: 40000. Length for fragment.
    --help|-h                 Display this help information.

Version: 1.0.0
USAGE
