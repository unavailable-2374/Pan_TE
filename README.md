# Zhou Lab @ AGIS Pan_TE 
Design for linear or graph genome TE detection and annotation

## Software Installation 

### 1.Download the latest Pipeline:

    git clone https://github.com/unavailable-2374/Pan_TE.git

### 2.Install:

    cd Pan_TE
    chomp 750 bin/*
    echo 'export PATH="/YOUR/PATH/TO/bin:$PATH"' >> ~/.bashrc //* like: echo 'export PATH=/public/home/soft/Pan_TE/bin:$PATH"' >> ~/.bashrc
    mamba env create -f env/pgta.yml
    ln -s /PATH/TO/miniconda3/envs/PGTA/bin/x86_64-conda-linux-gnu-g++ /PATH/TO/miniconda3/envs/PGTA/bin/g++
    export PERL5LIB=/home/tool/Pan_TE/share:/home/tool/miniconda3/envs/PGTA/share/RepeatMasker:$PERL5LIB
    conda activate PGTA

### [Look4LTRs](https://github.com/BioinformaticsToolsmith/Look4LTRs) intallation section:
    git clone https://github.com/BioinformaticsToolsmith/Look4LTRs.git
    cd Look4LTRs
    mkdir bin
    cd bin
    cmake ..
    vim ../src/Util.h add [#include <cstdint>] on the top
    vim ../src/KmerHistogram.h add [#include <cstdint>] on the top
    make look4ltrs
    export PATH=/YOUR/PATH/TO/Look4LTRs/bin >> ~/.bashrc

### [ClassifyTE](https://github.com/manisa/ClassifyTE/tree/master) intallation section:
    git clone https://github.com/manisa/ClassifyTE.git
go to [this link](https://drive.google.com/file/d/1CuDciG0Ru5zRBhffjQmgJdqSMQB89mfh/view?usp=sharing)
    
- Click on **ClassifyTE_Models.zip**. This will automatically download models built with TE datasets.
- Unzip and copy all the models from "ClassifyTE_Models" directory into the folder **model** inside the root folder **ClassifyTE**
 
## Usage
    Usage:
    perl $0 [options]

    Example:
        perl $0 --genome genome.fasta --cpu 80 --ClassifyTE_dir Your_Path_To_ClassifyTE 
    
    Parameters:
    [General]
        --genome <string>         Required. Genome file in FASTA format.
        --ClassifyTE_dir <string> Provide path to ClassifyTE for classification.
    
    [Other]
        --D <string>           Default: NA. Path for VCF, see gfa.list for format.
        --out <string>            Default: current directory. The work directory.
        -M <int>                  Memory limit (in MB), default: 0 (unlimited).
        --cpu <int>               Default: 4. Number of threads, preferably in multiples of 4.
        --fragment_size <int>     Default: 40000. Length for fragment.
        --help|-h                 Display this help information.
    
    Version: 1.0.0
    USAGE
