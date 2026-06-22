# Zhou Lab @ AGIS Pan_TE 
Design for linear or graph genome TE detection and annotation

## Install via conda (bioconda)

Once the recipes in [`recipes/`](recipes/) are published to bioconda, the whole pipeline
(plus every tool dependency) installs in one command:

    conda create -n pan_te -c conda-forge -c bioconda pan_te
    conda activate pan_te
    # Fetch the non-redistributable data assets (Dfam libraries + ClassifyTE models):
    pan_te-setup-data --dfam-list
    pan_te-setup-data --classifyte-dir ~/Pan_TE_data/ClassifyTE \
        --dfam-url https://www.dfam.org/releases/Dfam_3.9/families/dfam39_full.0.h5.gz

The bioconda packaging (4 recipes: `look4ltrs`, `mdl-repeat`, `te-looker`, `pan_te`) and the
submission steps still pending (release tags, LICENSE files, sha256, PRs) are documented in
[`recipes/BIOCONDA_SUBMISSION.md`](recipes/BIOCONDA_SUBMISSION.md). To verify the recipes
build from local source before submitting, run `recipes/build_local.sh build`.

## Software Installation (from source)

### 1.Download the latest Pipeline:

    git clone --recursive https://github.com/unavailable-2374/Pan_TE.git

### 2.Install:

    cd Pan_TE
    chmod 750 bin/*
    echo 'export PATH=$HOME/tools/Pan_TE/bin:$PATH' >> ~/.bashrc
    mamba env create -f env/pgta.yml
    echo 'export PERL5LIB=$HOME/tools/Pan_TE/share:$HOME/tools/miniconda3/envs/PGTA/share/RepeatMasker:$PERL5LIB' >> ~/.bashrc
    conda activate PGTA
    cd submodule/ClassifyTE
    conda env create -f environment.yml
    cd submodule/LookLTRs
    mkdir bin
    cd bin
    export CC=$HOME/tools/miniconda3/envs/PGTA/bin/x86_64-conda-linux-gnu-gcc
    export CXX=$HOME/tools/miniconda3/envs/PGTA/bin/x86_64-conda-linux-gnu-c++
    cmake ..
    make look4ltrs
    echo 'export PATH=$HOME/tools/Pan_TE/submodule/Look4LTRs/bin:$PATH' >> ~/.bashrc


### [ClassifyTE](https://github.com/manisa/ClassifyTE/tree/master) database download:
go to [this link](https://drive.google.com/file/d/1CuDciG0Ru5zRBhffjQmgJdqSMQB89mfh/view?usp=sharing)
    
- Click on **ClassifyTE_Models.zip**. This will automatically download models built with TE datasets.
- Unzip and copy all the models from "ClassifyTE_Models" directory into the folder **models** inside the root folder **ClassifyTE**
 
## Usage
    Usage:
    perl $0 [options]

    Example:
        perl $0 --genome genome.fasta --cpu 80 --model-dir Your_Path_To_ClassifyTE 
    
    Parameters:
    [General]
        --genome <string>         Required. Genome file in FASTA format.
        --model-dir <string>      Provide path to ClassifyTE for classification.
    
    [Other]
        --vcf-dir <string>        Default: NA. Path for VCF, see gfa.list for format.
        --out <string>            Default: current directory. The work directory.
        -M <int>                  Memory limit (in MB), default: 0 (unlimited).
        --threads <int>           Default: 4. Number of threads, preferably in multiples of 4.
        --fragment_size <int>     Default: 40000. Length for fragment.
        --help|-h                 Display this help information.
    
    Version: 1.0.0
    USAGE
