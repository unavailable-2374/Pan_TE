# Zhou Lab @ AGIS Pan_TE 
Design for linear or graph genome TE detection and annotation

![The Pan_TE workflow](https://github.com/unavailable-2374/Pan_TE/tree/main/image/PGTA.png)

## Software Installation 

### 1.Download the latest Pipeline:

    git clone https://github.com/unavailable-2374/Pan_TE.git

### 2.Install:

    cd Pan_TE
    export PATH=/PATH/TO/bin >> ~/.bashrc
    mamba env create -f pgta.yml
    pip install sklearn
    conda activate pgta

### Look4LTRs intallation section:
    git clone https://github.com/BioinformaticsToolsmith/Look4LTRs.git
    mkdir bin
    cd bin
    (If your default compiler meets the version requirement) 
    cmake ..
    (Or if you would like to specify a different compiler that meets the requirement)
    cmake .. -DCMAKE_CXX_COMPILER=your_compiler_name_for_example_g++-7
    (Or if this fails, try using this to set a different compiler. Replace the paths with your own.)
    cmake .. -DCMAKE_CXX_COMPILER=$HOME/C++/GCC/bin/g++ -DCMAKE_C_COMPILER=$HOME/C++/GCC/bin/gcc -DCMAKE_PREFIX_PATH=$HOME/C++/GCC
    make look4ltrs
 
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
          
              --help|-h Display this help info
          
          Version: 1.0.0
          USAGE
