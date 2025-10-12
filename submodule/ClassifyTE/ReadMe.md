# ClassifyTE: Stacking-based Machine Learning Framework for Hierarchical Classification of transposable 

- [Paper](https://academic.oup.com/bioinformatics/advance-article-abstract/doi/10.1093/bioinformatics/btab146/6158037?redirectedFrom=fulltext#)
- The basic framework has been adapted from the state-of-the-art method [1].

### Table of Content
- [Setup](#getting-started)
	- [Download and install code](#download-and-install-code)
	- [Download Models](#download-models)
	- [Demo](#demo)
- [Training](#training)
	- [Download Datasets](#download-datasets)
- [References](#references)


# Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

## Prerequisites

You would need to install the following software before replicating this framework in your local or server machine.

```
Java JDK
Python version 3.5+
Aanaconda version 3+
```

## Download and install code
- Retrieve the code
```sh
git clone https://github.com/manisa/ClassifyTE.git
cd ClassifyTE
```

- Create and activate the virtual environment with python dependendencies. 
```
conda env create -f environment.yml python==3.7
conda activate ClassifyTE_env
```


## Download Models
- Go to [this link](https://drive.google.com/file/d/1CuDciG0Ru5zRBhffjQmgJdqSMQB89mfh/view?usp=sharing).
- Click on **ClassifyTE_Models.zip**. This will automatically download models built with TE datasets.
- Unzip and copy all the models from "ClassifyTE_Models" directory into the folder **model** inside the root folder **ClassifyTE**.
- Your directory structure should look like this:

```
ClassifyTE/
	models/
		ClassifyTE_combined.pkl
		ClassifyTE_repbase.pkl
		ClassifyTE_pgsb.pkl
```

## Demo
To run the program on test TE sequence:
- Your directory structure should look like this with the demo.fasta file.

```
ClassifyTE/
	data/
		demo.fasta
```
- Firstly, Run following python command from the root directory to generate feature file.
```
python generate_feature_file.py -f demo.fasta -d demo_features -o demo_features.csv
```

#### *Parameters*
For **generate_feature_file.py**, the user has to provide two parameters:
- -f for fasta filename from  **data** directory.
- -o for resulting feature file name with **.csv** extension [**Optional**]  [By default the feature filename is **feature_file.csv**.]
- -d for feature folder name which is by default *features* unless user provide a feature folder name.
- Then run following python command from the root directory to get the prediction on new TE sequences. Prior following command, please make sure that all the model files have already been added to **models** directory. 

```
python evaluate.py -f demo_features.csv -n node.txt -d demo_features -m ClassifyTE_combined.pkl -a lcpnb
```
#### *Parameters*
For **evaluate.py**, the user has to provide following parameters:
- -f for feature file name which is by default *feature_file.csv* unless user have provided a feature filename in earlier step.
-d for feature folder name which is by default *features* unless user have provided a feature folder name while generating features.
- -n for node filename which is by default *node.txt*. Node file consists of numerical representation of taxonomy of the dataset. Please check nodes folder for other node files for each dataset.
- -m for model filename which has **.pkl** as file extension. All the model files must have been added in **models** directory.
- -a for algorithm choice (lcpnb or nllcpn)

- Finally, check files inside **output** folder for predicted label of the TE sequence/s.

## Deployment
To run the program on **new** TE sequence:
- Place your fasta file (with single or multiple fasta sequences) inside **data** folder.
- Your directory structure should look like this

```
ClassifyTE/
	data/
		[your_fasta_file]
```
- Firstly, Run following python command from the root directory to generate feature file.
```
python generate_feature_file.py -f your_fasta_file_name -o your_feature_file_name -d your_feature_directory
```

#### *Parameters*
For **generate_feature_file.py**, the user has to provide two parameters:
- -f for fasta filename from  **data** directory.
- -d for feature folder name if you want to replace the name of the feature directory so as to generate features for multiple fasta sequences.
- -o for resulting feature file name with **.csv** extension [**Optional**]  [By default the feature filename is **feature_file.csv**.]

- Then run following python command from the root directory to get the prediction on new TE sequences. Prior following command, please make sure that all the model files have already been added to **models** directory. 

```
python evaluate.py -f your_feature_file_name -d your_feature_directory -n node_file -m model_name
```
#### *Parameters*
For **evaluate.py**, the user has to provide following parameters:
- -f for feature file name which is by default *feature_file.csv* unless user have provided a feature filename in earlier step.
- -n for node filename which is by default *node.txt*. Node file consists of numerical representation of taxonomy of the dataset. Each node file is associated with the respective models trained on respective datasets. Please check below under **nodes** section for details.
- -m for model filename which has **.pkl** as file extension. All the model files must have been added in **models** directory.

- Finally, check **predicted_result.csv** file inside **output** folder for predicted label of the TE sequence/s.

## Download datasets
- Go to [this link](https://drive.google.com/file/d/1a18Kcv6PEJiShWm1Fm3aKdbhpznwXcKO/view?usp=sharing).
- Click on **ClassifyTE_Datasets.zip**. This will automatically download the datasets used to generate models.
- Unzip and copy all the datasets from **ClassifyTE_Datasets** directory into the folder **data** inside the root folder **ClassifyTE**.
- Your directory structure should look like this:

```
ClassifyTE/
	data/
		pgsb_feature_file.csv
		repbase_feature_file.csv
		combined.csv
```

## Training
- To replicate the training procedure, follow following command line
```
python train.py -f csv_file_name -n txt_node_file -m model_filename -c SVM_cost_parameter -g SVM_gamma_parameter
```

#### *Parameters*
For **train.py**, the user has to provide following parameters:
- -f for feature file name.
- -m for model file name
- -n for node filename. Node file consists of numerical representation of taxonomy of the dataset.
- -c for cost parameter of SVM with RBF kernel
- -g for gamma parameter of SVM with RBF kernel

- We have optimized cost and gamma parameters of SVM with RBF kernel for all three datasets. The cost and gamma parameters for training each datasets would be different. You will have to pass the hyper-parameters accordingly. 
```
For PGSB dataset : C=32, gamma=0.03125
For REPBASE dataset : C=128.0, gamma=0.0078125
For combined dataset : C=512.0, gamma=0.0078125
```

##### nodes 
- Under **node** directory you will find three files. These files consists of all the nodes in each of the corresponding datasets. These node files consist of numerical representation of taxonomy in the datasets.
	1.	node_pgsb.txt for pgsb dataset
	2.	node_repbase.txt for repase dataset
	3.	node.txt for combined dataset

- If you would like to train the model on your machine, the training example would look like as below:
```
python train.py -f combined.csv -n node.txt -m ClassifyTE_combined -c 512.0 -g 0.0078125
```


## Authors
Manisha Panta, Avdesh Mishra, Md Tamjidul Hoque, Joel Atallah
## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## References
[1] Nakano, F.K., et al. Top-down Strategies for Hierarchical Classification of Tranposable Elements with Neural Networks. In, IEEE. 2017.

