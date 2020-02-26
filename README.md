
# CRISPRIdentifier 

CRISPR-ML-Identifier is a tool to search for CRISPR arrays which utilises machine learning approach
 for distinguishing false candidates from true CRISPRS

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. 

### Prerequisites

First you need to install Miniconda
Then create an environment and install the required libraries in it


### Creating a Miniconda environment 

First we install Miniconda for python 3.
Miniconda can be downloaded from here:

https://docs.conda.io/en/latest/miniconda.html 

Then Miniconda should be installed. On a linux machine the command is similar to this one: 

```
bash Miniconda3-latest-Linux-x86_64.sh
```

Then we create an environment. The necessary setup is provided in the "environment.yml" file inside the "for_environment" directory

In order to install the corresponding environment one can execute the following command from the "for_environment" directory

```
conda env create -f environment.yml
```

### Additional preparations
HMM models for the hmmsearch tool might contain a lot of information.
The Cas gene models exceed the file size limit provided by github.
In order to include the Cas genes into the project we put that file as a zip archive.

Please unzip models_cas_genes.zip into the same folder.


### Activation of environment

Before running the tool one need to activate the corresponding environment

```
conda activate crispr_identifier_env
```

## Running the tool

We suggest running the command line interface

for example:

```
python CRISPRidentifier.py --input_folder TestInput
```

### Flags

#### Mandatory flags
The only mandatory parameter which has to be specifier is the input.
Our approach has two options to handle the input. User has to specify either the path to the folder with the input fasta files
or the full path to a single fasta input file.

--input_folder "Input folder path"

Specifies the folder with fasta files which will be used as input for the tool.

--file "Input file path"

#### Optional flags

##### Output

--result_folder "Result folder name"

Specifies the path and name of the folder with the output results. If not specified the results will appear in "Results" folder

--pickle_report "Folder name for pickle report"

Specifies if found CRISPR arrays should be stored also as python objects. Turned off by default.


##### Classification parameters

--model "Model to use"


Takes values: 8, 9, 10, ALL and specifies the classification model. The default value is 8.
If the value "ALL" is picked for the flag the certainty score will be calculated as average between all three available models.


##### Performance speed

--parallel  "Bool"

Specifies if multiprocessing is used or not. Default value is True.

--cpu "int"

--fast_run "Bool"

Specifies if the repeat enhancement step should be skipped which drastically speeds up the process but might decrease the recall quality.
The default value is False

--enhancement_max_min "Bool"

Specifies if the filter approximation based on the max. and min. elements should be built
The default value is True 

--enhancement_start_end

Specifies if the start/end omitting of the repeat candidates should be done to enrich the candidate set.
The default value is True




##### Filtering criteria 


--min_len_rep  "int"

Specifies the minimum length of repeats in a CRISPR array. The default value: 21

--max_len_rep "int"

Specifies the maximum length of repeats in a CRISPR array. The default value: 55

--min_len_spacer "int"

Specifies the minimum average length of spacers in a CRISPR array. The default value: 18

--max_len_spacer "int"

Specifies the maximum average length of spacers in a CRISPR array. The default value: 78

--min_repeats "int"

Specifies the minimum number of repeats in a CRISPR array. The default value: 3


##### Additional computations


--cas "Bool"

Specifies if Cas genes should be predicted. The default value is False.

--is_element "Bool"

Specifies if IS-Elements should be predicted. The default value is False.



## Running the tool

The simplest possible way to execute the tool is just to provide the only mandatory flag which is the input folder.

We prepared the test folder which can be used for the test run.

```
python CRISPRidentifier.py --input_folder TestInput
```




