
# CRISPRidentified: Identification of CRISPR arrays using machine learning approach

CRISPRIdentifier is a tool to search for CRISPR arrays which utilises 
machine learning approach for distinguishing false candidates from true CRISPRS.
CRISPRIdentifier, performs three steps: detection, feature extraction and 
classification based on manually curated sets of positive and negative examples of CRISPR arrays.
The identified CRISPR arrays are then reported to the user accompanied by detailed annotation.
We demonstrate that our approach identifies not only previously detected CRISPR arrays,
but also CRISPR array candidates not detected by other tools. Compared to other methods,
our tool has a drastically reduced false positive rate. In contrast to the existing tools, CRISPRIdentifier
approach not only provides the user with the basic statistics on the identified CRISPR arrays
but also produces a certainty score as an intuitive measure of the likelihood that a given
genomic region is a CRISPR array.

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

Before running the tool one need to activate the corresponding environment.

```
conda activate crispr_identifier_env
```

## Running the tool 

We suggest running the command line interface.
We prepared the test folder which can be used for the test run.

Example of running the tool over a folder of files:

```
python CRISPRidentifier.py --input_folder TestInput
```

Example of running the tool over a single multiline fasta input:
```
python CRISPRidentifier.py --file TestInputMultiline/MultilineFasta.fasta
```

### Flags

#### Mandatory flags
The only mandatory parameter which has to be specifier is the input.
Our approach has two options to handle the input. User has to specify either the path to the folder with the input fasta files
or the full path to a single fasta input file.

##### Input as a folder of fasta files

--input_folder "Input folder path"

Specifies the mode where a folder with fasta files which will be used as the input for the tool. The CRISPR array search will be
then conducted separately for each file in the corresponding input folder

```
python CRISPRidentifier.py --input_folder TestInput
```

##### Input as a single file

--file "Input file path"

Specifies the mode where a singe file is used as the tool input. The file might contain a single entry or multiple entries. 
The CRISPR array search will be done for each entry independently.

For example:

```
python CRISPRidentifier.py --file InputFile
```


#### Optional flags

##### Output

--result_folder "Result folder name"

Specifies the path and name of the folder with the output results. If not specified the results will appear in "Results" folder


For example:

```
python CRISPRidentifier.py --input_folder TestInput --result_folder Results
```

--pickle_report "Folder name for pickle report"

Specifies if found CRISPR arrays should be stored also as python objects. Turned off by default.


For example:

```
python CRISPRidentifier.py --input_folder TestInput --pickle_report PickleReportFolder
```


##### Classification parameters

--model "Model to use"


Takes values: 8, 9, 10, ALL and specifies the classification model. The default value is ALL.
If the value "ALL" is picked for the flag the certainty score will be calculated as average between all three available models.


For example:

```
python CRISPRidentifier.py --input_folder TestInput --model 8
```


```
python CRISPRidentifier.py --input_folder TestInput --model ALL
```


##### Performance speed
--fast_run "Bool"

Specifies if the repeat set enhancement step should be skipped which drastically speeds up the process but might decrease the recall quality.
Only matching pairs found with Vmatch will be used as repeat candidates. Automatically turns off filter approximation and start_end approximation (see enhancement_max_min and enhancement_start_end)
Turned off by default.

For example:

```
python CRISPRidentifier.py --input_folder TestInput --fast_run True
```

--enhancement_max_min "Bool"

Specifies if the filter approximation based on the max. and min. elements should be built
The default value is True 

--enhancement_start_end

Specifies if the start/end omitting of the repeat candidates should be done to enrich the candidate set.
The default value is True


For example:

```
python CRISPRidentifier.py --input_folder TestInput --enhancement_max_min True --enhancement_start_end False
```

##### Candidate filtering criteria 


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


For example:

```
python CRISPRidentifier.py --input_folder TestInput --min_len_rep 25 --max_len_rep 50 --min_repeats 2 
```


##### Additional computations


--cas "Bool"

Specifies if Cas genes should be predicted. The default value is False.

--is_element "Bool"

Specifies if IS-Elements should be predicted. The default value is False.


```
python CRISPRidentifier.py --input_folder TestInput --cas True --is_element True 
```

## Output files

The output folder for each input entries consist of the following files.

There are 5 files of text representations of the found CRISPR arrays what correspond to:
Best candidates

Alternative candidates

Possible candidates

Discarded possible candidates 

Bad candidates

Each text representation is complemented with a .bed format summary which specifies the CRISPR array start and end positions as well as the strand.

Finally the output contains a .csv summary of all the identified CRISPR arrays.







