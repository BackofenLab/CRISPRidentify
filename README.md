
# CRISPRidentify: Identification of CRISPR arrays using machine learning approach

CRISPRidentify is a tool to search for CRISPR arrays which utilises 
machine learning approach for distinguishing false candidates from true CRISPRS.
CRISPRidentify, performs three steps: detection, feature extraction and 
classification based on manually curated sets of positive and negative examples of CRISPR arrays.
The identified CRISPR arrays are then reported to the user accompanied by detailed annotation.
We demonstrate that our approach identifies not only previously detected CRISPR arrays,
but also CRISPR array candidates not detected by other tools. Compared to other methods,
our tool has a drastically reduced false positive rate. In contrast to the existing tools, CRISPRidentify
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

Then we create an environment. The necessary setup is provided in the "environment.yml" file.

In order to install the corresponding environment one can execute the following command.

```
conda env create -f environment.yml
```

We recommend to install mamba package manager which is a faster alternative to conda.

```
conda install -c conda-forge mamba
```

Then we can create the environment using mamba.
```
mamba env create -f environment.yml
```

<sub><sub>We want to acknowledge Richard Stöckl @richardstoeckl for his contribution to the environment.yml file.</sub></sub>


### Additional preparations

CRISPRidentify utilizes CRISPRcasIdentifier for the detection of the cas genes.
If you are interested in cas gene result please install CRISPRcasIdentifier.

Please make sure that after you downloaded CRISPRcasIdentifier its relative path is:

```
tools/CRISPRcasIdentifier/CRISPRcasIdentifier/CRISPRcasIdentifier.py
```

You can find the CRISPRcasIdentifier tool and its description [here](https://github.com/BackofenLab/CRISPRcasIdentifier)

You need to make two steps:

Firstly, you need to download the CRISPRcasIdentifier tool:
```
wget https://github.com/BackofenLab/CRISPRcasIdentifier/archive/v1.1.0.tar.gz
tar -xzf v1.1.0.tar.gz
```
Secondly, you need to download the models:

Due to GitHub's file size constraints, authors made their HMM and ML models available in Google Drive. You can download them [here](https://drive.google.com/file/d/1YbTxkn9KuJP2D7U1-6kL1Yimu_4RqSl1/view?usp=sharing) and [here](https://drive.google.com/file/d/1Nc5o6QVB6QxMxpQjmLQcbwQwkRLk-thM/view?usp=sharing). Save both tar.gz files inside CRISPRcasIdentifier's directory.


### Activation of the environment

Before running CRISPRidentify one need to activate the corresponding environment.

```
conda activate crispr_identify_env
```

## Running CRISPRidentify

We prepared the test folder which can be used for the test run.

Example of running CRISPRidentify over a folder of files:

```
python CRISPRidentify.py --input_folder TestInput
```

Example of running CRISPRidentify over a single multiline fasta input:
```
python CRISPRidentify.py --file TestInputMultiline/MultilineFasta.fasta
```

### Flags

You can see the help by using the `-h` option

```

python CRISPRidentify.py -h

```

#### Mandatory flags
The only mandatory parameter which has to be specified is the input.
Our approach has two options to handle the input. User has to specify either the path to the folder with the input fasta files
or the full path to a single fasta input file.

##### Input as a folder of fasta files

* `--input_folder <path_to_the_folder>`

Specifies the mode where a folder with fasta files which will be used as the input for CRISPRidentify. The CRISPR array search will be
then conducted separately for each file in the corresponding input folder.

```
python CRISPRidentify.py --input_folder TestInput
```

##### Input as a single file

* `--file <path_to_the_file>`

Specifies the mode where a singe file is used as the input for the algorithm. The file might contain a single entry or multiple entries. 
The CRISPR array search will be done for each entry independently.

For example:

```
python CRISPRidentify.py --file InputFile
```
##### Input as a folder of multiline fasta files

* `-- input_folder_multifasta <path_to_the_folder>`

Specifies the mode where a folder with fasta files which will be used as the input for CRISPRidentify. The CRISPR array search will be
then conducted separately for each file in the corresponding input folder. The difference between this mode and the previous one is that
in this mode the input files can contain multiple entries.

For example:

```
python CRISPRidentify.py --input_folder_multifasta TestFolderMultiline
```

#### Optional flags

##### Output

* `--result_folder [paht_to_the_result_folder]`

Specifies the path and name of the folder with the output results. If not specified the results will appear in "Results" folder


For example:

```
python CRISPRidentify.py --input_folder TestInput --result_folder Results
```

* `--pickle_report [folder_to_put_pickle_results]`

Specifies if found CRISPR arrays should be stored also as python objects. Turned off by default.


For example:

```
python CRISPRidentify.py --input_folder TestInput --pickle_report PickleReportFolder
```


##### Classification parameters

* `--model [8/9/10/ALL]`


Takes values: 8, 9, 10, ALL and specifies the classification model. The default value is `ALL`.
If the value `ALL` is picked for the flag the certainty score will be calculated as average between all three available models.


For example:

```
python CRISPRidentify.py --input_folder TestInput --model 8
```


```
python CRISPRidentify.py --input_folder TestInput --model ALL
```


##### Performance speed
* `--fast_run [True/False]`

Specifies if the repeat set enhancement step should be skipped which drastically speeds up the process but might decrease the recall quality.
Only matching pairs found with Vmatch will be used as repeat candidates. Automatically turns off filter approximation and start_end approximation (see enhancement_max_min and enhancement_start_end)
Turned off by default.

For example:

```
python CRISPRidentify.py --input_folder TestInput --fast_run True
```

* `--enhancement_max_min [True/False]`

Specifies if the filter approximation based on the max. and min. elements should be built
The default value is True 

* `--enhancement_start_end [True/False]`

Specifies if the start/end omitting of the repeat candidates should be done to enrich the candidate set.
The default value is True


For example:

```
python CRISPRidentify.py --input_folder TestInput --enhancement_max_min True --enhancement_start_end False
```

##### Candidate filtering criteria 


* `--min_len_rep  [integer]`

Specifies the minimum length of repeats in a CRISPR array. The default value: 21

* `--max_len_rep [integer]`

Specifies the maximum length of repeats in a CRISPR array. The default value: 55

* `--min_len_spacer [integer]`

Specifies the minimum average length of spacers in a CRISPR array. The default value: 18

* `--max_len_spacer [integer]`

Specifies the maximum average length of spacers in a CRISPR array. The default value: 78

* `--min_repeats [integer]`

Specifies the minimum number of repeats in a CRISPR array. The default value: 3


For example:

```
python CRISPRidentify.py --input_folder TestInput --min_len_rep 25 --max_len_rep 50 --min_repeats 2 
```

#####Candidate Enhancement 

* `--degenerated' [True/False]`

Allows search for degenerated repeat candidates on both ends of the CRISPR array candidate. The default value: True

* `--margin_degenerated [int]`

Specifies the maximum length difference between a new spacer sequence (obtained with the search of degenerated repeats) and the average value of spacer length in the array. The default value: 30

* `--max_edit_distance_enhanced [int]`

Specifies the number of editing operations for candidate enhancement. The default value: 6


##### Additional computations

* `--strand[True/False]`

Specifies if the array orientation should be predicted. The default value is True.

* `--cas [True/False]`

Specifies if cas genes should be predicted. The default value is False.

* `--is_element [True/False]`

Specifies if IS-Elements should be predicted. The default value is False.


```
python CRISPRidentify.py --input_folder TestInput --cas True --is_element True 
```

## Output files

The output folder for each input entries consist of the following files:

* Bona-Fide_Candidates. The file will contain the representation of the found CRISPR arrays complemented with the support information. 
For each candidate the output will contain the values for extracted features as well as the certainty score of the used classifier.
On top of that in the support information you can find the orientation for each array, leader and downstream regions, cas genes and IS-elements (if the corresponding flags were selected).

* Alternative_Candidates. In this file we demonstrate alternative representations of bona-fide arrays. These alternative representations also got a high score from the classifier but this score was lower than the corresponding score of the bona-fide representation.
Alternative representation of a CRISPR array usually corresponds to a slightly longer/shorter repeat sequence but represents the same genomic region.

The candidates with the certainty scores between 0.4 and 0.75 are stored in Possible_Candidates and Possible_Discarded_Candidates

* Possible_Candidates. In this file the algorithm stores the candidate with the highest certainty score.

* Possible_Discarded.  Here are collected all the other representations 


The algorithm also demonstrates CRISPR-looking structures which obtained certainty score lower than 0.4 from the classifier.

* Low_score_candidates. The user can find these structures in this file.


On top of that the algorithm builds a csv summary. 

* Summary.csv

Following information can be found in the summary:

1. Array index
2. Array start
3. Array end
4. Array length
5. Consensus repeat
6. Repeat length
7. Average length of the spacers
8. Number of spacers
9. Array orientation
10. Array category

## Metagenomic analysis

CRISPRidentify is suitable for easy and powerful metagenomic analysis
When `--file` or `--input_folder` flag is used the pipeline with automatically generate two complete summaries 
: 

1. For all the identified arrays
2. For all labeled Cas genes


On top of that the user might use the flag:

`--fasta_report True`

This option with create three fasta files:
1. All the array sequences with their origins in the header
2. All the repeat sequences with their origins and locations in the arrays
3. All the spacer sequences with their origins and locations in the arrays

## Improving CRISPRidentify

We are constantly working on the improvements of CRISPRidentify. If you found a bug or incorrect/missing CRISPR array representation please submit via github issue interface.




