import subprocess
import os

conda_activation_cmd = "export PATH=~/miniconda3/bin:$PATH && conda activate crispr_strand_env && python CRISPRstrand.py -r -i Example/Input3.fa && conda deactivate"
subprocess.run(conda_activation_cmd, shell=True)


#os.system("source activate crispr_strand_env")