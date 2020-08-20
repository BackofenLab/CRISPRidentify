import os
import subprocess


def rev_compliment(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', "K": "K"}
    try:
        compliment_seq = "".join([complement[nt] for nt in seq])
    except KeyError:
        compliment_seq = ""
        for char in seq:
            if char in complement:
                compliment_seq += complement[char]
            else:
                compliment_seq += char
    return compliment_seq[::-1]


def to_rna(seq):
    return seq.replace("T", "U")


def get_orientation(sequence):
    rev_comp = rev_compliment(sequence)

    orig_rna, rev_comp_rna = to_rna(sequence), to_rna(rev_comp)

    with open("forward.txt", "w") as f:
        f.write(orig_rna)

    with open("reversed.txt", "w") as f:
        f.write(rev_comp_rna)

    cmd = "tools/strand_prediction/EDeN -a TEST -i forward.txt -M 1 -r 3 -d 3 -f SEQUENCE -g" \
          " DIRECTED -m tools/strand_prediction/DR_Repeat_model"

    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    process.communicate()
    f = open("prediction", "r")
    val_in_plus = float((f.readline().split())[1])
    f.close()
    cmd = "tools/strand_prediction/EDeN -a TEST -i reversed.txt -M 1 -r 3 -d 3 -f SEQUENCE -g" \
          " DIRECTED -m tools/strand_prediction/DR_Repeat_model"

    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    process.communicate()
    f = open("prediction", "r")
    val_in_minus = float((f.readline().split())[1])
    f.close()

    os.remove("forward.txt")
    os.remove("reversed.txt")
    os.remove("prediction")
    val_in_minus *= -1

    val = (val_in_plus + val_in_minus) / 2
    if val > 0:
        return 1
    else:
        return 0


class StrandComputationSingleCrispr:
    def __init__(self, crispr_array):
        self.crispr_array = crispr_array
        self.strand = None

        self._compute_strand()

    def _compute_strand(self):
        consensus = self.crispr_array.consensus
        self.strand = get_orientation(consensus)
        self.strand = "Forward" if self.strand else "Reverse"

    def output(self):
        return self.strand


class StrandComputation:
    def __init__(self, list_of_crisprs):
        self.list_of_crisprs = list_of_crisprs
        self.dict_strands = {}

        self._compute_all_strands()

    def _compute_all_strands(self):
        for index, crispr in enumerate(self.list_of_crisprs):
            consensus = crispr.consensus
            strand = get_orientation(consensus)
            strand = "Forward" if strand else "Reverse"
            self.dict_strands[index] = strand

    def output(self):
        return self.dict_strands
