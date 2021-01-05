#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  3 20:44:47 2019

@author: ekrem
"""
import re
import os
import numpy as np
import pandas as pd
from collections import namedtuple
import sklearn.model_selection as skms

#from Models.repeat import RepeatTyper
#from Models.xgb import XGB

# =============================================================================
# DICTIONARIES
# =============================================================================

one_hot_encoding_dict    = {
                           'A': np.array([0, 0, 0, 1]).reshape((-1,1)), # A
                           'T': np.array([0, 0, 1, 0]).reshape((-1,1)), # T
                           'G': np.array([0, 1, 0, 0]).reshape((-1,1)), # G
                           'C': np.array([1, 0, 0, 0]).reshape((-1,1)), # C
                           'N': np.array([0, 0, 0, 0]).reshape((-1,1)), # N
                           '-': np.array([0, 0, 0, 0]).reshape((-1,1)), # N
                           } 
complement_encoding_dict = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 'N':'N', '-':'N'}

# =============================================================================
# Decoders/Encoders/Helpers
# =============================================================================

'''
Performs one-hot encoding on repeat sequences.
Reshapes each encoded nucleotide to a column vector, then concatenates them.
Returns representation with shape (4, sequence_length, 1). Last dimension is 
channels (always 1 in this case) to comply with 2D-CNN.
'''
def one_hot_encode_sequence(repeat):
    return np.concatenate(list(map(one_hot_encoding_dict.get, repeat)), axis = 1).reshape(4, -1, 1)

'''
Returns the reverse complement of given repeat sequence.
'''
def reverse_complement(repeat):
    complement = list(map(complement_encoding_dict.get, repeat))
    reverse_complement = ''.join(reversed(complement))
    return reverse_complement

'''
Returns the length of all repeat sequences in the dataframe.
'''
def get_seq_len(repeat):
    return repeat.shape[1]

'''
Pads all repeat sequences with zero from right-hand side to the length of 
maximum sequence.
'''
def pad_seq(repeat, max_seq_len):
    seq_len = get_seq_len(repeat)
    if max_seq_len>=seq_len:
        return np.pad(repeat, ((0,0),(0, max_seq_len-seq_len), (0,0)), 'constant')
    else:
        return repeat[:,:max_seq_len]

'''
Parses the ID and splits accession number, program and strand into different fields. Also implies start and end indices in case.
''' 
def parse_ID(raw_sample):
    ID, program, strand = raw_sample['ID'].split('-')
    program, start, end = program.split('_')
    raw_sample['Accession'] = ID
    raw_sample['Program'] = program
    raw_sample['Strand'] = strand
    raw_sample['Start'] = start
    raw_sample['End'] = end
    return raw_sample

# =============================================================================
# FUNCTIONS TO GET NEGATIVE DATASET
# =============================================================================
'''
Wrapper method.
Forms negative dataset.
'''
def seperate_to_pos_neg(df_pos, _all = False, keep_orig = False):
    
    if keep_orig:
        df_pos = df_pos.copy()
    # Form negative dataset
    df_neg = get_neg(df_pos, _all)
    
    df_pos['Label'] = 1
    df_neg['Label'] = 0

    return df_pos, df_neg

'''
Applies reverse complement to repeat sequences and to strand.
'''
def get_neg(df_pos, _all = False):
    
    df_neg = df_pos.copy()
    df_neg['Cons'] = df_neg['Cons'].apply(reverse_complement)
    if _all:
        df_neg['Conservation'] = df_neg['Conservation'].apply(np.flip)
        df_neg['Edge_Conservations'] = df_neg['Edge_Conservations'].apply(np.flip)
        df_neg['Up_Down_AT_content'] = df_neg['Up_Down_AT_content'].apply(np.flip)
        df_neg['Up_Down_AT_content'] = df_neg['Up_Down_AT_content'].apply(lambda x: 1-x)

    return df_neg

'''
Reads and formats the data.
'''
def read_xls(path = None):
    df = pd.read_excel(path)
    df.rename(columns={'Consensus repeat': 'Cons'}, inplace=True)
    df = df.dropna()
    df = df.apply(parse_ID, axis = 1)
    return df

# =============================================================================
# 
# =============================================================================

'''
Finds conflicting samples in less trusted dataset
'''
def compare_datasets(df_reliable_pos, df_corrupt_neg):
    conflict_dict = check_datasets(df_corrupt_neg, df_reliable_pos)
    keys = list(conflict_dict.keys())
    return keys

'''
Cleans conflicting samples
'''
def clean_datasets(df_pos, df_neg, conflict_dict = None, _all = False):
    if conflict_dict == None:
        conflict_dict = check_datasets(df_pos, df_neg)
    keys = list(conflict_dict.keys())
    print('Cleaning the dataset from conflicting samples...')
    df_pos = df_pos[~df_pos['ID'].isin(keys)]
    return seperate_to_pos_neg(df_pos, _all = _all, keep_orig = False)

'''
Finds conflicting samples
'''
def check_datasets(df_pos, df_neg):
    print('Checking the dataset for conflicting samples...')
    conflict_dict = dict()
    for idx, sample in df_pos.iterrows():
        consensus = sample['Cons']
        acc_id = sample['ID']
        df_conf = df_neg[df_neg['Cons'] == consensus]
        if (len(df_conf)>0):
            conflict_dict[acc_id] = (df_conf, consensus)
    return conflict_dict

'''
'''
def issue_conflicts(conflict_dict, path):
    f = open(path, 'w+')
    for conf_acc_id, conf_samples in conflict_dict.items():
        df_conf, consensus = conf_samples
        f.write('Sample '+conf_acc_id+'\t'+consensus+' conflicts with the negative samples down below.\n')
        for _, sample in df_conf.iterrows():
            f.write('\t'.join([str(s) for s in sample]))
            f.write('\n')
        f.write(45*'_'+'\n')
    print('Saved conflicting samples into issue file in path %s' % (path))
    f.close()

# =============================================================================
# FUNCTIONS BELOW PROCESSES DATAFRAME TO FORMAT INPUT
# =============================================================================
'''
Performs operations below to prepare input for neural network.
For input with only repeats.
'''
def process_dataset_for_training(df):
    df = df.copy()
    df['Cons'] = df['Cons'].apply(one_hot_encode_sequence)
    
    # Pad repeats to equal size
    seq_lens = df['Cons'].apply(get_seq_len)
    max_seq_len = max(seq_lens)
    df['Cons'] = df['Cons'].apply(pad_seq, args = [max_seq_len])
    
    # Split into features-labels
    df_X = df['Cons']
    df_y = df['Label']
    
    # Return train and also length of the sequence 
    # (width of the input)
    return df_X.values, df_y.values, max_seq_len

'''
''' 
def process_dataset_for_test(df, max_seq_len, clean = False):
    
    # To remove conflicting samples
    if clean:
        df_pos, df_neg = seperate_to_pos_neg(df)
        print('Number of conflicting samples:', len(check_datasets(df_pos, df_neg)))
        print('Cleaning...')
        df_pos, df_neg = clean_datasets(df_pos, df_neg, con = False)
        print('SUCCESS' if len(check_datasets(df_pos, df_neg)) == 0 else 'Failure')
    
        df = pd.concat([df_pos, df_neg]).reset_index(drop = True)
    
    df_encoded = df.copy()
    # One-hot encode
    df_encoded['Cons'] = df_encoded['Cons'].apply(one_hot_encode_sequence)
    
    # Pad repeats to model input shape
    df_encoded['Cons'] = df_encoded['Cons'].apply(pad_seq, args = [max_seq_len])
    
    # Split into features-labels
    df_X = df_encoded['Cons']
    df_y = df_encoded['Label']
    
    return df, df_X.values, df_y.values

'''
'''
def fasta_to_df(fasta_file_path):
    f= open(fasta_file_path, "r")
    contents = f.read()
    contents = re.split(r'>', contents) # Split every id
    samples = []
    for content in contents[1:]: # First row is an empty partition
        temp = re.split(r'\n', content) # Split from new line or tab
        samples.append(temp[:2])
    df = pd.DataFrame.from_records(data=samples, columns = ['ID', 'Cons'])
    return df

def txt_to_df(path):
    f = open(path, 'r')
    repeats = f.readlines()
    f.close()

    samples = []
    for idx, repeat in enumerate(repeats):
        ID = 'SAMPLE_'+str(idx)
        samples.append((ID, repeat.strip()))

    return pd.DataFrame.from_records(data=samples, columns = ['ID', 'Cons'])
'''
'''

def read_datasets_prepare_train_test(folder, filenames, test_size = 0.25, usecols = None):
    print('Reading the files...')
    df_list = []
    for filename in filenames:
        
        if 'xls' in filename:
            df = read_xls(filename, usecols = usecols, header = 0, names = ['ID', 'Cons'])
        elif '.txt' in filename:
            df = txt_to_df(filename)
        elif 'tsv' or 'tab' in filename:
            df = pd.read_csv(filename, sep = '\t', usecols = usecols, header = 0, names = ['ID', 'Cons'])
        elif '.fasta' in filename:
            df = fasta_to_df(filename)
        
        df_train, df_test = skms.train_test_split(df, test_size = test_size)
        df_list.append((df_train, df_test))

def read_datasets(filenames, test_size = 0.25, usecols = None):
    print('Reading the files...')
    df_list = []
    for filename in filenames:
        if '.xls' in filename:
            df = read_xls(filename, usecols = usecols, header = 0, names = ['ID', 'Cons'])
        elif '.txt' in filename:
            df = txt_to_df(filename)
        elif ('.tsv' in filename) or ('.tab' in filename):
            df = pd.read_csv(filename, sep = '\t', usecols = usecols, header = 0, names = ['ID', 'Cons'])
        elif ('.fasta' in filename) or ('.fa' in filename):
            df = fasta_to_df(filename)
        
        df_list.append(df)
    return df_list

'''
'''
def concat_datasets(df_list):
    if len(df_list[0]) == 2:
        df_train =pd.concat([dfs[0] for dfs in df_list], ignore_index = True)
        df_test =pd.concat([dfs[1] for dfs in df_list], ignore_index = True)
        return df_train, df_test
    else:
        df = pd.concat(df_list, ignore_index = True)
        return df
# =============================================================================


def prepare_input(filenames, output_folder, check_conflict = True, usecols = None, do_type = False):
    df_list = read_datasets(filenames, usecols = usecols)
    if len(df_list) == 1:
        df_pos = df_list[0]
    else:
        df_pos = concat_datasets(df_list)

    #if do_type:
        #df_pos = find_repeat_type(df_pos)

    df_pos, df_neg = seperate_to_pos_neg(df_pos)

    if check_conflict:
        conflict_dict = check_datasets(df_pos, df_neg)
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
        issue_conflicts(conflict_dict, os.path.join(output_folder, 'conflict_issues.txt'))
        df_pos, df_neg = clean_datasets(df_pos, df_neg, conflict_dict)

    df = concat_datasets((df_pos, df_neg))
    return df

if __name__ == '__main__':

    nov_2020_filenames = ['I-E_repeat_seqs.xls', 'VI-B_repeat_seqs.xls', 'VI-A_repeat_seqs.xls','I-F_repeat_seqs.xls', 'II-A_repeat_seqs.xls', 'II-B_repeat_seqs.xls', 'II-C_repeat_seqs.xls']

    df_list = read_datasets(folder, ['clean_repeat_seqs.tsv'], True)
    df_pos_train, df_pos_test = concat_datasets(df_list)

    print(df_pos_train)
    print(df_pos_test)
    df_pos_train.to_csv(os.path.join(folder, 'train_repeat_seqs.tsv'), sep = '\t', index = False)
    df_pos_test.to_csv(os.path.join(folder, 'test_repeat_seqs.tsv'), sep = '\t', index = False)
