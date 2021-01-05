#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 11 19:10:47 2019

@author: ekrem
"""
import os

import numpy as np
import pandas as pd

from sklearn.metrics import roc_auc_score

def test(classifier, df, X, y, output_folder = 'Results', batch_size = 64):
    
    X = np.stack(X)
    probs = classifier.predict(X, batch_size = batch_size)
    df['Probs'] = probs
    
    orig_len = len(df)//2
    samples = []
    for i in range(orig_len):
        pos_sample = df.iloc[i]
        neg_sample = df.iloc[orig_len+i]
    
        ID = pos_sample.ID
        assert ID == neg_sample.ID, 'Positive sample ID does not match the negative.'

        input_cons = pos_sample['Cons']
        sample = pos_sample if pos_sample['Probs'] > neg_sample['Probs'] else neg_sample
        strand = 'Forward' if pos_sample['Probs'] > neg_sample['Probs'] else 'Reverse'
        predicted_cons = sample['Cons']

        prob = sample['Probs']
        if prob >= 0.7:
            confidence = 'High'
        elif 0.5 < prob < 0.7:
            confidence = 'Medium'
        else:
            confidence = 'Low'
        
        if 'Type' in df.columns.values:
            samples.append((ID, input_cons, predicted_cons, strand, confidence, sample.Type))
            columns = ['ID', 'Input Sequence', 'Predicted Sequence', 'Strand', 'Confidence', 'Type']
        else:
            samples.append((ID, input_cons, predicted_cons, strand, confidence))
            columns = ['ID', 'Input Sequence', 'Predicted Sequence', 'Strand', 'Confidence']

    df_out = pd.DataFrame.from_records(data=samples, columns = columns)

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    df_out.to_csv(os.path.join(output_folder, 'CRISPRstrand_Summary.tsv'), sep='\t', index = False)

    return df_out, roc_auc_score(y, probs)
    

